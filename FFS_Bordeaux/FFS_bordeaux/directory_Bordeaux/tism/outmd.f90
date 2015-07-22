!BM----------------------
Module outMD

contains

!BS------------------------------
subroutine outputMD(imd,phasepoint,syst,pot,dyn,out_param)
use types;use assign_objects;use alloc_objects
use wannier
implicit none
integer, intent(in)::iMD
type(phasepoint_type), intent(in)::phasepoint
type(output_type),       intent(in)::out_param 
type(system_type),       intent(in)::syst
type(potential_type),    intent(in)::pot
type(dynamics_type),     intent(in)::dyn
integer::i,IU

call check_cross(imd,phasepoint,syst,out_param,pot)

if (modulo(imd,out_param%skipE)==0) then
  IU=out_param%IUEN
  call write_energies(imd,phasepoint,syst,pot,dyn,IU)
endif

if (modulo(imd,out_param%skipO)==0) then
  call write_orderparameters(imd,phasepoint%phasexv,syst,out_param,pot)
endif

if (modulo(imd,out_param%skipT)==0) then
  write(out_param%IUTRAJ,*) "***", imd
  do i=1,syst%Npart
    write(out_param%IUTRAJ,'(i8,a3,10f16.8)') i,syst%atom_types(i),phasepoint%phasexv%x(i,:), &
                                           phasepoint%phasexv%v(i,:)
  enddo
  if ((pot%POTENTIAL=="EXTERNAL").AND.(syst%NWANNIER>0)) &
      call out_wan(out_param%IUTRAJ,syst%Npart,phasepoint%electrons)
endif

end subroutine outputMD
!ES-------------------------------

!BS-----------------------------------------
subroutine check_cross(imd,phasepoint,syst,out_param,pot)
use types;use assign_objects;use alloc_objects
use orderparameter
implicit none
integer, intent(in)::iMD
type(phasepoint_type),   intent(in)        ::phasepoint
type(system_type),       intent(in)::syst
type(output_type),       intent(in)::out_param
type(potential_type)::pot
logical, save::saveflag 
integer, parameter::ncross=100
logical,save::leftside(ncross)
double precision::op
integer::i

if (ncross < out_param%NCROSSPLANES) then
  print *,"ERROR check_cross:ncross < out_param%NCROSSPLANE"
  stop
endif

op=orderp(phasepoint%phasexv%x,syst,phasepoint%phasexv%v,pot)
if (saveflag) then
  do i=1,out_param%NCROSSPLANES
    if ((op> out_param%crossing_plane(i)).and.(leftside(i))) then
      leftside(i)=.false.
      write(out_param%IUCR,*) imd, i, "+" 
      if (i==2) write(out_param%IUCROSSPOINTS,*) phasepoint%phasexv%x(1,1), phasepoint%phasexv%v(1,1)
    else if ((op< out_param%crossing_plane(i)).and.(.not.leftside(i))) then
      leftside(i)=.true.
      write(out_param%IUCR,*) imd, i, "-"
    endif
  enddo
else
  leftside(1:out_param%NCROSSPLANES)=.false.
  do i=1,out_param%NCROSSPLANES
    if (op<out_param%crossing_plane(i)) leftside(i)=.true.
  enddo
  saveflag=.true.
endif
  

end subroutine check_cross
!ES-----------------------------------------

!BS--------------------------------
subroutine write_energies(imd,phasepoint,syst,potential,dyn,iU)
use types;use assign_objects;use alloc_objects
use forcefield
use orderparameter !!hust for the test
implicit none
integer, intent(in)::iMD
type(phasepoint_type),   intent(in)        ::phasepoint
type(system_type),       intent(in)        ::syst
type(potential_type),    intent(in)        ::potential
type(dynamics_type),     intent(in)        ::dyn
integer,                 intent(in)        ::IU
type(NOSE_HOOVER_type)::NHparam
double precision::KIN,POT,ETOT,EHAM,TEMP_inst
!character::char
double precision::Eelectrons
integer::i

if (potential%potential=="EXTERNAL") then
  KIN=phasepoint%electrons%KIN
  POT=phasepoint%electrons%POT
  ETOT=phasepoint%electrons%ETOT
  EHAM=phasepoint%electrons%EHAM
  TEMP_inst=phasepoint%electrons%TEMP_inst
  Eelectrons=phasepoint%electrons%Eelectrons
  !test
  print *,"KIN-CPMD,KIN(v)",real(KIN),real(EKIN(phasepoint%phasexv%v,syst%masses_2,syst%Npart,syst%dim))
  print *,"TEM-CPMD,TEM(v)",real(TEMP_inst),real(EKIN(phasepoint%phasexv%v,syst%masses_2,syst%Npart,syst%dim)*syst%x2_kbNFtrans)
  

else
  KIN=EKIN(phasepoint%phasexv%v,syst%masses_2,syst%Npart,syst%dim)
  POT=Epot(phasepoint%phasexv%x,syst,potential)
  ETOT=KIN+POT
  Temp_inst=syst%x2_kbNF*KIN
  if (dyn%DYNAMICS=="NOSEHOOVER") then
    NHparam=dyn%Nose_Hoover
    EHAM=ETOT+E_NOSE_THERMOS(phasepoint%extcoord,NHparam,syst)
  else  
    EHAM=ETOT
  endif
  Eelectrons=0.d0
endif

write(IU,'(i12,4E25.10,f15.4,E25.10)') iMD,POT,KIN,ETOT,EHAM,Temp_inst,& 
                                       Eelectrons
end subroutine write_energies
!ES--------------------------------

!BS--------------------------------
subroutine write_orderparameters(imd,phasepoint,syst,out_param,pot) 
use types;use assign_objects;use alloc_objects
use orderparameter
use arimetric
use forcefieldwca
use forcefieldiot
implicit none
integer, intent(in)::iMD
type(phasexv_type), intent(in)::phasepoint
type(system_type),       intent(in)::syst
type(output_type),       intent(in)::out_param
type(potential_type),          intent(in)::pot
character(LEN=MSTR)::POTENTIAL
integer::IU
double precision::OP_dp_array(50)
integer::OP_int_array(10),NLdp,NLint
double precision::xi(syst%dim),xj(syst%dim),vi(syst%dim),vj(syst%dim)
double precision::dv(syst%dim),dr(syst%dim),absr,absv,dabsr_dt
integer::i,j,k,iSi,jO,O_array(8)
double precision::NN,NC,Fcoop
double precision::eps,sig,eps2,sig2,r0,r02,Vshift,Rij


POTENTIAL=pot%POTENTIAL
IU=out_param%IUOP
NLdp=0;NLint=0
select case (POTENTIAL)
  case("HARMOSC") 
    NLdp=2
    OP_dp_array(1)=phasepoint%x(1,1)
    OP_dp_array(2)=phasepoint%v(1,1)
  case("DOUBLEWELL")
    NLdp=2
    OP_dp_array(1)=phasepoint%x(1,1)
    OP_dp_array(2)=phasepoint%v(1,1)
  case("2DHARM")    
    NLdp=4    
    OP_dp_array(1)=phasepoint%x(1,1)    
    OP_dp_array(2)=phasepoint%x(1,2)
    OP_dp_array(3)=phasepoint%v(1,1)        
    OP_dp_array(4)=phasepoint%v(1,2)
  case("PBD")
    NLdp=1
    OP_dp_array(1)=minval(phasepoint%x(:,1))
  case("NICOLIS1D")
    NLdp=2
    OP_dp_array(1)=phasepoint%x(1,1)
    OP_dp_array(2)=phasepoint%v(1,1)
  case("NICOLIS2D")
    NLdp=4
    OP_dp_array(1)=phasepoint%x(1,1)
    OP_dp_array(2)=phasepoint%x(1,2)
    OP_dp_array(3)=phasepoint%v(1,1)
    OP_dp_array(4)=phasepoint%v(1,2)
  case("EXTERNAL")
    if (syst%REACTION_COORDINATE=="SILICWAT64") then
      NLdp=9
      k=0
      do i=1,2
        iSi=208+i
        do j=1,4
          k=k+1
          jO=145-k
          xi(1:3)=phasepoint%x(iSi,1:3)
          xj(1:3)=phasepoint%x(jO,1:3)
          OP_dp_array(k+1)=distance(xi,xj,3,syst%BOXLENGTH)
        enddo
      enddo
      OP_dp_array(1)=minval(OP_dp_array(2:9))
    else if (syst%REACTION_COORDINATE=="SILICWAT64ANION") then  
      NLdp=9
      k=0
      O_array=(/ 137,138,141,142,136,139,140,143 /)
      do i=1,2
        iSi=207+i
        do j=1,4
          k=k+1
          jO=O_array(k)
          xi(1:3)=phasepoint%x(iSi,1:3)
          xj(1:3)=phasepoint%x(jO,1:3)
          OP_dp_array(k+1)=distance(xi,xj,3,syst%BOXLENGTH)
          !!print *,"hhh",i,j
        enddo
      enddo
      OP_dp_array(1)=minval(OP_dp_array(2:9)) 
      print *,"!! check", OP_dp_array(1),orderp(phasepoint%x,syst,phasepoint%v,pot)
    else
      NLdp=1
      OP_dp_array(1)=orderp(phasepoint%x,syst,phasepoint%v,pot)
    endif
  case("NONE")
    NLdp=1
    OP_dp_array(1)=orderp(phasepoint%x,syst,phasepoint%v,pot)
  case("WCA")
    NLdp=9
    xi(:)=phasepoint%x(1,:)
    xj(:)=phasepoint%x(2,:)
    vi(:)=phasepoint%v(1,:)
    vj(:)=phasepoint%v(2,:)
    dr(:)=xi(:)-xj(:)
    dv(:)=vi(:)-vj(:)
    absr=sqrt(dot_product(dr,dr))
    absv=sqrt(dot_product(dv,dv))
    dabsr_dt=dot_product(dr,dv)/absr
    !!if (syst%icrash==457) then
      !!print *,"xi",xi
      !!print *,"xj",xj
     !! print *,"vi",vi
     !! print *,"vj",vj
     !! print *,"dr",dr
     !! print *,"dv",dv
     !! print *,"absr",absr
     !! print *,"absv",absv
     !! print *,"dabsr_dt",dabsr_dt

    !!  stop
    !!endif
    
    OP_dp_array(1)=orderp(phasepoint%x,syst,phasepoint%v,pot)
    OP_dp_array(2)=absr
    OP_dp_array(3)=absv
    OP_dp_array(4)=dabsr_dt 
    OP_dp_array(5)=0.25D0*syst%masses(1)*absv**2
    OP_dp_array(6)=0.25D0*syst%masses(1)*dabsr_dt**2
    OP_dp_array(7)=Vdw(absr,pot%potWCA%h,pot%potWCA%w,pot%potWCA%r0)
    OP_dp_array(8)=OP_dp_array(5)+OP_dp_array(7)
    OP_dp_array(9)=OP_dp_array(6)+OP_dp_array(7)
    !!write(111,*)  imd,syst%masses(1)*(vi(:)+vj(:))
  case("IONTRANS")
    NLdp=11
    OP_dp_array(1)=orderp(phasepoint%x,syst,phasepoint%v,pot)
    xi(:)=phasepoint%x(1,:)
    xj(:)=phasepoint%x(2,:)
    OP_dp_array(2)= distance(xi,xj,syst%dim,syst%boxlength)
    xj(:)=phasepoint%x(3,:)
    OP_dp_array(3)= distance(xi,xj,syst%dim,syst%boxlength)
    xi(:)=phasepoint%x(2,:)
    OP_dp_array(4)= distance(xi,xj,syst%dim,syst%boxlength)
    OP_dp_array(5)= Coordination_number(phasepoint%x,syst,pot%potiot)
    NN=pot%potIOT%NN
    NC=pot%potIOT%NC
    Fcoop=pot%potIOT%Fcoop
    OP_dp_array(6)=FCOORD(OP_dp_array(5),Fcoop,NC,NN)

    eps=pot%potIOT%eps
    sig=pot%potIOT%sig
    eps2=pot%potIOT%eps2
    sig2=pot%potIOT%sig2
    r0=pot%potIOT%r0
    r02=pot%potIOT%r02
    Vshift=pot%potIOT%Vshift
    OP_dp_array(7)=VpairIOT(OP_dp_array(2),1,2,eps,sig,eps2,sig2,r0,r02,Vshift)
    OP_dp_array(8)=VpairIOT(OP_dp_array(3),1,3,eps,sig,eps2,sig2,r0,r02,Vshift)
    OP_dp_array(9)=VpairIOT(OP_dp_array(4),2,3,eps,sig,eps2,sig2,r0,r02,Vshift)
    OP_dp_array(10)=OP_dp_array(7)+OP_dp_array(8)+OP_dp_array(9)
    OP_dp_array(11)=Epotiot(phasepoint%x,syst,pot%potiot)
  case default
    print *,"ERROR write_orderparameters POTENTIAL=",POTENTIAL
    stop
end select

write(IU,'(i12,50f15.8)') iMD,OP_dp_array(1:NLdp),OP_int_array(1:NLint)
end subroutine write_orderparameters
!ES--------------------------------

end Module outMD
!EM----------------------
