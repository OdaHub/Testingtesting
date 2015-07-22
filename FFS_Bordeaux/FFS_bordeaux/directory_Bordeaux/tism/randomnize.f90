!BM------------------------
module randomnize 

contains

!BS------------------------------------------------------------------------
! This subroutine makes a random modification to the timeslice. For now,
! this implies that only the velocities of the phasexv point are changed. 
! sigdp_sqrtm gives how strong the change is. If aimless=.true., the
! velocities are completely regenerated from a Maxwellian distribution.
! For dynamics="EXTERNAL", the electronicr-structure file and index are 
! renamed.
!BS------------------------------------------------------------------------
subroutine kick_timeslice(timeslice,syst,pot,dyn,sigdp_sqrtm,aimless,dK,Eold,Enew)
use types;use assign_objects;use alloc_objects
use mdstep
use forcefield
use shell
use cpmd_subr
use random
use orderparameter
implicit none
type(timeslice_type),  intent(inout)::timeslice
type(system_type),        intent(inout)   ::syst
type(potential_type),        intent(in)   ::pot
type(dynamics_type),      intent(in)   ::dyn
logical, intent(in)::aimless
double precision, intent(in)::sigdp_sqrtm(syst%Npart)
double precision, optional, intent(out)::dK
character(LEN=LSTR), optional, intent(out)::Eold,Enew
double precision:: vel(syst%Npart,syst%dim),sigma_v(syst%Npart)
double precision::kin_old, kin_new
integer::N,d,Lfilenm
character(LEN=LSTR)::Efile1,Efile2,command
integer::efile_index1, efile_index2
logical::noprint
double precision:: vel12(2,syst%dim),m12_2(2),velno12(syst%Npart-2,syst%dim)
double precision::mno12_2(syst%Npart-2),KIN12,KINno12,op,EP
double precision::mno12(syst%Npart-2),P12(syst%dim)

!extract some info from the timeslice
N=timeslice%N
d=timeslice%d
sigma_v=syst%sigma_v
!extract the velocities
vel=timeslice%phasepoint%phasexv%v

print *,"kick op before:", timeslice%OPS(1)

!get old kinetic energy
if ((dyn%DYNAMICS=="NVE").AND.(dyn%NVE%RESCALE_ENERGY)) then
  EP=Epot(timeslice%phasepoint%phasexv%x,syst,pot)
  KIN_old=dyn%NVE%NVE_energy-EP
  !slightly more accurate than kin_old=ekin(...)
  !velocity rescaling will be done to this energy 
else if (dyn%DYNAMICS=="NOSEHOOVER") then
  kin_old=timeslice%MDinout_param%kin
else if (dyn%DYNAMICS=="EXTERNAL") then
  kin_old=timeslice%phasepoint%electrons%KIN
else
  kin_old=ekin(vel,syst%masses_2,syst%Npart,syst%dim)
endif


!Make a velocity change
if (aimless) then  !velocities regenerated from Maxwellian distribution
    call set_maxwellian_velocities(vel,N,d,sigma_v)
else !Velocity change by adding small random velocity displacements
  call soft_velocity_change(vel,N,d,sigdp_sqrtm) 
endif


!take care of linear momentum and energy conservation if needed
if (dyn%DYNAMICS=="NVE") then
  !if linear momentum is conserved, we to shift the velocities
  if (dyn%NVE%SET_LINMOM_ZERO) then
    call set_linmom(vel,syst%masses,N,d)
  endif

  if (dyn%NVE%RESCALE_ENERGY) then
    !if energy is conserved, we rescale the velocities
    call rescale(vel,syst%masses_2,N,d,KIN_old)
  endif
endif


!update some v-dependent properties of the timeslice
if (syst%REACTION_COORDINATE=="WCAJCP1") timeslice%OPS(1)=orderp(timeslice%phasepoint%phasexv%x,syst,vel,pot)



if (dyn%DYNAMICS=="NOSEHOOVER") then
  timeslice%MDinout_param%kin= ekin(vel,syst%masses_2,syst%Npart,syst%dim)
else if (dyn%DYNAMICS=="EXTERNAL") then
  kin_new=ekin(vel,syst%masses_2,syst%Npart,syst%dim)
  timeslice%phasepoint%electrons%kin=kin_new
  timeslice%phasepoint%electrons%Temp_inst=timeslice%phasepoint%electrons%kin*syst%x2_kbNFtrans
  timeslice%phasepoint%electrons%ETOT=timeslice%phasepoint%electrons%ETOT+kin_new-kin_old
  timeslice%phasepoint%electrons%EHAM=timeslice%phasepoint%electrons%EHAM+kin_new-kin_old
endif
if (present(dK)) dK=kin_new-kin_old





!Save the new velocities back into the timeslice
timeslice%phasepoint%phasexv%v=vel

!For the "external" option, we have to do some more things
if (dyn%DYNAMICS=="EXTERNAL") then
  !Now take care of the electrons
  Efile1=timeslice%phasepoint%electrons%ESTRUCFILE
  Lfilenm=LEN(TRIM(Efile1))
  Efile2=Efile1
  Efile_index1=timeslice%phasepoint%electrons%EFILE_INDEX
  Efile_index2=0
  !Change the _a extension into _b and vice versa
  if (Efile1(Lfilenm:Lfilenm)=="a") then
      Efile2(Lfilenm:Lfilenm)="b"
  else
      Efile2(Lfilenm:Lfilenm)="a"
  endif
  !save this into timeslice
  timeslice%phasepoint%electrons%ESTRUCFILE=Efile2
  !set EFILE_INDEX to efile_index2 
  timeslice%phasepoint%electrons%EFILE_INDEX=efile_index2
  !cp previous EFILE using the new name
    print *,"RENAME ESTRUCFILE WHILE MAKING TRIAL"
    command="mv -f "//trim(efile1)//dyn%ext_dyn%ext(efile_index1)//" "// &
                 trim(efile2)//dyn%ext_dyn%ext(efile_index2)
    !cp is too expensive as CPMD RESTART file can be in the order og 1 GB
    !Therefore we rename it. 
    !Beware to rename it back if the MC is rejected!
    call do_shell_cpmd(command,noprint,syst%icrash,syst%NCPMD_UNSAVED)
  if (present(Eold)) then
    Eold=trim(efile1)//dyn%ext_dyn%ext(efile_index1)
    Enew=trim(efile2)//dyn%ext_dyn%ext(efile_index2)
  endif

endif

print *,"kick op after:", timeslice%OPS(1)

end subroutine kick_timeslice 
!ES----------------------------------------------------------------

!BS----------------------------------------------------------------
subroutine Metropolis_momenta_change(dK,beta,aimless,op,lb,ub,accmom)
use random
implicit none
double precision, intent(in)::dK,beta
logical, intent(in)::aimless
logical, intent(out)::accmom
double precision, intent(in)::op,lb,ub
double precision::prob,ran

if ((op < lb).OR.(op>ub)) then
  accmom=.false.
else if (aimless) then
  accmom=.true.
else
  if (dK < 0) then
    accmom=.true.
  else
    prob=exp(-beta*dK)
    !!call random_number(ran)
    ran=random01()
    if (ran<prob) then
      accmom=.true.
    else
      accmom=.false.
    endif
  endif
endif

end subroutine Metropolis_momenta_change 
!ES----------------------------------------------------------------

!BS------------------------------------------------------------------
subroutine soft_velocity_change(vel,N,d,sigdp_sqrtm)
use random
implicit none
integer, intent(in)::N,d
double precision, intent(in)::sigdp_sqrtm(N)
double precision, intent(inout)::vel(N,d)
integer::i,j
double precision::sdp,ran

do i=1,N
 sdp=sigdp_sqrtm(i)
 do j=1,d
   vel(i,j)=vel(i,j)+rangaussian(sdp) 
 enddo
enddo

end subroutine soft_velocity_change
!ES------------------------------------------------------------------

 !BS-------------------------------------------------------
  subroutine randomnize_v(v,syst,dyn,x,pot,potE)
  use types;use assign_objects;use alloc_objects
  use forcefield
  use random
  implicit none
  type(system_type),  intent(in   )::syst
  type(dynamics_type),intent(in   )::dyn
  double precision,intent(inout)::v(syst%Npart,syst%dim)
  double precision,   optional, intent(in)::x(syst%Npart,syst%dim)
  type(potential_type),optional, intent(in)::pot
  double precision,   optional, intent(in)::potE
  integer::Npart,dim
  double precision::sigma_v(syst%Npart)
  double precision::EP,EK

  Npart=syst%Npart
  dim=syst%dim
  sigma_v=syst%sigma_v

  call set_maxwellian_velocities(v,Npart,dim,sigma_v)
  if (dyn%DYNAMICS=="NVE") then
    if (dyn%NVE%SET_LINMOM_ZERO) then
     call set_linmom(v,syst%masses,syst%Npart,syst%dim)
    endif

    if (dyn%NVE%RESCALE_ENERGY) then
      if (present(potE)) EP=PotE
      if (present(x)) EP=Epot(x,syst,pot)
      EK=dyn%NVE%NVE_ENERGY-EP
      if (EK<0) then
        print *,"ERROR randomnize_v, EP TOO HIGH, VELOCITIES CAN'T BE RESCALED"
        stop
      endif
      call rescale(v,syst%masses_2,syst%Npart,syst%dim,EK)
    endif
  endif

  end subroutine randomnize_v
  !ES-------------------------------------------------------

  !BS----------------------------------------------------------
  subroutine rescale(v,m_2,N,d,EK)
  use forcefield
  implicit none
  integer, intent(in)::N,d
  double precision,intent(inout)::v(N,d)
  double precision,intent(in)::m_2(N) !masses/2
  double precision,intent(in)::EK
  double precision::EKo,factor

  EKo=EKIN(v,m_2,N,d)
  factor=sqrt(EK/EKo)
  v=v*factor
  end subroutine rescale
  !ES-----------------------------------------------------------

  !BS-----------------------------------------------------------
  SUBROUTINE set_linmom(v,masses,N,dim,Pin)
  implicit none
  integer, intent(in)::N,dim
  double precision, intent(inout)::V(N,dim)
  double precision, intent(in)::masses(N)
  double precision, optional, intent(in)::Pin(dim)
  double precision::totmass, totP(dim),dv(dim),Pset(dim)
  integer::i

  Pset=0.d0
  if (present(Pin)) Pset=Pin 
  totP=0.d0 
  totmass=0.d0
  do i=1,N
   totP(:)=totP(:)+masses(i)*V(i,:)
   totmass=totmass+masses(i)
  enddo
  dv(:)=(totP(:)-Pset(:))/totmass
  do i=1,N
   v(i,:)=v(i,:)-dv(:)
  enddo
    
  END SUBROUTINE set_linmom
  !ES------------------------------------------------------------



end module randomnize 
!EM------------------------
