!BM------------------------
module mdintegrators

contains
!BS----------------------------------------------------------------
! To start a MD loop with the velocity Verlet algorithm, the force
! of the start configuration has to be given once, before calling
! this subroutine. It makes first a half-timestep-update with
! the velocities, then a full timestep update positions, recalcualte
! new force and updates another half-timestep the velocities
! New positions, velocities, and forces are given as output
!------------------------------------------------------------------
subroutine velocity_Verlet(phasexv,syst,pot,F,timestep)
use types;use assign_objects;use alloc_objects
use forcefield
implicit none
type(phasexv_type),intent(inout)::phasexv
type(system_type),      intent(in)::syst
double precision,       intent(inout)::F(syst%Npart,syst%dim)
type(potential_type),   intent(in)::pot
type(timestep_type),    intent(in)   ::timestep
double precision::dt,dt_2m(syst%Npart)
double precision:: V_2(syst%Npart,syst%dim)
integer::i,j,dim,N


dt=timestep%dt
dt_2m=timestep%dt_2m
N=syst%Npart
dim=syst%dim
!dt_2m=0.5d0*dt/m
  do i=1,N
    do j=1,dim
      V_2(i,j)=phasexv%v(i,j) + dt_2m(i)*F(i,j)
    enddo
  enddo
 
  phasexv%x=phasexv%x+dt*V_2
  F=FORCE(phasexv%x,syst,pot)

  do i=1,N
    do j=1,dim
      phasexv%v(i,j)=V_2(i,j) + dt_2m(i)*F(i,j)
    enddo
  enddo
end subroutine velocity_Verlet
!ES----------------------------------------------------------------

!BS--------------------------------------------------------------
subroutine Andersen_velocity_Verlet(phasexv,syst,pot,F,Andersen_param,timestep)
use types;use assign_objects;use alloc_objects
use thermo
implicit none
type(phasexv_type),intent(inout)::phasexv
type(system_type),      intent(in)::syst
double precision,       intent(inout)::F(syst%Npart,syst%dim)
type(potential_type),   intent(in)::pot
type(Andersen_type),    intent(in)::Andersen_param
type(timestep_type),    intent(in)::timestep
double precision::sigma_v(syst%Npart),freqdt
integer::N,dim

N=syst%Npart
dim=syst%dim
sigma_v=syst%sigma_v
freqdt=Andersen_param%freqdt

call velocity_Verlet(phasexv,syst,pot,F,timestep)
call Andersen_velocity_change(phasexv%v,N,dim,sigma_v,freqdt)

end subroutine Andersen_velocity_Verlet
!ES--------------------------------------------------------------

!BS------------------------------------------------------------------------
subroutine Lange_Verlet(phasexv,syst,pot,F,Langevin_param)
use types;use assign_objects;use alloc_objects
use random
use forcefield
implicit none
type(phasexv_type),intent(inout)::phasexv
type(system_type),      intent(in)   ::syst
double precision,       intent(inout)::F(syst%Npart,syst%dim)
type(Langevin_type),     intent(in)::Langevin_param
type(potential_type),   intent(in)   ::pot

if (Langevin_param%high_friction_limit) then
  call Lange_Verlet_overdamped(phasexv,syst,pot,F,Langevin_param)
else
  call Lange_Verlet_inertia(phasexv,syst,pot,F,Langevin_param)
endif 

end subroutine Lange_Verlet
!ES------------------------------------------------------------------------


!BS------------------------------------------------------------------------
subroutine Lange_Verlet_overdamped(phasexv,syst,pot,F,Langevin_param)
use types;use assign_objects;use alloc_objects
use random
use forcefield
implicit none
type(phasexv_type),intent(inout)::phasexv
type(system_type),      intent(in)   ::syst
double precision,       intent(inout)::F(syst%Npart,syst%dim)
type(Langevin_type),     intent(in)::Langevin_param
type(potential_type),   intent(in)   ::pot
integer::i,j,N,dim
double precision::bDdt,sig,dt,ran

N=syst%Npart
dim=syst%dim

F=FORCE(phasexv%x,syst,pot)
do i=1,N
  bDdt=Langevin_param%bDdt(i)
  sig=Langevin_param%sigma_lange(i)
  do j=1,dim 
     ran=rangaussian(sig)
     phasexv%x(i,j)=phasexv%x(i,j)+bDdt*F(i,j)+ran
     phasexv%v(i,j)=ran
     !as velocities are not used, it can still be handy to keep track
     !of the "random displacement velocity"
  enddo
enddo

end subroutine Lange_Verlet_overdamped
!ES------------------------------------------------------------------------


!BS------------------------------------------------------------------------
subroutine Lange_Verlet_inertia(phasexv,syst,pot,F,Langevin_param) 
use types;use assign_objects;use alloc_objects
use random
use forcefield
implicit none
type(phasexv_type),intent(inout)::phasexv
type(system_type),      intent(in)   ::syst
double precision,       intent(inout)::F(syst%Npart,syst%dim)
type(Langevin_type),     intent(in)::Langevin_param
type(potential_type),   intent(in)   ::pot
double precision::randx,randv
integer::i,j,N,dim
double precision::s12os11,sqrts11,sqrtSos11,c0,a1,a2,b1,b2
double precision::V_2(syst%Npart,syst%dim)

N=syst%Npart
dim=syst%dim
c0=Langevin_param%c0
a1=Langevin_param%a1
do i=1,N
   s12os11=Langevin_param%s12os11(i)
   sqrts11=Langevin_param%sqrts11(i)
   sqrtSos11=Langevin_param%sqrtSos11(i)
   a2=Langevin_param%a2(i)
   b1=Langevin_param%b1(i)
   
   do j=1,dim
     if (Langevin_param%gamma > 0 ) then
       call gssbivar(randx,randv,s12os11,sqrts11,sqrtSos11)
     else
       randx=0.d0;randv=0.d0
     endif
     phasexv%x(i,j)=phasexv%x(i,j)+a1*phasexv%v(i,j)+a2*F(i,j)+randx
     V_2(i,j)=c0*phasexv%v(i,j)+b1*F(i,j)+randv
   enddo
enddo
F=FORCE(phasexv%x,syst,pot)
do i=1,N
  b2=Langevin_param%b2(i)
  do j=1,dim
     phasexv%v(i,j)=V_2(i,j)+b2*F(i,j)
  enddo
enddo

end subroutine Lange_Verlet_inertia
!ES------------------------------------------------------------------------

!BS--------------------------------------------------------------------------
subroutine Nose_Hoover_step(phasepoint,syst,pot,NH_param,KIN,timestep)
use types;use assign_objects;use alloc_objects
use thermo
use forcefield
implicit none
type(phasepoint_type),intent(inout)  ::phasepoint
type(system_type),      intent(in)   ::syst
type(Nose_Hoover_type),    intent(in)::NH_param
type(potential_type),   intent(in)   ::pot
double precision, intent(inout)      ::KIN
type(timestep_type),    intent(in)::timestep
double precision::masses_2(syst%Npart)
integer::N,dim
double precision::V(syst%Npart,syst%dim),X(syst%Npart,syst%dim)

N=syst%Npart
dim=syst%dim
masses_2=syst%masses_2
X=phasepoint%phasexv%x
V=phasepoint%phasexv%v
call chainNose(V,phasepoint%extcoord,N,dim,KIN,NH_param,syst,timestep)
call pos_vel_Nose(X,V,syst,pot,timestep)
KIN=EKIN(V,masses_2,N,dim)
call chainNose(V,phasepoint%extcoord,N,dim,KIN,NH_param,syst,timestep)
phasepoint%phasexv%x=X
phasepoint%phasexv%v=V
end subroutine Nose_Hoover_step
!ES--------------------------------------------------------------------------

!BS---------------------------------------------------------------------------
subroutine externalMD(phasep,syst,ext_dyn,dt)
use cpmd_subr
use stringlengths
use shell
use types;use assign_objects;use alloc_objects
implicit none
type(phasepoint_type),intent(inout)  ::phasep
type(system_type),      intent(inout)   ::syst
type(external_dynamics_type),      intent(in)::ext_dyn
double precision, intent(in)::dt
integer::twait
logical::wait,RSEXIST
integer::NOMORE,EFILE_INDEX,NX
character(LEN=LSTR)::RS,RE,FILEBASE
character(LEN=LSTR)::RS2,RE2
character(LEN=LSTR)::command
logical::noprint
double precision::timestep
integer::tot_waiting_time,maxwait

print *,"ffff",syst%icrash,syst%NCPMD_UNSAVED
!!if (syst%icrash < syst%NCPMD_UNSAVED) then
!!  !No need to call CPMD
!!  !Just read geometry, energies and wannier centers from CRASHFILE
!!  print *,"GET GEOMETRY, ENERGIES FROM CRASH-FILE"
!!  call read_geo_crashfile(syst,phasep)
!!else

  timestep=dt
  if (phasep%electrons%rev_velec) then
    !The electron velocities in the electronic
    !configuration file are flipped compared
    !to the correct ionic velocities
    !Instead of regenerating a new electronic binary file,
    !we apply another simple solution. We flip the ionic
    !velocities instead and supply a negative timestep to cpmd.
    !Then we reverse the velocities again.
    phasep%phasexv%v=-phasep%phasexv%v
    timestep=-timestep
  endif
  !write coordinates into GEOMETRY FILE
  call make_GEOFILE(phasep%phasexv%x,phasep%phasexv%v,syst%Npart, & 
                    syst%inv_a0,syst%inv_atomicv_redv)

  twait=syst%twait
  maxwait=syst%maxwait
  NOMORE=syst%NSUBCYCLES
  NX=syst%NX
  FILEBASE=phasep%electrons%ESTRUCFILE
  EFILE_INDEX=phasep%electrons%EFILE_INDEX
  RS=trim(FILEBASE)//ext_dyn%EXT(EFILE_INDEX)
  !!if (ext_dyn%RENAME_ESTRUCTURE_FILE) then
    if (dt>0) then 
      EFILE_INDEX=EFILE_INDEX+1
    else
      EFILE_INDEX=EFILE_INDEX-1
    endif
  !!endif
  if (EFILE_INDEX>=NX) EFILE_INDEX=EFILE_INDEX-NX
  RE=trim(FILEBASE)//ext_dyn%EXT(EFILE_INDEX)
  phasep%electrons%EFILE_INDEX=EFILE_INDEX

  RS2=RS
  RE2=RE
  if (FILEBASE(4:4)=="/") then
    !RESTART file is situated in a higher directory: 0../ESTRUC_...
    !CPMD cannot read filenames that include "/"
    RS2=RS(5:);RE2=RE(5:)
    command="mv -f "//TRIM(RS)//" "//TRIM(RS2)
    call do_shell_cpmd(command,noprint,syst%icrash,syst%NCPMD_UNSAVED)
  endif

  !!INQUIRE(FILE=RS2,EXIST=RSEXIST)
  !!if (.NOT.RSEXIST) THEN
  !!  print *,"FOLLOWING RESTART FILE NOT FOUND: ",RS2
  !!  print *,"REGENERATING NEW RESTART FILE BY OPTIMIZE WAVEFUNCTION"
  !!  print *,"RESTART FILE SHOULD BE REGENERATED FROM IONIC POSITIONS"
  !!  print *,"SORRY,THIS OPTION IS NOT YET IMPLEMENTED"
  !!  stop
  !!endif

  call writePSAMPLE(NOMORE,timestep,RS2,RE2,1,NOMORE,&
                       syst%ICRASH,syst%NCPMD_UNSAVED)
  wait=.true.
  !CPMD PROGRAM WILL READ INFORMATION FROM "PSAMPLE"
  !TO PERFORM A MD TRAJECTORY. WHEN IT IS FINISHED
  !"PSAMPLE" WILL BE DELETED

  tot_waiting_time=0
  do while(wait)
   call sleep(twait)
   !print *,"call sleep() blocked, please remove !!";stop
   INQUIRE(FILE="PSAMPLE",EXIST=wait)
   !IF EXIST:CPMD STILL BUSY
   tot_waiting_time=tot_waiting_time+twait
   if (tot_waiting_time>maxwait)  then
       print *,"**************************************************************"
       print *,"EXCEEDED MAX WAITING TIME"
       print *,"CPMD MIGHT HAVE CRASHED"
     if (syst%RESUBMIT) then
       !!call do_shell("cp -fr GEOMETRY GEOMETRY00")
       print *,"**************************************************************"
       print *,"**********  RESUBMIT TISMOL **********************************"
       print *,"**************************************************************"
       call do_shell("./RESUBMIT")
     endif
     stop
   endif
   
  enddo
  !CPMD HAS FINISHED ITS TRAJECTORY
  !Read new coordinates/velocities
  call read_cpmd_coordinates("GEOMETRY",syst,phasep=phasep)


  command="mv -f "//trim(RE2)//"1"//" "//trim(RE)     
  !CPMD patches a 1 to the RE file that we need to get rid off
  call do_shell_cpmd(command,noprint,syst%icrash,syst%NCPMD_UNSAVED)                      
  if (FILEBASE(4:4)=="/") then
     command="mv -f "//TRIM(RS2)//" "//TRIM(RS)
     call do_shell_cpmd(command,noprint,syst%icrash,syst%NCPMD_UNSAVED)
     !move the start file back as well
  endif


  if (phasep%electrons%rev_velec) then
    !Now we reverse the velocities again.
    phasep%phasexv%v=-phasep%phasexv%v
  endif

!!endif

!!syst%icrash=syst%icrash+1
!!if (syst%icrash > syst%NCPMD_UNSAVED) then
!!  print *,"ADD GEOMETRY, ENERGIES TO CRASH-FILE"
!!  call print_geo_crashfile(syst,phasep)
!!endif


end subroutine externalMD
!ES---------------------------------------------------------------------------

!BS--------------------------------------------------------------------------
subroutine make_randomwalk(x)
use random
implicit none
double precision, intent(inout)::x
double precision::ran
!!call random_number(ran)
ran=random01()
print *,"mdintegrators RANDOM_NUMBER:",ran
if (ran>.5d0) then
  x=x+1.d0
else
  x=x-1.d0
endif

end subroutine make_randomwalk
!ES--------------------------------------------------------------------------
end module mdintegrators
!EM------------------------
