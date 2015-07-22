!-----------------------------------------------------------------------
! Program TISMOL 
! 20-10-2006 Start of project. 
! This Program employs the TIS (Transition Interface Sampling) algorithm
! for rare events. In addition some other type of simulations can be 
! performed as well (standard MD, reactive flux approach)
! Several potential/systems have been implemented or will 
! be implemented in the near future.
! These are the 1D Peyrard-Bishop model for DNA, The 1D and 2D
! Nicolis potential, The 1D Harmonic oscilator. 
! Several types of dynamics can be used: Newtonian (NVE) MD, Langevin, 
! Brownian, Nose-hoover thermostatted MD, Andersen thermostatted MD.
! Parallel Path Swapping is included as well. 
! Not yet implemented are: PPTIS, FFS, Reaction Coordinate Free Sampling
! Implementation to combine TIS with CPMD is in development
!---------------------------------------------------------------------------
!BP-------------- TISMOL  -------------------------------------------
!--------------------------------------------------------------------------
! Beginning of Program
!--------------------------------------------------------------------------
Program TISMOL

!use statements
use parameters
use init
use exec
use close

!no implicit
implicit none

!definitions of variables
character*8::date
character*10::time
double precision::cpu

  open(22,file="cpu.dat",position="append")
  !File to keep track of the cpu cost

  call date_and_time(date,time)
  !!print *,"DATE:"//trim(date)//",TIME:"//trim(time)
  print *,"DATE:",date,",TIME:",time
  write(22,'(2a8,2a10)') "start:",trim(date)," time:",trim(time)

  print *,"STARTING TISMOL"
  print *,"----------------"
  !Reading input parameters and perform initializations
  print *, "READING INPUT/DEFAULT VARIABLES" 
  call setparam
  
  print *,"INITIALIZE"
  call initialize

  call date_and_time(date,time)
  write(22,'(2a8,2a10)') "after initialization:",trim(date)," time:",trim(time) 
  !Write the time after initialization to keep track of 
  !initialization/execution time


  print *, "EXECUTE"
  call execute

  print *,"----------------"
  print *,"FINISHED"
  call date_and_time(date,time)
  write(22,'(2a8,2a10)') "close:",trim(date)," time:",trim(time) 

  print *,"CLOSURE"
  call closure

  call date_and_time(date,time)

  print *,"DATE: ", date, ", TIME: ",time
  write(22,'(2a8,2a10)') "end:",trim(date)," time:",trim(time)

  call CPU_TIME(cpu)

  print *,"TOTAL COMPUTATION TIME:"
  print *,cpu," SECONDS"
  write(22,*) "cpu",cpu
  close(22)

end Program TISMOL

