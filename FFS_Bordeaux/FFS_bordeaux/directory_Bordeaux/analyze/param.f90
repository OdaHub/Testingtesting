Module param
implicit none

contains

!BS-------------------------------------------------
subroutine setpar
use var
use read_input_file
implicit none
character(LEN=LSTR)::inputfile,inputTIS
character(LEN=XLSTR)::dir_inputTIS
logical::FULL_LINE

print *," "
print *, "setting parameters for the analyze"
inputfile="input.analyze"
inputTIS="input.TISMOL"

call read_inputparameter(PSF_SIM,inputfile,"PSF_SIM")
call read_inputparameter(NGRID,inputfile,"NGRID")
call read_inputparameter(MAXBLOCKLENGTH,inputfile,"MAXBLOCKLENGTH")
call read_inputparameter(BLOCKSKIP,inputfile,"BLOCKSKIP")
call read_inputparameter(MSDLENGTH,inputfile,"MSDLENGTH")
call read_inputparameter(SKIPCROSS,inputfile,"SKIPCROSS")
call read_inputparameter(NMOVMAX,inputfile,"NMOVMAX")
!!call read_inputparameter(LBOX,inputfile,"LBOX")
call read_inputparameter(dxinteg,inputfile,"DXINTEG")

call read_inputparameter(dir,inputfile,"DIR",FULL_LINE)
if (LEN(TRIM(dir))+LEN(TRIM(inputTIS)) > XLSTR) then
  print *,"dir+TISMOL-inputfile name contain too many characters"
  stop
endif
dir_inputTIS=TRIM(dir)//TRIM(inputTIS)
call read_inputparameter (POTENTIAL,dir_inputTIS,"POTENTIAL" )
call read_inputparameter(dt,dir_inputTIS,"DT")
call read_inputparameter(TEMP,dir_inputTIS,"TEMP")
call read_inputparameter(NPART,dir_inputTIS,"NPART")
call read_inputparameter(NWANNIER,dir_inputTIS,"NWANNIER")
call read_inputparameter(DIM,dir_inputTIS,"DIM")
call read_inputparameter(gamma,dir_inputTIS,"GAMMA")
call read_inputparameter(mass,dir_inputTIS,"MASS")
call read_inputparameter(TASK,dir_inputTIS,"TASK" )
call read_inputparameter(DYNAMICS,dir_inputTIS,"DYNAMICS")
call read_inputparameter(HIGH_FRICTION_LIMIT,dir_inputTIS,"HIGH_FRICTION_LIMIT")
call read_inputparameter(TWO_point_method,dir_inputTIS,"TWO_POINT_METHOD")
call read_inputparameter(LBOX,dir_inputTIS,"BOXLENGTH")
call read_inputparameter(forcefieldmatching,dir_inputTIS,"FORCEFIELDMATCHING")



select case(POTENTIAL)
  case("PBD")
    call setparPBD(dir_inputTIS)
  case("WCA")
    call setparWCA(dir_inputTIS)
  case("HARMOSC")
    call read_inputparameter(kharm,dir_inputTIS,"KHARM")
  case("DOUBLEWELL") 
    call read_inputparameter(DOUBLEWELLK4,dir_inputTIS,"DOUBLEWELLK4")
    call read_inputparameter(DOUBLEWELLK2,dir_inputTIS,"DOUBLEWELLK2")
  case default
    print *, "ERROR setpar: POTENTIAL=",POTENTIAL
end select

end subroutine setpar
!ES--------------------------------------------------

!BS--------------------------------------------------
subroutine setparPBD(inputfile)
use var
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
 print *,"PBD parameters:"
  call     read_inputparameter(DAT,inputfile,"DAT")
  call     read_inputparameter(DGC,inputfile,"DGC")
  call     read_inputparameter(AAT,inputfile,"AAT")
  call     read_inputparameter(AGC,inputfile,"AGC")
  call     read_inputparameter(S,inputfile,"S")
  call     read_inputparameter(RHO,inputfile,"RHO")
  call     read_inputparameter(ALPHA,inputfile,"ALPHA")

end subroutine setparPBD
!ES--------------------------------------------------

!BS--------------------------------------------------
subroutine setparWCA(inputfile)
use var
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
 print *,"WCA parameters:"
  call     read_inputparameter(WCA_EPSILON,inputfile,"WCA_EPSILON")
  call     read_inputparameter(WCA_SIGMA,inputfile,"WCA_SIGMA")
  call     read_inputparameter(WCA_H,inputfile,"WCA_H")
  call     read_inputparameter(WCA_W,inputfile,"WCA_W ")
end subroutine setparWCA
!ES--------------------------------------------------


end Module param
