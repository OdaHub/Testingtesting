!BM--------------------------------------
module init

contains

!BS------------------------------
subroutine initialize
use var
implicit none
double precision::ev_J,u_kg,A_m,ns_s
character*3::dum
integer::i
character*50::ru

A_m=1.D-10               !Angstrom in m
ev_J=1.60217646D-19      !eV in J
u_kg=1.66053886D-27      !umass in kg
ns_s=1.d-9               !ns in  sec
                

open(11,file="fexc.sh")


ru="reduced units"
select case(POTENTIAL)
  case("PBD")
   chartu="ns"
   chareu="eV"
   charlu="A"
   chartempu="K"
   timeunit=A_m*sqrt(u_kg/eV_J)/(ns_s)
  case Default
   chartu=ru
   chareu=ru
   charlu=ru
   chartempu=ru
   timeunit=1.d0
end select  

dt_ut=timeunit*dt
gamma_invut=gamma/timeunit


pi=4*atan(1.d0)

 if ((POTENTIAL=="WCA").OR.(POTENTIAL=="IONTRANS").OR.(POTENTIAL=="DOUBLEWELL")) then
    print *,"USING REDUCED UNITS kb=1"
    kb=1.d0
    mass_eVns2_A2=mass
  else
    kb=1.3806503D-23     ! J/K
    kb=kb/ev_J
    print *,"kb=", kb, " eV/K"
    mass_eVns2_A2=(mass*u_kg/ev_J)*(A_m)**2/(ns_s**2)
  endif
kbT=kb*Temp
beta=1.d0/kbT                !inverse temperature


do i=0,999
  if (i <10) then
    write(dum,'(A2,i1)')"00",i
  else if (i <100) then
    write(dum,'(A1,i2)')"0",i
  else if (i <1000) then
    write(dum,'(i3)')i
  endif
  EXT(i)=dum
enddo

!check interfaces and pathstat.files
call initMD
call inittrans
call initanalyzeTIS
NPSFMAX=(numint+1)*PSF_SIM
NPSF=0
allocate(PSFILES(NPSFMAX))
Brownian=.false.
if ((dynamics=="LANGEVIN").AND.(high_friction_limit)) Brownian=.true.

end subroutine initialize
!ES-------------------------------

!BS--------------------------------------
subroutine initMD
use var
use read_input_file
use inquire
implicit none
character(LEN=LSTR):: crossfile
character(LEN=XLSTR)::fullcrossfile
character(LEN=XLSTR)::inputMD
integer::i

crossfile="CROSS.dat"
fullcrossfile=trim(dir)//"/md/"//trim(crossfile)

!if TASK=MD and CROSS.dat is in main directory take this one.
inputMD=trim(dir)//"/input.TISMOL"
if ((TASK/="MD").AND.(NONEMPTY(fullcrossfile))) inputMD=trim(dir)//"/md/input.TISMOL"

call read_inputparameter(MD_INTFL,inputMD,"INTERFACEL")
call read_inputparameter(MD_INTFM,inputMD,"INTERFACEM")
call read_inputparameter(MD_INTFR,inputMD,"INTERFACER")
call read_inputparameter(MDCYCLES,inputMD,"NCYCLES")

end subroutine initMD
!ES--------------------------------------

!BS--------------------------------------
subroutine inittrans
use var
use read_input_file
use inquire
implicit none
character(LEN=LSTR):: transfile
character(LEN=XLSTR)::fulltransfile
character(LEN=XLSTR)::inputtrans

transfile="TRANSMISSION.dat"
fulltransfile=trim(dir)//"/trans/"//trim(transfile)

!if TASK=TRANSMISSION and TRANSMISSION.dat is in main directory take this one.
inputtrans=trim(dir)//"/input.TISMOL"
if ((TASK/="TRANSMISSION").AND.(NONEMPTY(fulltransfile))) inputtrans=trim(dir)//"/trans/input.TISMOL"

call read_inputparameter(trans_INTFL,inputtrans,"INTERFACEL")
call read_inputparameter(trans_INTFM,inputtrans,"INTERFACEM")
call read_inputparameter(trans_INTFR,inputtrans,"INTERFACER")

end subroutine inittrans
!ES--------------------------------------


!BS---------------------------------
subroutine initanalyzeTIS
use var
use read_input_file
use inquire
use results
implicit none
character(LEN=LSTR)::pathstat
character(LEN=XLSTR)::fullpathstat
character(LEN=XLSTR)::inputTIS
integer::i,numint_pps
character(LEN=LSTR)::input_an
input_an="input.analyze"

pathstat="PATH.dat"
fullpathstat=trim(dir)//trim(pathstat)

!if TASK=TIS and PATH.dat is in main directory take this one.
!Hence, one single TIS analysis will be performed
!if TASK=PPS all interfaces must be read from central input file
if (TASK=="PPS") then
  inputTIS=trim(dir)//"/input.TISMOL"
  call read_inputparameter(NUMINT_PPS,inputTIS,"NUMINT_PPS")
  allocate(PPS_INTERFACES(NUMINT_PPS))
  numint=numint_pps-1
  allocate(INTFL(numint));allocate(INTFM(numint));allocate(INTFR(numint))
  call read_inputparameter(PPS_INTERFACES,inputTIS,"PPS_INTERFACES")
  do i=1,numint
    INTFL(i)=PPS_INTERFACES(1)
    INTFM(i)=PPS_INTERFACES(i)
    INTFR(i)=PPS_INTERFACES(numint_pps)    
  enddo

else if ((TASK=="TIS").AND.(NONEMPTY(fullpathstat))) then

  NUMINT=1
  inputTIS=trim(dir)//"/input.TISMOL"
  allocate(INTFL(1));allocate(INTFM(1));allocate(INTFR(1))
  call read_inputparameter(INTFL(1),inputTIS,"INTERFACEL")
  call read_inputparameter(INTFM(1),inputTIS,"INTERFACEM")
  call read_inputparameter(INTFR(1),inputTIS,"INTERFACER") 

else 

  numint=0
  do 
    fullpathstat=trim(dir)//"/"//trim(EXT(numint+1))//"/"//trim(pathstat)
    inputTIS=trim(dir)//"/"//trim(EXT(numint+1))//"/"//"input.TISMOL" 
    if (.not.NONEMPTY(fullpathstat)) exit
    numint=numint+1
  enddo

  allocate(INTFL(numint));allocate(INTFM(numint));allocate(INTFR(numint))
  do i=1,numint
    inputTIS=trim(dir)//"/"//trim(EXT(i))//"/"//"input.TISMOL"
    if (.not.nonempty(inputTIS)) then
      print *,"non-existing or empty ",inputTIS
      numint=0;exit
    endif
    call read_inputparameter(INTFL(i),inputTIS,"INTERFACEL")
    call read_inputparameter(INTFM(i),inputTIS,"INTERFACEM")
    call read_inputparameter(INTFR(i),inputTIS,"INTERFACER")
  enddo

endif

allocate(Crossprobab(numint))
allocate(relerrcr(numint))
allocate(corcr(numint))
allocate(Lacc(numint))
allocate(Ltr(numint))
allocate(TIScycles(numint))
allocate(accrate(numint))
allocate(PMAXBLOCKLENGTH(numint))
allocate(PBLOCKSKIP(numint))


call read_inputparameter(PMAXBLOCKLENGTH,input_an,"PMAXBLOCKLENGTH")
call read_inputparameter(PBLOCKSKIP,input_an,"PBLOCKSKIP")

end subroutine initanalyzeTIS
!ES---------------------------------


end module init
!EM--------------------------------------
