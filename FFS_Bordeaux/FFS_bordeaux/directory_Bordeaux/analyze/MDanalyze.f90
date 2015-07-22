Module MDanalyze 
contains
!BS---------------------------------------------------------
subroutine MDout
use var
use results
use inquire
use analyze_EN
use analyze_OP
use analyze_flux
implicit none
character(LEN=XLSTR)::fullfilename
character(LEN=LSTR)::ENFILE,MDdir,OPFILE,CROSSFILE
character(LEN=LSTR)::effcrossfile1,effcrossfile2,effcrossfile3
character(LEN=LSTR)::fluxfile1,fluxfile2,fluxfile3
integer::NL,nlflux

print *,"Analyze MD results"

MDdir="md"
ENFILE="ENERGIES.dat"
OPFILE="ORDERPARAM.dat"
CROSSFILE="CROSS.dat"

  call check_datafile(dir,MDdir,ENFILE,fullfilename,NL,TASK,"MD")
  if (NL>1) then
    print *,"  analyze energies"
    call  Evstime(fullfilename,NL)
    call  runavE(fullfilename,NL)
    call  Edistr(fullfilename,NL)
    call  E_block_analyze(fullfilename,NL)
  endif

  call check_datafile(dir,MDdir,OPFILE,fullfilename,NL,TASK,"MD")
  if (NL>1) then
     print *,"  analyze orderparameter"
     call OPvstime(fullfilename,NL)
     call runavOP(fullfilename,NL)
     call OPdistr(fullfilename,NL)
     call OP_block_analyze(fullfilename,NL)
     call mean_square_disp_OP(fullfilename,NL)
  endif

  call check_datafile(dir,MDdir,CROSSFILE,fullfilename,NL,TASK,"MD")
  if (NL>1) then
    print *,"  analyze crossings"
    !make eff-crossing files for 3 interfaces
    call make_effcrossings_file(fullfilename,nl,effcrossfile1,1,MDcycles, &
                  NSA,NSB,NOSA,NOSB,NCROSS(1),NEFFCROSS(1))
    call make_effcrossings_file(fullfilename,nl,effcrossfile2,2,MDcycles,& 
              NCROSSo=NCROSS(2),NEFFCROSSo=NEFFCROSS(2))
    call make_effcrossings_file(fullfilename,nl,effcrossfile3,3,MDcycles,&
              NCROSSo=NCROSS(3),NEFFCROSSo=NEFFCROSS(3))
    !make flux file
    call make_flux_file(effcrossfile1,1,NOSA,skipcross,dt,fluxfile1,nlflux)
    call make_flux_file(effcrossfile2,2,NOSA,skipcross,dt,fluxfile2)
    call make_flux_file(effcrossfile3,3,NOSA,skipcross,dt,fluxfile3)
    if (nlflux>1) then
      print *,"  analyze flux"
      call runavFlux(fluxfile1,fluxfile2,fluxfile3,nlflux,flux1,flux2,flux3)
      call flux_block_analyze(fluxfile1,fluxfile2,fluxfile3,nlflux,&
           relerrflux1,relerrflux2,relerrflux3) !,& 
           !corflux1,corflux2,corflux3)
    endif
        
  endif



end subroutine MDout
!ES----------------------------------------------------------

!BS-------------------------------------------
subroutine make_effcrossings_file(IFILE,NL,ofile,intf,MDcycles,&
                                  NSAo,NSBo,NOSAo,NOSBo,NCROSSo,NEFFCROSSo)
use stringlengths
implicit none
character(LEN=*), intent(in)::ifile
integer,intent(in)::nl,intf,MDcycles
character(len=LSTR), intent(out)::ofile
integer, optional, intent(out)::NSAo,NSBo,NOSAo,NOSBo
integer, optional, intent(out)::NCROSSo,NEFFCROSSo
integer::MDtime,interface
character::sign
integer::startOSA,endOSA,startOSB,endOSB,line,startSB,endSB,startSA,endSA
logical::overallstateA,overallstateB,firstcross
integer::NSA,NSB,NOSA,NOSB
integer::NCROSS,NEFFCROSS

! This routine writes down the first crossings with interface 
! intf and substracts the time that
! the system was in overall-state B. Optionally, also the the number of 
! timesteps that the system spends in stable-state A,B and overall-state A,B
! can be an output


firstcross=.true.
write(ofile,'(a6,i1,a4)')"effcr.",intf,".tmp"
startOSA=0;endOSA=0
startOSB=0;endOSB=0
startSB=0;endSB=0
startSA=0;endSA=0
NSA=0;NSB=0;NOSA=0;NOSB=0
NCROSS=0;NEFFCROSS=0

!first determine which overallstate at start
open(1,file=ifile)
  read(1,*) MDtime, interface, sign
close(1)  
overallstateA=.true.;overallstateB=.false.
if ((interface==3).AND.(sign=="-")) then
  overallstateA=.false.
  overallstateB=.true.
endif

open(1,file=ifile)
open(2,file=ofile)
  do line=1,nl
    read(1,*) MDtime, interface, sign
  
    if ((interface==1).AND.(sign=="-")) then
      startSA=MDtime
      firstcross=.true.
      if (overallstateB) then
        endOSB = MDtime
        startOSA=MDtime
        NOSB=NOSB+endOSB-startOSB
        overallstateA=.true.;overallstateB=.false.
      endif
    endif

    if ((interface==1).AND.(sign=="+")) then
      endSA=MDtime
      NSA=NSA+endSA-startSA
    endif

    if ((interface==3).AND.(sign=="+")) then
      startSB=MDtime
      if (overallstateA) then
        endOSA = MDtime
        startOSB=MDtime
        NOSA=NOSA+endOSA-startOSA
        overallstateA=.false.;overallstateB=.true.
      endif
    endif

    if ((interface==3).AND.(sign=="-")) then
      endSB=MDtime
      NSB=NSB+endSB-startSB
    endif

    if ((interface==intf).AND.(sign=="+")) then
      NCROSS=NCROSS+1
      if (firstcross) then
        write(2,*) MDtime-NOSB,MDtime
        firstcross=.false.
        NEFFCROSS=NEFFCROSS+1
      endif
    endif
    write(1111,*) line,NOSA
  enddo
close(2)
close(1)

if (overallstateA)  NOSA=NOSA+MDcycles-startOSA
if (overallstateB)  NOSB=NOSB+MDcycles-startOSB
if ((interface==3).AND.(sign=="+")) NSB=NSB+MDcycles-startSB
if ((interface==1).AND.(sign=="-")) NSA=NSA+MDcycles-startSA

if (present(NSAo)) NSAo=NSA
if (present(NSBo)) NSBo=NSB
if (present(NOSAo)) NOSAo=NOSA
if (present(NOSBo)) NOSBo=NOSB
if (present(NCROSSo)) NCROSSo=NCROSS
if (present(NEFFCROSSo)) NEFFCROSSo=NEFFCROSS

end subroutine make_effcrossings_file
!ES--------------------------------------------

!BS------------------------------------------------------
subroutine make_flux_file(ifile,interface,NOSA,skip,dt,ofile,nlflux)
implicit none
character(LEN=*), intent(in)::ifile
integer,intent(in)::interface,NOSA,skip
double precision, intent(in)::dt
character(LEN=*), intent(out)::ofile
integer, optional, intent(out)::nlflux
integer::NL,i,time,crtime,Ncross
integer::status

write(ofile,'(a5,i1,a4)') "FLUX.",interface,".tmp"
crtime=0
NL=int(1.d0*NOSA/skip)

open(1,file=ifile)
open(2,file=ofile)
read(1,*,iostat=status) crtime
if (status /=0 ) crtime=nosa

do i=1,NL

  time=i*skip
  Ncross=0
  do while (crtime < time) 
    Ncross=Ncross+1
    read(1,*,iostat=status) crtime
    if (status /=0 ) crtime=nosa 
  enddo
  write(2,*) time,Ncross,1.d0*Ncross/(skip*dt)
enddo
close(2)
close(1)

if (present(nlflux)) nlflux=nl

end subroutine make_flux_file
!ES------------------------------------------------------


end Module MDanalyze
!EM-------------------------------------------------------------
