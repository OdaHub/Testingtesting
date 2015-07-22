Module pathanalyze
implicit none

CONTAINS
!BS-----------------------------
subroutine pathout
use var
use inquire
use results
implicit none
character(LEN=XLSTR)::fullfilename
character(LEN=LSTR)::PATHFILE,PATHDIR
double precision::intfdetect
integer::i,NL

PATHFILE="PATH.dat"
fullfilename=trim(dir)//"/"//trim(PATHFILE)
if ((TASK=="TIS").AND.(nonempty(fullfilename))) then
  call pathstat(fullfilename,ext(1),intfl(1),intfm(1),intfr(1),intfr(1),&
                crossprobab(1),relerrcr(1),corcr(1),Lacc(1),Ltr(1),&  
                TISCYCLES(1),accrate(1),PMAXBLOCKLENGTH(1),PBLOCKSKIP(1),&
                forcefieldmatching)
else
  print *,"analyzing TIS directories. Number of ensmenbles ",numint
  do i=1,numint
    fullfilename=trim(dir)//"/"//ext(i)//"/"//trim(PATHFILE)
    intfdetect=intfr(i);if (i<numint) intfdetect=intfm(i+1) 
    if (intfdetect>intfr(i)) then
        print *,"intfdetect>intfr(i)"
        stop
    endif
    call pathstat(fullfilename,ext(i),intfl(i),intfm(i),intfr(i),intfdetect, &
                  crossprobab(i),relerrcr(i),corcr(i),Lacc(i),Ltr(i),&
                  TIScycles(i),accrate(i),PMAXBLOCKLENGTH(i),PBLOCKSKIP(i),&
                  forcefieldmatching)
  enddo 
endif

if (TASK=="PPS") then
  !these ones are for free
  T1=Lacc(1)-2
  NL1=TIScycles(1)
  LTRT1=Ltr(1)
  ACCRATET1=accrate(1)
  call PPSout(T0,T1,NL0,NL1,reT0,reT1,NCORT0,NCORT1,LTRT0,LTRT1,ACCRATET0,&
              ACCRATET1)
endif


end subroutine pathout
!ES--------------------------------

!BS-------------------------------------
subroutine PPSout(T0,T1,NL0,NL1,reT0,reT1,NCORT0,NCORT1,LTRT0,LTRT1,ACCRATET0,ACCRATET1)
use var
use stringlengths
use inquire
use analyze_pcross
!use results
use data_analysis
use gnuplot
implicit none
double precision, intent(out)::T0,reT0,reT1,NCORT0,NCORT1,LTRT0,ACCRATET0
double precision, intent(inout)::T1,LTRT1,ACCRATET1
integer, intent(out)::NL0
integer, intent(inout)::NL1
character(LEN=XLSTR)::ifile0,ifile1
character(LEN=LSTR)::PATHFILE
character(LEN=3)::ooo
character(len=LSTR)::accfileT0,accfileT1,binaryfile,TcycT0
double precision::L1,L2
character(LEN=LSTR)::ofile0,ofile1,name0,name1
character(LEN=LSTR)::commands(20)
integer::col


PATHFILE="PATH.dat"
ifile0=trim(dir)//"/000/"//trim(PATHFILE)
ifile1=trim(dir)//"/001/"//trim(PATHFILE)
ooo="000"



call countlines(ifile0,nl0)
call makeACCPATH(ifile0,nl0,ooo,accfileT0,accRateT0)
call makeTcyc(ifile0,nl0,ooo,TcycT0,forcefieldmatching)
accfileT1="ACCPATH.001.tmp"

call pathdistr(accfileT0,TcycT0,nl0,ooo,L1,L2)
T0=L1-2
LTRT0=L2


name0="errT0.000"
name1="errT1.001"
ofile0=trim(name0)//".tmp"
ofile1=trim(name1)//".tmp"
col=7
call block_error(accfileT0,NL0,col,ofile0,maxblocklength,blockskip,relative_e=reT0,&
                  Ncorrelation=NcorT0)
commands(1)="set xlab ""block length"""
write(commands(2),*) "set title ""T0 rel.err, Ncor",reT0,NcorT0,""""
commands(3)="pl """//trim(ofile0)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",reT0,"w l lw 5"
call make_plot(name0,4,commands)

call block_error(accfileT1,NL1,col,ofile1,maxblocklength,blockskip,relative_e=reT1,&
                  Ncorrelation=NcorT1)
write(commands(2),*) "set title ""T1 rel.err, Ncor",reT0,NcorT0,""""
commands(3)="pl """//trim(ofile1)//"""   u 1:3 w l lw 5"
write(commands(4),*) "repl",reT1,"w l lw 5"
call make_plot(name1,4,commands)

reT0=reT0*(T0+2)/T0
reT1=reT1*(T1+2)/T1


end subroutine PPSout
!ES-------------------------------------


!BS-------------------------------------
subroutine pathstat(ifile,extension,IL,IM,IR,Idetect,cprob,relerr,Ncor,L1,L2, &
                     NL,accRate,mbl,blsk,forcefieldmatching)
use stringlengths
use inquire
use analyze_pcross
implicit none
character(LEN=*), intent(in)::ifile,extension
double precision, intent(in)::IL,IM,IR,Idetect
integer, intent(in)::mbl,blsk
logical, intent(in)::forcefieldmatching
double precision, intent(out)::cprob,relerr,Ncor,L1,L2,accrate
integer, intent(out)::NL
character(len=LSTR)::accfile,binaryfile,Tcyc

call countlines(ifile,nl) 
call makeACCPATH(ifile,nl,extension,accfile,accRate)
call makeTcyc(ifile,nl,extension,Tcyc,forcefieldmatching)
call SUCCESSrate(accfile,nl,extension,binaryfile,il,im,idetect)

call CROSSPROB(accfile,nl,extension,im,ir,idetect)
call runavPCROSS(binaryfile,nl,extension,cprob)
call PCROSS_error(binaryfile,nl,extension,relerr,Ncor,mbl,blsk)
call pathdistr(accfile,Tcyc,nl,extension,L1,L2)

end subroutine pathstat
!ES-------------------------------------

!BS------------------------------------
subroutine makeACCPATH(ifile,nl,ext,ofile,accRate)
use stringlengths
implicit none
character(LEN=*), intent(in)::ifile,ext
integer, intent(in)::nl
character(len=LSTR), intent(out)::ofile
double precision, intent(out)::accrate
character(len=XLSTR)::line,nline
integer::i,NLacc

NLacc=0
ofile="ACCPATH."//trim(EXT)//".tmp"
nline=""
open(1,file=ifile)
open(2,file=ofile)
do i=1,nl
  read(1,'(A200)') line
  if (line(48:50)=="ACC") then
    nline=line
    NLacc=NLacc+1
  endif
  write(2,'(A200)') nline
enddo


close(2)
close(1)
accrate=1.d0*NLacc/NL
end subroutine makeACCPATH
!ES------------------------------------

!BS-------------------------------------------
subroutine SUCCESSrate(ifile,nl,ext,ofile,il,im,idetect)
use stringlengths
implicit none
character(LEN=*), intent(in)::ifile,ext
integer, intent(in)::nl
double precision, intent(in)::il,im,idetect
character(len=LSTR), intent(out)::ofile
character::dum(50)
double precision::lmin,lmax
integer::i,success

ofile="SUCCPATH."//trim(EXT)//".tmp"
open(1,file=ifile)
open(2,file=ofile)
do i=1,nl
  read(1,*) dum(1:9), lmin,lmax
  if (lmin > il) then
      print *,"lmin > il at line",i,ifile
      !stop
  endif
  if (lmax < im) then
    print *,"lmax < im at line",i,ifile
    stop
  endif
  success=0
  if (lmax>idetect) success=1
  write(2,*) i,success
enddo
close(2)
close(1)

end subroutine SUCCESSrate
!ES---------------------------------------------

!BS-------------------------------------------
subroutine MakeTcyc(ifile,nl,ext,ofile,forcefieldmatching)
use stringlengths
implicit none
character(LEN=*), intent(in)::ifile,ext
integer, intent(in)::nl
logical, intent(in)::forcefieldmatching
character(len=LSTR), intent(out)::ofile
character::dum(50)
integer::i,Tcyc,Lpath,stepscopied
character(LEN=2)::MCmove

stepscopied=2
if (forcefieldmatching) stepscopied=1

ofile="Tcyc."//trim(EXT)//".tmp"
open(1,file=ifile)
open(2,file=ofile)
do i=1,nl
  read(1,*) dum(1:6), Lpath, dum(1),MCmove

  select case (MCmove)
    case("sh")
      Tcyc=Lpath-1
    case("tr")
      Tcyc=0
    case("s+")
      Tcyc=0
      if (ext=="000") Tcyc=Lpath- stepscopied
    case("s-")
      Tcyc=0
      if (ext=="001") Tcyc=Lpath- stepscopied
    case("00")
      Tcyc=0
    case default
      print *,"MakeTcyc error: MCmove=",MCmove
  end select

  write(2,*) i,Tcyc
enddo
close(2)
close(1)

end subroutine makeTcyc 
!ES---------------------------------------------



end Module pathanalyze

