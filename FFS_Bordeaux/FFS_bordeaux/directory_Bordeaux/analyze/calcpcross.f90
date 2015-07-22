module calcpcross

CONTAINS

!BS--------------------------
subroutine calc_crossp(l1,l2,ngrid,fpath,Nlp,name,f_out,crossp)
use gnupf
implicit none
double precision, intent(in)::l1,l2
integer, intent(in)::ngrid,Nlp
character*100, intent(in)::fpath,name
character*100,optional::f_out
double precision, optional::crossp
double precision::dlam,arrl(0:ngrid),lam 
character::dum(10)
integer::i,il
character*100::ofile,xlab
character*100::title

dlam=(l2-l1)/ngrid
arrl(:)=0.d0

open(1,file=fpath)
  do i=1,Nlp
    read(1,*) dum(1:4),lam 
    il=int((lam-l1)/dlam)
    if (il > ngrid ) il=ngrid
    if (il < 0) then
      print *,"calc_crossp-error in file:",fpath,i
      print *,"l1,l2,lam",l1,l2,lam
      stop
    endif
    arrl(0:il)=arrl(0:il)+1.d0
  enddo
close(1)

arrl(:)=arrl(:)/Nlp
ofile=trim(name)//".tmp"

open(1,file=ofile)
  do i=0,ngrid
    write(1,*) i,l1+i*dlam,arrl(i)
  enddo
close(1)
xlab="lambda"
write(title,'(A21,A20,f8.4)')"crossing probability:",name,arrl(ngrid)
call picture(ofile,2,3,name,title,xlab)
if (present(f_out)) f_out=ofile
if (present(crossp)) crossp=arrl(ngrid)

end subroutine calc_crossp
!ES--------------------------

!BS---------------------------------------------
subroutine calc_actEN(l1,l2,ngrid,fpath,Nlp,name,f_out)
use efcol
use gnupf
implicit none
double precision, intent(in)::l1,l2
integer, intent(in)::ngrid,Nlp
character*100, intent(in)::fpath,name
character*100, optional::f_out
double precision::dlam,lam,arrActE(0:ngrid),ActEN 
character::dum(10)
integer::i,il,Ncount(0:ngrid),NAEd
character*100::ofile,xlab,fileAEd,name2
character*100::title

dlam=(l2-l1)/ngrid
arrActE(:)=0.d0
Ncount=0
fileAEd="d"//trim(name)//".tmp"
NAEd=0

open(1,file=fpath)
open(2,file=fileAEd)
  do i=1,Nlp
    read(1,*) dum(1:3),ActEN,lam
    il=int((lam-l1)/dlam)
    if (il > ngrid ) il=ngrid
    if (il < 0) then
      print *,"calc_actEN-error in file:",fpath,i
      stop
    endif
    arrActE(0:il)=arrActE(0:il)+ActEN
    Ncount(0:il)=Ncount(0:il)+1
    if (il==ngrid) then
      NAEd=NAEd+1
      write(2,*) NAEd, ActEN
    endif
  enddo
close(2)
close(1)

!avoid division by zero
do i=0,ngrid
  if (Ncount(i)==0) Ncount(i)=1
enddo
arrActE(:)=arrActE(:)/Ncount(:)

ofile=trim(name)//".tmp"
open(1,file=ofile)
  do i=0,ngrid
    write(1,*) i,l1+i*dlam,arrActE(i)
  enddo
close(1)

xlab="lambda"
write(title,'(A18,A20,f8.4)')"activation energy:",name,arrActE(ngrid)
call picture(ofile,2,3,name,title,xlab)
if (present(f_out)) f_out=ofile

if (NAEd == 0) return
xlab="(eV)"
name2="endEN"
call distribution(fileAEd,NAEd,2,name2,ngrid,xlab)

end subroutine calc_actEN
!ES---------------------------------------------

!BS----------------------------------------------------------------
subroutine mochrej(ifile,ofile,NL,ngrid)
use efcol
implicit none
character*100, intent(in)::ifile,ofile
integer, intent(in)::NL,ngrid
double precision::maxi,mini
double precision::dx,val,xx
integer::Ncount(0:ngrid),i,im 
double precision:: distr(0:ngrid),av,stdev
character::dum
character*3::AR

Ncount(:)=0
distr(:)=0.d0
call avminmax(ifile,NL,2,av,mini,maxi,stdev)
dx=( maxi-mini )/ngrid
distr(:)=0.d0
Ncount(:)=0
open(1,file=ifile)
  do i=1,NL
    read(1,*) dum,val,AR
    im=int((val-mini)/dx)
    if ((im>=0).AND.(im <= ngrid)) then
      Ncount(im)=Ncount(im)+1.d0
      if (AR=="MCR") distr(im)=distr(im)+1.d0
    else
      print *,"mochrej-error:file=",ifile,i;stop
    endif
  enddo
close(1)
!avoid division by zero
do i=0,ngrid
  if (Ncount(i)==0) Ncount(i)=1
enddo
distr(:)=distr(:)/Ncount(:)

open(1,file=ofile)
  do i=0,ngrid
    xx=mini+dx*i
    write(1,'(i8,20f16.8)') i,xx,distr(i)
  enddo
close(1)

end subroutine mochrej
!ES----------------------------------------------------------------

end module calcpcross

