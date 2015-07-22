Module efcol
contains

!BS-------------------------------------------------------------
subroutine runav(ifile,NL,col,name,xlab,dx,ngrid,pict)
use gnupf
implicit none
character*100, intent(in)::ifile,name
integer, intent(in)::NL,col
character*100, optional::xlab
double precision, optional::dx
integer, optional::ngrid
character,optional::pict
character*100::xxlab,name2,ofile
character*100::title
double precision::ddx,val,av,hl(1)
integer::NL2,skip,i
character::dum(col-1)

ddx=1.d0
xxlab=" "
if (present(dx)) ddx=dx
if (present(xlab)) xxlab=xlab

name2="runav"//trim(name)
ofile=trim(name2)//".tmp"

av=0.d0
NL2=0
skip=1
if (present(ngrid)) skip=max( int(1.d0*NL/(ngrid+1)) , 1 )

open(1,file=ifile)
open(2,file=ofile)
  do i=1,NL
     read(1,*) dum(1:col-1),val
     av=av+val
     if (modulo(i-1,skip)==0) then
        write(2,*) (i-1)*ddx,av/i
        NL2=NL2+1
     endif
  enddo
close(2)
close(1)
av=av/NL


if (present(pict)) then
  if (pict=="n") return
  write(title,'(A40,f20.8)') name2,av
  hl(1)=av
  call picture(ofile,1,2,name2,title,xxlab,nhl=1,hl=hl )  
endif

end subroutine runav
!ES-------------------------------------------------------------

!BS------------------------------------------
subroutine sparse(ifile,ofile,NL,ngrid,NL2)
implicit none
character*100, intent(in)::ifile,ofile
integer, intent(in)::NL,ngrid
integer, intent(out)::NL2
character*200::LINE
integer::skip,i

skip=max( int(1.d0*NL/(ngrid+1)) , 1 )
NL2=0
open(1,file=ifile)
open(2,file=ofile)
  do i=1,NL
     read(1,'(A200)') LINE
     if (modulo(i-1,skip)==0) then
       write(2,'(A200)') LINE
       NL2=NL2+1
     endif
  enddo
close(2)
close(1)
end subroutine sparse
!ES------------------------------------------

!BS-------------------------------------------------------------
subroutine distribution(ifile,NL,col,name,ngrid,xlab,nlx,ax,avo,minimum,maximum,spacing)
use gnupf
implicit none
character*100, intent(in)::ifile,name,xlab
integer, intent(in)::NL,col,ngrid
integer, optional::nlx
double precision, optional::avo,minimum,maximum
character*100, optional::ax(:)
double precision, optional::spacing
character*100::name2,ofile
character*100::title
double precision::av,mini,maxi,stdev,dx,distr(0:ngrid),val,xx,spac
character::dum(col-1)
integer::i,im
 
call avminmax(ifile,NL,col,av,mini,maxi,stdev)
dx=( maxi-mini )/ngrid
distr(:)=0.d0

open(1,file=ifile)
  do i=1,NL
    read(1,*) dum(1:col-1),val
    im=int((val-mini)/dx)
    if ((im>=0).AND.(im <= ngrid)) then
      distr(im)=distr(im)+1.d0
    else
      print *,"error-distribution: IMC=",i,"file=",ifile ;stop
    endif
  enddo
close(1)
distr(:)=distr(:)/(NL*dx)

name2="dist"//trim(name)
ofile=trim(name2)//".tmp"

open(1,file=ofile)
  do i=0,ngrid
    xx=mini+dx*i
    write(1,'(i8,2f16.8)') i,xx,distr(i)
  enddo
close(1)

spac=1.d0
if (present(spacing)) spac=spacing
write(title,'(A20,A5,f16.8,A8,f16.8)') name2,", av=",av,", stdev=",stdev
if (present(nlx)) then
  call picture(ofile,2,3,name2,title,xlab,nlx=nlx,ax=ax,dx=spac)
else
  call picture(ofile,2,3,name2,title,xlab,dx=spac)
endif
if (present(avo)) avo=av
if (present(minimum)) minimum=mini
if (present(maximum)) maximum=maxi

end subroutine distribution 
!ES-------------------------------------------------------------

!BS--------------------------------------
subroutine avminmax(ifile,NL,col,av,mini,maxi,stdev)
implicit none
character*100, intent(in)::ifile
integer, intent(in)::col,NL
double precision, intent(out)::av,mini,maxi,stdev
double precision:: val
integer::i 
character::dum(col-1)

av=0.d0
maxi=-9.d9
mini=+9.d9
stdev=0.d0
 
open(1,file=ifile)
  do i=1,NL
     read(1,*) dum(1:col-1),val
     av=av+val
     stdev=stdev+val**2
     if (val < mini  ) mini = val 
     if (val > maxi  ) maxi = val  
   enddo
close(1)
av=av/NL
stdev=stdev/NL
stdev=stdev-av**2
stdev=sqrt( (NL*stdev)/(NL-1) )

end subroutine avminmax
!ES---------------------------------------

!BS---------------------------------------------------------
subroutine block_error(ifile,NL,col,name,lbm,rerr,corr,stdev)
use gnupf
implicit none
character*100, intent(in)::ifile,name
integer, intent(in)::NL,col,lbm
double precision, optional::rerr,corr,stdev
double precision::block(lbm),sumb2(lbm),sumavb(lbm),avb,val,std
double precision::Neff,err,error,reler,err1,hl(1)
integer::nblock(lbm),i,lb,lbmax,k,jj,nhl
character::dum(col-1)
character*100::name2,ofile,xlab
character*100::title

block(:)=0.d0;sumb2(:)=0.d0;nblock(:)=0;sumavb(:)=0.d0
lbmax=min(lbm, int(NL/2.d0) )

open(1,file=ifile)
  do i=1,NL
    read(1,*) dum(1:col-1),val
     do lb=1,lbmax
       block(lb)=block(lb)+val
       if ( modulo(i,lb)==0) then
         avb=block(lb)/lb
         sumavb(lb)=sumavb(lb)+avb
         sumb2(lb)=sumb2(lb)+avb**2
         block(lb)=0.d0
         nblock(lb)=nblock(lb)+1
       endif
     enddo
  enddo
close(1)

std=sqrt(  ( sumb2(1)-(sumavb(1)**2/nblock(1)) )/ (nblock(1)-1)  )
error=0.d0;jj=0
name2="err"//trim(name)
ofile=trim(name2)//".tmp"
open(1,file=ofile)
  do lb=1,lbmax
    k=nblock(lb)
    sumavb(lb)=sumavb(lb)/k
    sumb2(lb)=sumb2(lb)/k
    err=sqrt( ( sumb2(lb)-sumavb(lb)**2 )/(k-1)  )
    if (lb==1) err1=err
    Neff=(err/err1)**2 
    if (lb > int(1.d0*lbmax/2)) then
      error=error+ err
      jj=jj+1 
    endif
    write(1,'(i8,30f15.8)') lb,err,err/sumavb(1),Neff
  enddo
close(1)
error=error/jj
reler=error/sumavb(1)
Neff=(error/err1)**2

xlab="block length"
hl(1)=reler
write(title,'(A20,A8,f8.4,A10,f8.4,A7,f8.4)') name2,", error=",error,", rel.err=",reler,", Neff=",Neff
call picture(ofile,1,3,name2,title,xlab,nhl=1,hl=hl)
if (present(rerr)) rerr=reler
if (present(corr)) corr=Neff
if (present(stdev)) stdev=std

end subroutine block_error
!ES---------------------------------------------------------

end Module efcol
