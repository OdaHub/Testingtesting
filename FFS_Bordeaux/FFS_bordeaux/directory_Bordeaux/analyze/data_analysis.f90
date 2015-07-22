Module data_analysis 
contains

!BS------------------------------------------
subroutine sparse(ifile,NL,ofile,ngrid)
use stringlengths
implicit none
character(LEN=*), intent(in)::ifile,ofile
integer, intent(in)::ngrid,NL
character(LEN=XLSTR)::LINE
integer::skip,i

skip=max( int((NL-1.d0)/ngrid) , 1 )
open(1,file=ifile)
  open(2,file=ofile)
    do i=1,NL
      read(1,'(A200)') LINE
      if (LINE(200:200)/=" ") then
        print *, "ERROR sparse: lines too long (>200) in", ifile
        print *, "LINE nr",i
        stop
      endif
      if (modulo(i-1,skip)==0) then
        write(2,'(A200)') LINE
      endif
    enddo
  close(2)
close(1)

end subroutine sparse
!ES-----------------------------------------------------------


!BS-------------------------------------------------------------
subroutine runav(ifile,NL,col,ofile,av,ngrid)
!use gnupf
implicit none
character(LEN=*), intent(in)::ifile,ofile
integer, intent(in)::NL,col
double precision, intent(out)::av
integer, optional::ngrid
character::dum(col-1)
integer::skip,i
double precision::xval,yval

av=0.d0
skip=1
if (present(ngrid)) skip=max( int(NL-1.d0)/ngrid, 1 )

open(1,file=ifile)
open(2,file=ofile)
  do i=1,NL
     read(1,*) xval,dum(1:col-2),yval
     av=av+yval
     if (modulo(i-1,skip)==0) then
        write(2,*) xval,av/i
     endif
  enddo
close(2)
close(1)
av=av/NL


end subroutine runav
!ES-------------------------------------------------------------

!BS-------------------------------------------------------------
subroutine distribution(ifile,NL,col,ofile,ngrid,av,minimum,maximum,stdev)
implicit none
character(LEN=*), intent(in)::ifile,ofile
integer, intent(in)::NL,col,ngrid
double precision, optional::av,minimum,maximum,stdev
double precision::avi,mini,maxi,sd,dx,distr(0:ngrid),val,xx,spac
character::dum(col-1)
integer::i,im
 
call avminmax(ifile,NL,col,avi,mini,maxi,sd)
if (present(av)) av=avi
if (present(minimum)) minimum=mini
if (present(maximum)) maximum=maxi
if (present(stdev)) stdev=sd

dx=( maxi-mini )/ngrid
if (dx==0.) dx=1.
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

open(1,file=ofile)
  do i=0,ngrid
    xx=mini+dx*i
    write(1,'(f16.8,e16.8)') xx,distr(i)
  enddo
close(1)

end subroutine distribution 
!ES-------------------------------------------------------------

!BS--------------------------------------
subroutine avminmax(ifile,NL,col,av,mini,maxi,stdev)
implicit none
character(LEN=*), intent(in)::ifile
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
if (stdev<0) then
  print *, "standar deviation ",trim(ifile)," col=",col
  print *,"could not be determined"
  stdev=0.d0
else
  stdev=sqrt( (NL*stdev)/(NL-1) )
endif
end subroutine avminmax
!ES---------------------------------------

!BS---------------------------------------------------------
subroutine block_error(ifile,NL,col,ofile,mbl,blockskip,absolute_e,relative_e,& 
                       Ncorrelation,stdev) 
implicit none
character(LEN=*), intent(in)::ifile,ofile
integer, intent(in)::NL,col,mbl,blockskip
double precision, optional::absolute_e,relative_e,Ncorrelation,stdev
double precision::block(mbl),sumb2(mbl),sumavb(mbl),avb,val,std,std_sqrtN
double precision::denom,Ncor,abs_err,rel_err,av_abs_err,tot_av
double precision::av_rel_err, av_Ncor
integer::nblock(mbl),i,lb,maxbl,countspb,i2ndhalve
character::dum(col-1)

block(:)=0.d0;sumb2(:)=0.d0;nblock(:)=0;sumavb(:)=0.d0
maxbl=min(mbl, int(NL/2.d0) )

open(1,file=ifile)
  do i=1,NL
    read(1,*) dum(1:col-1),val
     do lb=1,maxbl,blockskip
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

std=0.d0
denom=1.d0*(nblock(1)-1)
if (denom > 0 ) then
  std=  ( sumb2(1)-(sumavb(1)**2/nblock(1)) )/ denom  
  if (std<0) then
     print *, "standar deviation ",trim(ifile)," col=",col
     print *,"could not be determined"
     std=0.d0 
  else 
    std=sqrt(std)
  endif
endif
std_sqrtN=std/sqrt(1.d0*NL)
tot_av=sumavb(1)/NL

av_abs_err=0.d0;i2ndhalve=0
open(1,file=ofile)
  do lb=1,maxbl,blockskip
    countspb=nblock(lb)
    sumavb(lb)=sumavb(lb)/countspb
    sumb2(lb)=sumb2(lb)/countspb
    abs_err=0.d0;denom=1.d0*(countspb-1)
    if (denom>0) then
       abs_err= ( sumb2(lb)-sumavb(lb)**2 )/denom  
       if (abs_err>0) then
         abs_err=sqrt(abs_err)
       else
         abs_err=0.d0
       endif
    endif
    rel_err=0.d0;if (abs(tot_av) > 0.d0) rel_err=abs_err/abs(tot_av)
    Ncor=1.d0;if (std_sqrtN>0) Ncor=(abs_err/std_sqrtN)**2 
    if (lb > int(1.d0*maxbl/2)) then
      av_abs_err=av_abs_err+ abs_err
      i2ndhalve=i2ndhalve+1 
    endif
    write(1,'(i8,30f15.8)') lb,abs_err,rel_err,Ncor
  enddo
close(1)
av_abs_err=av_abs_err/i2ndhalve
av_rel_err=0.d0;if (abs(tot_av) > 0.d0) av_rel_err=av_abs_err/abs(tot_av)
av_Ncor=1.d0;if (std_sqrtN>0) av_Ncor=(av_abs_err/std_sqrtN)**2
if (present(absolute_e)) absolute_e=av_abs_err
if (present(relative_e)) relative_e=av_rel_err
if (present(Ncorrelation)) Ncorrelation=av_Ncor 
if (present(stdev)) stdev=std

end subroutine block_error
!ES---------------------------------------------------------

!BS---------------------------------------------------------
subroutine msd(ifile,NL,col,ofile,msdl)
implicit none
character(LEN=*), intent(in)::ifile,ofile
integer, intent(in)::NL,col,msdl
double precision::sumdr2(msdl),r(msdl),shift(msdl)
integer::countspp(msdl)
double precision::time0,time1,dtime,val
character::dum(col-1)
integer::maxdl,i,j

sumdr2(:)=0.d0
countspp(:)=0
maxdl=min(msdl, NL-1 )
open(1,file=ifile)
  read(1,*) time0
  read(1,*) time1
close(1)
dtime=time1-time0

open(1,file=ifile)
  do i=1,maxdl
    read(1,*) dum(1:col-1),val
    r(i)=val
    do j=1,i-1
      sumdr2(j)=sumdr2(j)+(val-r(i-j))**2
      countspp(j)=countspp(j)+1
    enddo
  enddo  
  
  do i=maxdl+1,NL
    read(1,*) dum(1:col-1),val
    do j=1,maxdl
      sumdr2(j)=sumdr2(j)+(r(maxdl+1-j)-val)**2
      countspp(j)=countspp(j)+1
    enddo
    shift(1:maxdl-1)=r(2:maxdl);shift(maxdl)=val
    r=shift
  enddo
close(1)
 
open(1,file=ofile)
write(1,*) 0.d0, 0.d0
do i=1,maxdl
  write(1,*) i*dtime,sumdr2(i)/countspp(i)
enddo
close(1)


end subroutine msd 
!ES---------------------------------------------------------

!BS--------------------------------------------------------
subroutine cross(ifile,NL,col,ofile,ngrid,l1,l2)
implicit none
character(LEN=*), intent(in)::ifile,ofile
integer, intent(in)::NL,col,ngrid
double precision, intent(in)::l1,l2
double precision::avi,mini,maxi,sd,dl,lambda,val,l3
character::dum(col-1)
double precision::arrayl(0:ngrid)
integer::i,im 

call avminmax(ifile,NL,col,avi,mini,maxi,sd)
l3=min(maxi,l2)

dl=( l3-l1 )/ngrid
arrayl(0:ngrid)=0.d0

open(1,file=ifile)
  do i=1,NL
    read(1,*) dum(1:col-1),val
    im=int((val-l1)/dl)
    if (im>ngrid) im=ngrid
    if ((im>=0).AND.(im <= ngrid)) then
      arrayl(0:im)=arrayl(0:im)+1.d0
    else
      print *,"error-cross: IMC=",i,"file=",ifile ;stop
    endif
  enddo
close(1)
arrayl(0:ngrid)=arrayl(0:ngrid)/NL

open(1,file=ofile)
  do i=0,ngrid
    lambda=l1+dl*i
    write(1,'(2f16.8)') lambda,arrayl(i)
  enddo
close(1)

end subroutine cross 
!ES--------------------------------------------------------

!BS---------------------------------------------------------
subroutine vline(val,ofile)
implicit none
double precision, intent(in)::val
character(LEN=*), intent(in)::ofile

open(1,file=ofile)
 write(1,*) val, 0
 write(1,*) val,1
close(1)

end subroutine vline
!ES---------------------------------------------------------

end Module data_analysis
!EM----------------------------------------------------------
