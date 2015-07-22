Module movies 
implicit none

CONTAINS
!BS-----------------------------
subroutine makemovies 
use var
use inquire
implicit none
character(LEN=XLSTR)::fullfilename
character(LEN=LSTR)::MOVIEFILE,MOVIEDIR
double precision::intfdetect
integer::i,NLMOV,NskiP
print *,"Making Movie file"

MOVIEFILE="TRAJECTORY.dat"
fullfilename=trim(dir)//"/"//trim(MOVIEFILE)
if (nonempty(fullfilename)) then
  call countlines(fullfilename,NLMOV)
  Nskip=int(1.d0*NLMOV/(NMOVMAX*(Npart+Nwannier+2)))
  if (nskip==0) nskip=1
  print *,"Make single moviefile"
  call movmak(fullfilename,ext(1),Npart+Nwannier,Nskip,LBOX,dim)
else if (numint>0) then
  print *,"Make TIS moviefiles"
  print *," TIS directories. Number of ensmenbles ",numint
    do i=0,numint
      fullfilename=trim(dir)//"/"//ext(i)//"/"//trim(MOVIEFILE)
      if (.NOT.nonempty(fullfilename)) cycle 
      call countlines(fullfilename,NLMOV)
      Nskip=int(1.d0*NLMOV/(NMOVMAX*(Npart+Nwannier+2)))
      if (Nskip==0) nskip=1
      call movmak(fullfilename,ext(i),Npart+Nwannier,nskip,LBOX,dim)
    enddo 
else
  fullfilename=trim(dir)//"/"//"md"//"/"//trim(MOVIEFILE)
  if (nonempty(fullfilename)) then
    call countlines(fullfilename,NLMOV)
    Nskip=int(1.d0*NLMOV/(NMOVMAX*(Npart+Nwannier+2)))
    if (nskip==0) nskip=1
    print *,"Make single moviefile"
    call movmak(fullfilename,ext(1),Npart+Nwannier,Nskip,LBOX,dim)
  endif
endif

end subroutine makemovies 
!ES--------------------------------

!BS-------------------------------------
subroutine movmak(ifile,ext,N,nskip,LBOX,dim)
use stringlengths
!use inquire
implicit none
character(LEN=*), intent(in)::ifile,ext
integer, intent(in)::N,nskip,dim
double precision, intent(in)::LBOX
character(len=LSTR)::ofile
character*100::line
integer::status,count,i
character::char
double precision::X(N,dim),X3D(3)
character*2::type(N)
ofile="MOVIE"//trim(EXT)//".xyz"
open(1,file=ifile,status="old")
open(2,file=ofile)
count=0

do 
  read(1,'(a100)',iostat=status) line 
  if (status /=0 ) exit
  count=count+1
  do i=1,N
    read(1,*) char,type(i),X(i,1:dim)
  enddo

  if (modulo(count-1,nskip)==0) then
    write(2,*) N
    write(2,*) line(1:60) 
    do i=1,N
      if (trim(type(i))=="X") type(i)="Au"
      X3D=0.d0;X3D(1:dim)=X(i,1:dim)
      X3D=X3D-LBOX*NINT(X3D/LBOX)
      write(2,*) type(i),X3D(1:3)
    enddo
  endif   
enddo
close(2)
close(1)
print *,ifile,ofile
end subroutine movmak 
!ES-------------------------------------


end Module movies 

