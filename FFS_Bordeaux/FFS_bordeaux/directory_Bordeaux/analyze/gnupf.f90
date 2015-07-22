module gnupf
contains

!BS---------------------------
subroutine picture(dfile,col1,col2,name,title,xlab,dx,nhl,hl,nlx,ax,multiplot)
implicit none
character*100,intent(in)::dfile,name
character*100,intent(in)::xlab
character*100,intent(in)::title
integer, intent(in):: col1,col2
integer, optional:: nhl,nlx
double precision, optional::hl(:)
character*100, optional::ax(:)
double precision, optional:: dx 
character, optional::multiplot
character*100, allocatable::axx(:)
character*6::plot
integer::Nxx,Nt,i
double precision::mu

Nt=4
if (present(nhl)) Nt=Nt+nhl
if (present(nlx)) Nt=Nt+nlx
allocate(axx(Nt))

axx(1)="set xlab """//trim(xlab)//""""
axx(2)="set title """//trim(title)//""""

mu=1.
if (present(dx)) mu=dx
write(axx(3),'(A3,f20.8)') "mu=",mu
plot="pl   """
if (present(multiplot)) then
   if (multiplot/="s") plot="repl """
endif

write(axx(4),'(A40,i2,A5,i2,A10)') plot//trim(dfile)//""" u ($",col1,"*mu):",col2," w l lw 5"
Nxx=4

if (present (nhl)) then
  do i=1,nhl
      write(axx(Nxx+i), '(A5,f20.8,A8)') "repl ",hl(i), "w l lw 5"
  enddo
  Nxx=Nxx+nhl
endif

if (present (nlx)) then
  do i=1,nlx
     axx(Nxx+i)= ax(i)
  enddo
  Nxx=Nxx+nlx
endif

if  (present(multiplot)) then
  call make_plot(name,nxx,axx,multiplot=multiplot)
else
  call make_plot(name,nxx,axx)
endif
deallocate(axx)

end subroutine picture 
!ES---------------------------


!BS---------------------------------------
subroutine make_plot(name,nlx,ax,multiplot)
use var_analyze
use system
implicit none
character*100,intent(in)::name
integer, intent(in):: nlx
character*100, intent(in)::ax(nlx)
character, optional::multiplot
character*100::string
character*100::gfile,psfile
integer::i
logical::mp

mp=.true. !make postscript (y/n)
gfile=trim(name)//".gnu"
psfile=trim(name)//".eps"

if (present(multiplot)) then
  if (multiplot=="s")  open(1,file=gfile)
  if (multiplot/="s")  open(1,file=gfile,position="append")
  if (multiplot/="f") mp=.false.
else
  open(1,file=gfile)
endif

  do i=1,nlx
    write(1,'(A80)')  ax(i)
  enddo

  if (mp) then
    write(1,*) "set ter pos  enhance col sol"
    write(1,*) "set ou """//trim(psfile)//""""
    write(1,*) "repl"
    write(1,*) "set ter x11"

    string =trim("gnuplot "//gfile)
    print *, trim(string)
    call SYSTEMx(string)
    iff=iff+1
    if (iff > NFMAX) then
      print *,"make_plot-error:max number of files reached",iff
      stop
    endif
    PSFILES(iff)=psfile
  endif
close(1)

end subroutine make_plot
!ES----------------------------------------

!BS---------------------------
subroutine barfig(dfile,col1,col2,name,xlab,title)
implicit none
character*100,intent(in)::dfile,name
character*100,intent(in)::xlab
character*100,intent(in)::title
integer, intent(in):: col1,col2
character*100::axx(3)

axx(1)="set xlab """//trim(xlab)//""""
axx(2)="set title """//trim(title)//""""

write(axx(3),'(A40,i2,A5,i2,A10)') "pl """//trim(dfile)//""" u ",col1,":",col2," w i lw 20"
call make_plot(name,3,axx)

end subroutine barfig
!ES---------------------------


end module gnupf

