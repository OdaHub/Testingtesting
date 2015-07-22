module order
implicit none

contains

!BS-----------------------------------------------
subroutine make_order(arr1,arr2,N)
implicit none
integer, intent(in)::N
double precision, intent(inout)::arr1(N),arr2(N)
double precision::big,arr1s(N),arr2s(N)
integer::i,minpos(1),mm
big=1.d10
do i=1,N
  minpos=minloc(arr1(:))
  mm=minpos(1)
  arr1s(i)=arr1(mm)
  arr2s(i)=arr2(mm)
  arr1(mm)=big
enddo
arr1(:)=arr1s(:)
arr2(:)=arr2s(:)

end subroutine make_order
!ES-----------------------------------------------

end module order
