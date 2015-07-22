!BM--------------------------------------------------
module arimetric

CONTAINS

!BF-----------------------------------------------
 function distance(x1,x2,d,Lbox)
  implicit none
  integer, intent(in)::d
  double precision, intent(in)::x1(d),x2(d)
  double precision, optional, intent(in)::Lbox
  double precision::distance
  double precision::vec(d)
     vec=x1-x2
     if (present(Lbox)) vec=vec-LBOX*NINT(vec/LBOX)
     distance=sqrt( dot_product(vec,vec) )
end function distance
!EF-----------------------------------------------

end module arimetric
!EM--------------------------------------------------
