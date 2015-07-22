!BM-----------------------------------------------------------------
module erfc

contains

!*******************************************************************************
! Complementary error-function  erfc(x) = 1 - erf(x)
! Reference: W.J. Cody, Mathematics of Computation 22 (1969), 631-637
! Taken from CERN library.
!Adapted by Titus
!*******************************************************************************

      function erfc01(x) 
      implicit none
      double precision, intent(in)  :: x
      double precision              :: erfc01
      integer               :: j
      double precision              :: a,b,v,y
      double precision, parameter   :: isqrtpi = 0.56418958354775629d0   ! 1/sqrt(pi)
      double precision, parameter   :: p1(0:3) = (/  2.426679552305318d+2 , &
                                             2.197926161829415d+1 , &
                                             6.996383488619136d+0 , &
                                            -3.560984370181538d-2   /)
      double precision, parameter   :: q1(0:3) = (/  2.150588758698612d+2 , &
                                             9.116490540451490d+1 , &
                                             1.508279763040779d+1 , &
                                             1.000000000000000d+0   /)
      double precision, parameter   :: p2(0:7) = (/  3.004592610201616d+2 , &
                                             4.519189537118729d+2 , &
                                             3.393208167343437d+2 , &
                                             1.529892850469404d+2 , &
                                             4.316222722205674d+1 , &
                                             7.211758250883094d+0 , &
                                             5.641955174789740d-1 , &
                                            -1.368648573827167d-7   /)
      double precision, parameter   :: q2(0:7) = (/  3.004592609569833d+2 , &
                                             7.909509253278980d+2 , &
                                             9.313540948506096d+2 , &
                                             6.389802644656312d+2 , &
                                             2.775854447439876d+2 , &
                                             7.700015293522947d+1 , &
                                             1.278272731962942d+1 , &
                                             1.000000000000000d+0   /)
      double precision, parameter   :: p3(0:4) = (/ -2.996107077035422d-3 , &
                                            -4.947309106232507d-2 , &
                                            -2.269565935396869d-1 , &
                                            -2.786613086096478d-1 , &
                                            -2.231924597341847d-2   /)
      double precision, parameter   :: q3(0:4) = (/  1.062092305284679d-2 , &
                                             1.913089261078298d-1 , &
                                             1.051675107067932d+0 , &
                                             1.987332018171353d+0 , &
                                             1.000000000000000d+0   /)


      v = abs(x)
      if (v <= 0.46875d0) then
          y = v**2
          a = p1(3)
          b = q1(3)
          do j=2,0,-1
              a = a*y + p1(j)
              b = b*y + q1(j)
          enddo
          erfc01 = 1.0d0 - x*a/b
          return
      elseif (v <= 4.0d0) then
          a = p2(7)
          b = q2(7)
          do j=6,0,-1
              a = a*v + p2(j)
              b = b*v + q2(j)
          enddo
          erfc01 = exp(-v**2)*a/b
      elseif (v <= 10.0d0) then
          y = 1.0d0/(v**2)
          a = p3(4)
          b = q3(4)
          do j=3,0,-1
              a = a*y + p3(j)
              b = b*y + q3(j)
          enddo
          erfc01 = exp(-v**2)*(isqrtpi+y*a/b)/v
      else
          erfc01 = 0.0d0
      endif
      if (x <= 0.0d0) erfc01 = 2.0d0 - erfc01

      end function erfc01

end module erfc
!EM----------------------------------------------------------------
