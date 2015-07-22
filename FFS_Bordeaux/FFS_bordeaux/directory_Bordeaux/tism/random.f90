!BM----------------------
module random
implicit none

CONTAINS

  !BF--------------------------------------------------
  function random01()
  implicit none
  double precision::random01
  double precision::ran
   call random_number(ran)
   random01=ran
   !!print *,"ran",ran
  end function random01
  !EF--------------------------------------------------

  !-------------------------------------
  ! This function returns a random number from
  ! a Gaussian distribution. 
  ! Technical Note: In fact, 2 random numbers 
  ! are automatically generated, but one is not used (waist).
  ! The a previous version was tested using a saveflag option
  ! so that the 2nd random number was given in a second call.
  ! No speed increment was observed while it required more
  ! difficulty for RESTART options etc.
  !BF---------------------------
  function rangaussian(sigma)
  implicit none
  double precision, optional::sigma
  double precision::rangaussian
  double precision::gx,s,r,sqrtr,gy

  
      s=2.
      do while (s>1.)
        !!call random_number(gx)
        gx=random01()
        !!call random_number(gy)
        gy=random01()
        gx=2*gx-1.
        gy=2*gy-1.
        s=gx**2+gy**2
      enddo
      r=-2.*log(s)/s
      sqrtr=sqrt(r)
      gx = gx * sqrtr

      rangaussian=sigma*gx

  end function rangaussian
  !EF--------------------------

  !BS----------------------------------------------------------------
  subroutine gssbivar(rdr,rdv,s12os11,sqrts11,sqrtSos11)
  implicit none
  double precision, intent(in)::s12os11,sqrts11,sqrtSos11
  double precision, intent(out)::rdr,rdv
  double precision::v1,v2,sv,r,fac

  sv=2.
  do while (sv>1.)
    !!call random_number(v1)
    v1=random01()
    !!call random_number(v2)
    v2=random01()
    !!print *,"gssbivar random_number",v1,v2
    v1=2*v1-1.
    v2=2*v2-1.
    sv=v1**2+v2**2
  enddo

  r=-2.*log(sv)/sv
  fac=sqrt(r)

  rdr=v1*fac*sqrts11
  rdv=v2*fac*sqrtsos11+rdr*s12os11

  end subroutine gssbivar
  !ES---------------------------------------------------

  !BS-----------------------------------------------------
  subroutine set_maxwellian_velocities(v,Npart,dim,sigma_v)
  implicit none
  integer, intent(in)::Npart,dim
  double precision, intent(in)::sigma_v(Npart)
  double precision, intent(out)::v(Npart,dim)
  integer::i,j
  double precision::sig

    do i=1,Npart
      sig=sigma_v(i)
      do j=1,dim
        v(i,j)=rangaussian(sig)
      enddo
    enddo

  end subroutine set_maxwellian_velocities
  !ES-----------------------------------------------------

end module random
!EM----------------------

