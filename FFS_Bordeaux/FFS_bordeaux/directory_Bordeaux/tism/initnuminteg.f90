!BM--------------------------
Module initNUMINTEG

Contains

  !BS------------------------
  subroutine init_numericinteg
  use inputpar 
  use integration_module
  implicit none
  integer::maxpoints,npoints(dim),i

  allocate(Numinteg_param%RANGE_LEFT(dim))
  allocate(Numinteg_param%RANGE_RIGHT(dim))
  allocate(Numinteg_param%INTSTEP(dim))
  Numinteg_param%RANGE_LEFT=RANGE_LEFT
  Numinteg_param%RANGE_RIGHT=RANGE_RIGHT
  Numinteg_param%INTSTEP=INTSTEP
  npoints=int((RANGE_RIGHT-RANGE_LEFT)/INTSTEP)
  maxpoints=maxval(npoints)
  allocate(Numinteg_param%npoints(dim))
  Numinteg_param%npoints=npoints
  allocate(Numinteg_param%FREE_EN_ARRAY(dim,0:maxpoints))
  allocate(Numinteg_param%prob_ARRAY(dim,0:maxpoints))
  Numinteg_param%prob_ARRAY=0.d0
  Numinteg_param%FREE_EN_ARRAY=0.d0
  Numinteg_param%dim=dim
  Numinteg_param%maxpoints=maxpoints
  Numinteg_param%phaseVol=1.d0
  do i=1,dim
    Numinteg_param%phaseVol=Numinteg_param%phaseVol*INTSTEP(i)
  enddo
   
  end subroutine init_numericinteg
  !ES------------------------

end Module initNUMINTEG
!EM--------------------------
