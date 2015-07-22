!BM-----------------------------------------------------
Module alloc_objects
use types
implicit none

!INTERFACE ASSIGMENTS/PROCEDURES: define operators on the objects(types)
  !BI------------------------------------------------
  interface alloc
    module procedure alloc_phasexv, alloc_extc, alloc_phase, &
                     alloc_MDinout_param, alloc_timeslice, alloc_path, &
                     alloc_TISparam, alloc_ensemble, alloc_output, &
                     alloc_elec
  end interface
  !EI-----------------------------------------------

  !BI------------------------------------------------
  interface dealloc
    module procedure dealloc_phasexv, dealloc_extc, dealloc_phase, &
                     dealloc_MDinout_param, dealloc_timeslice, dealloc_path, &
                     dealloc_TISparam, dealloc_ensemble, dealloc_output, &
                     dealloc_elec
  end interface
  !EI-----------------------------------------------

contains
  !ALLOCATIONS

  !BS----------------------------------------------------
  subroutine alloc_phasexv(p,N,d)
  implicit none
  integer, intent(in)::N,d
  type(phasexv_type), intent(inout)::p
   allocate(p%x(N,d))
   allocate(p%v(N,d))
   p%N=N
   p%d=d
  end subroutine alloc_phasexv
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine alloc_extc(p,NT)
  implicit none
  integer, intent(in)::NT
  type(extended_coordinates_type), intent(inout)::p
   allocate(p%xi(NT))
   allocate(p%vxi(NT))
   p%NT=NT
  end subroutine alloc_extc
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine alloc_elec(p,NWANNIER,dwc)
  implicit none
  integer, intent(in)::NWANNIER,dwc
  type(electrons_type), intent(inout)::p
   allocate(p%WCENT(NWANNIER,dwc))
   p%NWANNIER=NWANNIER
   p%dwc=dwc
  end subroutine alloc_elec
  !ES---------------------------------------------------


  !BS----------------------------------------------------
  subroutine alloc_phase(p,N,d,NT,NWANNIER,dwc)
  implicit none
  integer, intent(in)::N,d,NT,NWANNIER,dwc
  type(phasepoint_type), intent(inout)::p
   call alloc(p%phasexv,N,d)
   call alloc(p%extcoord,NT)
   call alloc(p%electrons,NWANNIER,dwc)
   p%N=N
   p%d=d
   p%NT=NT
   p%NWANNIER=NWANNIER
   p%dwc=dwc
   end subroutine alloc_phase
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine alloc_MDinout_param(p,N,d)
  implicit none
  integer, intent(in)::N,d
  type(MDinout_param_type), intent(inout)::p
    allocate(p%F(N,d))
    p%N=N
    p%d=d
  end subroutine alloc_MDinout_param
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine alloc_timeslice(p,N,d,NT,NOPS,NWANNIER,dwc)
  implicit none
  integer, intent(in)::N,d,NT,NOPS,NWANNIER,dwc
  type(timeslice_type), intent(inout)::p
   call alloc(p%phasepoint,N,d,NT,NWANNIER,dwc)
   call alloc(p%MDinout_param,N,d)
   allocate(p%OPS(Nops))
   p%N=N
   p%d=d
   p%NT=NT
   p%NOPS=NOPS
   p%NWANNIER=NWANNIER
   p%dwc=dwc
  end subroutine alloc_timeslice
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine alloc_path(p,N,d,NT,NOPS,NX,NWANNIER,dwc)
  implicit none
  integer, intent(in)::N,d,NT,NOPS,NX,NWANNIER,dwc
  type(path_type), intent(inout)::p
  integer::i
  allocate(p%timeslices(NX))
  do i=1,NX
    call alloc(p%timeslices(i),N,d,NT,NOPS,NWANNIER,dwc)
  enddo
  p%N=N
  p%d=d
  p%NT=NT
  p%NOPS=NOPS
  p%NX=NX
  p%NWANNIER=NWANNIER
  p%dwc=dwc
  end subroutine alloc_path
  !ES---------------------------------------------------

  !BS-------------------------------------------------------
  subroutine alloc_output(p,Ncrossplanes)
  implicit none
  integer, intent(in)::Ncrossplanes
  type(output_type), intent(inout)::p
  allocate(p%crossing_plane(Ncrossplanes) )
  p%ncrossplanes=Ncrossplanes
  end subroutine alloc_output
  !ES--------------------------------------------------------

  !BS-----------------------------------------------------
  subroutine alloc_TISparam(p,N)
  implicit none
  integer, intent(in)::N
  type(TIS_type), intent(inout)::p
  allocate(p%sigdp_sqrtm(N))
  p%N=N
  end subroutine alloc_TISparam
  !ES-----------------------------------------------------

  !BS----------------------------------------------------
  subroutine alloc_ensemble(p,N,d,NT,NOPS,NX,NUMINT,NWANNIER,dwc)
  implicit none
  integer, intent(in)::N,d,NT,NOPS,NX,NUMINT,NWANNIER,dwc
  type(path_ensemble), intent(inout)::p
  integer::i
  allocate(p%PATHS(NUMINT))
  do i=1,NUMINT
    call alloc(p%PATHS(i),N,d,NT,NOPS,NX,NWANNIER,dwc) 
  enddo
  allocate(p%PPS_interfaces(NUMINT))
  allocate(p%PPS_sigdp_sqrtm(NUMINT,N))
  allocate(p%PPS_aimless(NUMINT))
  allocate(p%relative_shootfreq(NUMINT))
  allocate(p%IUPATH(NUMINT));allocate(p%IUEN(NUMINT));allocate(p%IUTRAJ(NUMINT))
  allocate(p%IUOP(NUMINT));allocate(p%IUCR(NUMINT))
  allocate(p%IUCROSSPOINTS(NUMINT))

  p%N=N
  p%d=d
  p%NT=NT
  p%NOPS=NOPS
  p%NX=NX
  p%NUMINT=NUMINT
  p%NWANNIER=NWANNIER
  p%dwc=dwc
  end subroutine alloc_ensemble
  !ES---------------------------------------------------

  !DEALLOCATE
  !BS----------------------------------------------------
  subroutine dealloc_phasexv(p)
  implicit none
  type(phasexv_type), intent(inout)::p
   deallocate(p%x)
   deallocate(p%v)
  end subroutine dealloc_phasexv
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine dealloc_extc(p)
  implicit none
  type(extended_coordinates_type), intent(inout)::p
   deallocate(p%xi)
   deallocate(p%vxi)
  end subroutine dealloc_extc
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine dealloc_elec(p)
  implicit none
  type(electrons_type), intent(inout)::p
   deallocate(p%WCENT)
  end subroutine dealloc_elec
  !ES---------------------------------------------------


 !BS----------------------------------------------------
  subroutine dealloc_phase(p)
  implicit none
  type(phasepoint_type), intent(inout)::p
   call dealloc(p%phasexv)
   call dealloc(p%extcoord)
  end subroutine dealloc_phase
  !ES---------------------------------------------------

 !BS----------------------------------------------------
  subroutine dealloc_MDinout_param(p)
  implicit none
  type(MDinout_param_type), intent(inout)::p
    deallocate(p%F)
  end subroutine dealloc_MDinout_param
  !ES---------------------------------------------------

 !BS----------------------------------------------------
  subroutine dealloc_timeslice(p)
  implicit none
  type(timeslice_type), intent(inout)::p
   call dealloc(p%phasepoint)
   call dealloc(p%MDinout_param)
   deallocate(p%OPS)
  end subroutine dealloc_timeslice
  !ES---------------------------------------------------

  !BS----------------------------------------------------
  subroutine dealloc_path(p)
  implicit none
  type(path_type), intent(inout)::p
  integer::i
  do i=1,p%NX
    call dealloc(p%timeslices(i))
  enddo
  deallocate(p%timeslices)
  end subroutine dealloc_path
  !ES---------------------------------------------------

  !BS-------------------------------------------------------
  subroutine dealloc_output(p)
  implicit none
  type(output_type), intent(inout)::p
    deallocate(p%crossing_plane)
  end subroutine dealloc_output
  !ES--------------------------------------------------------

  !BS-----------------------------------------------------
  subroutine dealloc_TISparam(p)
  implicit none
  type(TIS_type), intent(inout)::p
    deallocate(p%sigdp_sqrtm)
  end subroutine dealloc_TISparam
  !ES-----------------------------------------------------

  !BS----------------------------------------------------
  subroutine dealloc_ensemble(p)
  implicit none
  type(path_ensemble), intent(inout)::p
  integer::i,NUMINT
   
  NUMINT=p%NUMINT
  do i=1,NUMINT
    call dealloc(p%PATHS(i))
  enddo
  deallocate(p%paths)
  deallocate(p%PPS_interfaces)
  deallocate(p%PPS_sigdp_sqrtm)
  deallocate(p%PPS_aimless)
  deallocate(p%relative_shootfreq)
  deallocate(p%IUPATH);deallocate(p%IUEN);deallocate(p%IUTRAJ)
  deallocate(p%IUOP);deallocate(p%IUCR);deallocate(p%IUCROSSPOINTS)

  end subroutine dealloc_ensemble
  !ES---------------------------------------------------

end module alloc_objects 
!EM----------------------------------------------
