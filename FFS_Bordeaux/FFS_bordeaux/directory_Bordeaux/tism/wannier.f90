!BM----------------------
Module wannier 

contains

!BS--------------------------------
subroutine out_wan(IU,Npart,electrons)
use types;use assign_objects;use alloc_objects
type(electrons_type), intent(in)::electrons
integer, intent(in)::iu,Npart
integer::Nwan,dwc,i

Nwan=electrons%NWANNIER
dwc=electrons%dwc
do i=1, Nwan
 write(IU,'(i8,a3,10f16.8)') i+Npart,"X",electrons%WCENT(i,1:dwc)
enddo

end subroutine out_wan 
!ES--------------------------------

end Module wannier 
!EM----------------------
