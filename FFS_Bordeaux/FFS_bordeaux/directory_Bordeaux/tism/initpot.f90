!BM-----------------------------------------------------
Module initpot
implicit none

CONTAINS

!BS--------------------------------------------------------
subroutine initpotential
use inputpar
use pot_module
use initPBD
use initWCA
use initiot
use initHARM
use initEXTERNAL
implicit none

  select case(POTENTIAL)
    case("PBD")
      call initializePBD
    case("HARMOSC")
      call initializeHO
    case("DOUBLEWELL")
      call initializeDW
    case("2DHARM")
      call initializeHO2D
    case("NICOLIS1D")
      call initializeNICOLIS
    case("NICOLIS2D")
      call initializeNICOLIS
    case("EXTERNAL")
      call initializeEXTERNAL
    case("NONE")
      !just continue
    case("WCA")
      call initializeWCA
    case("IONTRANS")
      call initializeIONTRANS
    case default
      print *,"INITIALIZE: NO POTENTIAL ",POTENTIAL
      stop
  end select


end subroutine initpotential
!ES---------------------------------------------------------


end module initpot
!EM-----------------------------------------------------
