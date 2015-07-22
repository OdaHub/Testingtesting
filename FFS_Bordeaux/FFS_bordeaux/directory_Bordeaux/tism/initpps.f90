!BM--------------------------
Module initPPS

Contains

  !BS------------------------
  subroutine init_PPS
  use inputpar 
  use TIS_module
  use path_module
  use pathensemble_module
  use initTIS
  use output_module
  use pot_module
  use types;use assign_objects;use alloc_objects 
  implicit none
  integer::i,j,dwc
  character*3, allocatable::EXT(:)
  character*3::dum
  type(potential_type)::pot_store

  print *,"INITIALIZE PPS"
  print *,"**************"

  dwc=dim+1
  allocate(EXT(Numint_pps))
  do i=1,numint_pps
    if (i-1 <10) then
      write(dum,'(A2,i1)')"00",i-1
    else if (i-1 <100) then
      write(dum,'(A1,i2)')"0",i-1
    else if (i-1 <1000) then
      write(dum,'(i3)')i-1
    endif
    EXT(i)=dum
  enddo
  


  allocate(pps_set%ext(numint_pps))
  pps_set%ext=ext

  pps_set%FORCEFIELDMATCHING=FORCEFIELDMATCHING
  !two different potentials are used: default potential for [0+],[1+] etc
  !2nd potential for [0-]
  if (FORCEFIELDMATCHING) then
    call init2ndpot
  endif
  


  print *,"THIS REQUIRES ",NUMINT_PPS," TIS INITIALIZATIONS"
  call alloc(PPS_SET,Npart,dim,NThermo,NOPS,NX,NUMINT_PPS,NWANNIER,dwc)
  do i=1,numint_pps
    print *,"INITIALIZE PPS PATH NUMBER",i
    print *,"------------------"
    !perform initialization of all TIS pathways
    PPS_SET%PPS_interfaces(i) = PPS_interfaces(i)
    !Set startcondition,sigdp,interfaceL,interfaceM,interfaceR
    !initialize TIS and save PATH,PPS_sigdp_sqrtm/PPS_aimless   
    
    if (i==1) then 
      startcondition="RIGHT"
      INTERFACEL=-99999999.d0
      INTERFACEM=PPS_interfaces(1)
      INTERFACER=PPS_interfaces(1)
      if (FORCEFIELDMATCHING) then 
        pot_store=pot
        pot=secpot
      endif
    else 
      startcondition="LEFT"
      INTERFACEL=PPS_interfaces(1)
      INTERFACEM=PPS_interfaces(i-1)
      INTERFACER=PPS_interfaces(NUMINT_PPS)
    endif 
    SIGDP=PPS_SIGDP(i)
    restartfile=EXT(i)//"/TISMOL.RESTART"
    output%input_CPMD_optwav=EXT(i)//"/input.CPMD_optwav"
    output%output_CPMD_optwav=EXT(i)//"/output.CPMD_optwav"
    output%input_CPMD_MD     =EXT(i)//"/input.CPMD_MD"
    output%output_CPMD_MD    =EXT(i)//"/output.CPMD_MD"
    output%input_CPMD_TIS=EXT(i)//"/input.CPMD_TIS"
    output%output_CPMD_TIS=EXT(i)//"/output.CPMD_TIS"
    output%GEOMETRY00=EXT(i)//"/GEOMETRY00"
    output%ESTRUC_a=EXT(i)//"/ESTRUC_a"


    call init_TIS
    if ((i==1).AND.(FORCEFIELDMATCHING)) pot=pot_store 

   

    PPS_SET%PPS_aimless(i)=TIS_param%aimless
    do j=1,Npart
      PPS_SET%PPS_sigdp_sqrtm(i,j)=TIS_param%sigdp_sqrtm(j)
    enddo
    PPS_SET%PATHS(i)=startpath
   
    !deallocate some arrays and close files that were openend by init_TIS 
    call dealloc(startpath);call dealloc(TIS_param)
    close(output%IUPATH);close(output%IUEN)
    close(output%IUTRAJ);close(output%IUOP);close(output%IUCR)
    close(output%IUCROSSPOINTS)

  enddo
  
  PPS_SET%timerevfreq=timerevfreq
  PPS_SET%swapfreq=swapfreq
  PPS_SET%RELATIVESHOOTS=RELATIVESHOOTS
  PPS_SET%NULLMOVES=NULLMOVES 
  PPS_SET%SWAPSIMUL=SWAPSIMUL
  if (RELATIVESHOOTS) PPS_SET%RELATIVE_SHOOTFREQ= &
       RELATIVE_SHOOTFREQ/SUM(RELATIVE_SHOOTFREQ)

  allocate(PPS_SET%IUPATH(numint_pps));allocate(PPS_SET%IUEN(numint_pps))
  allocate(PPS_SET%IUTRAJ(numint_pps));allocate(PPS_SET%IUOP(numint_pps))
  allocate(PPS_SET%IUCR(numint_pps));allocate(PPS_SET%IUCROSSPOINTS(numint_pps))
  do i=1,numint_pps
    PPS_SET%IUPATH(i)=1000+i
    PPS_SET%IUEN(i)=  2000+i
    PPS_SET%IUTRAJ(i)=3000+i
    PPS_SET%IUOP(i)  =4000+i
    PPS_SET%IUCR(i)  =5000+i
    PPS_SET%IUCROSSPOINTS(i) =6000+i
    open(PPS_SET%IUPATH(i),file=EXT(i)//"/"//trim(output%PATH_FILE),position="append")
    open(PPS_SET%IUEN(i),file=EXT(i)//"/"//trim(output%ENERGY_FILE),position="append")
    open(PPS_SET%IUTRAJ(i),file=EXT(i)//"/"//trim(output%TRAJECTORY_FILE),position="append")
    open(PPS_SET%IUOP(i),file=EXT(i)//"/"//trim(output%ORDER_PARAM_FILE),position="append")
    open(PPS_SET%IUCR(i),file=EXT(i)//"/"//trim(output%CROSSFILE),position="append")
    open(PPS_SET%IUCROSSPOINTS(i),file=EXT(i)//"/"//trim(output%CROSSPOINTS),position="append")
  enddo
  
  end subroutine init_PPS
  !ES------------------------

  !BS------------------------------------------------------------------------
  subroutine init2ndpot
  use inputpar
  use pot_module
  use parameters
  use initpot
  use types;use assign_objects;use alloc_objects
  implicit none
  type(potential_type)::pot_store
  character(LEN=LSTR)::inputfile
  character(LEN=LSTR)::PSTORE

  if (POTENTIAL=="PBD") then
    print *,"ERROR init2ndpot: OPTION PBD AND FORCEFIELDMATCHING REQUIRES ALLOCATIONS"
    print *,"NOT YET IMPLEMENTED"
    stop 
  else if (POTENTIAL=="EXTERNAL") then
    print *,"ERROR init2ndpot: OPTION EXTERNAL AND FORCEFIELDMATCHING NOT YET AVAILLABLE"
    stop
  endif
  !store old values
  pot_store=pot
  PSTORE=POTENTIAL

  inputfile="input2ndpot.TISMOL"
  !the input file for the 2nd potential
  call read_potential_par(inputfile)
  !init second potential which will overwrite pot in pot_module
  call initpotential
  !now move this new pot to secpot
  secpot=pot
  !reset the old values
  pot=pot_store
  POTENTIAL=PSTORE
  end subroutine init2ndpot
  !ES------------------------------------------------------------------------


end Module initPPS
!EM--------------------------
