!BM--------------------------------------------------------------------
Module parameters
    

Contains

!BS-------------------------------------------------------------------
subroutine setparam
use inputpar
use read_input_file
implicit none
character(LEN=LSTR)::inputfile

inputfile="input.TISMOL"

call read_inputparameter(TASK           ,inputfile,"TASK"        )
call read_inputparameter(DYNAMICS       ,inputfile,"DYNAMICS"    )
call read_inputparameter(Npart          ,inputfile,"NPART"       )
call read_inputparameter(DIM            ,inputfile,"DIM"         )
call read_inputparameter(Ncyc           ,inputfile,"NCYCLES"     )
call read_inputparameter(Temp           ,inputfile,"TEMP"        )
call read_inputparameter(dt             ,inputfile,"DT"          )
call read_inputparameter(mass           ,inputfile,"MASS"        )
call read_inputparameter(Restart        ,inputfile,"RESTART"     )
call read_inputparameter(NCYCLE_RESTART ,inputfile,"NCYCLE_RESTART"     )
call read_inputparameter(RENUM          ,inputfile,"RENUM"       )
call read_inputparameter(PBC            ,inputfile,"PBC"         )
call read_inputparameter(skipE          ,inputfile,"SKIPE"       )
call read_inputparameter(skipO          ,inputfile,"SKIPO"       )
call read_inputparameter(skipT          ,inputfile,"SKIPT"       )
call read_inputparameter(skipP          ,inputfile,"SKIPP"       )
call read_inputparameter(STARTRANDOM    ,inputfile,"STARTRANDOM" )
call read_inputparameter(INTERFACEL     ,inputfile,"INTERFACEL"  )
call read_inputparameter(INTERFACEM     ,inputfile,"INTERFACEM"  )
call read_inputparameter(INTERFACER     ,inputfile,"INTERFACER"  )
call read_inputparameter(NOPRINT        ,inputfile,"NOPRINT"     )
call read_inputparameter(SKIPINITOPTWAVE,inputfile,"SKIPINITOPTWAVE")
call read_inputparameter(SKIPINITMDSTEP ,inputfile,"SKIPINITMDSTEP" )
call read_inputparameter(NX            ,inputfile,"NX"            )
call read_inputparameter(BOXLENGTH,inputfile,"BOXLENGTH")

call read_inputparameter(REACTION_COORDINATE,inputfile,"REACTION_COORDINATE")
if (REACTION_COORDINATE=="DISTANCE") &
         call read_inputparameter(REACDIST_ATNR ,inputfile,"REACDIST_ATNR"  )

call read_inputparameter(READINITIAL_POSITIONS, &
                                     inputfile,"READINITIAL_POSITIONS")
if (RESTART) call read_inputparameter(RESTARTFILE,inputfile,"RESTARTFILE" )
if (READINITIAL_POSITIONS) then
  allocate(initpos(Npart,dim))
  call read_inputparameter(initpos,inputfile,"POSITIONS")
endif

allocate(ATOM_TYPES(Npart))
call read_inputparameter(ATOM_TYPES,inputfile,"ATOM_TYPES")

call read_inputparameter(masspolydisp,inputfile,"MASSPOLYDISP")
if (masspolydisp) then
  call read_inputparameter(MASS_INPUT_SPECIFICATION,inputfile,"MASS_INPUT_SPECIFICATION")
  select case(MASS_INPUT_SPECIFICATION)
    case("AMU")
      allocate(masses(Npart))
      call read_inputparameter(masses,inputfile,"MASSES")
    case("ATOM_TYPES")
      !!allocate(ATOM_TYPES(Npart))
      !!call read_inputparameter(ATOM_TYPES,inputfile,"ATOM_TYPES")
    case DEFAULT
      print *,"ERROR setpar MASS_INPUT_SPECIFICATION=" ,MASS_INPUT_SPECIFICATION
      stop
  end select
endif

call read_potential_par(inputfile)

select case(DYNAMICS)
  case("LANGEVIN")
    call setparLANGEVIN(inputfile)
  case("NVE")
    call read_inputparameter(RESCALE_ENERGY,inputfile,"RESCALE_ENERGY")
    call read_inputparameter(NVE_ENERGY,inputfile,"NVE_ENERGY")
    call read_inputparameter(SET_LINMOM_ZERO,inputfile,"SET_LINMOM_ZERO")
  case("ANDERSEN")
    call setparAndersenNVT(inputfile)
  case("NOSEHOOVER")
    call setparNHNVT(inputfile)
  case("EXTERNAL")
    !nothing more 
  case("RANDOMWALK")
    if (POTENTIAL/="NONE") then
       print *,"RANDOMWALK REQUIRES POTENTIAL=NONE"
       stop
    endif 
  case default
    print *,"ERROR: SETPARAM DYNAMICS=",DYNAMICS
    stop
end select

select case (TASK)
  case("TIS")
    call setparTIS(inputfile)
  case("PPTIS")
    call setparTIS(inputfile)
  case("MD")
    continue
  case("PPS")
    call setparPPS(inputfile)
  case("TRANSMISSION")
    call setparTrans(inputfile) 
  case("FFS")
    call setparFFS(inputfile)
  case("NUMERICINTEG")
    call setparNumInt(inputfile)
  case DEFAULT
    print *,"SET PARAM ERROR TASK=",TASK
    stop
end select

end subroutine setparam
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine read_potential_par(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile

call read_inputparameter(POTENTIAL      ,inputfile,"POTENTIAL"   )

select case(POTENTIAL)
  case("PBD")
    call setparPBD(inputfile)
  case("HARMOSC")
    call read_inputparameter(kharm,inputfile,"KHARM")
  case("DOUBLEWELL")
    call read_inputparameter(doublewellk4,inputfile,"DOUBLEWELLK4")
    call read_inputparameter(doublewellk2,inputfile,"DOUBLEWELLK2")
  case("2DHARM")
    call read_inputparameter(kxharm,inputfile,"KXHARM")
    call read_inputparameter(kyharm,inputfile,"KYHARM")
    call read_inputparameter(xharm,inputfile,"XHARM")
    call read_inputparameter(yharm,inputfile,"YHARM")
    call read_inputparameter(ycut,inputfile,"YCUT")
  case("NICOLIS1D")
    call setparNICOLIS(inputfile)
  case ("NICOLIS2D")
    call setparNICOLIS(inputfile)
  case ("EXTERNAL")
    call setparEXTPROG(inputfile)
  case ("NONE")
    !just continue
  case("WCA")
    call setparWCA(inputfile)
  case("IONTRANS")
    call setpariontrans(inputfile)
  case default
    print *,"ERROR parameters: POTENTIAL=",POTENTIAL
    stop
end select


end subroutine read_potential_par
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparPBD(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"PBD PARAMETERS:"
  call read_inputparameter(DAT  ,inputfile,"DAT"  )
  call read_inputparameter(DGC  ,inputfile,"DGC"  )
  call read_inputparameter(AAT  ,inputfile,"AAT"  )
  call read_inputparameter(AGC  ,inputfile,"AGC"  )
  call read_inputparameter(S    ,inputfile,"S"    )
  call read_inputparameter(RHO  ,inputfile,"RHO"  )
  call read_inputparameter(ALPHA,inputfile,"ALPHA")
  
  call read_inputparameter(BIAS      ,inputfile,"BIAS"      )
  if (BIAS) then
    call read_inputparameter(BIASEXP   ,inputfile,"BIASEXP"   )
    call read_inputparameter(BIASPREFAC,inputfile,"BIASPREFAC")
    call read_inputparameter(BIASCUTOFF,inputfile,"BIASCUTOFF")
  endif
  

  allocate(seq(Npart))
  call read_inputparameter(seq,inputfile,"SEQUENCE")
  
end subroutine setparPBD
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparWCA(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"WCA PARAMETERS:"
  call read_inputparameter(WCA_EPSILON,inputfile,"WCA_EPSILON")
  call read_inputparameter(WCA_SIGMA  ,inputfile,"WCA_SIGMA"  )
  call read_inputparameter(WCA_H      ,inputfile,"WCA_H"      )
  call read_inputparameter(WCA_W      ,inputfile,"WCA_W"      )

end subroutine setparWCA
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setpariontrans(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"IONTRANS PARAMETERS:"
  call read_inputparameter(IOT_EPSILON,inputfile,"IOT_EPSILON")
  call read_inputparameter(IOT_SIGMA,inputfile,"IOT_SIGMA")
  call read_inputparameter(IOT_EPS2,inputfile,"IOT_EPS2")
  call read_inputparameter(IOT_SIG2,inputfile,"IOT_SIG2")
  call read_inputparameter(IOT_COOPERATIVITY,inputfile,"IOT_COOPERATIVITY")
  call read_inputparameter(IOT_FCOOP,inputfile,"IOT_FCOOP")
  call read_inputparameter(IOT_ND,inputfile,"IOT_ND")
  call read_inputparameter(IOT_NN,inputfile,"IOT_NN")
  call read_inputparameter(IOT_NC,inputfile,"IOT_NC")
  call read_inputparameter(IOT_RCOOP,inputfile,"IOT_RCOOP")


end subroutine setpariontrans
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparNicolis(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"NICOLIS-POTENTIAL PARAMETERS:"
  call read_inputparameter(LAM_NIC,inputfile,"LAM_NIC")
  call read_inputparameter(MU_NIC ,inputfile,"MU_NIC" )
  call read_inputparameter(GAM_NIC,inputfile,"GAM_NIC")

  if (Dynamics/="LANGEVIN") call read_inputparameter(GAMMA,inputfile,"GAMMA")
  !GAMMA is needed for Nicolis potential, therefore it should be read
  !even if Langevin is not applied 


end subroutine setparNicolis
!ES--------------------------------------------------------------------


!BS--------------------------------------------------------------------
subroutine setparLANGEVIN(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"LANGEVIN PARAMETERS:"
  call read_inputparameter(GAMMA                        ,inputfile,"GAMMA"    )
  call read_inputparameter(HIGH_FRICTION_LIMIT,inputfile,"HIGH_FRICTION_LIMIT")

end subroutine setparLANGEVIN
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparAndersenNVT(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"ANDERSEN-NVT PARAMETERS:"
  call read_inputparameter(AN_FREQ,inputfile,"AN_FREQ")
end subroutine setparAndersenNVT
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparNHNVT(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"NOSE HOOVER-NVT PARAMETERS:"
  call read_inputparameter(NTHERMO,inputfile,"NTHERMO")
end subroutine setparNHNVT
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparEXTPROG(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"EXTERNAL PROGRAM PARAMETERS:"
  call read_inputparameter(EXTERNAL_PROGRAM,inputfile,"EXTERNAL_PROGRAM")
  call read_inputparameter(CPMD_UNITS,inputfile,"CPMD_UNITS")
  call read_inputparameter(NCPU,inputfile,"NCPU")
  call read_inputparameter(TWAIT,inputfile,"TWAIT")
  call read_inputparameter(MAXWAIT,inputfile,"MAXWAIT")
  call read_inputparameter(RESUBMIT,inputfile,"RESUBMIT")
  call read_inputparameter(NSUBCYCLES,inputfile,"NSUBCYCLES")
 !! call read_inputparameter(WANNIER_CENTERS,inputfile,"WANNIER_CENTERS")
  call read_inputparameter(NWANNIER,inputfile,"NWANNIER")
  !!call read_inputparameter(BOXLENGTH,inputfile,"BOXLENGTH")
  
end subroutine setparEXTPROG
!ES--------------------------------------------------------------------


!BS--------------------------------------------------------------------
subroutine setparTIS(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"(PP)TIS PARAMETERS:"
  !INTERFACEL,INTERFACEM,INTERFACER are already set in main setparam routine
  call read_inputparameter(NOPS          ,inputfile,"NOPS"          )
  call read_inputparameter(timerevfreq   ,inputfile,"TIMEREVFREQ"   )
  call read_inputparameter(sigdp         ,inputfile,"SIGDP"         ) 
  call read_inputparameter(PATHINFO      ,inputfile,"PATHINFO"      )
  call read_inputparameter(STARTCONDITION,inputfile,"STARTCONDITION")
  call read_inputparameter(BIASV,inputfile,"BIASV")
  
end subroutine setparTIS
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparTrans(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"TRANS PARAMETERS:"
  call read_inputparameter(NOPS    ,inputfile,"NOPS"    )
  call read_inputparameter(NX      ,inputfile,"NX"      )
  call read_inputparameter(PATHINFO,inputfile,"PATHINFO")
  !(Same as in setparTIS )

  call read_inputparameter(NTRANSRUN       ,inputfile,"NTRANSRUN"       )
  call read_inputparameter(TRANSALGORITHM  ,inputfile,"TRANSALGORITHM"  )
  call read_inputparameter(TWO_POINT_METHOD,inputfile,"TWO_POINT_METHOD")

end subroutine setparTrans
!ES--------------------------------------------------------------------

!BS--------------------------------------------------------------------
subroutine setparFFS(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
  print *,"FFS PARAMETERS:"
  call read_inputparameter(NOPS    ,inputfile,"NOPS"    )
  call read_inputparameter(NX      ,inputfile,"NX"      )
  call read_inputparameter(PATHINFO,inputfile,"PATHINFO")
  !(Same as in setparTIS )
end subroutine setparFFS
!ES--------------------------------------------------------------------


!BS--------------------------------------------------------------------
subroutine setparPPS(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
integer::i

print *,"PPS PARAMETERS"
print *,"INCLUDE ALL TIS PARAMETERS"
call setparTIS(inputfile)

call read_inputparameter(NUMINT_PPS    ,inputfile,"NUMINT_PPS"    )
call read_inputparameter(SWAPFREQ      ,inputfile,"SWAPFREQ"      )
call read_inputparameter(NULLMOVES     ,inputfile,"NULLMOVES"     )
call read_inputparameter(SWAPSIMUL     ,inputfile,"SWAPSIMUL"     )
call read_inputparameter(RELATIVESHOOTS,inputfile,"RELATIVESHOOTS")

if (RELATIVESHOOTS) then
  allocate(relative_shootfreq(NUMINT_PPS))
  call read_inputparameter(RELATIVE_SHOOTFREQ,inputfile,"RELATIVE_SHOOTFREQ")
endif

allocate(PPS_INTERFACES(NUMINT_PPS))
allocate(PPS_SIGDP(NUMINT_PPS))
call read_inputparameter(PPS_INTERFACES,inputfile,"PPS_INTERFACES")
call read_inputparameter(PPS_SIGDP     ,inputfile,"PPS_SIGDP"     )
call read_inputparameter(FORCEFIELDMATCHING,inputfile,"FORCEFIELDMATCHING") 
end subroutine setparPPS
!ES--------------------------------------------------------------------

!BS---------------------------------------------------------------------
subroutine setparNumInt(inputfile)
use inputpar
use read_input_file
implicit none
character(LEN=*), intent(in)::inputfile
allocate(range_left(dim));allocate(range_right(dim))
allocate(intstep(dim))
call read_inputparameter(RANGE_LEFT ,inputfile,"RANGE_LEFT" )
call read_inputparameter(RANGE_RIGHT,inputfile,"RANGE_RIGHT")
call read_inputparameter(INTSTEP    ,inputfile,"INTSTEP"    )

end subroutine setparNumInt
!ES---------------------------------------------------------------------

end Module parameters
!EM--------------------------------------------------------------------
