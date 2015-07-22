!BM--------------------------
Module initEXTERNAL

Contains

  !BS------------------------
  subroutine initializeEXTERNAL
  use stringlengths
  use inputpar 
  use pot_module
  use system_module
  use cpmd_subr
  use random
  use shell
  use output_module
  implicit none
  character(LEN=LSTR)::command
  double precision::a0,atomicv_redv,inv_a0,inv_atomicv_redv,hartree_eV
  integer::IUENCPMD,IUWANCPMD !!,NTOT
  double precision::x(Npart,dim),v(Npart,dim)
  
  print *, "INITIALIZE EXTERNAL PROGRAM PARAMETERS"
  if (DIM/=3) then
    print *,"dim=",dim
    print *,"INCOMPATIBLE WITH POTENTIAL EXTERNAL"
    stop
  endif
  Nthermo=0

  syst%EXTERNAL_PROGRAM=EXTERNAL_PROGRAM
  syst%NCPU=NCPU
  syst%NSUBCYCLES=NSUBCYCLES
  syst%twait=TWAIT
  syst%MAXWAIT=MAXWAIT
  syst%RESUBMIT=RESUBMIT
  select case (CPMD_UNITS)
    case("BOHR")
      a0=0.5291772108        !Bohr in Angstrom  
    case("ANGSTROM") 
      a0=1.d0
    case default
      print *,"ERROR initexternal.f90: CPMD_UNITS=",CPMD_UNITS;stop
  end select
  !atomic unit of time is 2.418884326505D-17 sec =0.02418884326505 fs 
  !unit time in [eV,A,amu]: sqrt(amu/eV)*A=10.1805056 fs
  atomicv_redv=(10.1805056/0.02418884326505)*a0
  !conversion factor to go from atomic velocities (as used in CPMD)
  !to reduced unit velocity (this program)
  hartree_eV=27.21138386

  inv_a0=1.d0/a0;inv_atomicv_redv=1.d0/atomicv_redv
  pot%external%a0=a0
  pot%external%atomicv_redv=atomicv_redv
  pot%external%inv_a0=inv_a0
  pot%external%inv_atomicv_redv=inv_atomicv_redv
  pot%external%hartree_eV=hartree_eV

  !put them is system as well
  syst%a0=a0
  syst%atomicv_redv=atomicv_redv
  syst%inv_a0=inv_a0
  syst%inv_atomicv_redv=inv_atomicv_redv
  syst%hartree_eV=hartree_eV

  output%input_CPMD_optwav ="input.CPMD_optwav"
  output%output_CPMD_optwav="output.CPMD_optwav"
  output%input_CPMD_MD     ="input.CPMD_MD"
  output%output_CPMD_MD    ="output.CPMD_MD"
  output%input_CPMD_TIS="input.CPMD_TIS"
  output%output_CPMD_TIS="output.CPMD_TIS"
  output%input_CPMD_optwav_HEAD   ="input.CPMD_optwav_HEAD"
  output%input_CPMD_MD_HEAD="input.CPMD_MD_HEAD"
  output%input_CPMD_TIS_HEAD="input.CPMD_TIS_HEAD"
  output%input_CPMD_TAIL="input.CPMD_TAIL"
  output%GEOMETRY00="GEOMETRY00"
  output%ESTRUC_a="ESTRUC_a"

  !kill possible runs of CPMD that are unintentionally running in the background
  !due to an incorrect closures of a previous run
  call kill_external_program(EXTERNAL_PROGRAM)
  call delete_CPMD_files
   if (READINITIAL_POSITIONS) then
     print *,"SORRY, READINITIAL_POSITIONS=",READINITIAL_POSITIONS
     print *,"BUT THIS OPTION IS INVALID WITH CPMD"
     !!print *,"COORDINATES SHOULD BE GIVEN IN input.CPMD_TAIL"
     print *,"COORDINATES SHOULD BE GIVEN IN GEOMETRY00"
     stop
   endif


  if (.NOT.RESTART) then 
     print *,"CREATE PRINCIPAL TIS-CPMD INPUT FILE"
     command="cat "//trim(output%input_CPMD_TIS_HEAD)//" "//trim(output%input_CPMD_TAIL)// &
            " > "//trim(output%input_CPMD_TIS)
      call do_shell(command)


  endif
  IUENCPMD=9000
  open(IUENCPMD,file="ENERGIES")
  syst%IUENCPMD=IUENCPMD
  !STRUCTURE OF CPMD FILE "ENERGIES":
  !NFI  EKINC  TEMPP     EKS      ECLASSIC          EHAM         DIS    TCPU
  !NFI:	Step number (number of finite iterations)
  !EKINC:(fictitious) kinetic energy of the electronic (sub-)system
  !TEMPP:Temperature (= kinetic energy / degrees of freedom) for atoms (ions)
  !EKS:	Kohn-Sham Energy, equivalent to the potential energy in classical MD
  !ECLASSIC:Equivalent to the total energy in a classical MD 
  !(ECLASSIC = EHAM - EKINC)
  !EHAM:total energy, should be conserved
  !DIS:	mean squared displacement of the atoms from the initial coordinates.
  !TCPU:(CPU) time needed for this step. 
  IUWANCPMD=9001
  open(IUWANCPMD,file="WANNIER_CENTER")
  syst%IUWANCPMD=IUWANCPMD
  syst%NWANNIER=NWANNIER
  syst%dwc=DIM+1 !xyz-WCENTERS + SPREADING VALUE
  !!syst%BOXLENGTH=BOXLENGTH 

  syst%IUCRASH=9002
  syst%CRASHFILE="CRASHFILE"
  syst%icrash=0
  !!NTOT=NPART+NWANNIER
  call get_n_crash(syst%CRASHFILE,syst%NCPMD_UNSAVED)
  
  end subroutine initializeEXTERNAL
  !ES------------------------

end Module initEXTERNAL
!EM--------------------------
