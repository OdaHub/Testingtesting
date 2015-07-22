!BM--------------------------
Module initMD

Contains

  !BS------------------------
  subroutine init_MD(skip_create_initial_point)
  use inputpar 
  use phase_module
  use system_module
  use dyn_module
  use shell
  use output_module
  use randomnize
  use random
  use restart_module
  use cpmd_subr
  use pot_module
  use types;use assign_objects;use alloc_objects 
  implicit none
  !!double precision::sig
  logical, optional, intent(in)::skip_create_initial_point
  double precision::x(Npart,dim), v(Npart,dim)
  integer::i,j,NCROSSPLANES
  character(len=LSTR)::ENERGY_FILE,TRAJECTORY_FILE,ORDER_PARAM_FILE,CROSSFILE,&
                       CROSSPOINTS
  character(len=LSTR)::command,comrun
  integer:: IUEN, IUTRAJ, IUOP,IUCR,IUCROSSPOINTS
  character(len=LSTR)::input_CPMD_optwav, output_CPMD_optwav,&
                       input_CPMD_MD,output_CPMD_MD
  logical::R00EXIST  !!,skip
  !!integer::NCRASH


  print *,"INITIALIZE MD PARAMETERS"
 
  input_CPMD_optwav =output%input_CPMD_optwav
  output_CPMD_optwav=output%output_CPMD_optwav
  input_CPMD_MD     =output%input_CPMD_MD
  output_CPMD_MD    =output%output_CPMD_MD

  NCROSSPLANES=3
  call alloc(output,NCROSSPLANES)
  call alloc(startpoint,Npart,dim,Nthermo,Nwannier,dim+1)

  ENERGY_FILE         ="ENERGIES.dat"  ;IUEN  =1000
  TRAJECTORY_FILE     ="TRAJECTORY.dat";IUTRAJ=1001
  ORDER_PARAM_FILE="ORDERPARAM.dat";IUOP  =1002
  CROSSFILE="CROSS.dat";IUCR=1003
  CROSSPOINTS="CROSSPOINTS.dat";IUCROSSPOINTS=1005
  open(IUEN,  file=ENERGY_FILE,position="append")
  open(IUTRAJ,file=TRAJECTORY_FILE,position="append")
  open(IUOP,  file=ORDER_PARAM_FILE,position="append")
  open(IUCR,file=CROSSFILE,position="append")
  open(IUCROSSPOINTS,file=CROSSPOINTS,position="append")
  output%ENERGY_FILE    =ENERGY_FILE;    output%IUEN  =IUEN
  output%TRAJECTORY_FILE=TRAJECTORY_FILE;output%IUTRAJ=IUTRAJ
  output%ORDER_PARAM_FILE=ORDER_PARAM_FILE;output%IUOP=IUOP
  output%CROSSFILE=CROSSFILE;output%IUCR=IUCR
  output%IUCROSSPOINTS=IUCROSSPOINTS
  output%CROSSPOINTS=CROSSPOINTS;output%IUCROSSPOINTS=IUCROSSPOINTS

  output%skipE=skipE
  output%skipO=skipO
  output%skipT=skipT
  
  output%crossing_plane(1)=INTERFACEL
  output%crossing_plane(2)=INTERFACEM
  output%crossing_plane(3)=INTERFACER


  if (present(skip_create_initial_point)) return
  !If this routine is calle inside a loop, we use the return-option 
  !!maybe not so neat. should be changed in the future


  if (RESTART) then
    output%icyclestart=0
    call readRestartMD(RESTARTFILE,startpoint,output%icyclestart)
    if (RENUM) output%icyclestart=0
    if (dyn%DYNAMICS=="EXTERNAL") call make_GEOFILE(startpoint%phasexv%x,&
                                                    startpoint%phasexv%v,&
                                  Npart,syst%inv_a0,syst%inv_atomicv_redv)
  else 

    if (potential=="EXTERNAL") then
      print *,"CREATE START CPMD INPUT FILES"
      command="cat "//trim(output%input_CPMD_optwav_HEAD)//" "//trim(output%input_CPMD_TAIL)// &
            " > "//trim(output%input_CPMD_optwav)
      call do_shell(command)
      command="cat "//trim(output%input_CPMD_MD_HEAD)//" "//trim(output%input_CPMD_TAIL)// &
            " > "//trim(output%input_CPMD_MD)
      call do_shell(command) 

      !Check if RESTART00 already exist. If not, create it
      inquire(file="RESTART00",EXIST=R00EXIST)
      if (.NOT.R00EXIST) then
        print *,"FILE RESTART00 WAS NOT YET CREATED."
        print *,"CREATE RESTART00 NOW"
        command="cat input.CPMD_R00_HEAD "//trim(output%input_CPMD_TAIL)//" > input.CPMD_R00"
        call do_shell(command)
        write(command,'(A11,i3,A80)') &
        "mpirun -np ",NCPU," ./"//trim(EXTERNAL_PROGRAM)//" input.CPMD_R00 > output.CPMD_R00"
        call do_shell(command)
        command="mv RESTART.1 RESTART00" 
        call do_shell(command)
      endif

      !!command="cp -f RESTART00 RESTART"
      !!call do_shell(command)
      open(1,file="LATEST")
        write(1,*) "RESTART00"
        write(1,*) 0
      close(1)

      command="cp -f "//trim(output%GEOMETRY00)//" GEOMETRY"
      call do_shell(command)


      !This RESTART File is not really been used. Still it needs
      !to be provided as the "RESTART GEOFILE" option requires the
      !existence of a RESTART file
      
      write(comrun,'(A11,i3,A20)') "mpirun -np ",NCPU," ./"//trim(EXTERNAL_PROGRAM)
        print *, "START EXTERNAL PROGRAM FOR WAVEFUNCTION OPTIMALIZATION"
        call mpirun_cpmd(comrun,input_cpmd_optwav,output_CPMD_optwav,syst%icrash,&
                         syst%NCPMD_UNSAVED,syst%RESUBMIT)
      print *,"COORDINATES ARE READ FROM GEOMETRY FILE"
      call read_cpmd_coordinates("GEOMETRY",syst,x=x,v=v,com="OPTWAV")

      call do_shell_cpmd("mv RESTART.1 RESTART.OPTWAV",.false.,syst%icrash,syst%NCPMD_UNSAVED)
      open(1,file="LATEST")
        write(1,*) "RESTART.OPTWAV"
        write(1,*) 0
      close(1)


      print *,"NEW VELOCITIES ARE TAKEN FROM A MAXWELLIAN DISTRIBUTION"
      call set_maxwellian_velocities(v,Npart,dim,syst%sigma_v)
      !!print *,"!! set v=0!!!!!!!!!!!!; testing fase"
      !!v=0.d0
      print *,"REWRITE NEW VELOCITIES IN GEOMETRY FILE"
      call make_GEOFILE(x,v,Npart,syst%inv_a0,syst%inv_atomicv_redv)

      startpoint%electrons%ESTRUCFILE=output%ESTRUC_a
      startpoint%electrons%EFILE_INDEX=0
      startpoint%electrons%rev_velec=.false.
        print *, "MAKING 1 VERY SMALL MD STEP"
        !THIS is neccessary to let cpmd know the new velocities
        call mpirun_cpmd(comrun,input_CPMD_MD,output_CPMD_MD,syst%icrash,&
                         syst%NCPMD_UNSAVED,syst%RESUBMIT)
        !read coordinates/velocities from TRAJECTORY file
        call read_cpmd_coordinates("GEOMETRY",syst,phasep=startpoint,com="MDSTEP")
        !rename and save the filename of the electronic structure file
        command="mv -f RESTART.1 "//trim(startpoint%electrons%ESTRUCFILE)//"+000"
        call do_shell_cpmd(command,.false.,syst%icrash,syst%NCPMD_UNSAVED)


      print *, "INITIALIZATION BY EXTERNAL PROGRAM FINISHED"


    else
      startpoint%electrons%ESTRUCFILE="xxx"  !dummy name
      startpoint%electrons%EFILE_INDEX=0     !dummy number
      if (readinitial_positions) then
        x=initpos
      else
        x=0.d0
      endif
      call randomnize_v(v,syst,dyn,x=x,pot=pot)
      startpoint%phasexv%x=x
      startpoint%phasexv%v=v
      output%icyclestart=0
    endif
    output%icyclestart=0
  endif

  !startpoint, output are still allocated
  !following files are still open:
  !IUEN:=ENERGY_FILE, IUTRAJ:=TRAJECTORY_FILE
  !IUOP:=ORDER_PARAM_FILE, IUCR:=CROSSFILE
  !Thus beware when init_MD is called inside a loop!
  print *, "END INIT MD" 
   
  end subroutine init_MD
  !ES------------------------

  !BS-------------------------------------------
  subroutine init_trans
  use inputpar
  use orderparameter
  use phase_module
  use output_module
  use system_module
  use trans_module
  use pot_module
  implicit none
  character(len=LSTR)::TRANSFILE
  integer:: IUTRANS
  double precision::op

  print *,"INITIALIZE TRANSMISSION PARAMETERS"
  TRANSFILE="TRANSMISSION.dat"  ;IUTRANS  =1004
  output%TRANSFILE=TRANSFILE;output%IUTRANS=IUTRANS
  call init_MD
  open(IUTRANS,file=TRANSFILE,position="append")

  op= orderp(startpoint%phasexv%x,syst,startpoint%phasexv%v,pot)
  startpoint%phasexv%x(:,1)=startpoint%phasexv%x(:,1)-op+interfacem
  ! assume translation is possible
  !not generic for all possible op's!

  trans_param%interfaceL=interfaceL
  trans_param%interfaceM=interfaceM
  trans_param%interfaceR=interfaceR
  trans_param%NOPS=NOPS
  trans_param%NX=NX
  trans_param%NTRANSRUN=NTRANSRUN
  trans_param%TRANSALGORITHM=TRANSALGORITHM
  trans_param%TWO_POINT_METHOD=TWO_POINT_METHOD

  output%PATH_FILE="PATH.dat"
  output%IUPATH=2000
  open(output%IUPATH,file=output%PATH_FILE,position="append")
  output%skipp=skipp
  output%PATHINFO=PATHINFO

  if ((INTERFACEL>INTERFACER).OR.(INTERFACEL>INTERFACEM).OR. &
      (INTERFACEM>INTERFACER)) then
    print *,"ERROR INIT_TRANS: POSITION INTERFACES"
    stop
  endif
  

  !startpoint, output are still allocated (via init_MD)
  !following files are still open:
  !IUEN:=ENERGY_FILE, IUTRAJ:=TRAJECTORY_FILE
  !IUOP:=ORDER_PARAM_FILE, IUCR:=CROSSFILE (all via init_MD)
  !IUTRANS:=TRANSFILE (direct via this subroutine)
  !output%IUPATH:=output%PATH_FILE)
 
  end subroutine init_trans
  !ES-------------------------------------------

!BS-------------------------------------------
  subroutine init_ffs
  use inputpar
  use orderparameter
  use phase_module
  use output_module
  use system_module
  use ffs_module
  use pot_module
  implicit none
  character(len=LSTR)::FFSFILE,CROSSPOINTS
  integer:: IUFFS
  double precision::op,x,v
  integer::nstartpoints,status

  
  CROSSPOINTS="CROSSPOINTS.DATA"
  print *,"INITIALIZE FFS PARAMETERS"
  FFSFILE="FFS.dat"  ;IUFFS  =1004
  output%FFSFILE=FFSFILE;output%IUFFS=IUFFS
  call init_MD
  open(IUFFS,file=FFSFILE,position="append")

  op= orderp(startpoint%phasexv%x,syst,startpoint%phasexv%v,pot)
  startpoint%phasexv%x(:,1)=startpoint%phasexv%x(:,1)-op+interfacem
  ffs_param%interfaceL=interfaceL
  ffs_param%interfaceM=interfaceM
  ffs_param%interfaceR=interfaceR
  ffs_param%NOPS=NOPS
  ffs_param%NX=NX


  output%PATH_FILE="PATH.dat"
  output%IUPATH=2000
  open(output%IUPATH,file=output%PATH_FILE,position="append")
  output%skipp=skipp
  output%PATHINFO=PATHINFO

  if ((INTERFACEL>INTERFACER).OR.(INTERFACEL>INTERFACEM).OR. &
      (INTERFACEM>INTERFACER)) then
    print *,"ERROR INIT_FFS: POSITION INTERFACES"
    stop
  endif


  if (syst%dim>1) then
    print *, "FFS-test implementation works only for 1D"
    stop
  endif
  close(output%IUCROSSPOINTS)
  open(output%IUCROSSPOINTS,file=CROSSPOINTS,status="old")
  nstartpoints=0
  do
    read(output%IUCROSSPOINTS,*,iostat=status) x,v 
    if (status /=0 ) exit 
    nstartpoints=nstartpoints+1
    ffs_param%startpoints(nstartpoints,1)=x
    ffs_param%startpoints(nstartpoints,2)=v
    if (nstartpoints==100000) exit
  enddo
  print *,"nstartpoints=",nstartpoints
  ffs_param%nstartpoints=nstartpoints

  end subroutine init_ffs
  !ES-------------------------------------------

end Module initMD
!EM--------------------------
