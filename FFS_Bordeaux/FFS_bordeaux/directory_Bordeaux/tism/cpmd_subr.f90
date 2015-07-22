!BM--------------------------
Module cpmd_subr 

Contains
  
  !BS-------------------------------------------------------
  subroutine kill_external_program(prog)
  use stringlengths
  use shell
  implicit none
  character(LEN=*), intent(in)::prog
  character(LEN=LSTR)::command

    call do_shell("killall "//trim(prog))

  end subroutine kill_external_program
  !ES-------------------------------------------------------

  !BS-------------------------------------------------------
  subroutine delete_CPMD_files
  use stringlengths
  use shell
  implicit none

    call do_shell("rm -f ENERGIES GEOMETRY TRAJECTORY TRAJEC.xyz PSAMPLE EXIT")
  
  end subroutine delete_CPMD_files
  !ES-------------------------------------------------------

  !BS-------------------------------------------------------
  subroutine read_cpmd_coordinates(GFILE,syst,x,v,phasep,com) 
  use types;use assign_objects;use alloc_objects
  implicit none
  character(LEN=*), intent(in)::GFILE
  type(system_type), intent(inout)::syst
  double precision, optional, intent(out)::x(syst%Npart,3),v(syst%Npart,3)
  type(phasepoint_type), optional, intent(inout):: phasep
  character(LEN=*), optional, intent(in)::com
  double precision::xx(syst%Npart,3),vv(syst%Npart,3)
  double precision::WCENT(syst%Npart,syst%dwc)
  integer::i,status,IU_EFILE,IUWANCPMD,NWANNIER,dwc
  character::char
  character*100::geoprint
  !!character*10::char1,char2
  double precision::Eelectrons
  double precision::KIN,POT,ETOT,EHAM,TEMP_inst
    print *,"READ GEOMETRY FILE"
    open(1,file=GFILE)
      do i=1,syst%Npart
        read(1,*) xx(i,1:3),vv(i,1:3)
        write(*,'(6f20.12)') xx(i,1:3),vv(i,1:3) 
        !!xx(i,1:3)=xx(i,1:3)*syst%a0
        !!vv(i,1:3)=vv(i,1:3)*syst%atomicv_redv
      enddo
    close(1)

    if (present(x)) then !simple case, only save coordinates
      x=xx*syst%a0;v=vv*syst%atomicv_redv
      Eelectrons=0.d0;TEMP_inst=0.d0;POT=0.d0;ETOT=0.d0;EHAM=0.d0
      WCENT=0.d0
    else !phasepoint case: save coordinates + energies + WC centers
      phasep%phasexv%x=xx*syst%a0;phasep%phasexv%v=vv*syst%atomicv_redv
       print *,"READ TIS_ENERGIES FILE" 
      open(1,file="TIS_ENERGIES")
        read(1,*) Eelectrons           ;print *,"Eelectrons=",Eelectrons
        read(1,*) TEMP_inst            ;print *,"TEMP_inst=",TEMP_inst
        read(1,*) POT                  ;print *,"POT=",POT
        read(1,*) ETOT                 ;print *,"ETOT=",ETOT
        read(1,*) EHAM                 ;print *,"EHAM=",EHAM
      close(1)
      
      phasep%electrons%KIN=TEMP_inst*syst%kbNFtrans_2
      phasep%electrons%POT=POT*syst%hartree_eV
      phasep%electrons%ETOT=phasep%electrons%KIN+phasep%electrons%POT
      phasep%electrons%EHAM=EHAM*syst%hartree_eV
      phasep%electrons%TEMP_inst=TEMP_inst
      phasep%electrons%Eelectrons=Eelectrons*syst%hartree_eV

      IUWANCPMD=syst%IUWANCPMD
      NWANNIER=syst%NWANNIER
      dwc=syst%dwc
        print *,"READ TIS_WCS FILE"
        open(1,file="TIS_WCS")
          do i=1,NWANNIER
            read(1,*) char, WCENT(i,1:dwc)
            write(*,'(4f20.12)') WCENT(i,1:4)
            phasep%electrons%WCENT(i,1:dwc)=WCENT(i,1:dwc)*syst%a0
            phasep%electrons%WCENT(i,1:3)=&
            phasep%electrons%WCENT(i,1:3)+syst%BOXLENGTH*syst%a0/2.d0
          enddo
        close(1)
    endif  

    syst%icrash=syst%icrash+1
    if (syst%icrash > syst%NCPMD_UNSAVED) then
      open(1,file=syst%crashfile,position="append")
        geoprint=" GEOPRINT CRASHFILE: PSAMPLE"
        if (present(com)) geoprint="GEOPRINT CRASHFILE:"//trim(com)
        write(1,'(i10,A120)') syst%icrash,geoprint
        write(1,'(3i10)') syst%icrash,syst%NPART,syst%NWANNIER
        !write positions and velocities
        do i=1,syst%Npart
          write(1,'(i10,6f20.12)') syst%icrash,xx(i,1:3),vv(i,1:3)
        enddo
        !print energies
        write(1,'(i10,f20.12)') syst%icrash,Eelectrons
        write(1,'(i10,f20.12)') syst%icrash,TEMP_inst
        write(1,'(i10,f20.12)') syst%icrash,POT
        write(1,'(i10,f20.12)') syst%icrash,ETOT
        write(1,'(i10,f20.12)') syst%icrash,EHAM
        !write wannier
        do i=1,syst%NWANNIER
          write(1,'(i10,4f20.12)') syst%icrash,WCENT(i,1:4)
        enddo
      close(1)
    endif
    
    
  end subroutine read_cpmd_coordinates 
  !ES-------------------------------------------------------

 !BS-------------------------------------------------------
  subroutine geo_crashfile(linenr) !syst,phasep)
  !!use types;use assign_objects;use alloc_objects
  implicit none
  integer, intent(in)::linenr
  !!type(system_type), intent(in)::syst
  !!type(phasepoint_type), optional, intent(inout):: phasep
  integer::NWANNIER,Npart,crashnr,i
  double precision::vector(6)
  !character::char
  double precision:: Eelectrons,TEMP_inst,POT,ETOT,EHAM


   !!icrash=syst%icrash
   !!NWANNIER=syst%NWANNIER
   !!dwc=syst%dwc
   !!Npart=syst%Npart
   !!NE=5

   !!NSKIP=icrash*(Npart+NWANNIER+NE+2)+1
   print *,"READING GEOMETRY/ENERGIES/ETC FROM CRASHLINES",linenr

   open(1,file="CRASHFILE")

   !skip first lines
   !!do i=1,Nskip
   !!     read(1,*) char
   !!enddo
   do
     read(1,*) crashnr
     if (crashnr==linenr) exit
   enddo

   read(1,*) crashnr,Npart,NWANNIER
   !read/write positions and velocities    
   open(2,file="GEOMETRY")
   do i=1,Npart
     read(1,*) crashnr,vector(1:6)
     write(2,'(6f20.12)') vector(1:6)
   enddo
   close(2)
 
   !read energies
   read(1,*) crashnr,Eelectrons
   read(1,*) crashnr,TEMP_inst
   read(1,*) crashnr,POT
   read(1,*) crashnr,ETOT
   read(1,*) crashnr,EHAM
   open(2,file="TIS_ENERGIES")
     write(2,*) Eelectrons
     write(2,*) TEMP_inst
     write(2,*) POT
     write(2,*) ETOT
     write(2,*) EHAM
   close(2)

   !read wannier
   open(2,file="TIS_WCS")
   do i=1,NWANNIER
     read(1,*) crashnr,vector(1:4)
     write(2,'(4f20.12)') vector(1:4)
   enddo
   close(2)
   close(1)

  end subroutine geo_crashfile
  !ES-------------------------------------------------------

  !!!BS-------------------------------------------------------
  !!subroutine print_geo_crashfile(syst,phasep)
  !!use types;use assign_objects;use alloc_objects
  !!implicit none
  !!type(system_type), intent(in)::syst
  !!type(phasepoint_type), optional, intent(inout):: phasep
  !!integer::i,NWANNIER,dwc,Npart,icrash
  !!character::char
  !!
  !!icrash=syst%icrash
  !!NWANNIER=syst%NWANNIER
  !!dwc=syst%dwc
  !!Npart=syst%Npart
!!
!!
   !!open(1,file=syst%crashfile,position="append")
   !!write(1,*) "*********** GEOPRINT CRASH INDEX:",icrash
   !!write(1,*) icrash,trim(phasep%electrons%ESTRUCFILE),&
 !!             phasep%electrons%EFILE_INDEX
!!
!!   !write positions and velocities 
!!   do i=1,Npart
!!     write(1,'(i10,6f20.12)') icrash,phasep%phasexv%x(i,1:3),phasep%phasexv%v(i,1:3)
!!   enddo
!!
!!   !print energies
!!   write(1,*) icrash,phasep%electrons%KIN
!!   write(1,*) icrash,phasep%electrons%POT
!!   write(1,*) icrash,phasep%electrons%ETOT
!!   write(1,*) icrash,phasep%electrons%EHAM
!!   write(1,*) icrash,phasep%electrons%TEMP_inst
!!
!!   !write wannier
!!   do i=1,NWANNIER
!!     write(1,'(i10,4f20.12)') icrash,phasep%electrons%WCENT(i,1:dwc)
!!   enddo
!!
!!   close(1)
!!
!!  end subroutine print_geo_crashfile
!!  !ES-------------------------------------------------------

 !BS----------------------------------------------------------
  subroutine mpirun_cpmd(comrun,ifile,ofile,icrash,Nsaved,resub)
  use shell
  implicit none
  character(LEN=*), intent(in)::comrun,ifile,ofile
  integer, intent(in)::icrash,Nsaved
  logical, intent(in)::resub
  character(LEN=100)::command
  !print *,ifile,ofile,"fffffddd"
  write(command,'(a100)') trim(comrun)//" "//trim(ifile)//" > "//trim(ofile)
  if (icrash < Nsaved) then
    print *,"ICRASH=",ICRASH,"<NSAVED=",NSAVED
    write(*,'(A100)')"SKIP:"//trim(command)
    print *,"EXTRACT GEOMETRY FILES FROM CRASHFILE"
    call geo_crashfile(icrash+1)
  else
    call do_shell(command)
    call checkfinish(ofile,resub)
  endif
  end subroutine mpirun_cpmd
  !ES--------------------------------------------------------------------
  
  !BS-----------------------------------------
  subroutine checkfinish(ofile,resub)
  use shell
  implicit none
  character(LEN=*), intent(in)::ofile
  logical, intent(in)::resub
  integer::status
  character::char

  open(1,file=ofile,status="old")
  do
    READ(1,*,iostat=status) char
    if (status /= 0 ) exit
  enddo
  close(1)
  if ((char/="=").AND.(char/="P")) then
   print *,"CPMD NOT PROPERLY FINISHED"
   if (RESUB) then
       print *,"**************************************************************"
       print *,"**********  RESUBMIT TISMOL **********************************"
       print *,"**************************************************************"
       call do_shell("./RESUBMIT")
   endif
   stop
  else
    print *,"CPMD RUN PROPERLY FINISHED"
  endif

  end subroutine checkfinish
  !ES-----------------------------------------


  !BS----------------------------------------------------------
  subroutine writePSAMPLE(NOMORE,dt,Rstart,Rend,i,j,icrash,Nsaved)
  use shell
  implicit none
  integer, intent(in)::NOMORE,i,j,NSAVED,icrash
  character*80, intent(in)::Rstart,Rend
  double precision, intent(in)::dt
  logical::noprint
 
  if (icrash < Nsaved) then
    print *,"ICRASH=",ICRASH,"<NSAVED=",NSAVED
    print *,"SKIP MAKING PSAMPLE TO RUN CPMD"
    print *,"EXTRACT GEOMETRY FILES FROM CRASHFILE"
    call geo_crashfile(icrash+1)
  else
    print *,"WRITING PSAMPLE FILE TO RUN CPMD"
    open(1,file="PSAMPLE")
      write(1,*) NOMORE  ,"NOMORE: number of MD steps"
      write(1,*) dt      , "CP time step (+ or -)"
      write(1,*) trim(Rstart)  , " RESTART FILE"
      write(1,*) trim(Rend)    , " NEW NAME RESTART FILES"
      write(1,*) i       , "NRESTF number of new restart files"
      write(1,*) j       , "position at which RESTART files need to be written"
    close(1)
    call do_shell("cp -f PSAMPLE PSAMPLE_LATEST",noprint)
  endif
  !Beware!:
  !icrash is not incremented. It will be incremented when GEOMETRY files
  !are read. i.e. when CPMD is actuallyu finished with this move 

  end subroutine writePSAMPLE
  !ES----------------------------------------------------------

  !BS-----------------------------------------------------------------
  subroutine make_GEOFILE(x,v,N,A_AU,redv_AU)
  implicit none
  integer, intent(in)::N
  double precision, intent(in)::x(N,3),v(N,3)
  double precision, intent(in)::A_AU,redv_AU
  integer::i
    print *,"MAKING GEOMETRY FILE:"
    open(1,file="GEOMETRY")
      do i=1,N
        write(1,'(6f20.12)') x(i,1:3)*A_AU,v(i,1:3)*redv_AU
        write(*,'(6f20.12)') x(i,1:3)*A_AU,v(i,1:3)*redv_AU
      enddo
    close(1)

  end subroutine make_GEOFILE
  !ES-----------------------------------------------------------------

  !BS-----------------------------------------------------------------
  subroutine delete_Efiles(Efile,revab,icrash,nsaved)
  use shell
  use stringlengths
  implicit none
  character(LEN=*), intent(in)::Efile
  logical, optional, intent(in):: revab 
  integer, intent(inout)::icrash
  integer, intent(in)::nsaved
  character(LEN=LSTR)::Efile2,command
  integer::L

  print *,"DELETE EFILES"
  Efile2=Efile
  if (present(revab)) then
    L=LEN(TRIM(Efile2))
    if (Efile2(L:L)=="a") then
      Efile2(L:L)="b"
    else
      Efile2(L:L)="a"
    endif
  endif
  command="rm -f "//TRIM(Efile2)//"*"
  call do_shell_cpmd(command,.false.,icrash,nsaved)

 
  
  
  end subroutine delete_Efiles 
  !ES-----------------------------------------------------------------

  !BS----------------------------------------------------------------
  Subroutine restartcpmd(prog,NCPU) 
  use stringlengths
  use shell
  implicit none
  character(LEN=*), intent(in)::prog
  integer, intent(in)::NCPU
  character(LEN=LSTR)::command

  call kill_external_program(prog)
  command="cp -fr GEOMETRY GEOMETRY00"
  call do_shell(command)
  write(command,'(A11,i3,A50)') "mpirun -np ",NCPU," ./"//trim(prog)// &
                                " input.CPMD_TIS > output.CPMD_TIS &"
  call do_shell(command)

  end Subroutine restartcpmd
  !ES----------------------------------------------------------------

  !BS------------------------------------------------------------------
  Subroutine get_n_crash(CRASHFILE,NCPMD_UNSAVED)
  use inquire
  implicit none
  character(Len=*), intent(in)::CRASHFILE
  !!integer, intent(in)::NTOT
  integer, intent(out)::NCPMD_UNSAVED
  integer::status !NL

  if (nonempty(crashfile) ) then
    open(1,file=CRASHFILE,status="old")
    do
      READ(1,*,iostat=status) NCPMD_UNSAVED 
      if (status /= 0 ) exit
    enddo
    close(1) 
  else
    NCPMD_UNSAVED=0
  endif
  print *,"NCPMD_UNSAVED=",NCPMD_UNSAVED
  !!f (modulo(NL,(NTOT+7))/=0) then
  !!  print *,"ERROR: CHECK CRASHFILE"
  !!  stop
  !!endif
  !!NCPMD_UNSAVED= NL/(NTOT+7)

  end Subroutine get_n_crash
  !ES-------------------------------------------------------------------
 
  !BS--------------------------------------------------------------------
  subroutine do_shell_cpmd(command,noprint,icrash,Nsaved)
  use shell
  implicit none
  character(LEN=*), intent(in)::command
  logical, intent(in)::noprint
  integer, intent(inout)::icrash
  integer, intent(in)::Nsaved
   
  if (icrash < Nsaved) then
    print *,"ICRASH=",ICRASH,"<NSAVED=",NSAVED
    print *,"SKIP:",trim(command)
  else
    call do_shell(command,noprint)
  endif
  ICRASH=ICRASH+1
  if (icrash > Nsaved) then
    open(1,file="CRASHFILE",position="append")
      write(1,'(i10,A100)') ICRASH, trim(command)
    close(1)
  endif 
  
  end subroutine do_shell_cpmd
  !ES--------------------------------------------------------------------




end Module cpmd_subr 
!EM--------------------------
