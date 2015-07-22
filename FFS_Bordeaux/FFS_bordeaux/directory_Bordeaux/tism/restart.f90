!BM----------------------
Module restart_module 

contains

!BS------------------------------
subroutine makeRestartMD(iMD,phasep,restartfile,iu,syst)
use types;use assign_objects;use alloc_objects
use shell
implicit none
integer, intent(in)::iMD,iu
type(phasepoint_type), intent(in)::phasep
character(len=*), intent(in)::restartfile
type(system_type),   intent(inout)::syst

print *,"MAKE RESTART MD"
open(iu,file=restartfile,form="unformatted")
write(iu) iMD
call writephasep(phasep,iu)
call writerandomseed(iu)
close(iu)

!reset crashfile (only needed for external dynamics)
call do_shell("mv -f "//trim(syst%CRASHFILE)//" CRASHFILE_OLD")
syst%icrash=0
syst%NCPMD_UNSAVED=0

end subroutine MakeRestartMD
!ES-------------------------------

!BS------------------------------
subroutine makeRestartPATH(icyc,path,restartfile,inforestartfile,iu,syst)
use types;use assign_objects;use alloc_objects
use shell
implicit none
integer, intent(in)::icyc,iu
type(path_type), intent(in)::path
character(len=*), intent(in)::restartfile,inforestartfile
type(system_type),   intent(inout)::syst
integer::i,j

open(iu,file=restartfile,form="unformatted")
write(iu) icyc
call writepath(path,iu)
call writerandomseed(iu)
close(iu)

open(1,file=inforestartfile)
write(1,*) "icycle:",icyc
write(1,*) "Length of Path:",  path%Lpath
write(1,*) "ORDERPARAMETER"
do i=1,path%Lpath
  write(1,*) i,path%timeslices(i)%OPS(1)
enddo
write(1,*) "ELECTRONIC STRUCTURE FILES"
do i=1,path%Lpath
  write(1,*) i,trim(path%timeslices(i)%phasepoint%electrons%ESTRUCFILE),&
               path%timeslices(i)%phasepoint%electrons%EFILE_INDEX
enddo
write(1,*) "KIN,POT,ETOT,EHAM,TEMP"
do i=1,path%Lpath
  write(1,'(i4,4E25.10,f15.4)') i,path%timeslices(i)%phasepoint%electrons%KIN,&
               path%timeslices(i)%phasepoint%electrons%POT,&
               path%timeslices(i)%phasepoint%electrons%ETOT,&
               path%timeslices(i)%phasepoint%electrons%EHAM,&
               path%timeslices(i)%phasepoint%electrons%TEMP_inst
enddo
write(1,*) "GEOMETRY.xyz files"
do i=1,path%Lpath
  write(1,*) "timeslice:",i
  do j=1,path%N
    write(1,'(6f20.12)') path%timeslices(i)%phasepoint%phasexv%x(j,1:3),&
                         path%timeslices(i)%phasepoint%phasexv%v(j,1:3)
  enddo
enddo
close(1)

!reset crashfile (only needed for external dynamics)
!!call do_shell("rm -rf "//trim(syst%CRASHFILE))
call do_shell("mv -f "//trim(syst%CRASHFILE)//" CRASHFILE_OLD")
syst%icrash=0
syst%NCPMD_UNSAVED=0


end subroutine MakeRestartPATH
!ES-------------------------------

!BS------------------------------
subroutine makeRestartPPS(icyc,pps_ens,wp1,syst)
use types;use assign_objects;use alloc_objects
use shell
implicit none
integer, intent(in)::icyc
type(path_ensemble), intent(in)::pps_ens
type(system_type),   intent(inout)::syst
type(path_type), intent(inout)::wp1
character(len=LSTR)::rfile,infofile
integer::i

do i=1,pps_ens%numint
  wp1=pps_ens%paths(i)
  rfile=pps_ens%EXT(i)//"/TISMOL.RESTART"
  infofile=pps_ens%EXT(i)//"/INFO.RESTART"
  call makeRestartPATH(icyc,wp1,rfile,infofile,1,syst)
enddo

end subroutine MakeRestartPPS
!ES-------------------------------


!BS------------------------------
subroutine readRestartPATH(icyc,path,restartfile,syst,dynamics)
use types;use assign_objects;use alloc_objects
use shell
implicit none
integer, intent(out)::icyc
type(path_type), intent(out)::path
character(len=*), intent(in)::restartfile
type(system_type),intent(in)::syst
character(LEN=*), intent(in)::DYNAMICS !,EXTERNAL_PROGRAM
character(LEN=LSTR)::command
logical, save::first_visit
data first_visit /.true./


open(1,file=restartfile,form="unformatted")
read(1) icyc
call readpath(path,1)
call readrandomseed(1)
close(1)


 if ((DYNAMICS=="EXTERNAL").AND.(first_visit)) then
    !!maybe not so elegant !!!!
    print *, "STARTING CPMD PATH SAMPLING OPTION"
    write(command,'(A11,i3,A50)') &
      "mpirun -np ",syst%NCPU," ./"//trim(syst%EXTERNAL_PROGRAM)// &
      " input.CPMD_TIS > output.CPMD_TIS &"
      call do_shell(command)
     first_visit=.false.
 endif


end subroutine readRestartPATH
!ES-------------------------------

!BS------------------------------
subroutine ReadRestartMD(restartfile,phasep,iMD)
use types;use assign_objects;use alloc_objects
implicit none
integer, intent(out)::iMD
type(phasepoint_type), intent(out)::phasep
character(len=*), intent(in)::restartfile

open(1,file=restartfile,form="unformatted")
read(1) iMD
call readphasep(phasep,1)
call readrandomseed(1)
close(1)

end subroutine ReadRestartMD
!ES-------------------------------

!BS-----------------------------------------
subroutine writephasep(phasep,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(phasepoint_type), intent(in)::phasep
integer, intent(in)::iu
integer::N,d,NT
integer::i,j
N=phasep%N
d=phasep%d
NT=phasep%NT

do i=1,N
  do j=1,d
    write(iu) phasep%phasexv%x(i,j)
    write(iu) phasep%phasexv%v(i,j) 
  enddo
enddo
do i=1,NT
  write(iu) phasep%extcoord%xi(i)
  write(iu) phasep%extcoord%vxi(i)
enddo
write(iu) phasep%electrons%ESTRUCFILE
write(iu) phasep%electrons%EFILE_INDEX
write(iu) phasep%electrons%KIN
write(iu) phasep%electrons%POT
write(iu) phasep%electrons%ETOT
write(iu) phasep%electrons%EHAM
write(iu) phasep%electrons%TEMP_inst
write(iu) phasep%electrons%rev_velec

end subroutine writephasep
!ES-----------------------------------------

!BS---------------------------------------
subroutine writeMDinout(MDio,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(MDinout_param_type), intent(in)::MDio
integer, intent(in)::iu
integer::N,d,i,j
N=MDio%N
d=MDio%d
do i=1,N
  do j=1,d
    write(iu) MDio%F(i,j)
  enddo
enddo
write(iu) MDio%KIN

end subroutine writeMDinout
!ES---------------------------------------

!BS---------------------------------------
subroutine readMDinout(MDio,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(MDinout_param_type), intent(out)::MDio
integer, intent(in)::iu
integer::N,d,i,j
N=MDio%N
d=MDio%d
do i=1,N
  do j=1,d
    read(1) MDio%F(i,j)
  enddo
enddo
read(1) MDio%KIN

end subroutine readMDinout
!ES---------------------------------------


!BS-----------------------------------------
subroutine writetimeslice(timeslice,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(timeslice_type), intent(in)::timeslice
integer, intent(in)::iu
integer::NOPS
integer::i
NOPS=timeslice%NOPS

call writephasep(timeslice%phasepoint,iu)
call writeMDinout(timeslice%MDinout_param,iu)
do i=1,NOPS
  write(iu) timeslice%OPS(i)
enddo

end subroutine writetimeslice
!ES-----------------------------------------

!BS-----------------------------------------
subroutine readtimeslice(timeslice,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(timeslice_type), intent(out)::timeslice
integer, intent(in)::iu
integer::NOPS
integer::i
NOPS=timeslice%NOPS

call readphasep(timeslice%phasepoint,iu)
call readMDinout(timeslice%MDinout_param,iu)
do i=1,NOPS
  read(1) timeslice%OPS(i)
enddo

end subroutine readtimeslice
!ES-----------------------------------------


!BS-----------------------------------------
subroutine writepath(path,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(path_type), intent(in)::path
integer, intent(in)::iu
integer::NX
integer::i
NX=path%NX

do i=1,NX
  call writetimeslice(path%timeslices(i),iu)
enddo
write(iu) path%Lpath
write(iu) path%opmin
write(iu) path%opmax
write(iu) path%iopmin
write(iu) path%iopmax
write(iu) path%index_acc
write(iu) path%index_shoot
write(iu) path%start
write(iu) path%end
write(iu) path%cross
write(iu) path%ACCREJ
write(iu) path%MCmove
write(iu) path%OPshoot
write(iu) path%iOPshoot_old
write(iu) path%iOPshoot_new

end subroutine writepath
!ES-----------------------------------------

!BS-----------------------------------------
subroutine readpath(path,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(path_type), intent(out)::path
integer, intent(in)::iu
integer::NX
integer::i
NX=path%NX

do i=1,NX
  call readtimeslice(path%timeslices(i),iu)
enddo
read(iu) path%Lpath
read(iu) path%opmin
read(iu) path%opmax
read(iu) path%iopmin
read(iu) path%iopmax
read(iu) path%index_acc
read(iu) path%index_shoot
read(iu) path%start
read(iu) path%end
read(iu) path%cross
read(iu) path%ACCREJ
read(iu) path%MCmove
read(iu) path%OPshoot
read(iu) path%iOPshoot_old
read(iu) path%iOPshoot_new

end subroutine readpath
!ES-----------------------------------------


!BS-----------------------------------------
subroutine readphasep(phasep,iu)
use types;use assign_objects;use alloc_objects
implicit none
type(phasepoint_type), intent(out)::phasep
integer, intent(in)::iu
integer::N,d,NT
integer::i,j
N=phasep%N
d=phasep%d
NT=phasep%NT
do i=1,N
  do j=1,d
    read(iu) phasep%phasexv%x(i,j)
    read(iu) phasep%phasexv%v(i,j)
  enddo
enddo
do i=1,NT
  read(iu) phasep%extcoord%xi(i)
  read(iu) phasep%extcoord%vxi(i)
enddo

read(iu) phasep%electrons%ESTRUCFILE
read(iu) phasep%electrons%EFILE_INDEX
read(iu) phasep%electrons%KIN
read(iu) phasep%electrons%POT
read(iu) phasep%electrons%ETOT
read(iu) phasep%electrons%EHAM
read(iu) phasep%electrons%TEMP_inst
read(iu) phasep%electrons%rev_velec

end subroutine readphasep
!ES-----------------------------------------

!BS-------------------------------------------
subroutine writerandomseed(iu)
implicit none
integer, intent(in)::iu
integer::isze,i
integer, allocatable::iseedarr(:)

call random_seed(size=isze) ! determine length of seed-array
                            ! random number generator
allocate(iseedarr(isze))
write(iu) isze
call random_seed(get=iseedarr(1:isze))
do i=1,isze
  write(iu) iseedarr(i)
enddo

deallocate(iseedarr)
end subroutine writerandomseed
!ES-------------------------------------------

!BS-------------------------------------------
subroutine readrandomseed(iu)
implicit none
integer, intent(in)::iu
integer::isze,i
integer, allocatable::iseedarr(:)

call random_seed(size=isze) ! determine length of seed-array
                            ! random number generator
allocate(iseedarr(isze))
read(iu) isze
do i=1,isze
  read(iu) iseedarr(i)
enddo
call random_seed(put=iseedarr(1:isze))

deallocate(iseedarr)
end subroutine readrandomseed
!ES-------------------------------------------


end Module restart_module
!EM----------------------
