!BM-----------------------------------------------------------
module modswap
implicit none

CONTAINS
!BS-------------------------------------------------------------------
subroutine make_swap(PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)
use types;use assign_objects;use alloc_objects
use propagator
use TISstep
use outTIS
use null
use random
implicit none
type(path_ensemble), intent(inout)::PPS_ensemble
type(TIS_type),       intent(inout)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(inout):: output
type(timestep_type),      intent(in)   ::timestep
integer,                  intent(in)   ::icyc 
type(path_type), intent(inout)::wp1,wp2,wp3,wp4   !workspace
type(timeslice_type), intent(inout)::wt1,wt2      !workspace
double precision::ran
integer::iswap,Lmax,j,startpos,iswap0,iswapend
logical::reverseTF,overlap



if  (PPS_ensemble%SWAPSIMUL) then
  !Make more than one swap simultaneously:
  !Either the swaps [0^-]<->[0^+],[1^+]<->[2^+],... are performed (iswap0=0)
  !or [0^+]<->[1^+],[2^+]<->[3^+],...  (iswap0=1)
  iswap0=0
  if (PPS_ensemble%numint>2) then
    !!call random_number(ran)
    ran=random01()
    if (ran<.5d0) iswap0=1 
  endif
  do iswap=iswap0,PPS_ensemble%numint-2,2
    call swap(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)
    iswapend=iswap
  enddo
  if (PPS_ensemble%nullmoves) then
    !Simply repeat the last accepted path in the datafiles of 
    !the ensembles
    if (iswap0==1) call make_nullmove(1,PPS_ensemble,syst,&
                                     pot,dyn,output,icyc,wp1,secpot) 
    if (iswapend==PPS_ensemble%numint-3) &
      call make_nullmove(PPS_ensemble%numint,PPS_ensemble,syst,&
                                     pot,dyn,output,icyc,wp1,secpot)
  endif
else 
  !!call random_number(ran)
  ran=random01()
  iswap=int(ran*(PPS_ensemble%numint-1))
  !iswap=[0,1,..,numint-2]
  !iswap=0: swap 000/PATH.dat with 001/PATH.dat etc
  !or equivalently PPS_ensemble%Paths(1) with PPS_ensemble%Paths(2)  etc
  call swap(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)
  if (PPS_ensemble%nullmoves) then
    !write outputs for the non-swapped pathways
    do j=1,pps_ensemble%numint
      if ((j==iswap+1).OR.(j==iswap+2)) cycle
      !set path and output for j
      call make_nullmove(j,PPS_ensemble,syst,pot,dyn,output,icyc,wp1,secpot) 
     enddo
   endif
endif

end subroutine make_swap
!ES----------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine swap(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)
use types;use assign_objects;use alloc_objects
use propagator
use TISstep
use outTIS
use convert
implicit none
integer, intent(in)::iswap
type(path_ensemble), intent(inout)::PPS_ensemble
type(TIS_type),       intent(inout)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(inout):: output
type(timestep_type),      intent(in)   ::timestep
integer,                  intent(in)   ::icyc
type(path_type), intent(inout)::wp1,wp2,wp3,wp4   !workspace
type(timeslice_type), intent(inout)::wt1,wt2      !workspace
double precision::ran
integer::Lmax,j,startpos
logical::reverseTF,overlap

if (iswap==0) then
  print *,"SWAPPING 1<->2, [0-] <-> [0+], INTEGRATION NEEDED"
  call swap0(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)

else

  !direct swapping
  print *,"SWAPPING ",iswap+1,"<->",iswap+2,"[",iswap-1,"+] <-> [",iswap,"+]"
  wp1=PPS_ensemble%PATHS(iswap+2)
  wp2=PPS_ensemble%PATHS(iswap+1)
  wp1%mcmove="s+";wp2%mcmove="s-"
  wp1%OPshoot=0.d0;wp1%iOPshoot_old=0;wp1%iOPshoot_new=0
  wp2%OPshoot=0.d0;wp2%iOPshoot_old=0;wp2%iOPshoot_new=0
  wp1%index_shoot=PPS_ensemble%PATHS(iswap+1)%index_shoot !keep the correct indices
  wp2%index_shoot=PPS_ensemble%PATHS(iswap+2)%index_shoot


  call set_TIS(iswap+2,pps_ensemble,output=output,tisparam=tisparam)
  !let's see if wp2 is suited for this ensemble iswap+2
  if (wp2%opmax > tisparam%interfaceM) then
    wp2%ACCREJ="ACC";wp1%ACCREJ="ACC"
    wp2%cross="M"
    wp1%index_acc=PPS_ensemble%PATHS(iswap+1)%index_acc+1
    wp2%index_acc=PPS_ensemble%PATHS(iswap+2)%index_acc+1
    PPS_ensemble%PATHS(iswap+1)=wp1
    PPS_ensemble%PATHS(iswap+2)=wp2

    
    !if wp2 is accepted than wp1 also automatically
    if (pot%potential=="EXTERNAL") call swap_estruc_files(PPS_ensemble,& 
                         iswap+1,iswap+2,syst%icrash,syst%NCPMD_UNSAVED)

  else
    wp2%ACCREJ="NCR";wp1%ACCREJ="NCR"
    wp2%cross="*"
    wp1%index_acc=PPS_ensemble%PATHS(iswap+1)%index_acc
    wp2%index_acc=PPS_ensemble%PATHS(iswap+2)%index_acc
  endif 
  call outputTIS(icyc,wp2,syst,pot,dyn,output)
  call set_TIS(iswap+1,pps_ensemble,output=output)
  call outputTIS(icyc,wp1,syst,pot,dyn,output)
endif

end subroutine swap
!ES--------------------------------------------------------------------------

!BS-------------------------------------------------------------------
subroutine swap0(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)
use types;use assign_objects;use alloc_objects
use propagator
use TISstep
use outTIS
use convert
use forcefield
use modswapfm
implicit none
integer, intent(in)::iswap
type(path_ensemble), intent(inout)::PPS_ensemble
type(TIS_type),       intent(inout)::TISparam
type(system_type),        intent(inout)   ::syst
type(potential_type),     intent(in)   ::pot,secpot
type(dynamics_type),      intent(in)   ::dyn
type(output_type), intent(inout):: output
type(timestep_type),      intent(in)   ::timestep
integer,                  intent(in)   ::icyc
type(path_type), intent(inout)::wp1,wp2,wp3,wp4   !workspace
type(timeslice_type), intent(inout)::wt1,wt2      !workspace
double precision::ran
integer::Lmax,j,startpos
logical::reverseTF,overlap
!double precision::ET0m,ET0p,ET0sm,ET0sp
integer::NX

NX=tisparam%NX

if (PPS_ensemble%forcefieldmatching) then
  call swapfm(iswap,PPS_ensemble,tisparam,syst,pot,dyn,output,timestep,&
                     icyc,wp1,wp2,wp3,wp4,wt1,wt2,secpot)

else

  !Some MD-integration is needed
  call set_TIS(1,pps_ensemble,output=output,tisparam=tisparam)
  !prepare TIS parameters for the [0^-] ensemble
  if (pot%potential=="EXTERNAL")  call swap_estruc_files(PPS_ensemble,1,2,&
                      syst%icrash,syst%NCPMD_UNSAVED,E_EXT=dyn%ext_dyn%EXT)
  
  LMAX=1
  wt1=PPS_ensemble%Paths(2)%timeslices(2)
  !second timeslice of [0+] path
  reverseTF=.false.  !forward
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,tisparam%INTERFACEL,&
                  tisparam%INTERFACER,LMAX,wp1,wt2)
  !no integration is here performed as wt1 is already at the 
  !right side of lambda_0. The timeslice wt1 is simply transformed into
  !a trajectory wp1 of one timeslice 

  reverseTF=.true. !backward
  LMAX=NX-1
  wt1=PPS_ensemble%Paths(2)%timeslices(1)
  !first timeslice of [0+] path
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,tisparam%INTERFACEL,&
                  tisparam%INTERFACER,LMAX,wp3,wt2)
  !backward traj. stored in wp3
  wp4=wp1
  call paste(wp1,wp3,wp4,tisparam%interfaceL,tisparam%interfaceM,&
                         tisparam%interfaceR, overlap=.false.)
  !prepare wp1 for output
  wp1%mcmove="s+"
  wp1%OPshoot=wt1%OPS(1);wp1%iOPshoot_old=1;wp1%iOPshoot_new=wp1%Lpath-1
  wp1%index_acc=pps_ensemble%Paths(1)%index_acc+1
  wp1%index_shoot=pps_ensemble%Paths(1)%index_shoot+1 
  !we increase the shooting index  as 
  ![0^-]<->[0^+] swapping requires integration
  wp1%ACCREJ="ACC"
  if (wp1%LPATH==NX) wp1%ACCREJ="BTX"

  call set_TIS(2,pps_ensemble,output=output,tisparam=tisparam)
  LMAX=1
  startpos=PPS_ensemble%Paths(1)%Lpath-1      !the one but last
  wt1=PPS_ensemble%Paths(1)%timeslices(startpos)
  reverseTF=.true.    !bacward
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,tisparam%INTERFACEL,&
                  tisparam%INTERFACER,LMAX,wp2,wt2)

  !Just like before: No integration occurs here, only timeslice wt1 is
  !transformed in 1-timeslice trajectory wp2
  reverseTF=.false. !forward
  LMAX=NX-1
  startpos=PPS_ensemble%Paths(1)%Lpath  !the last one
  wt1=PPS_ensemble%Paths(1)%timeslices(startpos)
  call onewaytraj(wt1,syst,pot,dyn,timestep,reverseTF,tisparam%INTERFACEL,&
                  tisparam%INTERFACER,LMAX,wp4,wt2)
  !forward traj stored in wp4
  wp3=wp2
  call paste(wp2,wp3,wp4,tisparam%interfaceL,tisparam%interfaceM,&
                         tisparam%interfaceR, overlap=.false.)
  !prepare wp2 for output
  wp2%mcmove="s-"
  wp2%OPshoot=wt1%OPS(1);wp2%iOPshoot_old=startpos;wp2%iOPshoot_new=2
  wp2%index_acc=pps_ensemble%Paths(2)%index_acc+1
  wp2%index_shoot=pps_ensemble%Paths(2)%index_shoot+1 !see previous comment
  wp2%ACCREJ="ACC"
  if (wp2%LPATH==NX) wp1%ACCREJ="FTX"

  !output for wp1 and wp2
  call set_TIS(1,pps_ensemble,output=output)
  call outputTIS(icyc,wp1,syst,pot,dyn,output)
  call set_TIS(2,pps_ensemble,output=output)
  call outputTIS(icyc,wp2,syst,pot,dyn,output)
  if ((wp1%ACCREJ=="ACC").AND.(wp2%ACCREJ=="ACC")) then
    pps_ensemble%Paths(1)=wp1
    pps_ensemble%Paths(2)=wp2
  endif
  
endif

end subroutine swap0
!ES----------------------------------------------------------------------


!BS-----------------------------------------------------------------------
!This routine swaps all the electronic structure files from one directory
!to the other. It also adaps the directory name for all electron-objects
!-------------------------------------------------------------------------
subroutine swap_estruc_files(PPS_set,i,j,icrash,nsaved,E_EXT) !reverse,E_EXT)
use types;use assign_objects;use alloc_objects
use stringlengths
!!use shell
use cpmd_subr
implicit none
type(path_ensemble), intent(inout)::PPS_set
integer, intent(in)::i,j,nsaved
!logical, optional, intent(in):: reverse
integer, intent(inout)::icrash
character*4,optional, intent(in)::E_EXT(-PPS_set%NX:PPS_set%NX)
character(LEN=3)::dir1,dir2
character(LEN=LSTR)::command,Ebase,Ej1st,Ej2nd,Eilast,Eiobl
integer::k,L
character*4::Eindex


print *,"SWAP ESTRUC FILES ENSMS",i,j

dir1=PPS_SET%EXT(i)
dir2=PPS_SET%EXT(j)

if (i==1) then
  !only the first two  ESTRUCFILES of j and the last two of i
  !need to be swapped. Other ESTRUCFILES can be deleted
  Ebase=PPS_SET%PATHS(j)%timeslices(1)%phasepoint%electrons%ESTRUCFILE 
  Eindex=E_EXT(PPS_SET%PATHS(j)%timeslices(1)%phasepoint%electrons%EFILE_INDEX)
  Ej1st=trim(Ebase)//trim(Eindex)
  Ebase=PPS_SET%PATHS(j)%timeslices(2)%phasepoint%electrons%ESTRUCFILE
  Eindex=E_EXT(PPS_SET%PATHS(j)%timeslices(2)%phasepoint%electrons%EFILE_INDEX)
  Ej2nd=trim(Ebase)//trim(Eindex)
  L=PPS_SET%PATHS(1)%Lpath
  Ebase=PPS_SET%PATHS(1)%timeslices(L)%phasepoint%electrons%ESTRUCFILE
  Eindex=E_EXT(PPS_SET%PATHS(1)%timeslices(L)%phasepoint%electrons%EFILE_INDEX)
  Eilast=trim(Ebase)//trim(Eindex)
  Ebase=PPS_SET%PATHS(1)%timeslices(L-1)%phasepoint%electrons%ESTRUCFILE
  Eindex=E_EXT(PPS_SET%PATHS(1)%timeslices(L-1)%phasepoint%electrons%EFILE_INDEX)
  Eiobl=trim(Ebase)//trim(Eindex)


  !mv these two ESTRUC FILES to parent directory
  command="mv -f "//trim(Ej1st)//" ."
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  command="mv -f "//trim(Ej2nd)//" ."
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  !remove all other ESTRUCFILES from this directory
  command="rm -f "//trim(dir2)//"/ESTRUC_*"
  call do_shell_cpmd(command,.false.,icrash,nsaved) 
  !mv the two ESTRUC FILES of dir1 to dir2 
  command="mv -f "//trim(Eilast)//" "//trim(dir2)//"/."
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  command="mv -f "//trim(Eiobl)//" "//trim(dir2)//"/."
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  !remove all other ESTRUCFILES from this directory
  command="rm -f "//trim(dir1)//"/ESTRUC_*"
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  !mv ESTRUCFILES PARENT DIRECTORY to 000/.
  command="mv -f ESTRUC_* "//trim(dir1)//"/."
  call do_shell_cpmd(command,.false.,icrash,nsaved)

  !Change first 3 characters of ESTRUCFILE so that it contains
  !the correct new directory
  PPS_SET%PATHS(1)%timeslices(L  )%phasepoint%electrons% &
                          ESTRUCFILE(1:3)=PPS_SET%EXT(j)
  PPS_SET%PATHS(1)%timeslices(L-1)%phasepoint%electrons% &
                            ESTRUCFILE(1:3)=PPS_SET%EXT(j)
  PPS_SET%PATHS(j)%timeslices(1)%phasepoint%electrons% &
                          ESTRUCFILE(1:3)=PPS_SET%EXT(1)
  PPS_SET%PATHS(j)%timeslices(2)%phasepoint%electrons%&
                         ESTRUCFILE(1:3)=PPS_SET%EXT(1)
else
  command="mv -f "//trim(dir1)//"/ESTRUC_* ."
  !mv all ESTRUC-files from dir1 to parent directory
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  command="mv -f "//trim(dir2)//"/ESTRUC_* "//trim(dir1)//"/."
  !mv all ESTRUC-files from dir2 to dir1
  call do_shell_cpmd(command,.false.,icrash,nsaved)
  command="mv -f ESTRUC_* "//trim(dir2)//"/."
  call do_shell_cpmd(command,.false.,icrash,nsaved)

  do k=1,PPS_SET%NX
    !ESTRUCFILE is of the form 000/ESTRUC_a, 002/ESTRUC_b etc
    !First 3 characters should correspond to new directory
    PPS_SET%PATHS(i)%timeslices(k)%phasepoint%electrons% &
                     ESTRUCFILE(1:3)=PPS_SET%EXT(i)
    PPS_SET%PATHS(j)%timeslices(k)%phasepoint%electrons% &
                     ESTRUCFILE(1:3)=PPS_SET%EXT(j)
  enddo
endif

end subroutine swap_estruc_files
!ES-----------------------------------------------------------------

end module modswap
!EM-----------------------------------------------------------

