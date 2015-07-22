!BM----------------------
Module outTIS

contains

!BS------------------------------
subroutine outputTIS(icyc,path,syst,pot,dyn,out_param)
use types;use assign_objects;use alloc_objects
use outMD
use wannier
implicit none
integer, intent(in)::icyc
type(path_type), intent(in)::path
type(output_type),       intent(in)::out_param 
type(system_type),       intent(in)::syst  
type(potential_type),    intent(in)::pot
type(dynamics_type),     intent(in)::dyn
integer::ip,IU,imd
logical::WREN,WROP,WRTRAJ

if (modulo(icyc,out_param%skipP)==0) then
  iu=out_param%IUPATH
  write(IU,'(3i11,3a2,i7,2a4,2f15.6,2i7,1f15.6,2i7)') icyc,path%index_acc,&
             path%index_shoot,path%start,path%cross,path%end,path%Lpath,&
             path%ACCREJ,path%MCMOVE,&
             path%opmin,path%opmax,path%iopmin,path%iopmax,&
             path%OPshoot,path%iOPshoot_old,path%iOPshoot_new
endif

WREN=.false.;WROP=.false.;WRTRAJ=.false.
select case (out_param%PATHINFO)
  case ("TRIALS") 
    if (modulo(icyc,out_param%skipE)==0) WREN=.true.
    if (modulo(icyc,out_param%skipO)==0) WROP=.true.
    if (modulo(icyc,out_param%skipT)==0) WRTRAJ=.true.
  case("SHOOTS")
    if ((path%MCmove=="sh").AND.(path%ACCREJ=="ACC")) then
      if (modulo(path%index_shoot,out_param%skipE)==0) WREN=.true.
      if (modulo(path%index_shoot,out_param%skipO)==0) WROP=.true.
      if (modulo(path%index_shoot,out_param%skipT)==0) WRTRAJ=.true.
    endif
  case("ACCEPTED")
    if (path%ACCREJ=="ACC") then
      if (modulo(path%index_shoot,out_param%skipE)==0) WREN=.true.
      if (modulo(path%index_shoot,out_param%skipO)==0) WROP=.true.
      if (modulo(path%index_shoot,out_param%skipT)==0) WRTRAJ=.true.
    endif
  case default
     print *,"ERROR outTIS, out_param%PATHINFO=",out_param%PATHINFO
     stop
end select

if (WREN) then
  IU=out_param%IUEN
  write(iu,*) "#indexes:",icyc,path%index_acc,path%index_shoot
  do imd=1,path%Lpath
    call write_energies(imd,path%timeslices(imd)%phasepoint,syst,pot,dyn,IU)
  enddo
endif

if (WROP) then
  IU=out_param%IUOP
  write(iu,*) "#indexes:",icyc,path%index_acc,path%index_shoot
  do imd=1,path%Lpath
    call write_orderparameters(imd,path%timeslices(imd)%phasepoint%phasexv,& 
                             syst,out_param,pot)
  enddo
endif


if (WRTRAJ) then
  IU=out_param%IUTRAJ
  do imd=1,path%Lpath
    write(IU,*) "#index+time:",icyc,path%index_acc,path%index_shoot,imd
    do ip=1,syst%Npart
      write(IU,'(i8,a3,10f16.8)') ip,syst%atom_types(ip),&
      path%timeslices(imd)%phasepoint%phasexv%x(ip,:), &
      path%timeslices(imd)%phasepoint%phasexv%v(ip,:)
    enddo
    if ((pot%POTENTIAL=="EXTERNAL").AND.(syst%NWANNIER>0)) &
       call out_wan(iu,syst%Npart,path%timeslices(imd)%phasepoint%electrons)
  enddo
endif


end subroutine outputTIS
!ES-------------------------------


end Module outTIS
!EM----------------------
