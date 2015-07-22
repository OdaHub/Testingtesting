!BM--------------------------------------------------------
Module biasvacc
implicit none

CONTAINS

!BS-------------------------------------------------------
subroutine biasv_acceptance(tshoot,path,ishoot,syst,timestep,pot,startcon,accmom)
use types;use assign_objects;use alloc_objects
use orderparameter
implicit none
type(timeslice_type),  intent(in)::tshoot
type(path_type),  intent(in)::path
integer, intent(in)::ishoot
type(system_type),intent(in)::syst
type(timestep_type),      intent(in)   ::timestep
type(potential_type),      intent(in)   ::pot
character, intent(in)   ::startcon
logical, intent(inout)::accmom
double precision::Bigdt,invBigdt
double precision::f1_m(syst%Npart,syst%dim),f2_m(syst%Npart,syst%dim),f_m(syst%Npart,syst%dim)
double precision::x1(syst%Npart,syst%dim),v1(syst%Npart,syst%dim)
double precision::x2(syst%Npart,syst%dim)
double precision::Bigdtf1_m(syst%Npart,syst%dim),Bigdtf2_m(syst%Npart,syst%dim),Bigdtf_m(syst%Npart,syst%dim)
double precision::op1,op2
double precision::opmin,opmax,opminnew,opmaxnew

  print *,"entering biasv: accmom=",accmom


  opmin=path%timeslices(ishoot)%OPS(1)
  if (path%timeslices(ishoot+1)%OPS(1) < opmin ) opmin=path%timeslices(ishoot+1)%OPS(1)
  if (path%timeslices(ishoot-1)%OPS(1) < opmin ) opmin=path%timeslices(ishoot-1)%OPS(1)

  opmax=path%timeslices(ishoot)%OPS(1)
  if (path%timeslices(ishoot+1)%OPS(1) > opmax ) opmax=path%timeslices(ishoot+1)%OPS(1)
  if (path%timeslices(ishoot-1)%OPS(1) > opmax ) opmax=path%timeslices(ishoot-1)%OPS(1)

  opminnew=tshoot%OPS(1)
  opmaxnew=tshoot%OPS(1)

  Bigdt=timestep%Bigdt_ru
  invBigdt=1.d0/Bigdt

  accmom=accmom
  Bigdtf1_m=(path%timeslices(ishoot+1)%phasepoint%phasexv%v-path%timeslices(ishoot)%phasepoint%phasexv%v) !*invBigdt
  Bigdtf2_m=(path%timeslices(ishoot)%phasepoint%phasexv%v-path%timeslices(ishoot-1)%phasepoint%phasexv%v) !*invBigdt
  Bigdtf_m=0.5d0*(Bigdtf1_m+Bigdtf2_m)
  v1=tshoot%phasepoint%phasexv%v+Bigdtf_m
  x1=tshoot%phasepoint%phasexv%x+tshoot%phasepoint%phasexv%v*Bigdt
  x2=x1+0.5d0*Bigdt*Bigdtf_m
  op1=orderp(x1,syst,v1,pot)
  op2=orderp(x2,syst,v1,pot)

  if (op2<opminnew) opminnew=op2
  if (op2>opmaxnew) opmaxnew=op2

  print *,"prediction forward: op1,op2",op1,op2
  v1=tshoot%phasepoint%phasexv%v-Bigdtf_m
  x1=tshoot%phasepoint%phasexv%x-tshoot%phasepoint%phasexv%v*Bigdt
  x2=x1+0.5d0*Bigdt*Bigdtf_m
  op1=orderp(x1,syst,v1,pot)
  op2=orderp(x2,syst,v1,pot)
  print *,"prediction backward: op1,op2",op1,op2
  if (op2<opminnew) opminnew=op2
  if (op2>opmaxnew) opmaxnew=op2
  if (startcon=="L") then
    print *,"opmax,opmaxnew",opmax,opmaxnew
    if (opmaxnew<=opmax) accmom=.false.
  else if (startcon=="R") then 
    print *,"opmin,opminnew",opmin,opminnew
    if (opminnew>=opmin) accmom=.false.
  endif
  print *,"leaving biasv: accmom=",accmom

end subroutine biasv_acceptance
!ES-------------------------------------------------------

end Module biasvacc
!EM--------------------------------------------------------
