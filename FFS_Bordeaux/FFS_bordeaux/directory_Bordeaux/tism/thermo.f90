!BM-----------------------------------------
module thermo

contains

!BS------------------------------------------
subroutine Andersen_velocity_change(v,N,dim,sigma_v,freqdt)
use random
implicit none
integer, intent(in)::N,dim
double precision, intent(inout)::v(N,dim)
double precision, intent(in)::sigma_v(N),freqdt
integer::i,j
double precision::ran,sig

do i=1,N
  sig=sigma_v(i)
  do j=1,dim
    !!call random_number(ran)
    ran=random01()
    if (ran< Freqdt) then
       v(i,j)= rangaussian(sig)
    endif
  enddo
enddo

end subroutine Andersen_velocity_change
!ES------------------------------------------

!BS---------------------------------------------------------
subroutine chainNose(V,extcoord,N,dim,KIN,NH_param,syst,timestep)
use types;use assign_objects;use alloc_objects
implicit none
integer,                         intent(in)   ::N,dim
double precision,                intent(inout)::V(N,dim)
type(extended_coordinates_type), intent(inout)::extcoord
double precision,                intent(inout)::KIN
type(Nose_Hoover_type),          intent(in)   ::NH_param
type(system_type),               intent(in)   ::syst
type(timestep_type),             intent(in)   ::timestep
double precision::G(NH_param%Nthermo),inv_Q_Nose(NH_param%Nthermo)
double precision::Q_Nose(NH_param%Nthermo),scale
double precision::vxi(NH_param%Nthermo),xi(NH_param%Nthermo)
double precision::dt_2,dt_4,dt_8
double precision::NFkbT,kBT
integer::NT,i

NT=NH_param%Nthermo
Q_Nose=NH_param%Q_Nose
inv_Q_Nose=NH_param%inv_Q_Nose
dt_2=timestep%dt_2
dt_4=timestep%dt_4
dt_8=timestep%dt_8
NFkbT=syst%NFkbT
kbT=syst%kbT

xi=extcoord%xi
vxi=extcoord%vxi

G(NT)=(Q_Nose(NT-1)*vxi(NT-1)**2-kbT)*inv_Q_Nose(NT)
vxi(NT)= vxi(NT)+  G(NT)*dt_4
do i=1,NT-2
   vxi(NT-i)=vxi(NT-i)*exp(-dt_8*vxi(NT-i+1))
   G(NT-i)=(Q_Nose(NT-i-1)*vxi(NT-i-1)**2-kbT)*inv_Q_Nose(NT-i)
   vxi(NT-i)= vxi(NT-i)+  G(NT-i)*dt_4
   vxi(NT-i)=vxi(NT-i)*exp(-dt_8*vxi(NT-i+1))
enddo
vxi(1)=vxi(1)*exp(-dt_8*vxi(2))
G(1)=(2*KIN-NFkbt)*inv_Q_Nose(1)
vxi(1)= vxi(1)+  G(1)*dt_4
vxi(1)=vxi(1)*exp(-dt_8*vxi(2))
xi(:)=xi(:)+dt_2*vxi(:)

 SCALE=exp(-vxi(1)*dt_2)
 v=SCALE*v
 KIN=KIN*SCALE**2

 vxi(1)=vxi(1)*exp(-dt_8*vxi(2))
 G(1)=(2*KIN-NFkbt)*inv_Q_Nose(1)
 vxi(1)= vxi(1)+  G(1)*dt_4
 vxi(1)=vxi(1)*exp(-dt_8*vxi(2))

 do i=2,NT-1
    vxi(i)=vxi(i)*exp(-dt_8*vxi(i+1))
    G(i)=(Q_Nose(i-1)*vxi(i-1)**2-kbT)*inv_Q_Nose(i)
    vxi(i)= vxi(i)+  G(i)*dt_4
    vxi(i)=vxi(i)*exp(-dt_8*vxi(i+1))
 enddo

 G(NT)=(Q_Nose(NT-1)*vxi(NT-1)**2-kbT)*inv_Q_Nose(NT)
 vxi(NT)= vxi(NT)+  G(NT)*dt_4

extcoord%xi=xi
extcoord%vxi=vxi

end subroutine chainNose
!ES--------------------------------------------------------------------------

!BS---------------------------------------------------------
subroutine pos_vel_Nose(X,V,syst,pot,timestep)
use types;use assign_objects;use alloc_objects
use forcefield
implicit none
type(system_type),    intent(in)   ::syst
double precision,     intent(inout)::X(syst%Npart,syst%dim)
double precision,     intent(inout)::V(syst%Npart,syst%dim)
type(potential_type), intent(in)   ::pot
type(timestep_type), intent(in)    ::timestep
double precision::F(syst%Npart,syst%dim)
double precision::dt_2,dt_m(syst%Npart)
integer::i,j,N,dim


  dt_2=timestep%dt_2
  dt_m=timestep%dt_m
  N=syst%Npart
  dim=syst%dim
 
  X=X+dt_2*V

  F=FORCE(X,syst,pot)
  do i=1,N
    do j=1,dim
      V(i,j)=V(i,j)+F(i,j)*dt_m(i)
    enddo
  enddo
  X=X+V*dt_2
end subroutine pos_vel_Nose
!ES---------------------------------------------------------

end module thermo
!EM------------------------------------------
