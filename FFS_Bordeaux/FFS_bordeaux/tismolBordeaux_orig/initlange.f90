!BM--------------------------
Module initLange

Contains

  !BS------------------------
  subroutine init_Langevin
  use inputpar 
  use system_module
  use dyn_module
  implicit none
  double precision::gdt,ex
  double precision::beta
  double precision::sig_r,sig_v,c_rv,s11,s22,s12,LS
  double precision::c0,c1,c2
  integer::i
  
  print *,"initialize Langevin parameters"

  Nthermo=0
  NWANNIER=0
  beta=syst%beta
 
  dyn%Langevin%GAMMA=GAMMA 
  dyn%Langevin%high_friction_limit=high_friction_limit
  if (high_friction_limit) then
    allocate(dyn%langevin%sigma_lange(Npart))
    allocate(dyn%langevin%bDdt(Npart))
    dyn%langevin%sigma_lange(1:Npart)=sqrt(2*dt/(masses(1:Npart)*beta*gamma))
    dyn%langevin%bDdt(1:Npart)=dt/(masses(1:Npart)*gamma)
  else
    allocate(dyn%langevin%s12os11(Npart))
    allocate(dyn%langevin%sqrts11(Npart))
    allocate(dyn%langevin%sqrtSos11(Npart))
    allocate(dyn%langevin%a2(Npart))
    allocate(dyn%langevin%b1(Npart))
    allocate(dyn%langevin%b2(Npart))

    gdt=dt*GAMMA
    ex=exp(-gdt)
    do i=1,Npart
      s11 = (dt/(beta*masses(i)*gamma)) * (2. - (3. - 4.*ex + ex**2)/gdt)
      sig_r = sqrt(s11)
      s22 = (1.0 - ex**2)/(beta*masses(i))
      sig_v = sqrt(s22) 
      c_rv  = (1.d0/(beta*masses(i)*gamma))* (1.0 - ex)**2/( sig_r*sig_v)
      s12 =c_rv * sig_r * sig_v
      LS = s11 * s22 - s12**2
      dyn%langevin%s12os11(i) = s12 / s11
      dyn%langevin%sqrts11(i)   = sig_r
      dyn%langevin%sqrtSos11(i) = sqrt(LS/s11) 
      if (gamma>0.d0) then
        c0=ex
        c1=(1.d0-c0)/gdt
        c2=(1.d0-c1)/gdt
      else
        c0=1.d0
        c1=1.d0
        c2=0.5d0
      endif  
      dyn%langevin%c0=c0
      dyn%langevin%a1=c1*dt    
      dyn%langevin%a2(i)=c2*dt**2/masses(i)
      dyn%langevin%b1(i)=(c1-c2)*dt/masses(i)
      dyn%langevin%b2(i)=c2*dt/masses(i)
    enddo
  endif

  end subroutine init_Langevin
  !ES------------------------


end Module initLange
!EM--------------------------
