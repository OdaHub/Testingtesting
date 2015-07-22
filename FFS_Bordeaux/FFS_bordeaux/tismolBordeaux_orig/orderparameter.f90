!BM---------------------------------------
module orderparameter 

contains

!BF------------------------------------------
  function orderp(x,syst,v,pot)
  use types;use assign_objects;use alloc_objects
  use arimetric
  implicit none
  type(system_type),       intent(in)::syst
  double precision, intent(in)::x(syst%Npart,syst%dim)
  double precision, optional, intent(in)::v(syst%Npart,syst%dim)
  type(potential_type), optional,intent(in)::pot
  double precision::orderp
  double precision::xi(syst%dim),xj(syst%dim),vi(syst%dim),vj(syst%dim)
  integer::i,j,k,iSi,jO,O_array(8),d,N
  double precision::distOSi(8),r,m,E
  
    
    d=syst%dim
    N=syst%Npart
    select case(syst%REACTION_COORDINATE)

      case("MINVAL") 
        orderp=minval(x(:,1))
      case("X11")
        orderp=X(1,1)
      case("DISTANCE") 
        i=syst%REACDIST_ATNR(1);j=syst%REACDIST_ATNR(2)
        xi(1:d)=x(i,1:d)
        xj(1:d)=x(j,1:d)
        orderp=distance(xi,xj,d,syst%BOXLENGTH)
      case("SILICWAT64")
      !The Si(OH)4+64 H2O cluster consists of 2 Si with 
      !atomnr 209 and 210. Si209 is connected by oxygens 
      !O137,O138,O139,O140. Si210 is connected by oxygens
      !O141, O142,O143,O144. Compare all distance with Si and
      !the opposite oxygens. The smallest distance is the RC  
        k=0
        do i=1,2 
          iSi=208+i
          do j=1,4
            k=k+1
            jO=145-k
            xi(1:d)=x(iSi,1:d)
            xj(1:d)=x(jO,1:d)
            distOSi(k)=distance(xi,xj,d,syst%BOXLENGTH)
          enddo
        enddo
        orderp=minval(distOSi)
      case("SILICWAT64ANION")
      !The SiO4H3- +64 H2O cluster consists of 2 Si with
      !atomnr 208 (negatively charged) and 209. Si208 is connected to oxygens
      !O136,O139,O140,O143. Si209 is connected by oxygens
      !O137,O138,O141,O142. Compare all distance with Si and
      !the opposite oxygens. The smallest distance is the RC
        k=0
        O_array=(/ 137,138,141,142,136,139,140,143 /)
        do i=1,2
          iSi=207+i
          do j=1,4
            k=k+1
            jO=O_array(k)
            xi(1:d)=x(iSi,1:d)
            xj(1:d)=x(jO,1:d)
            distOSi(k)=distance(xi,xj,d,syst%BOXLENGTH)
          enddo
        enddo
        orderp=minval(distOSi)

      case("WCAJCP1")
        xi(:)=x(1,:)
        xj(:)=x(2,:)
        vi(:)=v(1,:)
        vj(:)=v(2,:)
        m=syst%masses(1)
        r=distance(xi,xj,d,syst%BOXLENGTH)
        if (r<1.2) then
          E=Edw(xi,xj,vi,vj,r,d,m,pot) 
          orderp=1.19d0
          if (E<1.5) orderp=1.18d0-((1.5d0-E)/0.5d0)*0.02
        else if (r>1.42) then
          E=Edw(xi,xj,vi,vj,r,d,m,pot)
          orderp=1.43
          if (E<5.d0) orderp=1.44d0+((5.d0-E)/0.5d0)*0.02 
        else
          orderp=r
        endif
      case("IONTRANS")
        xi(:)=x(1,:)
        xj(:)=x(3,:)
        r=distance(xi,xj,d,syst%BOXLENGTH)
        if (r>.4D0) then
          orderp=-r
        else
          xj(:)=x(2,:)
          r=distance(xi,xj,d,syst%BOXLENGTH)
          orderp=+r
        endif
      case DEFAULT
        print *,"ERROR orderp: REACTION_COORDINATE=",syst%REACTION_COORDINATE
        stop
    end select
end function orderp 
!EF------------------------------------------

!BF------------------------------------------
  function v_orderp(x,v,syst)
  use types;use assign_objects;use alloc_objects
  implicit none
  type(system_type),       intent(in)::syst
  double precision, intent(in)::x(syst%Npart,syst%dim)
  double precision, intent(in)::v(syst%Npart,syst%dim)
  double precision::v_orderp
  integer::i(1)
    i=minloc(x(:,1))
    v_orderp=v(i(1),1)
end function v_orderp
!EF------------------------------------------

!BS-------------------------------------------------
  function Edw(xi,xj,vi,vj,r,d,m,pot)
  use types;use assign_objects;use alloc_objects
  use forcefieldwca
  implicit none
  integer, intent(in)::d
  double precision, intent(in)::xi(d),xj(d),vi(d),vj(d),r,m
  type(potential_type),intent(in)::pot
  double precision:: Edw
  double precision:: dx(d),dv(d),dxdt
    dx(:)=xi(:)-xj(:)
    dv(:)=vi(:)-vj(:)
    dxdt=dot_product(dx,dv)/r
    Edw=0.25d0*m*dxdt**2+Vdw(r,pot%potWCA%h,pot%potWCA%w,pot%potWCA%r0)
  end function Edw  
!EF-------------------------------------------------


end module orderparameter
!EM---------------------------------------
