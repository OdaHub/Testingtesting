!BM--------------------------------------
module numericinteg 

contains

  !BS------------------------------------
  subroutine runNI(param,syst,pot)
  use types;use assign_objects;use alloc_objects
  use forcefield
  implicit none
  type(numinteg_type), intent(inout)::param
  type(system_type), intent(in)::syst
  type(potential_type), intent(in)::pot
  integer::i,j,idx,jdx
  double precision::X(1,param%dim)
  double precision::POT_ENERGY,lowest_energy_point

  !for the time being this algorithm works only for dim=2
  if (param%dim>2) then
    print *,"ERROR runNI:" 
    print *,"for the time being this algorithm works only for dim=1,2"
    stop
  endif

  lowest_energy_point=1.d10
  if (param%dim==1) then 
    do idx=0,param%npoints(1)
      X(1,1)=param%range_left(1)+idx*param%intstep(1)
      pot_energy=Epot(x,syst,pot)
      if (pot_energy<lowest_energy_point) lowest_energy_point=pot_energy
    enddo
  else
    do idx=0,param%npoints(1)   
      X(1,1)=param%range_left(1)+idx*param%intstep(1) 
      do jdx=0,param%npoints(2)
        X(1,2)=param%range_left(2)+jdx*param%intstep(2)    
        pot_energy=Epot(x,syst,pot)
        if (pot_energy<lowest_energy_point) lowest_energy_point=pot_energy
      enddo
    enddo
  endif
  print *,"lowest energy",lowest_energy_point
 
  param%prob_array=0.d0
  do i=1,param%dim 
    do idx=0,param%npoints(i)   !range of RC values
      X(1,i)=param%range_left(i)+idx*param%intstep(i)    !so this one is fixed
      if (param%dim>1) then
        do j=1,param%dim            !run over all over dimensions
          if (j==i) cycle           !except the RC which is fixed
          do jdx=0,param%npoints(j) !for each dimension run over all values
            X(1,j)=param%range_left(j)+jdx*param%intstep(j)
            pot_energy=Epot(x,syst,pot)
            pot_energy=pot_energy-lowest_energy_point !to avoid overflow
            param%prob_array(i,idx)=param%prob_array(i,idx)+&
                                    exp(-syst%beta*pot_energy)*param%phaseVol
          enddo
        enddo
      else
         pot_energy=Epot(x,syst,pot)
         pot_energy=pot_energy-lowest_energy_point !to avoid overflow
         param%prob_array(i,idx)=exp(-syst%beta*pot_energy)*param%phaseVol
      endif
    enddo
  enddo 
  param%prob_array(1,:)=param%prob_array(1,:)/sum(param%prob_array(1,:))
  param%prob_array(1,:)=param%prob_array(1,:)/param%intstep(1)
  if (param%dim>1) then
    param%prob_array(2,:)=param%prob_array(2,:)/sum(param%prob_array(2,:))
    param%prob_array(2,:)=param%prob_array(2,:)/param%intstep(2)
  endif

  open(1,file="NUMINT.dat") 
  if (param%dim>1) then
      do i=0,param%maxpoints
        write(1,'(6E16.4)') param%range_left(1)+i*param%intstep(1), &
          param%prob_array(1,i),-log(param%prob_array(1,i))/syst%beta,&
          param%range_left(2)+i*param%intstep(2), &
          param%prob_array(2,i),-log(param%prob_array(2,i))/syst%beta
      enddo
  else
     do i=0,param%maxpoints
        write(1,'(3E16.4)') param%range_left(1)+i*param%intstep(1), &
          param%prob_array(1,i),-log(param%prob_array(1,i))/syst%beta
      enddo
  endif
  close(1)
   
  end subroutine runNI
  !ES------------------------------------

end module numericinteg 
!EM--------------------------------------
