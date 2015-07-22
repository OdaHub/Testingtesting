!BM----------------------------------------------------
module system_module
use types;use assign_objects;use alloc_objects
 type(system_type)::syst 
end module system_module 
!EM----------------------------------------------------

!BM----------------------------------------------------
module timestep_module
use types;use assign_objects;use alloc_objects
 type(timestep_type)::timestep
end module timestep_module
!EM----------------------------------------------------

!BM----------------------------------------------------
module phase_module
use types;use assign_objects;use alloc_objects
 type(phasepoint_type)::startpoint
end module phase_module
!EM----------------------------------------------------

!BM----------------------------------------------------
module path_module
use types;use assign_objects;use alloc_objects
 type(path_type)::startpath
end module path_module
!EM----------------------------------------------------

!BM----------------------------------------------------
module pot_module
use types;use assign_objects;use alloc_objects
 type(potential_type)::pot
 type(potential_type)::secpot
end module pot_module
!EM----------------------------------------------------



!BM----------------------------------------------------
module dyn_module
use types;use assign_objects;use alloc_objects
 type(dynamics_type)::dyn
end module dyn_module
!EM----------------------------------------------------

!BM--------------------------------------------------
module output_module
use types;use assign_objects;use alloc_objects
  type(output_type)::output
end module output_module
!EM--------------------------------------------------

!BM---------------------------------------------------
module TIS_module
use types;use assign_objects;use alloc_objects
  type(TIS_type)::TIS_param
end module TIS_module
!EM---------------------------------------------------

!BM---------------------------------------------------
module trans_module
use types;use assign_objects;use alloc_objects
  type(trans_type)::trans_param
end module trans_module
!EM---------------------------------------------------

!BM---------------------------------------------------
module ffs_module
use types;use assign_objects;use alloc_objects
  type(ffs_type)::ffs_param
end module ffs_module
!EM---------------------------------------------------

!BM-----------------------------------------------------
module pathensemble_module
use types;use assign_objects;use alloc_objects
  type(path_ensemble)::PPS_SET
end module pathensemble_module
!EM-----------------------------------------------------

!BM---------------------------------------------------------
module integration_module
use types;use assign_objects;use alloc_objects
  type(numinteg_type)::Numinteg_param
end module integration_module
!EM---------------------------------------------------------
