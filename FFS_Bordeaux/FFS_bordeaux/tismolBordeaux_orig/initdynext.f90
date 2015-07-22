!BM--------------------------
Module initdynext

Contains

  !BS------------------------
  subroutine init_dynamics_external
  use inputpar 
  !!use system_module
  use dyn_module
  use charext
  use timestep_module
  implicit none
  
  print *,"INITIALIZE EXTERNAL DYNAMICS"
  Nthermo=0
 
  !!dyn%ext_dyn%ESTRUCFILE="ELECSTRUC" 
  
  !!if (TASK=="MD") then 
  !!  !!dyn%ext_dyn%RENAME_ESTRUCTURE_FILE=.false.
  !!  !new electronic structure file will overwrite the old one
  !!  allocate(dyn%ext_dyn%EXT(0:0))
  !!  call set_char_extensions(dyn%ext_dyn%EXT,0)

  !!else
  !!  !!dyn%ext_dyn%RENAME_ESTRUCTURE_FILE=.true. 
  !!  !new electronic structure file will not overwrite the old one
  !!  !The new file name has the same basis as the previous one,
  !!  !but with a different extension
  !!  if (NX>999) then !Too many electron structure files!
  !!                   !will take too much disk space
  !!    print *,"ERROR init_dynamics_external NX > 999:NX=",NX
  !!    stop
  !!  endif
    
    allocate(dyn%ext_dyn%EXT(-NX:NX))
    call set_char_extensions(dyn%ext_dyn%EXT,NX)
    !!startpos0=.true.
    !Create extension character array EXT (..,"-001","+000","+001",..,)
    
  !!endif
  timestep%Bigdt_ru=dt*nsubcycles*0.02418884326505/10.1805056


  end subroutine init_dynamics_external
  !ES------------------------

end Module initdynext
!EM--------------------------
