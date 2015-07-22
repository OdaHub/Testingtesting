! Program Analyse 
! This program is developed to analyse the output files of the TISMOL.exe
! program 
!---------------------------------------------------------------------------
!BP-------------- PROGRAM ANALYZE.exe -------------------------------------------
!---------------------------------------------------------------------------

  program ANALYSE 
  use param
  use init
  use funcplot
  use MDanalyze 
  use transmission
  use pathanalyze
  use numintanalyze
  use summary
  use match
  use output
  use movies
  use analyze_shoot
  implicit none

  print *," "
  print *," "
  print *,"This is program ANALYSE 27-10 2006:"
  print *,"A program developed to analyses the output files of TISMOL.exe"
  print *," "
    
  print *," "
  print *,"set input parameters"
  call setpar
  call initialize
  
  print *," "
  print *, "making output"
  
  call plot_functions
  call MDout 
  call transout
  call numintout
  call pathout
  call shootout
  call match_hist
  call make_summary
  call makemovies
  call makeoutput
  end program ANALYSE 
!EP--------------------------------------------------------------------------
