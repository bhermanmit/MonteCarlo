program main

  use execute,  only: run_problem,print_tallies
  use global
  use geometry, only: read_geometry
  use material, only: read_material

  implicit none

  ! read input from user
  write(*,*) 'How many particles would you like to simulate?'
  read(*,*) nhist
  call read_geometry(geo) 
  call read_material(mat)

  ! allocate the problem
  call allocate_problem()

  ! run problem
  call run_problem()

  ! print results to user
  call print_tallies()

  ! free memory
  call free_memory()

  ! terminate the program
  stop

end program main
