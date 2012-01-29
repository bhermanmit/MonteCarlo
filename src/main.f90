!-----------------------------------------------------------------------------!
!> @mainpage A Simple Monte Carlo Code for Slabs
!>
!> @section Overview
!>
!> This is a simple program for performing Monte Carlo particle
!> transport in a homogeneous slab for monoenergetic particles and isotropic
!> scattering in the lab system.
!>
!> @section Compiling
!>
!> Compiling is straightforward with the Makefile.  Just do the following:
!>
!> @verbatim
!>   you@yourcomputer:~/yourfolder/$ make monte
!> @endverbatim
!>
!> The defaults are for debugging.  Just change the fairly obvious variables
!> in the Makefile to suit your needs.
!>
!> @section Running
!>
!> Just as easy!  Do the following and follow the prompts. 
!>
!> @verbatim
!>   you@yourcomputer:~/yourfolder/$ ./monte
!> @endverbatim
!>
!> @section Documentation
!>
!> This is an example of Doxygen documentation.  It's all automated;
!> everything you see here is right in the source code.  Sometimes, the
!> the look and feel of the Doxygen syntax is slightly awkward, but it's
!> worth it when at the end of the day you get automated html and LaTeX 
!> output for all your wonderful modules.  Plus, it can produce call and
!> caller graphs, which helps you keep track of how modules and functions
!> depend on one another---a very important feature if you're given someone
!> else's code to use and modify!
!-----------------------------------------------------------------------------!
program main

  use execute,  only: run_problem, print_tallies
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
