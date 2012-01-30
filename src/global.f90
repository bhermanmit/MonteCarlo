!==============================================================================!
! MODULE: global
!
!> @author Bryan Herman
!>
!> @brief Global data and problem allocation.
!==============================================================================!
module global

  use geometry, only: geometry_type
  use material, only: material_type
  use particle, only: particle_type
  use tally,    only: tally_type
  use timing,   only: Timer

  implicit none
  save 

  !> define geometry
  type(geometry_type) :: geo

  !> define material
  type(material_type) :: mat

  !> define tally array
  type(tally_type), allocatable :: tal(:), local_tal(:)

  !> the neutron that we follow
  type(particle_type) :: neutron

  !> execution timer
  type(Timer) :: timer_run

!$omp threadprivate ( neutron, local_tal)

  !> number of particles to simulate
  integer :: nhist

  integer, parameter :: I8 = selected_int_kind(18)

contains

  !=============================================================================
  !> @brief Allocate the tallies.
  !=============================================================================
  subroutine allocate_problem()

    if (.not. allocated(tal)) allocate(tal(geo%n_slabs))

  end subroutine allocate_problem

  !=============================================================================
  !> @brief Free the tallies.
  !=============================================================================
  subroutine free_memory()

    if (allocated(tal)) deallocate(tal)

  end subroutine free_memory

end module global
