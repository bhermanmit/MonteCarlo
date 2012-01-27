module global

  use geometry, only: geometry_type
  use material, only: material_type
  use particle, only: particle_type
  use tally,    only: tally_type

  implicit none
  save

  ! define geometry
  type(geometry_type) :: geo

  ! define material
  type(material_type) :: mat

  ! define tally
  type(tally_type),allocatable :: tal(:)

  ! define particle
  type(particle_type) :: neutron

  ! number of particles to simulate
  integer :: nhist

contains

  subroutine allocate_problem()

    if (.not. allocated(tal)) allocate(tal(geo%n_slabs))

  end subroutine allocate_problem

  subroutine free_memory()

    if (allocated(tal)) deallocate(tal)

  end subroutine free_memory

end module global
