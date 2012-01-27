module particle

  use pdfs

  implicit none
  private
  public :: particle_init

  type, public :: particle_type

    integer :: slab   ! the slab id number
    real    :: xloc   ! location of particle in slab
    real    :: mu     ! cosine of angle
    logical :: alive  ! is the particle alive

  end type particle_type

contains

  subroutine particle_init(this,L)

    type(particle_type) :: this
    real :: L   ! length of slab

    ! get the x location
    this%xloc = get_particle_pos(L)

    ! get the angle
    this%mu = get_particle_mu()

    ! set to alive
    this%alive = .true.

  end subroutine particle_init

end module particle
