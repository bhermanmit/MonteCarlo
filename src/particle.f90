!==============================================================================!
! MODULE: particle
!
!> @author Bryan Herman
!>
!> @brief Routines for defining and manipulating a particle.
!==============================================================================!
module particle

  use pdfs

  implicit none
  private
  public :: particle_init

  !=============================================================================
  !> @brief Defines a particle.
  !=============================================================================
  type, public :: particle_type
    integer :: slab           ! the slab id number
    double precision :: xloc  ! location of particle in slab
    double precision :: mu    ! cosine of angle
    logical :: alive          ! is the particle alive
  end type particle_type

contains

  !=============================================================================
  !> @brief Initialize a particle.
  !>
  !> @param[in]     this      A reference to a particle object
  !> @param[in]     L         Total width of domain
  !=============================================================================
  subroutine particle_init(this, L)

    type(particle_type), intent(inout) :: this
    double precision, intent(in) :: L   ! length of slab

    ! get the x location
    this%xloc = get_particle_pos(L)

    ! get the angle
    this%mu = get_particle_mu()

    ! set to alive
    this%alive = .true.

  end subroutine particle_init

end module particle
