!==============================================================================!
! MODULE: geometry
!
!> @author Bryan Herman
!>
!> @brief Routines for defining and manipulating the problem geometry.
!==============================================================================!
module geometry

  implicit none

  private

  public :: read_geometry

  !=============================================================================
  !> @brief Defines the geometry.
  !=============================================================================
  type, public :: geometry_type
    integer :: n_slabs          ! number of sub-slabs
    double precision :: length  ! length of slab
    double precision :: dx      ! length of a sub-slab 
  end type geometry_type

contains

  !=============================================================================
  !> @brief Read user input to define the geometry.
  !>
  !> @param[in]     this      A reference to a geometry object
  !=============================================================================
  subroutine read_geometry(this)

    type(geometry_type), intent(inout) :: this

    ! ask user for length of slab
    write(*,*) 'Enter the length of the slab:'
    read(*,*) this%length

    ! ask user for number of sub-slabs
    write(*,*) 'How many tally regions?'
    read(*,*) this%n_slabs

    ! calculate dx
    this%dx = this%length/this%n_slabs

  end subroutine read_geometry

end module geometry
