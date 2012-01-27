module geometry

  implicit none
  private
  public :: read_geometry

  type, public :: geometry_type

    integer :: n_slabs ! number of sub-slabs
    real    :: length  ! length of slab
    real    :: dx      ! length of a sub-slab 

  end type geometry_type

contains

  subroutine read_geometry(this)

    type(geometry_type) :: this

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
