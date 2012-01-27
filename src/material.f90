module material

  implicit none
  private
  public :: read_material

  type, public :: material_type

    real :: totalxs  ! macroscopic total cross section
    real :: absxs    ! macroscopic absorption cross section
    real :: scattxs  ! macroscopic scattering cross section

  end type material_type

contains

  subroutine read_material(this)

    type(material_type) :: this

    ! read in total xs
    write(*,*) 'Enter macroscopic total cross section:'
    read(*,*) this%totalxs

    ! read in absorption xs
    write(*,*) 'Enter macroscopic absorption cross section:'
    read(*,*) this%absxs

    ! read in scattering xs
    write(*,*) 'Enter macroscopic scattering cross section:'
    read(*,*) this%scattxs

  end subroutine read_material

end module material
