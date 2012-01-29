module material

  implicit none
  private
  public :: read_material, default_material

  type, public :: material_type

    double precision :: totalxs  ! macroscopic total cross section
    double precision :: absxs    ! macroscopic absorption cross section
    double precision :: scattxs  ! macroscopic scattering cross section

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

  subroutine default_material(this)

    type(material_type) :: this

    ! set total xs
    this%totalxs = 1.0

    ! set absorption xs
    this%absxs = 0.5

    ! set scattering xs
    this%scattxs = 0.5

  end subroutine default_material


end module material
