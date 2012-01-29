module tally

  implicit none
  private
  public :: tally_init,tally_reset,bank_tally,perform_statistics

  type, public :: tally_type

    real :: c1    ! collision accumulator
    real :: c2    ! square of collision accumulator
    real :: s1    ! path accumulator
    real :: s2    ! square of path accumulator 
    real :: smean ! mean for tracklength est
    real :: cmean ! mean for collision est
    real :: svar  ! variance of tracklength est
    real :: cvar  ! variance of collision est
    real :: track ! the temp track var
    real :: coll  ! the temp coll var

  end type tally_type

contains

!===============================================================================
!
!===============================================================================

  subroutine tally_init(this)

    type(tally_type) :: this

    this%c1 = 0.0
    this%c2 = 0.0
    this%s1 = 0.0
    this%s2 = 0.0
    this%smean = 0.0
    this%cmean = 0.0
    this%svar = 0.0
    this%cvar = 0.0
    this%track = 0.0
    this%coll = 0.0

  end subroutine tally_init

!===============================================================================
!
!===============================================================================

  subroutine tally_reset(this)

    type(tally_type) :: this

    this%track = 0.0
    this%coll = 0.0

  end subroutine tally_reset

!===============================================================================
!
!===============================================================================

  subroutine bank_tally(this)

    type(tally_type) :: this

    ! bank tally
    this%c1 = this%c1 + this%coll
    this%c2 = this%c2 + this%coll**2
    this%s1 = this%s1 + this%track
    this%s2 = this%s2 + this%track**2


  end subroutine bank_tally

!===============================================================================
!
!===============================================================================

  subroutine perform_statistics(this,n_hist,dx)

    type(tally_type) :: this
    real    :: dx
    integer :: n_hist

    ! compute mean
    this%smean = this%s1/(dx * n_hist)
    this%cmean = this%c1/(dx * n_hist)

  end subroutine perform_statistics

end module tally
