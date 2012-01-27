module tally

  implicit none
  private
  public :: tally_reset,bank_tally,perform_statistics

  type, public :: tally_type

    real :: s1    ! s1 accumulator
    real :: s2    ! s2 accumulator 
    real :: mean  ! mean of tally
    real :: var   ! variance of tally
    real :: track ! the temp track var

  end type tally_type

contains

!===============================================================================
!
!===============================================================================

  subroutine tally_reset(this)

    type(tally_type) :: this

    this%track = 0.0

  end subroutine tally_reset

!===============================================================================
!
!===============================================================================

  subroutine bank_tally(this)

    type(tally_type) :: this

    ! bank tally
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
    this%mean = this%s1/(dx * n_hist)

  end subroutine perform_statistics

end module tally
