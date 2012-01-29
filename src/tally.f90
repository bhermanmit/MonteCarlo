!==============================================================================!
! MODULE: tally
!
!> @author Bryan Herman
!>
!> @brief Provides routines for tallying particle interactions.
!==============================================================================!
module tally

  implicit none

  private

  public :: tally_reset, bank_tally, perform_statistics

  !=============================================================================
  !> A tally for a single slab region.
  !>
  !> @note If the members are not initialized, valgrind complains.  It's
  !>       arguably good practice to initialize variables.
  !=============================================================================
  type, public :: tally_type

    double precision :: s1    = 0.0  ! s1 accumulator
    double precision :: s2    = 0.0  ! s2 accumulator 
    double precision :: mean  = 0.0  ! mean of tally
    double precision :: var   = 0.0  ! variance of tally
    double precision :: track = 0.0  ! the temp track var

  end type tally_type

contains

  !=============================================================================
  !> @brief Reset a tally.
  !>
  !> @param[in]     this      A reference to a tally object
  !=============================================================================
  subroutine tally_reset(this)

    type(tally_type), intent(inout) :: this

    this%track = 0.0

  end subroutine tally_reset

  !=============================================================================
  !> @brief Reset a tally.
  !>
  !> @param[in]     this      A reference to a tally object
  !=============================================================================
  subroutine bank_tally(this)

    type(tally_type) :: this

    ! bank tally
    this%s1 = this%s1 + this%track
    this%s2 = this%s2 + this%track**2

  end subroutine bank_tally

  !=============================================================================
  !> @brief Compute the mean.
  !>
  !> @param[in]     this      A reference to a tally object
  !> @param[out]    this      A reference to a tally object
  !> @param[in]     n_hist    Number of histories
  !> @param[in]     dx        Slab width
  !=============================================================================
  subroutine perform_statistics(this, n_hist, dx)

    type(tally_type), intent(inout) :: this
    double precision, intent(in)                :: dx
    integer, intent(in)             :: n_hist

    ! compute mean
    this%mean = this%s1/(dx * n_hist)

  end subroutine perform_statistics

end module tally
