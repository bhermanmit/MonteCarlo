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

    double precision :: c1    = 0.0_8 ! collision accumulator
    double precision :: c2    = 0.0_8 ! square of collision accumulator
    double precision :: s1    = 0.0_8 ! path accumulator
    double precision :: s2    = 0.0_8 ! square of path accumulator
    double precision :: smean = 0.0_8 ! mean for tracklength est
    double precision :: cmean = 0.0_8 ! mean for collision est
    double precision :: svar  = 0.0_8 ! variance of tracklength est
    double precision :: cvar  = 0.0_8 ! variance of collision est
    double precision :: track = 0.0_8 ! the temp track var
    double precision :: coll  = 0.0_8 ! the temp coll var

  end type tally_type

contains

  !=============================================================================
  !> @brief Reset a tally.
  !>
  !> @param[in]     this      A reference to a tally object
  !=============================================================================
  subroutine tally_reset(this)

    type(tally_type), intent(inout) :: this

    this%track = 0.0_8
    this%coll  = 0.0_8

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
    this%c1 = this%c1 + this%coll
    this%c2 = this%c2 + this%coll**2

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
    double precision, intent(in)    :: dx
    integer, intent(in)             :: n_hist

    ! compute mean
    this%smean = this%s1/(dx * dble(n_hist))
    this%cmean = this%c1/(dx * dble(n_hist))

    ! compute variance
    this%svar = 1.0_8/dble(n_hist - 1)*(this%s2/dble(n_hist) -                 &
                (this%s1/dble(n_hist))**2)
    this%cvar = 1.0_8/dble(n_hist - 1)*(this%c2/dble(n_hist) -                 &
                (this%c1/dble(n_hist))**2)

  end subroutine perform_statistics

end module tally
