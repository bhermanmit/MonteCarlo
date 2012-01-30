!==============================================================================!
! MODULE: random_number_generator
!
!> @author Jeremy Roberts
!>
!> @brief Provides tools for sampling pseudo-random numbers.
!>
!> Note, this uses the Fortran intrinsic functions "random_seed" and
!> "random_number".  The implementation of these is compiler-specific.  In
!> particular, the gfortran implementation of random_number uses 
!> Marsaglia's "Keep It Simple Stupid" generator.  While the implementation 
!> <em>is</em> thread safe, it is done in such a way that its parallel 
!> performance is really bad when heavily used (as we do in Monte Carlo).  
!> For parallel applications, a different generator is essential.
!==============================================================================!
module random_number_generator

  implicit none

  private

  public :: initialize_rng,rand

contains

  !============================================================================! 
  !> @brief Initialize the random number generator. 
  !============================================================================!
  subroutine initialize_rng()
  
    integer                            :: i, seed_length
    integer, dimension(:), allocatable :: seed

    ! Get the minimum size of the array we can use for seeding the generator.
    call random_seed(size = seed_length)
    
    ! Then allocate the seed array.
    allocate(seed(seed_length))

    ! Hard code a seed.  One could also use the system clock, via
    !   call system_clock(count = clock)
    ! for an integer "clock", and then put that in place of 11235.
    seed = 11235 + 81321 * (/ (i - 1, i = 1, seed_length) /)

    ! Actually set the seed and deallocate it.
    call random_seed(put = seed)
    deallocate(seed)

  end subroutine initialize_rng

  !============================================================================! 
  !> @brief Initialize the random number generator for a new history.
  !>
  !> For this implementation of the PRNG library, this function is not 
  !> needed, but the interface must be supported anyway.
  !>
  !> @param[in] 	n 		Current neutron history
  !============================================================================!
  subroutine initialize_rng_history(n)
  
    integer, intent(in) :: n
    ! do nothing.

  end subroutine initialize_rng_history

  !============================================================================!
  !> @brief Obtain a single random number.
  !>
  !> @return random number
  !============================================================================!
  double precision function rand()

    call random_number(rand)

  end function rand


end module random_number_generator
