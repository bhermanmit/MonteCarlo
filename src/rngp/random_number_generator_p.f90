!==============================================================================!
! MODULE: random_number_generator
!
!> @author Jeremy Roberts
!>
!> @brief Provides a simple implementation of the MCNP standard PRNG
!>
!> MCNP and most other Monte Carlo codes employ LCGs, defined by the relation
!> @f[
!>      S(k+1) = \textrm{modulo}[ (g \times S(k) + c), 2^M ] \, ,      
!> @f]
!> with the random number being defined
!> @f[
!>            r(k+1) = S(k+1)/2^M \, .
!> @f]
!> For the standard MCNP LCG,
!> @f[
!>            g      = 5^19 \,
!> @f]
!> @f[
!>            c      = 0 \,
!> @f]
!> @f[
!>            S(0)   = 5^19 \,
!> @f]
!> and 
!> @f[
!>       M      = 48 .  
!> @f]
!> (Recall that @f$ \textrm{modulo}[a,b] = a\, \textrm{mod} \, b @f$).
!>
!> For parallel jobs, we have several options for RNs. We could let each process 
!> have its own (random?) starting seed, and then perform its batch of particle
!> histories using its own sequence.  This is not the best approach, since then
!> the same problem run using different numbers of processes would yield 
!> different results, and we like reproducing results for debugging.
!>
!> The solution: let each history have its own sequence of 
!> a certain length called the "stride".  This works since any subset of the 
!> infinite sequence defined by the first two equations above
!> should be random.  To use this 
!> approach, we need a fast way to get any @f$ S(k) @f$ given @f$ S(0) @f$.  
!> This is done by
!> @f[
!>            S(k) = \textrm{modulo}[ 
!>                   (g^k \times S(0) + c \times (g^k-1)/(g-1), 2^M ] \, ,
!> @f]
!> which must use exact integer arithmetic.  For the case that @f$ c=0 @f$, we 
!> need only compute 
!. @f$ S(k) = G(k) \times S(0) = \textrm{modulo}[g^k \times S(0), 2^M] @f$ 
!> exactly.  Brown gives the
!> following algorithm:
!>
!> @verbatim
!>     G = 1; h = g; i = modulo[k+2^m,2^m];
!>     do while (i .gt. 0 )
!>         if ( i is odd ) G = modulo[ G*h, 2^m ]
!>         h = modulo[ h^2 , 2^m ]
!>         i = floor(i/2)
!>     end 
!> @endverbatim
!>
!> A similar algorithm exists for the other term.  See the "Fundamentals..." 
!> notes for more.  
!>
!> It's worth pointing out use of the Fortran intrinsics "iand" and "btest"
!> below.  These are "bitwise and" and "bitwise test", respectively.  To give a 
!> flavor...
!> @verbatim
!>   iand(1,2) compares [1000...] to [0100...] = 0
!>   iand(1,3) compares [1000...] to [1100...] = 0*2^0 + 1*2^1 ... = 2
!> @endverbatim
!> and so on.  A good exercise would be to show why
!> @verbatim
!>   S(k) = modulo[ (g^k*S(0)+c*(g^k-1)/(g-1), 2^M ]
!>        = iand( iand( g*S(0), 2^M-1 ) + c,  2^M-1 ) 
!> @endverbatim
!> Hint: Show 
!> @verbatim 
!>   modulo[a+b,c] = modulo[modulo[a,c]+b,c]
!> @endverbatim
!> Then show 
!> @verbatim
!>   iand(a,c-1) = modulo(a,c) 
!> @endverbatim
!> when c is a power of 2.
!>
!> References: 
!> - "THE MCNP5 RANDOM NUMBER GENERATOR", LA-UR-02-3782
!>   http://mcnp-green.lanl.gov/publication/pdf/LA-UR-02-3782.pdf
!> - "FUNDAMENTALS OF MONTE CARLO TRANSPORT", LA-UR-05-4983
!>   http://...gov/publication/pdf/LA-UR-05-4983_Monte_Carlo_Lectures.pdf
!-------------------------------------------------------------------------------
module random_number_generator


	private ! default all module variables to be local within module routines

	! 64-bit real and integers (could make available to transport routines)
	integer, parameter :: R8 = selected_real_kind(15,307)
	integer, parameter :: I8 = selected_int_kind(18)

	! Public functions and subroutines for this module
	public ::  rand                   ! Gives a uniform random number
	public ::  initialize_rng         ! Initialize generator
	public ::  initialize_particle    ! New subsequence for a history

	! Default LCG values
  integer(I8), save :: g            =  19073486328125_I8 ! g, 5^19
  integer(I8), save :: c            =               0_I8 ! c = 0
  integer(I8), save :: Stride       =          152917_I8 ! default stride
  integer(I8), save :: S0           =  19073486328125_I8 ! S0, 5^19
  integer(I8), save :: TwoToMMinus1 = 281474976710655_I8 ! 2^48 - 1
  real(R8),    save :: OneOver2toM  = 1._R8 / 281474976710656._R8 ! 1/2^48

  ! Current seed, needs to be private to a thread
  integer(I8), save :: S   = 19073486328125_I8 
  !$omp threadprivate ( S )
    
contains


  !============================================================================!
  !> @brief Obtain a single 64 bit real pseudo random number.
  !>
  !> @return random number
  !============================================================================!
  real(R8) function rand()

      implicit none

      S    = iand( iand( g*S, TwoToMMinus1) + c,  TwoToMMinus1 )
      rand = S * OneOver2toM 

  end function rand

  !============================================================================!
  !> @brief Initiates the random number generator for a new problem.
  !>
  !> See "Fundamentals..." slides 2-18 and 2-19 to understand the algorithm.
  !>
  !> @param   Seed    Initial seed (i.e. S(0))
  !> @param   Skip    The value "k" for S(k) that we'd like
  !> @return          The first seed for a new history's sequence, S(k)
  !============================================================================!
  integer(I8) function RN_skip_ahead( Seed, Skip )

      implicit none
      integer(I8), intent(in)  :: Seed, Skip

      ! Local Variables:
      integer(I8)              :: nskip, BigG, tmpg, BigC, tmpc, gp, rn, SeedOld

      SeedOld = Seed
      nskip = skip
      nskip = iand( nskip, TwoToMMinus1 ) ! modulo(nskip,2^M)
      BigG  = 1           
      tmpg  = g
      BigC  = 0
      tmpc  = c
      do while( nskip .gt. 0_I8 )
        if( btest(nskip,0) )  then      ! Checks 0th bit of nskip: 0=even, 1=odd
          BigG = iand( BigG*tmpg, TwoToMMinus1 )  
          BigC = iand( BigC*tmpg, TwoToMMinus1 )
          BigC = iand( BigC+tmpc, TwoToMMinus1 )
        end if
        gp    = iand( tmpg+1,  TwoToMMinus1 )
        tmpg  = iand( tmpg*tmpg,  TwoToMMinus1 )
        tmpc  = iand( gp*tmpc, TwoToMMinus1 )
        nskip = ishft( nskip, -1 )
      end do
      rn = iand( BigG*SeedOld, TwoToMMinus1 )
      rn = iand( rn + BigC, TwoToMMinus1 )
      RN_skip_ahead = rn

  end function RN_skip_ahead

  !============================================================================!
  !> @brief Initiates the random number generator for a new problem.
  !>
  !> Warning: ONLY CALL THIS FROM THE MASTER THREAD!
  !============================================================================!
  subroutine initialize_rng()

      implicit none

      integer(I8) :: UserS
      integer(I8) :: UserStride

      UserS = 1234567_I8
      UserStride = 0_I8

      ! Update defaults if user gives nonzero seed and/or stride
      if( UserS .gt. 0_I8 ) then
        S0  = UserS
      endif
      if( UserStride .gt. 0_I8 ) then
        Stride = Stride
      endif

      ! set the initial particle seed
      S  = S0

  end subroutine initialize_rng

  !============================================================================!
  !> @brief Sample from a new location within the LCG sequence.
  !> 
  !> Begin sampling from the LCG sequence at an index defined by the 
  !> given history number, stride, and initial seed.
  !>
  !> @param[in]   History       Current particle history
  !============================================================================!
  subroutine initialize_rng_history(History)

      implicit none
      integer(I8), intent(in) :: History
      S  = RN_skip_ahead( S0, History*Stride )
      return

  end subroutine initialize_rng_history

end module random_number_generator
