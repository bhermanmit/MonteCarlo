!==============================================================================!
! MODULE: pdfs
!
!> @author Bryan Herman
!>
!> @brief Routines for sampling various distributions.
!>
!> Note, all values are double precision, i.e. about 16 digits.  We could use
!> real's (about 8 digits), but what happens when we sample a billion times?
!==============================================================================!
module pdfs

  use random_number_generator, only: rand

  implicit none

contains

  !=============================================================================
  !> @brief Get a particle position in the slab
  !>
  !> @param[in]     L         Total width of domain.
  !> @return                  A random location in the domain
  !=============================================================================
  function get_particle_pos(L)

    double precision              :: get_particle_pos
    double precision, intent(in)  :: L

    get_particle_pos = rand()*L

  end function get_particle_pos

  !=============================================================================
  !> @brief Sample the particle's direction (for the source).
  !>
  !> This assumes an <em>isotropic</em> source in the LAB system.
  !>
  !> @return                  The particle's direction at birth.
  !=============================================================================
  function get_particle_mu()

    double precision :: get_particle_mu

    get_particle_mu = 2.0*rand() - 1.0

  end function get_particle_mu

  !=============================================================================
  !> @brief Sample a particle's direction following a scattering event.
  !>
  !> This assumes <em>isotropic</em> scattering in the LAB system. 
  !>
  !> @return                  The particle's new direction.
  !=============================================================================
  function get_scatter_mu()

    double precision :: get_scatter_mu

    get_scatter_mu = 2.0*rand() - 1.0

  end function get_scatter_mu

  !=============================================================================
  !> @brief Sample a particle's distance to next collision.
  !>
  !> This only applies to the slab with the given cross section.  If the 
  !> particles gets to a boundary, a new distance must be sampled with the
  !> appropriate cross section.
  !>
  !> @param[in]     xs        The cross section of the slab traversed.
  !> @return                  The distance to be traveled.
  !=============================================================================
  function get_collision_distance(xs)

    double precision              :: get_collision_distance
    double precision, intent(in)  :: xs

    get_collision_distance = -log(rand()) /xs

  end function get_collision_distance

  !=============================================================================
  !> @brief Sample the type of collision a particle undergoes.
  !>
  !> Cross sections are nothing but probabilities.  If a particle can undergo
  !> either absorption or scattering, the probability of each is simply
  !> @f[
  !>     P(\textrm{absorption}) = \Sigma_a / \Sigma_t 
  !> @f]
  !> and
  !> @f[
  !>     P(\textrm{scattering}) = \Sigma_s / \Sigma_t \, ,
  !> @f]
  !> respectively.  Things get way more complicated in multigroup transport
  !> with multiple collision types.  It's even worse for continuous energy
  !> treatments.  Just ask Paul.
  !>
  !> @param[in]     absxs     The cross section of the slab traversed.
  !> @param[in]     scattxs   The cross section of the slab traversed.
  !> @param[in]     totxs     The cross section of the slab traversed.
  !> @return                  The type of collision (1 for abs, 2 for scatter)
  !=============================================================================
  function get_collision_type(absxs, scattxs, totxs)

    integer                       :: get_collision_type
    double precision, intent(in)  :: absxs, scattxs, totxs
    double precision              :: rn

    ! get rn
    rn = rand()

    ! absorption in first bin, scattering in second
    if (rn < absxs/totxs) then
      get_collision_type = 1
    else
      get_collision_type = 2
    end if

  end function get_collision_type

end module pdfs
