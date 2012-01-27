module pdfs

  implicit none

contains

!===============================================================================
!
!===============================================================================

  subroutine initialize_rng()

    real :: rn
    rn = rand(234231)

  end subroutine initialize_rng

!===============================================================================
!
!===============================================================================

  function get_particle_pos(L)

    real :: get_particle_pos
    real :: L

    get_particle_pos = rand(0)*L

  end function get_particle_pos

!===============================================================================
!
!===============================================================================

  function get_particle_mu()

    real :: get_particle_mu

    get_particle_mu = 2.0*rand(0) - 1.0

  end function get_particle_mu

!===============================================================================
!
!===============================================================================

  function get_scatter_mu()

    real :: get_scatter_mu

    get_scatter_mu = 2.0*rand(0) - 1.0

  end function get_scatter_mu

!===============================================================================
!
!===============================================================================

  function get_collision_distance(xs)

    real :: get_collision_distance
    real :: xs

    get_collision_distance = -log(rand(0))/xs

  end function get_collision_distance

!===============================================================================
!
!===============================================================================

  function get_collision_type(absxs,scattxs,totxs)

    integer :: get_collision_type
    real :: absxs,scattxs,totxs
    real :: rn
    
    ! get rn
    rn = rand(0)

    ! absorption in first bin, scattering in second
    if (rn < absxs/totxs) then
      get_collision_type = 1
    else
      get_collision_type = 2
    end if

  end function get_collision_type

end module pdfs
