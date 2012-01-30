!==============================================================================!
! MODULE: execute
!
!> @author Bryan Herman
!>
!> @brief Main routines for running the problem.
!==============================================================================!
module execute

  use global
  use pdfs
  use random_number_generator, only: initialize_rng

  implicit none
  private
  public :: run_problem,print_tallies

! Uncomment this if diagnostic output is needed.  Note, this requires -cpp.
!#define DEBUG

contains

  !=============================================================================
  !> @brief Run the simulation.
  !=============================================================================
  subroutine run_problem()

    use particle,  only: particle_init
    use tally, only: tally_init


    integer :: n  ! history loop counter
    integer :: i  ! counter for tallies

    ! initialize rng
    call initialize_rng()

    ! initialize tallies
    do i = 1,geo%n_slabs
      call tally_init(tal(i))
    end do

    ! initialize the random number generator
    call initialize_rng()

    ! begin history loop
    HISTORY: do n = 1, nhist

      ! initialize particle
      call particle_init(neutron, geo%length)

#ifdef DEBUG
     print *,'Neutron BORN'
     print *,'Neutron is at x:', neutron%xloc
     print *,'Neutron mu is:', neutron%mu
#endif

      ! reset track tallies
      call reset_tallies()

      ! determine which slab the particle is in
      neutron%slab = get_slab_id()

#ifdef DEBUG
      print *,'neutron is slab:', neutron%slab
#endif

      ! begin loop around particle's life
      LIFE: do while (neutron%alive)

        ! tranpsort neutron
        call transport()

        ! get interaction type
        if(neutron%alive) call interaction()

      end do LIFE

      ! bank tallies
      call bank_tallies()

      ! update user
      if ( mod(n, 100000) == 0 ) then
        write(*,'("Successfully transported: ",I0," particles...")') n
      end if

    end do HISTORY

  end subroutine run_problem

  !=============================================================================
  !> @brief Perform transport of a single particle
  !=============================================================================
  subroutine transport()

    double precision :: s     ! free flight distance
    double precision :: newx  ! the temp newx location
    double precision :: neig  ! nearest neighbor surface in traveling direction
    logical :: resample ! resample the distance

    ! set resample
    resample = .true.

    ! begin while loop until collide
    do while (resample)

      ! get the distance to next collision
      s = get_collision_distance(mat%totalxs)

      ! compute x component
      newx = neutron%xloc + s*neutron%mu 

#ifdef DEBUG
      print *,'s =',s
#endif

      ! compute x component
      newx = neutron%xloc + neutron%mu * s

#ifdef DEBUG
      print *,'newx =',newx
#endif

      ! get nearest neigbor
      if (neutron%mu > 0.0) then
        neig = float(neutron%slab)*geo%dx
      else
        neig = float(neutron%slab - 1)*geo%dx
      end if


#ifdef DEBUG
      print *,'Neighbor is:',neig
#endif

      ! check for surface crossing
      if ( (neutron%mu < 0.0 .and. newx < neig) .or.                           &
     &     (neutron%mu > 0.0 .and. newx > neig) ) then

        ! check for global boundary crossing
        if ( (newx < 0.0 .and. neutron%slab == 1) .or.                         &
       &     (newx > geo%length .and. neutron%slab == geo%n_slabs)) then

          ! kill particle
          neutron%alive = .false.

          ! no resample needed
          resample = .false.

        end if

        ! record tally
        tal(neutron%slab)%track = tal(neutron%slab)%track +                    &
                                  (neig - neutron%xloc) / neutron%mu

        ! move particle to surface and resample
        neutron%xloc = neig

        ! change slab number
        if (neutron%mu > 0.0) then
          neutron%slab = neutron%slab + 1
        else
          neutron%slab = neutron%slab - 1
        end if

      else ! collision occurred

        ! record distance in tally
        tal(neutron%slab)%track = tal(neutron%slab)%track + s 

        ! move neutron
        neutron%xloc = newx

        ! set resample to false
        resample = .false.

#ifdef DEBUG
        print *,'Collision occurred move neutron to:', neutron%xloc
#endif

      end if

    end do 

  end subroutine transport

  !=============================================================================
  !> @brief Determine and handle an interaction.
  !=============================================================================
  subroutine interaction()

    integer :: id

    ! record tally
    tal(neutron%slab)%coll = tal(neutron%slab)%coll + 1.0/mat%totalxs

    ! get reaction type
    id = get_collision_type(mat%absxs, mat%scattxs, mat%totalxs)

    if ( id == 1 ) then

#ifdef DEBUG
      print *,'Neutron absorbed'
#endif

      ! kill particle
      neutron%alive = .false.

    else

      ! sample new angle
      neutron%mu = get_scatter_mu()

#ifdef DEBUG
      print *,'Neutron scattered'
#endif

      ! sample new angle
      neutron%mu = get_scatter_mu()

#ifdef DEBUG
      print *,'New angle:',neutron%mu
#endif

    end if

  end subroutine interaction

  !=============================================================================
  !> @brief Determine the slab in which a particle resides.
  !> @return                  Slab ID
  !=============================================================================
  function get_slab_id()

    integer :: get_slab_id

    get_slab_id = ceiling(neutron%xloc / geo%dx)

  end function get_slab_id

  !=============================================================================
  !> @brief Reset the tally array.
  !=============================================================================
  subroutine reset_tallies()

    use tally, only: tally_reset

    integer :: i

    ! loop around and reset
    do i = 1,geo%n_slabs

      call tally_reset(tal(i))

    end do 

  end subroutine reset_tallies

  !=============================================================================
  !> @brief Add the current tallies to the mean and variance.
  !=============================================================================
  subroutine bank_tallies()

    use tally, only: bank_tally

    integer :: i  ! counter

    ! begin loop around tallies
    do i = 1, geo%n_slabs

      ! bank tally
      call bank_tally(tal(i))

    end do

  end subroutine bank_tallies

  !=============================================================================
  !> @brief Print the tally of each slab.
  !=============================================================================
  subroutine print_tallies()

    use tally, only: perform_statistics

    integer :: i

    ! set results
    write(*,'(///,"Results",/,"=======",/)')

    do i = 1, geo%n_slabs

      ! compute stat
      call perform_statistics(tal(i), nhist, geo%dx)

      ! print mean
      write(*,'("Slab ",I0,T10," Flux: ",F0.4,T30,F0.4)')                      &
     &     i,tal(i)%smean,tal(i)%cmean

    end do

  end subroutine print_tallies

end module execute
