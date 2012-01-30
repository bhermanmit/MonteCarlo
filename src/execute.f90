!==============================================================================!
! MODULE: execute
!
!> @author Bryan Herman
!>
!> @brief Main routines for running the problem.
!==============================================================================!
module execute

  use global
  use particle,  only: particle_init
  use pdfs
  use random_number_generator, only: initialize_rng, initialize_rng_history

  implicit none

! Uncomment this if diagnostic output is needed.  Note, this requires -cpp.
!#define DEBUG

contains

  !=============================================================================
  !> @brief Run the simulation.
  !=============================================================================
  subroutine run_problem()

    integer(I8) :: n  ! history loop counter
    integer     :: i,j  ! loop index for tally update
    double precision :: omp_get_wtime, t
    integer    :: omp_get_thread_num

!$ t = omp_get_wtime()

    ! initialize the random number generator (master only!!)
    call initialize_rng()

    do i = 1, geo%n_slabs

      ! zero global tallies
      tal(i)%s1 = 0.0_8
      tal(i)%s2 = 0.0_8

    end do

!$omp parallel private(n, i) shared(tal, nhist, mat, geo)

    ! local tally allocation 
    if (.not. allocated(local_tal)) allocate(local_tal(geo%n_slabs))
    !print *, "i am ", omp_get_thread_num(), " num slabs = ",geo%n_slabs
    

    do i = 1, geo%n_slabs

      ! zero local tallies
      local_tal(i)%s1 = 0.0_8
      local_tal(i)%s2 = 0.0_8
      local_tal(i)%c1 = 0.0_8
      local_tal(i)%c2 = 0.0_8

    end do

    ! begin history loop

!$omp do 
    HISTORY: do n = 1, nhist

      ! initialize the random number generator for this history
      call initialize_rng_history(n)

      ! initialize particle
      call particle_init(neutron, geo%length)

      ! reset track tallies
      call reset_tallies()

      ! determine which slab the particle is in
      neutron%slab = get_slab_id()

      ! begin loop around particle's life
      LIFE: do while (neutron%alive)

        ! tranpsort neutron
        call transport(n)

        ! get interaction type
        if(neutron%alive) call interaction()

      end do LIFE

      ! bank tallies
      call bank_tallies()

      ! update user
!      if ( mod(n, 100000_I8) == 0 ) then
!        write(*,'("Successfully transported: ",I0," particles...")') n
!      end if

    end do HISTORY
!$omp end do

    ! Now combine the results in the master tally.  The "(tally)" denotation 
    ! is optional and must be unique in a code.  Here, it just helps clarify.

!$omp critical (tally)
    do i = 1, geo%n_slabs
      !print *, " i am ", omp_get_thread_num(), " j = ", i, " tal(j)s1 = ", local_tal(i)%s1
      tal(i)%s1 = tal(i)%s1 + local_tal(i)%s1
      tal(i)%s2 = tal(i)%s2 + local_tal(i)%s2
      tal(i)%c1 = tal(i)%c1 + local_tal(i)%c1
      tal(i)%c2 = tal(i)%c2 + local_tal(i)%c2
    end do
!$omp end critical (tally)

!$omp end parallel

  end subroutine run_problem

  !=============================================================================
  !> @brief Perform transport of a single particle
  !=============================================================================
  subroutine transport(n)

    integer(I8) :: n          ! history index
    double precision :: s     ! free flight distance
    double precision :: newx  ! the temp newx location
    double precision :: neig  ! nearest neighbor surface in traveling direction
    logical :: resample       ! resample the distance

    ! set resample
    resample = .true.

    ! begin while loop until collide
    do while (resample)

      ! get the distance to next collision
      s = get_collision_distance(mat%totalxs)

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
           (neutron%mu > 0.0 .and. newx > neig) ) then

        ! check for global boundary crossing
        if ( (newx < 0.0 .and. neutron%slab == 1) .or.                         &
              newx > geo%length .and. neutron%slab == geo%n_slabs) then

          ! kill particle
          neutron%alive = .false.

          ! no resample needed
          resample = .false.

        end if

        ! record tally
        local_tal(neutron%slab)%track = local_tal(neutron%slab)%track +                    &
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
        local_tal(neutron%slab)%track = local_tal(neutron%slab)%track + s 

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

    ! tally collision
    local_tal(neutron%slab)%coll = local_tal(neutron%slab)%coll +              &
   &                               1.0_8/mat%totalxs

    ! get reaction type
    id = get_collision_type(mat%absxs, mat%scattxs, mat%totalxs)

    if ( id == 1 ) then

#ifdef DEBUG
      print *,'Neutron absorbed'
#endif

      ! kill particle
      neutron%alive = .false.

    else

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
    do i = 1, geo%n_slabs

      call tally_reset(local_tal(i))

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
      call bank_tally(local_tal(i))

    end do

  end subroutine bank_tallies

  !=============================================================================
  !> @brief Print the tally of each slab.
  !=============================================================================
  subroutine print_tallies()

    use global, only: timer_run
    use tally,  only: perform_statistics

    integer :: i

    ! set results
    write(*,'(///,"Results",/,"=======",/)')

    ! print time
    write(*,'("Execution Time:",2X,F0.4," seconds",/)') timer_run%elapsed

    ! print tally header
    write(*,'("Slab #",T10,"Flux - Tracklength",T30,"Flux - Collision")')
    write(*,'("------",T10,"------------------",T30,"----------------")')

    do i = 1, geo%n_slabs

      ! compute stat
      call perform_statistics(tal(i), nhist, geo%dx)

      ! print mean
      write(*,'(I0,T10,F0.4,T30,F0.4)') i,tal(i)%smean,tal(i)%cmean

    end do

  end subroutine print_tallies

end module execute
