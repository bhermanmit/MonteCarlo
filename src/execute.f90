module execute

  use global
  use pdfs

  implicit none
  private
  public :: run_problem,print_tallies

contains

!===============================================================================
!
!===============================================================================

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

    ! begin history loop
    HISTORY: do n = 1,nhist

      ! initialize particle
      call particle_init(neutron,geo%length)

      ! reset track tallies
      call reset_tallies()

      ! determine which slab the particle is in
      neutron%slab = get_slab_id()

      ! begin loop around particles life
      LIFE: do while (neutron%alive)

        ! tranpsort neutron
        call transport()

        ! get interaction type
        if(neutron%alive) call interaction()

      end do LIFE

      ! bank tallies
      call bank_tallies()

      ! update user
      if ( mod(n,1000) == 0 ) then
        write(*,'("Successfully transported: ",I0," particles...")') n
      end if

    end do HISTORY

  end subroutine run_problem

!===============================================================================
!
!===============================================================================

  subroutine transport()

    real :: s     ! free flight distance
    real :: newx  ! the temp newx location
    real :: neig  ! nearest neighbor surface in traveling direction
    logical :: resample ! resample the distance

    ! set resample
    resample = .true.

    ! begin while loop until collide
    do while ( resample )

      ! get the distance to next collision
      s = get_collision_distance(mat%totalxs)

      ! compute x component
      newx = neutron%xloc + s*neutron%mu 

      ! get nearest neigbor
      if (neutron%mu > 0.0) then
        neig = float(neutron%slab)*geo%dx
      else
        neig = float(neutron%slab - 1)*geo%dx
      end if

      ! check for surface crossing
      if ( (neutron%mu < 0.0 .and. newx < neig) .or.                           &
     &     (neutron%mu > 0.0 .and. newx > neig) ) then

        ! check for global boundary crossing
        if ( (newx < 0.0 .and. neutron%slab == 1) .or.                         &
       &      newx > geo%length .and. neutron%slab == geo%n_slabs) then

          ! kill particle
          neutron%alive = .false.

          ! no resample needed
          resample = .false.

        end if

        ! record tally
        tal(neutron%slab)%track = tal(neutron%slab)%track +                    &
       &                          (neig - neutron%xloc) / neutron%mu

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

      end if

    end do 

  end subroutine transport

!===============================================================================
!
!===============================================================================

  subroutine interaction()

    integer :: id

    ! record tally
    tal(neutron%slab)%coll = tal(neutron%slab)%coll + 1.0/mat%totalxs

    ! get reaction type
    id = get_collision_type(mat%absxs,mat%scattxs,mat%totalxs)

    if ( id == 1 ) then

      ! kill particle
      neutron%alive = .false.

    else

      ! sample new angle
      neutron%mu = get_scatter_mu()

    end if

  end subroutine interaction

!===============================================================================
!
!===============================================================================

  function get_slab_id()

    integer :: get_slab_id

    get_slab_id = ceiling(neutron%xloc / geo%dx)

  end function get_slab_id

!===============================================================================
!
!===============================================================================

  subroutine reset_tallies()

    use tally, only: tally_reset

    integer :: i

    ! loop around and reset
    do i = 1,geo%n_slabs

      call tally_reset(tal(i))

    end do 

  end subroutine reset_tallies

!===============================================================================
!
!===============================================================================

  subroutine bank_tallies()

    use tally, only: bank_tally

    integer :: i  ! counter

    ! begin loop around tallies
    do i = 1,geo%n_slabs

      ! bank tally
      call bank_tally(tal(i))

    end do

  end subroutine bank_tallies

!===============================================================================
!
!===============================================================================

  subroutine print_tallies()

    use tally, only: perform_statistics

    integer :: i

    ! set results
    write(*,'(///,"Results",/,"=======",/)')

    do i = 1,geo%n_slabs

      ! compute stat
      call perform_statistics(tal(i),nhist,geo%dx)

      ! print mean
      write(*,'("Slab ",I0,T10," Flux: ",F0.4,T30,F0.4)') i,tal(i)%smean,tal(i)%cmean

    end do

  end subroutine print_tallies

end module execute
