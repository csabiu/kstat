!------------------------------------------------------------
! Test program for jackknife module
!
! This program:
! 1. Generates random test data (galaxy and random points)
! 2. Runs the Fortran jackknife_assign subroutine
! 3. Outputs results that can be compared with jk.sh shell script
!
! Usage:
!   ./test_jackknife <Ngal> <Nran> <Ndiv> <Ndim> <seed>
!
! The program writes:
!   - test_gal.dat: Galaxy coordinates
!   - test_ran.dat: Random coordinates
!   - test_fortran_gal.jk: Fortran jackknife results for galaxies
!   - test_fortran_ran.jk: Fortran jackknife results for randoms
!------------------------------------------------------------
program test_jackknife
  use jackknife_module
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  integer :: Ngal, Nran, Ndiv, Ndim, seed
  real(dp), allocatable :: coords(:,:)
  integer, allocatable :: region_id(:)
  integer :: i, Ntotal
  character(len=256) :: arg
  real(dp) :: harvest

  ! Parse command line arguments
  if (command_argument_count() < 5) then
    print *, 'Usage: test_jackknife <Ngal> <Nran> <Ndiv> <Ndim> <seed>'
    print *, '  Ngal  - Number of galaxy points'
    print *, '  Nran  - Number of random points'
    print *, '  Ndiv  - Number of divisions per dimension'
    print *, '  Ndim  - Dimensionality (2 or 3)'
    print *, '  seed  - Random seed for reproducibility'
    stop 1
  endif

  call get_command_argument(1, arg)
  read(arg, *) Ngal
  call get_command_argument(2, arg)
  read(arg, *) Nran
  call get_command_argument(3, arg)
  read(arg, *) Ndiv
  call get_command_argument(4, arg)
  read(arg, *) Ndim
  call get_command_argument(5, arg)
  read(arg, *) seed

  ! Validate inputs
  if (Ndim /= 2 .and. Ndim /= 3) then
    print *, 'Error: Ndim must be 2 or 3'
    stop 1
  endif

  print *, '================================================'
  print *, 'Jackknife Test Program'
  print *, '================================================'
  print *, 'Parameters:'
  print *, '  Ngal  =', Ngal
  print *, '  Nran  =', Nran
  print *, '  Ndiv  =', Ndiv
  print *, '  Ndim  =', Ndim
  print *, '  seed  =', seed
  print *, '  Expected regions:', Ndiv**Ndim
  print *, '================================================'

  Ntotal = Ngal + Nran

  ! Allocate arrays
  allocate(coords(Ndim, Ntotal))
  allocate(region_id(Ntotal))

  ! Initialize random number generator with seed
  call init_random_seed(seed)

  ! Generate random coordinates
  print *, 'Generating test data...'
  do i = 1, Ntotal
    call random_number(harvest)
    coords(1, i) = harvest * 1000.0_dp  ! X in [0, 1000]
    call random_number(harvest)
    coords(2, i) = harvest * 1000.0_dp  ! Y in [0, 1000]
    if (Ndim == 3) then
      call random_number(harvest)
      coords(3, i) = harvest * 1000.0_dp  ! Z in [0, 1000]
    endif
  enddo

  ! Write galaxy coordinates to file (for shell script comparison)
  print *, 'Writing galaxy data to test_gal.dat...'
  open(unit=10, file='test_gal.dat', status='replace')
  if (Ndim == 2) then
    do i = 1, Ngal
      write(10, '(2(ES20.12,1X))') coords(1, i), coords(2, i)
    enddo
  else
    do i = 1, Ngal
      write(10, '(3(ES20.12,1X),F5.1)') coords(1, i), coords(2, i), coords(3, i), 1.0
    enddo
  endif
  close(10)

  ! Write random coordinates to file
  print *, 'Writing random data to test_ran.dat...'
  open(unit=10, file='test_ran.dat', status='replace')
  if (Ndim == 2) then
    do i = Ngal + 1, Ntotal
      write(10, '(2(ES20.12,1X))') coords(1, i), coords(2, i)
    enddo
  else
    do i = Ngal + 1, Ntotal
      write(10, '(3(ES20.12,1X),F5.1)') coords(1, i), coords(2, i), coords(3, i), 1.0
    enddo
  endif
  close(10)

  ! Run Fortran jackknife assignment
  print *, 'Running Fortran jackknife assignment...'
  call jackknife_assign(coords, Ntotal, Ndiv, Ndim, region_id)

  ! Write Fortran results for galaxies
  print *, 'Writing Fortran results for galaxies to test_fortran_gal.jk...'
  open(unit=10, file='test_fortran_gal.jk', status='replace')
  if (Ndim == 2) then
    do i = 1, Ngal
      write(10, '(2(ES20.12,1X),I6)') coords(1, i), coords(2, i), region_id(i)
    enddo
  else
    do i = 1, Ngal
      write(10, '(3(ES20.12,1X),F5.1,1X,I6)') coords(1, i), coords(2, i), coords(3, i), 1.0, region_id(i)
    enddo
  endif
  close(10)

  ! Write Fortran results for randoms
  print *, 'Writing Fortran results for randoms to test_fortran_ran.jk...'
  open(unit=10, file='test_fortran_ran.jk', status='replace')
  if (Ndim == 2) then
    do i = Ngal + 1, Ntotal
      write(10, '(2(ES20.12,1X),I6)') coords(1, i), coords(2, i), region_id(i)
    enddo
  else
    do i = Ngal + 1, Ntotal
      write(10, '(3(ES20.12,1X),F5.1,1X,I6)') coords(1, i), coords(2, i), coords(3, i), 1.0, region_id(i)
    enddo
  endif
  close(10)

  ! Print summary statistics
  print *, ''
  print *, 'Summary of region assignments:'
  call print_region_stats(region_id, Ngal, Nran, Ndiv, Ndim)

  print *, ''
  print *, 'Test data files created:'
  print *, '  test_gal.dat        - Galaxy coordinates'
  print *, '  test_ran.dat        - Random coordinates'
  print *, '  test_fortran_gal.jk - Fortran JK results (galaxies)'
  print *, '  test_fortran_ran.jk - Fortran JK results (randoms)'
  print *, ''
  print *, 'To compare with shell script, run:'
  print *, '  ../bin/jk.sh test_gal.dat test_ran.dat', Ndiv, Ndim
  print *, '  Then compare test_gal.dat.jk with test_fortran_gal.jk'
  print *, ''

  ! Cleanup
  deallocate(coords)
  deallocate(region_id)

  print *, 'Test program completed successfully.'

contains

  subroutine init_random_seed(seed_val)
    implicit none
    integer, intent(in) :: seed_val
    integer :: n
    integer, allocatable :: seed_array(:)

    call random_seed(size=n)
    allocate(seed_array(n))
    seed_array = seed_val
    call random_seed(put=seed_array)
    deallocate(seed_array)
  end subroutine init_random_seed

  subroutine print_region_stats(region_id, Ngal, Nran, Ndiv, Ndim)
    implicit none
    integer, intent(in) :: region_id(:), Ngal, Nran, Ndiv, Ndim

    integer :: Nregions, i, count_gal, count_ran
    integer, allocatable :: gal_counts(:), ran_counts(:)

    Nregions = Ndiv**Ndim
    allocate(gal_counts(Nregions))
    allocate(ran_counts(Nregions))

    gal_counts = 0
    ran_counts = 0

    ! Count galaxies in each region
    do i = 1, Ngal
      if (region_id(i) >= 1 .and. region_id(i) <= Nregions) then
        gal_counts(region_id(i)) = gal_counts(region_id(i)) + 1
      endif
    enddo

    ! Count randoms in each region
    do i = Ngal + 1, Ngal + Nran
      if (region_id(i) >= 1 .and. region_id(i) <= Nregions) then
        ran_counts(region_id(i)) = ran_counts(region_id(i)) + 1
      endif
    enddo

    print *, ''
    print *, 'Region distribution:'
    print '(A8, A12, A12)', 'Region', 'Galaxies', 'Randoms'
    print '(A8, A12, A12)', '------', '--------', '-------'
    do i = 1, Nregions
      print '(I8, I12, I12)', i, gal_counts(i), ran_counts(i)
    enddo
    print *, ''
    print '(A, I8)', '  Total galaxies: ', sum(gal_counts)
    print '(A, I8)', '  Total randoms:  ', sum(ran_counts)
    print '(A, I8)', '  Min gal/region: ', minval(gal_counts)
    print '(A, I8)', '  Max gal/region: ', maxval(gal_counts)
    print '(A, I8)', '  Min ran/region: ', minval(ran_counts)
    print '(A, I8)', '  Max ran/region: ', maxval(ran_counts)

    deallocate(gal_counts)
    deallocate(ran_counts)

  end subroutine print_region_stats

end program test_jackknife
