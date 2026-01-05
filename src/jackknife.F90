!------------------------------------------------------------
! Jackknife Module
!
! This module provides subroutines for partitioning spatial data points
! into jackknife subsamples based on equal point density partitioning.
! This is a Fortran implementation of the jk.sh shell script.
!
! The partitioning works by:
! 1. Sorting all points by X coordinate, splitting into N equal groups
! 2. Within each X group, sorting by Y, splitting into N equal groups
! 3. For 3D: Within each (X,Y) group, sorting by Z, splitting into N equal groups
! 4. Assigning region IDs: 1 + X_region + N*Y_region + N^2*Z_region
!
!------------------------------------------------------------
module jackknife_module
  implicit none

  integer, parameter :: jk_dp = selected_real_kind(15, 307)

contains

!------------------------------------------------------------
! Main jackknife assignment subroutine
!
! Arguments:
!   coords    - Input coordinates array (D x Npoints)
!   Npoints   - Number of points
!   Ndiv      - Number of divisions in each dimension
!   Ndim      - Dimensionality (2 or 3)
!   region_id - Output array of region IDs (1 to Ndiv^Ndim)
!------------------------------------------------------------
subroutine jackknife_assign(coords, Npoints, Ndiv, Ndim, region_id)
  implicit none

  integer, intent(in) :: Npoints, Ndiv, Ndim
  real(jk_dp), intent(in) :: coords(Ndim, Npoints)
  integer, intent(out) :: region_id(Npoints)

  ! Local variables
  integer, allocatable :: indices(:), x_region(:), y_region(:), z_region(:)
  integer, allocatable :: sorted_idx(:), group_start(:), group_end(:)
  integer :: i, j, k, n, group_size, remainder, start_idx, end_idx
  integer :: x_grp, y_grp

  ! Allocate working arrays
  allocate(indices(Npoints))
  allocate(x_region(Npoints))
  allocate(y_region(Npoints))
  if (Ndim == 3) allocate(z_region(Npoints))
  allocate(sorted_idx(Npoints))
  allocate(group_start(Ndiv))
  allocate(group_end(Ndiv))

  ! Initialize indices
  do i = 1, Npoints
    indices(i) = i
  enddo

  ! Step 1: Sort all points by X coordinate and assign X regions
  call index_sort_by_coord(coords, Ndim, Npoints, indices, 1, sorted_idx)

  ! Calculate group sizes for X partition
  call calculate_group_bounds(Npoints, Ndiv, group_start, group_end)

  ! Assign X region IDs (0 to Ndiv-1)
  do i = 1, Ndiv
    do j = group_start(i), group_end(i)
      x_region(sorted_idx(j)) = i - 1
    enddo
  enddo

  ! Step 2: For each X region, sort by Y and assign Y regions
  do x_grp = 0, Ndiv - 1
    ! Collect indices of points in this X region
    n = 0
    do i = 1, Npoints
      if (x_region(i) == x_grp) then
        n = n + 1
        indices(n) = i
      endif
    enddo

    if (n > 0) then
      ! Sort these points by Y coordinate
      call index_sort_by_coord(coords, Ndim, n, indices(1:n), 2, sorted_idx(1:n))

      ! Calculate group sizes for Y partition within this X region
      call calculate_group_bounds(n, Ndiv, group_start, group_end)

      ! Assign Y region IDs
      do i = 1, Ndiv
        do j = group_start(i), group_end(i)
          y_region(sorted_idx(j)) = i - 1
        enddo
      enddo
    endif
  enddo

  ! Step 3: For 3D, for each (X,Y) region, sort by Z and assign Z regions
  if (Ndim == 3) then
    do x_grp = 0, Ndiv - 1
      do y_grp = 0, Ndiv - 1
        ! Collect indices of points in this (X,Y) region
        n = 0
        do i = 1, Npoints
          if (x_region(i) == x_grp .and. y_region(i) == y_grp) then
            n = n + 1
            indices(n) = i
          endif
        enddo

        if (n > 0) then
          ! Sort these points by Z coordinate
          call index_sort_by_coord(coords, Ndim, n, indices(1:n), 3, sorted_idx(1:n))

          ! Calculate group sizes for Z partition
          call calculate_group_bounds(n, Ndiv, group_start, group_end)

          ! Assign Z region IDs
          do i = 1, Ndiv
            do j = group_start(i), group_end(i)
              z_region(sorted_idx(j)) = i - 1
            enddo
          enddo
        endif
      enddo
    enddo
  endif

  ! Step 4: Calculate final region IDs
  if (Ndim == 2) then
    do i = 1, Npoints
      region_id(i) = 1 + x_region(i) + Ndiv * y_region(i)
    enddo
  else ! Ndim == 3
    do i = 1, Npoints
      region_id(i) = 1 + x_region(i) + Ndiv * y_region(i) + Ndiv * Ndiv * z_region(i)
    enddo
  endif

  ! Cleanup
  deallocate(indices)
  deallocate(x_region)
  deallocate(y_region)
  if (allocated(z_region)) deallocate(z_region)
  deallocate(sorted_idx)
  deallocate(group_start)
  deallocate(group_end)

end subroutine jackknife_assign

!------------------------------------------------------------
! Calculate group boundaries for equal-size partitioning
! This matches the shell script behavior: (total + N - 1) / N lines per group
!------------------------------------------------------------
subroutine calculate_group_bounds(Npoints, Ndiv, group_start, group_end)
  implicit none

  integer, intent(in) :: Npoints, Ndiv
  integer, intent(out) :: group_start(Ndiv), group_end(Ndiv)

  integer :: lines_per_group, i, current_pos

  ! Match shell script: lines_per_file = (total_lines + num_files - 1) / num_files
  lines_per_group = (Npoints + Ndiv - 1) / Ndiv

  current_pos = 1
  do i = 1, Ndiv
    group_start(i) = current_pos
    group_end(i) = min(current_pos + lines_per_group - 1, Npoints)
    current_pos = group_end(i) + 1
  enddo

end subroutine calculate_group_bounds

!------------------------------------------------------------
! Index sort by a specific coordinate dimension
! Uses simple insertion sort for stability (matches 'sort -g' behavior)
!------------------------------------------------------------
subroutine index_sort_by_coord(coords, Ndim, N, indices, coord_dim, sorted_indices)
  implicit none

  integer, intent(in) :: Ndim, N, coord_dim
  real(jk_dp), intent(in) :: coords(Ndim, *)
  integer, intent(in) :: indices(N)
  integer, intent(out) :: sorted_indices(N)

  integer :: i, j, key_idx
  real(jk_dp) :: key_val
  real(jk_dp), allocatable :: values(:)

  allocate(values(N))

  ! Copy indices and extract coordinate values
  do i = 1, N
    sorted_indices(i) = indices(i)
    values(i) = coords(coord_dim, indices(i))
  enddo

  ! Stable insertion sort
  do i = 2, N
    key_idx = sorted_indices(i)
    key_val = values(i)
    j = i - 1

    do while (j >= 1)
      if (values(j) <= key_val) exit
      sorted_indices(j + 1) = sorted_indices(j)
      values(j + 1) = values(j)
      j = j - 1
    enddo

    sorted_indices(j + 1) = key_idx
    values(j + 1) = key_val
  enddo

  deallocate(values)

end subroutine index_sort_by_coord

!------------------------------------------------------------
! Wrapper subroutine for use with the 2pcf program
! This assigns jackknife regions to both galaxy and random points
!
! Arguments:
!   my_array  - Coordinate array from 2pcf (d x (Ndata+Nrand))
!   Ndata     - Number of galaxy points
!   Nrand     - Number of random points
!   Ndiv      - Number of divisions per dimension
!   Ndim      - Dimensionality (2 or 3)
!   boot      - Output array for region assignments
!------------------------------------------------------------
subroutine jackknife_2pcf(my_array, Ndata, Nrand, Ndiv, Ndim, boot)
  implicit none

  integer, intent(in) :: Ndata, Nrand, Ndiv, Ndim
  real(jk_dp), intent(in) :: my_array(Ndim, Ndata + Nrand)
  integer, intent(out) :: boot(Ndata + Nrand)

  integer :: Ntotal

  Ntotal = Ndata + Nrand

  ! Assign jackknife regions to all points (galaxies and randoms together)
  call jackknife_assign(my_array, Ntotal, Ndiv, Ndim, boot)

end subroutine jackknife_2pcf

!------------------------------------------------------------
! Utility: Quick sort for larger arrays (more efficient than insertion sort)
! Uses a hybrid approach with insertion sort for small subarrays
!------------------------------------------------------------
recursive subroutine quicksort_indexed(values, indices, left, right)
  implicit none

  real(jk_dp), intent(inout) :: values(:)
  integer, intent(inout) :: indices(:)
  integer, intent(in) :: left, right

  integer :: i, j, pivot_idx
  real(jk_dp) :: pivot_val, temp_val
  integer :: temp_idx

  ! Use insertion sort for small arrays
  if (right - left < 10) then
    call insertion_sort_indexed(values, indices, left, right)
    return
  endif

  ! Choose pivot (median of three)
  pivot_idx = left + (right - left) / 2
  pivot_val = values(pivot_idx)

  i = left
  j = right

  do while (i <= j)
    do while (values(i) < pivot_val)
      i = i + 1
    enddo
    do while (values(j) > pivot_val)
      j = j - 1
    enddo

    if (i <= j) then
      ! Swap
      temp_val = values(i)
      values(i) = values(j)
      values(j) = temp_val

      temp_idx = indices(i)
      indices(i) = indices(j)
      indices(j) = temp_idx

      i = i + 1
      j = j - 1
    endif
  enddo

  if (left < j) call quicksort_indexed(values, indices, left, j)
  if (i < right) call quicksort_indexed(values, indices, i, right)

end subroutine quicksort_indexed

subroutine insertion_sort_indexed(values, indices, left, right)
  implicit none

  real(jk_dp), intent(inout) :: values(:)
  integer, intent(inout) :: indices(:)
  integer, intent(in) :: left, right

  integer :: i, j, key_idx
  real(jk_dp) :: key_val

  do i = left + 1, right
    key_val = values(i)
    key_idx = indices(i)
    j = i - 1

    do while (j >= left)
      if (values(j) <= key_val) exit
      values(j + 1) = values(j)
      indices(j + 1) = indices(j)
      j = j - 1
    enddo

    values(j + 1) = key_val
    indices(j + 1) = key_idx
  enddo

end subroutine insertion_sort_indexed

end module jackknife_module
