!*******************************************************************************
! Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
!
! This file is part of the SOD package.
!
! SOD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SOD is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!*******************************************************************************
!
! Modern shared-module update of the classic configuration enumeration logic.
!
!*******************************************************************************
module configurations
    use consts, only: error_unit, ip
    use utils, only: binomial_int64, next_combination, sort_int_ascending
    implicit none
    private

    public :: canonicalize_subset
    public :: collect_subset_stabilizer_operators
    public :: enumerate_unique_subsets
    public :: find_subset_index
    public :: read_outsod_file
    public :: read_outsod_header
    public :: write_outsod_file
    public :: write_outsod_unit

contains

    subroutine enumerate_unique_subsets(total_sites, level, eqmatrix, nop, unique_subsets, unique_deg, unique_count, &
        neighbor_paths, used_recursive, candidate_count, chosen_neighbor)
        integer, intent(in) :: total_sites, level, nop
        integer, intent(in) :: eqmatrix(:,:)
        integer, allocatable, intent(out) :: unique_subsets(:,:)
        integer, allocatable, intent(out) :: unique_deg(:)
        integer, intent(out) :: unique_count
        character(len=*), intent(in), optional :: neighbor_paths(:)
        logical, intent(out), optional :: used_recursive
        integer(ip), intent(out), optional :: candidate_count
        character(len=*), intent(out), optional :: chosen_neighbor
        integer :: neighbor_level, neighbor_unique_count
        integer(ip) :: recursive_candidate_count
        logical :: recursive_ok
        character(len=1024) :: recursive_neighbor

        if (level < 0) then
            write(error_unit,'(A)') 'Error: enumerate_unique_subsets received a negative level.'
            stop 1
        end if

        if (present(used_recursive)) used_recursive = .false.
        if (present(candidate_count)) candidate_count = 0_ip
        if (present(chosen_neighbor)) chosen_neighbor = ''

        if (level == 0) then
            allocate(unique_subsets(0,1))
            allocate(unique_deg(1))
            unique_deg(1) = 1
            unique_count = 1
            return
        end if

        recursive_ok = .false.
        recursive_neighbor = ''
        recursive_candidate_count = 0_ip
        if (present(neighbor_paths)) then
            call enumerate_unique_subsets_recursive(total_sites, level, eqmatrix, nop, neighbor_paths, unique_subsets, &
                unique_deg, unique_count, recursive_ok, recursive_candidate_count, recursive_neighbor, neighbor_level, &
                neighbor_unique_count)
        end if

        if (recursive_ok) then
            if (present(used_recursive)) used_recursive = .true.
            if (present(candidate_count)) candidate_count = recursive_candidate_count
            if (present(chosen_neighbor)) chosen_neighbor = trim(recursive_neighbor)
            return
        end if

        call enumerate_unique_subsets_direct(total_sites, level, eqmatrix, nop, unique_subsets, unique_deg, unique_count)
    end subroutine enumerate_unique_subsets

    subroutine enumerate_unique_subsets_direct(total_sites, level, eqmatrix, nop, unique_subsets, unique_deg, unique_count)
        integer, intent(in) :: total_sites, level, nop
        integer, intent(in) :: eqmatrix(:,:)
        integer, allocatable, intent(out) :: unique_subsets(:,:)
        integer, allocatable, intent(out) :: unique_deg(:)
        integer, intent(out) :: unique_count
        integer :: ntc
        integer :: unique_capacity
        integer :: current_rank, transformed_rank
        integer :: op
        integer(ip), allocatable :: binom(:,:)
        logical, allocatable :: visited(:)
        integer, allocatable :: subset(:), transformed(:)
        integer, allocatable :: work_subsets(:,:), work_deg(:)

        ntc = validated_total_combination_count(total_sites, level)
        call build_binomial_table(level, total_sites, binom)
        allocate(visited(ntc))
        visited = .false.

        unique_capacity = max(16, level)
        allocate(work_subsets(level, unique_capacity))
        allocate(work_deg(unique_capacity))
        work_subsets = 0
        work_deg = 0
        unique_count = 0

        allocate(subset(level))
        allocate(transformed(level))
        do current_rank = 1, level
            subset(current_rank) = current_rank
        end do

        do
            current_rank = rank_subset_index(subset, level, total_sites, binom)
            if (.not. visited(current_rank)) then
                unique_count = unique_count + 1
                call ensure_unique_capacity(level, work_subsets, work_deg, unique_capacity, unique_count)
                work_subsets(1:level, unique_count) = subset(1:level)
                work_deg(unique_count) = 1
                visited(current_rank) = .true.

                do op = 1, nop
                    transformed(1:level) = eqmatrix(op, subset(1:level))
                    call sort_int_ascending(transformed, level)
                    transformed_rank = rank_subset_index(transformed, level, total_sites, binom)
                    if (transformed_rank /= current_rank .and. .not. visited(transformed_rank)) then
                        visited(transformed_rank) = .true.
                        work_deg(unique_count) = work_deg(unique_count) + 1
                    end if
                end do
            end if
            if (.not. next_combination(subset, total_sites)) exit
        end do

        deallocate(subset)
        deallocate(transformed)
        deallocate(visited)
        deallocate(binom)

        allocate(unique_subsets(level, unique_count))
        allocate(unique_deg(unique_count))
        if (unique_count > 0) then
            unique_subsets(1:level,1:unique_count) = work_subsets(1:level,1:unique_count)
            unique_deg(1:unique_count) = work_deg(1:unique_count)
        end if

        deallocate(work_subsets)
        deallocate(work_deg)
    end subroutine enumerate_unique_subsets_direct

    subroutine enumerate_unique_subsets_recursive(total_sites, level, eqmatrix, nop, neighbor_paths, unique_subsets, &
        unique_deg, unique_count, success, candidate_count, chosen_neighbor, neighbor_level, neighbor_unique_count)
        integer, intent(in) :: total_sites, level, nop
        integer, intent(in) :: eqmatrix(:,:)
        character(len=*), intent(in) :: neighbor_paths(:)
        integer, allocatable, intent(out) :: unique_subsets(:,:)
        integer, allocatable, intent(out) :: unique_deg(:)
        integer, intent(out) :: unique_count
        logical, intent(out) :: success
        integer(ip), intent(out) :: candidate_count
        character(len=*), intent(out) :: chosen_neighbor
        integer, intent(out) :: neighbor_level, neighbor_unique_count
        integer :: npos_check
        integer, allocatable :: parent_subsets(:,:), parent_deg(:)
        logical :: found_neighbor

        success = .false.
        candidate_count = 0_ip
        chosen_neighbor = ''
        neighbor_level = -1
        neighbor_unique_count = 0

        call select_best_neighbor_outsod(total_sites, level, neighbor_paths, chosen_neighbor, neighbor_level, &
            neighbor_unique_count, found_neighbor)
        if (.not. found_neighbor) return

        call read_outsod_file(trim(chosen_neighbor), neighbor_level, npos_check, parent_unique_count=neighbor_unique_count, &
            parent_subsets=parent_subsets, parent_deg=parent_deg, success=found_neighbor)
        if (.not. found_neighbor) return
        if (npos_check /= total_sites) then
            deallocate(parent_subsets)
            deallocate(parent_deg)
            return
        end if

        call enumerate_from_parent_candidates(total_sites, level, eqmatrix, nop, neighbor_level, parent_subsets, &
            neighbor_unique_count, unique_subsets, unique_deg, unique_count, candidate_count)
        deallocate(parent_subsets)
        deallocate(parent_deg)
        success = .true.
    end subroutine enumerate_unique_subsets_recursive

    subroutine canonicalize_subset(subset, level, eqmatrix, nop, canonical)
        integer, intent(in) :: subset(:)
        integer, intent(in) :: level, nop
        integer, intent(in) :: eqmatrix(:,:)
        integer, intent(out) :: canonical(:)
        integer :: temp(level)
        integer :: op

        canonical(1:level) = subset(1:level)
        call sort_int_ascending(canonical, level)

        do op = 1, nop
            temp(1:level) = eqmatrix(op, subset(1:level))
            call sort_int_ascending(temp, level)
            if (lexicographically_less(temp, canonical, level)) then
                canonical(1:level) = temp(1:level)
            end if
        end do
    end subroutine canonicalize_subset

    subroutine collect_subset_stabilizer_operators(subset, level, eqmatrix, nop, stabilizer_ops, stabilizer_count)
        integer, intent(in) :: subset(:)
        integer, intent(in) :: level, nop
        integer, intent(in) :: eqmatrix(:,:)
        integer, intent(out) :: stabilizer_ops(:)
        integer, intent(out) :: stabilizer_count
        integer :: canonical(level)
        integer :: temp(level)
        integer :: op

        if (size(stabilizer_ops) < nop) then
            write(error_unit,'(A)') 'Error: stabilizer operator buffer is too small.'
            stop 1
        end if

        canonical(1:level) = subset(1:level)
        call sort_int_ascending(canonical, level)

        stabilizer_count = 0
        do op = 1, nop
            temp(1:level) = eqmatrix(op, subset(1:level))
            call sort_int_ascending(temp, level)
            if (all(temp(1:level) == canonical(1:level))) then
                stabilizer_count = stabilizer_count + 1
                stabilizer_ops(stabilizer_count) = op
            end if
        end do
    end subroutine collect_subset_stabilizer_operators

    pure integer function find_subset_index(subset, level, stored_subsets, stored_count)
        integer, intent(in) :: subset(:)
        integer, intent(in) :: level, stored_count
        integer, intent(in) :: stored_subsets(:,:)
        integer :: idx, j
        logical :: match

        find_subset_index = 0
        do idx = 1, stored_count
            match = .true.
            do j = 1, level
                if (stored_subsets(j, idx) /= subset(j)) then
                    match = .false.
                    exit
                end if
            end do
            if (match) then
                find_subset_index = idx
                return
            end if
        end do
    end function find_subset_index

    logical function write_outsod_file(filename, level, total_sites, unique_count, unique_deg, unique_subsets, &
        legacy_format) result(ok)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_deg(:)
        integer, intent(in), optional :: unique_subsets(:,:)
        logical, intent(in), optional :: legacy_format
        integer :: unit_outsod
        integer :: ios

        ok = .false.
        open(newunit=unit_outsod, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) return
        call write_outsod_unit(unit_outsod, level, total_sites, unique_count, unique_deg, unique_subsets, legacy_format)
        close(unit_outsod)
        ok = .true.
    end function write_outsod_file

    subroutine write_outsod_unit(unit_outsod, level, total_sites, unique_count, unique_deg, unique_subsets, &
        legacy_format)
        integer, intent(in) :: unit_outsod
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_deg(:)
        integer, intent(in), optional :: unique_subsets(:,:)
        logical, intent(in), optional :: legacy_format
        logical :: use_legacy
        integer :: idx, site

        use_legacy = .false.
        if (present(legacy_format)) use_legacy = legacy_format

        if (use_legacy) then
            write(unit_outsod,*) level, ' substitutions in ', total_sites, 'sites'
            write(unit_outsod,*) int(unique_count, kind=8), ' configurations'
        else
            write(unit_outsod,'(I12,"  substitutions in",I12," sites")') level, total_sites
            write(unit_outsod,'(I12,"  configurations")') unique_count
        end if

        do idx = 1, unique_count
            if (use_legacy) then
                write(unit_outsod,'(I6,1X,I6)', advance='no') idx, unique_deg(idx)
            else
                write(unit_outsod,'(I6,1X,I12)', advance='no') idx, unique_deg(idx)
            end if
            if (level > 0 .and. present(unique_subsets)) then
                do site = 1, level
                    if (use_legacy) then
                        write(unit_outsod,'(1X,I4)', advance='no') unique_subsets(site, idx)
                    else
                        write(unit_outsod,'(1X,I6)', advance='no') unique_subsets(site, idx)
                    end if
                end do
            end if
            write(unit_outsod,*)
        end do
    end subroutine write_outsod_unit

    pure logical function lexicographically_less(a, b, level)
        integer, intent(in) :: a(:), b(:)
        integer, intent(in) :: level
        integer :: i

        do i = 1, level
            if (a(i) < b(i)) then
                lexicographically_less = .true.
                return
            else if (a(i) > b(i)) then
                lexicographically_less = .false.
                return
            end if
        end do
        lexicographically_less = .false.
    end function lexicographically_less

    subroutine ensure_unique_capacity(level, subsets, deg, capacity, required)
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: subsets(:,:)
        integer, allocatable, intent(inout) :: deg(:)
        integer :: new_capacity
        integer, allocatable :: tmp_subsets(:,:)
        integer, allocatable :: tmp_deg(:)

        if (required <= capacity) return

        new_capacity = max(required, capacity + max(4, capacity / 2))
        if (new_capacity <= capacity) new_capacity = capacity + 4

        allocate(tmp_subsets(level, new_capacity))
        allocate(tmp_deg(new_capacity))
        tmp_subsets = 0
        tmp_deg = 0
        if (capacity > 0) then
            tmp_subsets(:,1:capacity) = subsets(:,1:capacity)
            tmp_deg(1:capacity) = deg(1:capacity)
        end if

        call move_alloc(tmp_subsets, subsets)
        call move_alloc(tmp_deg, deg)
        capacity = new_capacity
    end subroutine ensure_unique_capacity

    subroutine select_best_neighbor_outsod(total_sites, level, neighbor_paths, chosen_neighbor, neighbor_level, &
        neighbor_unique_count, found_neighbor)
        integer, intent(in) :: total_sites, level
        character(len=*), intent(in) :: neighbor_paths(:)
        character(len=*), intent(out) :: chosen_neighbor
        integer, intent(out) :: neighbor_level, neighbor_unique_count
        logical, intent(out) :: found_neighbor
        integer :: idx, candidate_level, candidate_unique_count, npos_check
        integer(ip) :: best_candidate_count, current_candidate_count
        logical :: ok

        chosen_neighbor = ''
        neighbor_level = -1
        neighbor_unique_count = 0
        found_neighbor = .false.
        best_candidate_count = 0_ip

        do idx = 1, size(neighbor_paths)
            if (len_trim(neighbor_paths(idx)) == 0) cycle
            call read_outsod_header(trim(neighbor_paths(idx)), candidate_level, npos_check, candidate_unique_count, ok)
            if (.not. ok) cycle
            if (npos_check /= total_sites) cycle
            if (abs(candidate_level - level) /= 1) cycle
            if (candidate_level < level) then
                current_candidate_count = int(candidate_unique_count, ip) * int(total_sites - candidate_level, ip)
            else
                current_candidate_count = int(candidate_unique_count, ip) * int(candidate_level, ip)
            end if
            if (.not. found_neighbor .or. current_candidate_count < best_candidate_count) then
                found_neighbor = .true.
                best_candidate_count = current_candidate_count
                chosen_neighbor = trim(neighbor_paths(idx))
                neighbor_level = candidate_level
                neighbor_unique_count = candidate_unique_count
            end if
        end do
    end subroutine select_best_neighbor_outsod

    subroutine read_outsod_header(filename, level, total_sites, unique_count, success)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: level, total_sites, unique_count
        logical, intent(out) :: success
        integer :: unit_in, ios, marker
        character(len=1024) :: line

        success = .false.
        level = -1
        total_sites = -1
        unique_count = 0

        open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        read(unit_in,'(A)',iostat=ios) line
        if (ios /= 0) then
            close(unit_in)
            return
        end if
        read(line,*,iostat=ios) level
        if (ios /= 0) then
            close(unit_in)
            return
        end if
        marker = index(line, ' in ')
        if (marker <= 0) then
            close(unit_in)
            return
        end if
        read(line(marker+4:),*,iostat=ios) total_sites
        if (ios /= 0) then
            close(unit_in)
            return
        end if

        read(unit_in,'(A)',iostat=ios) line
        if (ios /= 0) then
            close(unit_in)
            return
        end if
        read(line,*,iostat=ios) unique_count
        close(unit_in)
        if (ios /= 0) return
        success = .true.
    end subroutine read_outsod_header

    subroutine read_outsod_file(filename, parent_level, total_sites, parent_unique_count, parent_subsets, parent_deg, &
        success)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: parent_level, total_sites, parent_unique_count
        integer, allocatable, intent(out) :: parent_subsets(:,:)
        integer, allocatable, intent(out) :: parent_deg(:)
        logical, intent(out) :: success
        integer :: unit_in, ios, idx, dummy
        character(len=1024) :: line

        success = .false.
        call read_outsod_header(filename, parent_level, total_sites, parent_unique_count, success)
        if (.not. success) return

        open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            success = .false.
            return
        end if
        read(unit_in,'(A)',iostat=ios) line
        read(unit_in,'(A)',iostat=ios) line
        if (ios /= 0) then
            close(unit_in)
            success = .false.
            return
        end if

        allocate(parent_subsets(parent_level, parent_unique_count))
        allocate(parent_deg(parent_unique_count))
        if (parent_level > 0) parent_subsets = 0
        parent_deg = 0

        if (parent_level == 0) then
            do idx = 1, parent_unique_count
                read(unit_in,*,iostat=ios) dummy, parent_deg(idx)
                if (ios /= 0) exit
            end do
        else
            do idx = 1, parent_unique_count
                read(unit_in,*,iostat=ios) dummy, parent_deg(idx), parent_subsets(1:parent_level, idx)
                if (ios /= 0) exit
            end do
        end if
        close(unit_in)

        if (ios /= 0) then
            deallocate(parent_subsets)
            deallocate(parent_deg)
            success = .false.
            return
        end if
        success = .true.
    end subroutine read_outsod_file

    subroutine enumerate_from_parent_candidates(total_sites, level, eqmatrix, nop, parent_level, parent_subsets, &
        parent_unique_count, unique_subsets, unique_deg, unique_count, candidate_count)
        integer, intent(in) :: total_sites, level, nop, parent_level, parent_unique_count
        integer, intent(in) :: eqmatrix(:,:)
        integer, intent(in) :: parent_subsets(:,:)
        integer, allocatable, intent(out) :: unique_subsets(:,:)
        integer, allocatable, intent(out) :: unique_deg(:)
        integer, intent(out) :: unique_count
        integer(ip), intent(out) :: candidate_count
        integer :: ncand_max, ncand
        integer :: ic, pos, jpar, j_remove, jj, k
        integer :: ntc, unique_capacity
        integer :: current_rank, transformed_rank, op, cand
        integer(ip), allocatable :: binom(:,:)
        logical, allocatable :: visited(:)
        integer, allocatable :: candidates(:,:), transformed(:)
        integer, allocatable :: work_subsets(:,:), work_deg(:)

        if (parent_level == level) then
            write(error_unit,'(A)') 'Error: recursive enumeration requires an adjacent OUTSOD level.'
            stop 1
        end if

        if (parent_level < level) then
            ncand_max = parent_unique_count * (total_sites - parent_level)
            allocate(candidates(level, ncand_max))
            ncand = 0
            do ic = 1, parent_unique_count
                jpar = 1
                do pos = 1, total_sites
                    if (jpar <= parent_level .and. parent_subsets(jpar, ic) == pos) then
                        jpar = jpar + 1
                    else
                        ncand = ncand + 1
                        if (jpar > 1) candidates(1:jpar-1, ncand) = parent_subsets(1:jpar-1, ic)
                        candidates(jpar, ncand) = pos
                        if (jpar <= parent_level) candidates(jpar+1:level, ncand) = parent_subsets(jpar:parent_level, ic)
                    end if
                end do
            end do
        else
            ncand_max = parent_unique_count * parent_level
            allocate(candidates(level, ncand_max))
            ncand = 0
            do ic = 1, parent_unique_count
                do j_remove = 1, parent_level
                    ncand = ncand + 1
                    k = 0
                    do jj = 1, parent_level
                        if (jj /= j_remove) then
                            k = k + 1
                            candidates(k, ncand) = parent_subsets(jj, ic)
                        end if
                    end do
                end do
            end do
        end if
        candidate_count = int(ncand, ip)

        ntc = validated_total_combination_count(total_sites, level)
        call build_binomial_table(level, total_sites, binom)
        allocate(visited(ntc))
        visited = .false.

        unique_capacity = max(16, level)
        allocate(work_subsets(level, unique_capacity))
        allocate(work_deg(unique_capacity))
        work_subsets = 0
        work_deg = 0
        unique_count = 0

        allocate(transformed(level))
        do cand = 1, ncand
            current_rank = rank_subset_index(candidates(1:level, cand), level, total_sites, binom)
            if (visited(current_rank)) cycle

            unique_count = unique_count + 1
            call ensure_unique_capacity(level, work_subsets, work_deg, unique_capacity, unique_count)
            work_subsets(1:level, unique_count) = candidates(1:level, cand)
            work_deg(unique_count) = 1
            visited(current_rank) = .true.

            do op = 1, nop
                transformed(1:level) = eqmatrix(op, candidates(1:level, cand))
                call sort_int_ascending(transformed, level)
                transformed_rank = rank_subset_index(transformed, level, total_sites, binom)
                if (transformed_rank /= current_rank .and. .not. visited(transformed_rank)) then
                    visited(transformed_rank) = .true.
                    work_deg(unique_count) = work_deg(unique_count) + 1
                end if
            end do
        end do

        allocate(unique_subsets(level, unique_count))
        allocate(unique_deg(unique_count))
        if (unique_count > 0) then
            unique_subsets(1:level,1:unique_count) = work_subsets(1:level,1:unique_count)
            unique_deg(1:unique_count) = work_deg(1:unique_count)
        end if

        deallocate(transformed)
        deallocate(visited)
        deallocate(binom)
        deallocate(candidates)
        deallocate(work_subsets)
        deallocate(work_deg)
    end subroutine enumerate_from_parent_candidates

    subroutine build_binomial_table(level, total_sites, binom)
        integer, intent(in) :: level, total_sites
        integer(ip), allocatable, intent(out) :: binom(:,:)
        integer :: k, n

        allocate(binom(0:level, 0:total_sites))
        binom = 0_ip
        binom(0, :) = 1_ip
        do k = 1, level
            do n = k, total_sites
                binom(k, n) = binom(k - 1, n - 1) + binom(k, n - 1)
            end do
        end do
    end subroutine build_binomial_table

    integer function rank_subset_index(subset, level, total_sites, binom)
        integer, intent(in) :: subset(:)
        integer, intent(in) :: level, total_sites
        integer(ip), intent(in) :: binom(0:,0:)
        integer(ip) :: rank_value
        integer :: i, remaining_n, remaining_k

        rank_value = binom(level, total_sites)
        do i = 1, level
            remaining_n = total_sites - subset(i)
            remaining_k = level - i + 1
            if (remaining_n >= remaining_k) then
                rank_value = rank_value - binom(remaining_k, remaining_n)
            end if
        end do

        if (rank_value < 1_ip .or. rank_value > int(huge(0), ip)) then
            write(error_unit,'(A)') 'Error: subset combinatorial rank is out of range.'
            stop 1
        end if
        rank_subset_index = int(rank_value)
    end function rank_subset_index

    integer function validated_total_combination_count(total_sites, level)
        integer, intent(in) :: total_sites, level
        integer(ip) :: total_comb

        total_comb = binomial_int64(total_sites, level)
        if (total_comb <= 0_ip) then
            write(error_unit,'(A)') 'Error: total configuration count overflowed in enumerate_unique_subsets.'
            stop 1
        end if
        if (total_comb > int(huge(0), ip)) then
            write(error_unit,'(A)') 'Error: total configuration count exceeds addressable array size.'
            stop 1
        end if
        validated_total_combination_count = int(total_comb)
    end function validated_total_combination_count

end module configurations
