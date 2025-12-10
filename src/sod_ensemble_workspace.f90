!*******************************************************************************
!    Copyright (c) 2025, Salvador R.G. Balestra
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!******************************************************************************

!*******************************************************************************
!  Workspace manager for ensemble Monte Carlo drivers.
!*******************************************************************************

module sod_ensemble_workspace
    use sod_ensemble_consts
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    implicit none
    private

    type :: mc_workspace_type
        integer :: max_level = 0
        integer :: max_capacity = 0
        integer :: max_sites = 0
        real(dp), pointer :: energies(:) => null()
        real(dp), pointer :: energies_low(:) => null()
        real(dp), pointer :: energies_high(:) => null()
        integer, pointer :: subset_buffer(:,:) => null()
        real(dp), pointer :: low_contribs(:,:) => null()
        real(dp), pointer :: high_contribs(:,:) => null()
        integer, pointer :: accept_attempt(:) => null()
        real(dp), pointer :: gulp_energies(:) => null()
        integer, pointer :: canonical_subset(:) => null()
    end type mc_workspace_type

    type(mc_workspace_type), allocatable, save, target :: workspaces(:)
    logical, save :: workspaces_ready = .false.
    public :: mc_workspace_prepare, mc_workspace_checkout, mc_workspace_finalize

    contains

    ! Prepares thread-local workspaces sized for the expected level of parallelism.
    subroutine mc_workspace_prepare(enable_parallel)
        logical, intent(in) :: enable_parallel
        integer :: desired_threads, i

        desired_threads = 1
        if (enable_parallel) desired_threads = max(1, omp_get_max_threads())

        call release_workspaces()

        allocate(workspaces(desired_threads))
        do i = 1, desired_threads
            call reset_workspace(workspaces(i))
        end do
        workspaces_ready = .true.
    end subroutine mc_workspace_prepare

    ! Provides callers with pointers to the cached buffers sized by ensure_capacity.
    subroutine mc_workspace_checkout(level, capacity, total_sites, subsets, energies, energies_low, energies_high, &
                                     low_contribs, high_contribs, accept_attempt, gulp_energies, canonical_subset)
        integer, intent(in) :: level, capacity, total_sites
        integer, pointer, intent(out) :: subsets(:,:)
        real(dp), pointer, intent(out) :: energies(:)
        real(dp), pointer, intent(out) :: energies_low(:)
        real(dp), pointer, intent(out) :: energies_high(:)
        real(dp), pointer, intent(out) :: low_contribs(:,:)
        real(dp), pointer, intent(out) :: high_contribs(:,:)
        integer, pointer, intent(out) :: accept_attempt(:)
        real(dp), pointer, intent(out) :: gulp_energies(:)
        integer, pointer, intent(out) :: canonical_subset(:)
        type(mc_workspace_type), pointer :: ws

        call ensure_capacity(level, capacity, total_sites, ws)

        subsets => ws%subset_buffer
        energies => ws%energies
        energies_low => ws%energies_low
        energies_high => ws%energies_high
        low_contribs => ws%low_contribs
        high_contribs => ws%high_contribs
        accept_attempt => ws%accept_attempt
        gulp_energies => ws%gulp_energies
        canonical_subset => ws%canonical_subset
    end subroutine mc_workspace_checkout

    subroutine mc_workspace_finalize()
        call release_workspaces()
    end subroutine mc_workspace_finalize

    ! Ensures reusable buffers are large enough for the requested level and sample count.
    subroutine ensure_capacity(level_req, capacity_req, sites_req, ws)
        integer, intent(in) :: level_req, capacity_req, sites_req
        type(mc_workspace_type), pointer, intent(out) :: ws
        integer :: new_level, new_capacity, new_sites
        integer :: prev_level, prev_capacity, prev_sites

        call get_workspace(ws)

        new_level = max(1, level_req)
        new_capacity = max(1, capacity_req)
        new_sites = max(1, sites_req)

        if (.not. associated(ws%subset_buffer)) then
            allocate(ws%subset_buffer(new_level, new_capacity))
        else if (size(ws%subset_buffer,1) < new_level .or. size(ws%subset_buffer,2) < new_capacity) then
            prev_level = size(ws%subset_buffer,1)
            prev_capacity = size(ws%subset_buffer,2)
            deallocate(ws%subset_buffer)
            allocate(ws%subset_buffer(max(new_level, prev_level), max(new_capacity, prev_capacity)))
        end if

        call ensure_real_capacity(ws%energies, new_capacity)
        call ensure_real_capacity(ws%energies_low, new_capacity)
        call ensure_real_capacity(ws%energies_high, new_capacity)
        call ensure_real_matrix(ws%low_contribs, 4, new_capacity)
        call ensure_real_matrix(ws%high_contribs, 4, new_capacity)
        call ensure_int_capacity(ws%accept_attempt, new_capacity)
        call ensure_real_capacity(ws%gulp_energies, new_capacity)

        if (.not. associated(ws%canonical_subset)) then
            allocate(ws%canonical_subset(new_sites))
        else if (size(ws%canonical_subset) < new_sites) then
            prev_sites = size(ws%canonical_subset)
            deallocate(ws%canonical_subset)
            allocate(ws%canonical_subset(max(new_sites, prev_sites)))
        end if

        ws%max_level = size(ws%subset_buffer, 1)
        ws%max_capacity = size(ws%subset_buffer, 2)
        ws%max_sites = size(ws%canonical_subset)
    end subroutine ensure_capacity

    ! Guarantees a real vector has at least the requested capacity.
    subroutine ensure_real_capacity(buffer, desired)
        real(dp), pointer, intent(inout) :: buffer(:)
        integer, intent(in) :: desired
        integer :: prev_capacity

        if (.not. associated(buffer)) then
            allocate(buffer(desired))
        else if (size(buffer) < desired) then
            prev_capacity = size(buffer)
            deallocate(buffer)
            allocate(buffer(max(desired, prev_capacity)))
        end if
    end subroutine ensure_real_capacity

    ! Guarantees an integer vector has at least the requested capacity.
    subroutine ensure_int_capacity(buffer, desired)
        integer, pointer, intent(inout) :: buffer(:)
        integer, intent(in) :: desired
        integer :: prev_capacity

        if (.not. associated(buffer)) then
            allocate(buffer(desired))
        else if (size(buffer) < desired) then
            prev_capacity = size(buffer)
            deallocate(buffer)
            allocate(buffer(max(desired, prev_capacity)))
        end if
    end subroutine ensure_int_capacity

    ! Guarantees a real matrix has at least the requested column capacity.
    subroutine ensure_real_matrix(buffer, rows, desired)
        real(dp), pointer, intent(inout) :: buffer(:,:)
        integer, intent(in) :: rows, desired
        integer :: prev_capacity

        if (.not. associated(buffer)) then
            allocate(buffer(rows, desired))
        else if (size(buffer,1) /= rows .or. size(buffer,2) < desired) then
            prev_capacity = size(buffer,2)
            deallocate(buffer)
            allocate(buffer(rows, max(desired, prev_capacity)))
        end if
    end subroutine ensure_real_matrix

    ! Retrieves the workspace assigned to the current OpenMP thread.
    subroutine get_workspace(ws)
        type(mc_workspace_type), pointer, intent(out) :: ws
        integer :: tid, idx

        if (.not. workspaces_ready) then
            call mc_workspace_prepare(.false.)
        end if

        tid = 0
        !$ tid = omp_get_thread_num()
        idx = tid + 1
        if (idx > size(workspaces)) then
            write(*,'(A,I0,A,I0)') 'Error: threads OpenMP activos (', idx, ') exceden el workspace asignado (', size(workspaces), ').'
            stop 2
        end if
        ws => workspaces(idx)
    end subroutine get_workspace

    ! Releases all allocated buffers from every thread workspace.
    subroutine release_workspaces()
        integer :: i

        if (.not. allocated(workspaces)) return

        do i = 1, size(workspaces)
            call dispose_workspace(workspaces(i))
        end do
        deallocate(workspaces)
        workspaces_ready = .false.
    end subroutine release_workspaces

    ! Resets a workspace, deallocating any backing buffers.
    subroutine reset_workspace(ws)
        type(mc_workspace_type), intent(inout) :: ws

        call dispose_workspace(ws)
    end subroutine reset_workspace

    ! Deallocates every pointer owned by a workspace instance.
    subroutine dispose_workspace(ws)
        type(mc_workspace_type), intent(inout) :: ws

        if (associated(ws%subset_buffer)) then
            deallocate(ws%subset_buffer)
            nullify(ws%subset_buffer)
        end if
        if (associated(ws%energies)) then
            deallocate(ws%energies)
            nullify(ws%energies)
        end if
        if (associated(ws%energies_low)) then
            deallocate(ws%energies_low)
            nullify(ws%energies_low)
        end if
        if (associated(ws%energies_high)) then
            deallocate(ws%energies_high)
            nullify(ws%energies_high)
        end if
        if (associated(ws%low_contribs)) then
            deallocate(ws%low_contribs)
            nullify(ws%low_contribs)
        end if
        if (associated(ws%high_contribs)) then
            deallocate(ws%high_contribs)
            nullify(ws%high_contribs)
        end if
        if (associated(ws%accept_attempt)) then
            deallocate(ws%accept_attempt)
            nullify(ws%accept_attempt)
        end if
        if (associated(ws%gulp_energies)) then
            deallocate(ws%gulp_energies)
            nullify(ws%gulp_energies)
        end if
        if (associated(ws%canonical_subset)) then
            deallocate(ws%canonical_subset)
            nullify(ws%canonical_subset)
        end if
        ws%max_level = 0
        ws%max_capacity = 0
        ws%max_sites = 0
    end subroutine dispose_workspace

end module sod_ensemble_workspace
