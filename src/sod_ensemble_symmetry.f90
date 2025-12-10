!*******************************************************************************
!  Symmetry helpers with cached EQMATRIX access for ensemble drivers.
!*******************************************************************************

module sod_ensemble_symmetry
    use sod_ensemble_consts
    use sod_ensemble_energy_calculations, only: get_eqmatrix
    implicit none
    private
    logical, save :: symmetry_ready = .false.
    integer, allocatable, target, save :: cached_eqmatrix(:,:)
    integer, save :: cached_nop = 0
    integer, save :: cached_npos = 0

    public :: symmetry_initialize, symmetry_get_matrix, symmetry_finalize

    contains

    subroutine symmetry_initialize()
        integer, allocatable :: tmp(:,:)
        integer :: nop, npos

        if (symmetry_ready) return

        call get_eqmatrix(tmp, nop, npos)
        if (.not. allocated(tmp) .or. nop <= 0 .or. npos <= 0) then
            if (allocated(tmp)) deallocate(tmp)
            symmetry_ready = .false.
            cached_nop = 0
            cached_npos = 0
            return
        end if

        call move_alloc(tmp, cached_eqmatrix)
        cached_nop = nop
        cached_npos = npos
        symmetry_ready = .true.
    end subroutine symmetry_initialize

    subroutine symmetry_get_matrix(eqmatrix, nop, npos)
        integer, pointer, intent(out) :: eqmatrix(:,:)
        integer, intent(out) :: nop, npos

        if (.not. symmetry_ready) call symmetry_initialize()
        if (.not. symmetry_ready) then
            nullify(eqmatrix)
            nop = 0
            npos = 0
            return
        end if

        eqmatrix => cached_eqmatrix
        nop = cached_nop
        npos = cached_npos
    end subroutine symmetry_get_matrix

    subroutine symmetry_finalize()
        if (allocated(cached_eqmatrix)) deallocate(cached_eqmatrix)
        cached_nop = 0
        cached_npos = 0
        symmetry_ready = .false.
    end subroutine symmetry_finalize

end module sod_ensemble_symmetry
