!*******************************************************************************
! Copyright (c) 2025, Salvador R.G. Balestra
!
! This file is part of the SOD package.
!
! SOD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!*******************************************************************************

!*******************************************************************************
! Symmetry helpers and EQMATRIX ownership for the ensemble drivers.
!*******************************************************************************

! Module `symmetry` builds, caches, validates, and applies symmetry mappings.
! El módulo `symmetry` construye, almacena, valida y aplica mapeos de simetría.
module symmetry
    use consts, only: dp, error_unit
    use inputs, only: sgo_file_data, insod_file_data, core_read_sgo_file => read_sgo_file, &
        core_read_insod_file => read_insod_file
    implicit none
    private

    logical, save :: symmetry_ready = .false.
    integer, allocatable, target, save :: cached_eqmatrix(:,:)
    integer, save :: cached_nop = 0
    integer, save :: cached_npos = 0

    character(len=*), parameter :: default_eqmatrix_filename = 'EQMATRIX'
    real(dp), parameter :: symop_tol = 1.0e-8_dp
    real(dp), parameter :: coord_tol = 1.0e-3_dp
    real(dp), parameter :: wrap_tol = 1.0e-4_dp

    public :: symmetry_finalize
    public :: symmetry_generate_matrix_from_files
    public :: symmetry_get_matrix
    public :: symmetry_initialize
    public :: symmetry_read_matrix_file
    public :: symmetry_restrict_supercell_operators
    public :: symmetry_validate_matrix
    public :: symmetry_wrap_fractional_vector
    public :: symmetry_write_matrix

contains

    subroutine symmetry_initialize(force_generate)
        logical, intent(in), optional :: force_generate
        logical :: must_generate
        logical :: found
        integer, allocatable :: tmp(:,:)
        integer :: nop, npos

        must_generate = .false.
        if (present(force_generate)) must_generate = force_generate

        if (symmetry_ready .and. .not. must_generate) return

        call symmetry_finalize()

        found = .false.
        if (.not. must_generate) then
            call symmetry_read_matrix_file(default_eqmatrix_filename, tmp, nop, npos, found=found, quiet_missing=.true.)
        end if
        if (.not. found) then
            call symmetry_generate_matrix_from_files(tmp, nop, npos)
        end if

        call cache_matrix(tmp, nop, npos)
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

    subroutine symmetry_write_matrix(filename)
        character(len=*), intent(in), optional :: filename
        character(len=512) :: output_name
        integer :: unit_id, ios, op

        if (.not. symmetry_ready) call symmetry_initialize()
        if (.not. symmetry_ready) then
            write(error_unit,'(A)') 'Error: EQMATRIX is not available for writing.'
            return
        end if

        output_name = default_eqmatrix_filename
        if (present(filename)) output_name = filename
        if (len_trim(output_name) == 0) output_name = default_eqmatrix_filename

        open(newunit=unit_id, file=trim(output_name), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: cannot open EQMATRIX output file', trim(output_name)
            return
        end if

        write(unit_id,*) cached_nop, cached_npos
        do op = 1, cached_nop
            write(unit_id,*) cached_eqmatrix(op,1:cached_npos)
        end do
        close(unit_id)
        write(*,'(A,1X,A)') 'symmetry_write_matrix: wrote', trim(output_name)
    end subroutine symmetry_write_matrix

    subroutine symmetry_finalize()
        if (allocated(cached_eqmatrix)) deallocate(cached_eqmatrix)
        cached_nop = 0
        cached_npos = 0
        symmetry_ready = .false.
    end subroutine symmetry_finalize

    subroutine symmetry_read_matrix_file(filename, out_mat, out_nop, out_npos, found, quiet_missing)
        character(len=*), intent(in) :: filename
        integer, allocatable, intent(out) :: out_mat(:,:)
        integer, intent(out) :: out_nop, out_npos
        logical, intent(out), optional :: found
        logical, intent(in), optional :: quiet_missing
        logical :: exists
        logical :: quiet
        integer :: unit_id, ios, op

        out_nop = 0
        out_npos = 0
        if (allocated(out_mat)) deallocate(out_mat)

        quiet = .false.
        if (present(quiet_missing)) quiet = quiet_missing

        inquire(file=trim(filename), exist=exists)
        if (present(found)) found = exists
        if (.not. exists) then
            if (.not. quiet) then
                write(error_unit,'(A,1X,A)') 'Error: cannot open EQMATRIX file', trim(filename)
            end if
            return
        end if

        open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: cannot open EQMATRIX file', trim(filename)
            stop 1
        end if

        read(unit_id,*,iostat=ios) out_nop, out_npos
        if (ios /= 0 .or. out_nop <= 0 .or. out_npos <= 0) then
            write(error_unit,'(A,1X,A)') 'Error: invalid EQMATRIX header in', trim(filename)
            close(unit_id)
            stop 1
        end if

        allocate(out_mat(out_nop, out_npos))
        do op = 1, out_nop
            read(unit_id,*,iostat=ios) out_mat(op,1:out_npos)
            if (ios /= 0) then
                write(error_unit,'(A,1X,A)') 'Error: invalid EQMATRIX row while reading', trim(filename)
                close(unit_id)
                stop 1
            end if
        end do
        close(unit_id)

        if (.not. symmetry_validate_matrix(out_mat, out_nop, out_npos, trim(filename), .true.)) then
            stop 1
        end if
    end subroutine symmetry_read_matrix_file

    subroutine symmetry_generate_matrix_from_files(out_mat, out_nop, out_npos)
        integer, allocatable, intent(out) :: out_mat(:,:)
        integer, intent(out) :: out_nop, out_npos
        integer :: nat1r, nat1_unit, nat, nop1, nop
        integer :: na, nb, nc, sptarget
        integer :: i, j, k, op, at, at0, at1, at1i, op1, t, pos, attmp
        logical :: found
        real(dp) :: prod
        real(dp) :: coordstemp(3)
        integer, allocatable :: natsp1(:), spat1r(:), spat_unit(:), spat(:)
        integer, allocatable :: pos2coord(:), coord2pos(:), fulleqmatrix(:,:)
        real(dp), allocatable :: coords1r(:,:), coords1(:,:), coords(:,:)
        real(dp), allocatable :: mgroup(:,:,:), vgroup(:,:), vt(:,:)
        type(sgo_file_data) :: sgo
        type(insod_file_data) :: insod

        out_nop = 0
        out_npos = 0
        if (allocated(out_mat)) deallocate(out_mat)

        call core_read_sgo_file(sgo)
        call core_read_insod_file(insod)
        nop1 = sgo%nop1
        na = insod%na
        nb = insod%nb
        nc = insod%nc
        sptarget = insod%sptarget

        nat1r = insod%nat0 * nop1
        allocate(coords1r(nat1r,3))
        allocate(spat1r(nat1r))
        at1 = 0
        do at0 = 1, insod%nat0
            do op1 = 1, nop1
                at1 = at1 + 1
                coords1r(at1,1:3) = matmul(sgo%mgroup1(op1,1:3,1:3), insod%coords0(at0,1:3)) + sgo%vgroup1(op1,1:3)
                call symmetry_wrap_fractional_vector(coords1r(at1,1:3))
                spat1r(at1) = insod%spat0(at0)
            end do
        end do

        allocate(coords1(nat1r,3))
        allocate(spat_unit(nat1r))
        coords1(1,1:3) = coords1r(1,1:3)
        spat_unit(1) = spat1r(1)
        at1 = 1
        do at0 = 2, nat1r
            found = .false.
            do at1i = 1, at1
                prod = dot_product(coords1r(at0,1:3) - coords1(at1i,1:3), coords1r(at0,1:3) - coords1(at1i,1:3))
                if (prod <= coord_tol) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                at1 = at1 + 1
                coords1(at1,1:3) = coords1r(at0,1:3)
                spat_unit(at1) = spat1r(at0)
            end if
        end do
        nat1_unit = at1

        allocate(natsp1(insod%nsp))
        natsp1 = 0
        do at1 = 1, nat1_unit
            natsp1(spat_unit(at1)) = natsp1(spat_unit(at1)) + 1
        end do

        call symmetry_restrict_supercell_operators(sgo%mgroup1, sgo%vgroup1, nop1, na, nb, nc)

        allocate(vt(na*nb*nc,3))
        t = 0
        do i = 0, na-1
            do j = 0, nb-1
                do k = 0, nc-1
                    t = t + 1
                    vt(t,1) = real(i,dp) / real(na,dp)
                    vt(t,2) = real(j,dp) / real(nb,dp)
                    vt(t,3) = real(k,dp) / real(nc,dp)
                end do
            end do
        end do

        nat = nat1_unit * na * nb * nc
        allocate(coords(nat,3))
        allocate(spat(nat))
        at = 0
        do at0 = 1, nat1_unit
            do t = 1, na*nb*nc
                at = at + 1
                coords(at,1) = vt(t,1) + coords1(at0,1) / real(na,dp)
                coords(at,2) = vt(t,2) + coords1(at0,2) / real(nb,dp)
                coords(at,3) = vt(t,3) + coords1(at0,3) / real(nc,dp)
                call symmetry_wrap_fractional_vector(coords(at,1:3))
                spat(at) = spat_unit(at0)
            end do
        end do

        out_npos = natsp1(sptarget) * na * nb * nc
        if (out_npos <= 0) then
            write(error_unit,'(A)') 'Error: no substitutable sites were found while generating EQMATRIX.'
            stop 1
        end if

        allocate(pos2coord(out_npos))
        pos = 0
        do at = 1, nat
            if (spat(at) == sptarget) then
                pos = pos + 1
                pos2coord(pos) = at
            end if
        end do

        out_nop = nop1 * na * nb * nc
        nop = out_nop
        allocate(mgroup(nop,3,3))
        allocate(vgroup(nop,3))
        op = 0
        do op1 = 1, nop1
            do t = 1, na*nb*nc
                op = op + 1
                mgroup(op,1:3,1:3) = sgo%mgroup1(op1,1:3,1:3)
                vgroup(op,1) = sgo%vgroup1(op1,1) / real(na,dp) + vt(t,1)
                vgroup(op,2) = sgo%vgroup1(op1,2) / real(nb,dp) + vt(t,2)
                vgroup(op,3) = sgo%vgroup1(op1,3) / real(nc,dp) + vt(t,3)
            end do
        end do

        allocate(fulleqmatrix(nop,nat))
        allocate(coord2pos(nat))
        coord2pos = 0
        do pos = 1, out_npos
            coord2pos(pos2coord(pos)) = pos
        end do

        do op = 1, nop
            do at = 1, nat
                coordstemp = matmul(mgroup(op,1:3,1:3), coords(at,1:3)) + vgroup(op,1:3)
                call symmetry_wrap_fractional_vector(coordstemp)
                found = .false.
                do attmp = 1, nat
                    if (abs(coordstemp(1) - coords(attmp,1)) < coord_tol .and. &
                        abs(coordstemp(2) - coords(attmp,2)) < coord_tol .and. &
                        abs(coordstemp(3) - coords(attmp,3)) < coord_tol) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    write(error_unit,'(A,I0,A,I0,A)') 'Error: operator ', op, ' applied to atom ', at, &
                        ' does not map onto a valid supercell atom.'
                    stop 1
                end if
                fulleqmatrix(op,at) = attmp
            end do
        end do

        allocate(out_mat(out_nop, out_npos))
        do op = 1, out_nop
            do pos = 1, out_npos
                out_mat(op,pos) = coord2pos(fulleqmatrix(op,pos2coord(pos)))
            end do
        end do

        if (.not. symmetry_validate_matrix(out_mat, out_nop, out_npos, 'generated EQMATRIX', .true.)) then
            stop 1
        end if
    end subroutine symmetry_generate_matrix_from_files

    logical function symmetry_validate_matrix(matrix, nop, npos, label, report) result(ok)
        integer, intent(in) :: matrix(:,:)
        integer, intent(in) :: nop, npos
        character(len=*), intent(in), optional :: label
        logical, intent(in), optional :: report
        logical :: verbose
        integer :: op, pos, value
        logical, allocatable :: seen(:)
        character(len=128) :: target_label

        verbose = .false.
        if (present(report)) verbose = report

        target_label = 'EQMATRIX'
        if (present(label)) target_label = label

        ok = .false.
        if (nop <= 0 .or. npos <= 0) then
            if (verbose) write(error_unit,'(A,1X,A)') 'Error: invalid dimensions while validating', trim(target_label)
            return
        end if
        if (size(matrix,1) < nop .or. size(matrix,2) < npos) then
            if (verbose) then
                write(error_unit,'(A,1X,A)') 'Error: array bounds are smaller than EQMATRIX dimensions for', &
                    trim(target_label)
            end if
            return
        end if

        allocate(seen(npos))
        do op = 1, nop
            seen = .false.
            do pos = 1, npos
                value = matrix(op,pos)
                if (value < 1 .or. value > npos) then
                    if (verbose) then
                        write(error_unit,'(A,1X,A)') 'Error: out-of-range EQMATRIX entry in', trim(target_label)
                        write(error_unit,'(A,I0,A,I0,A,I0)') '  operator=', op, ' position=', pos, ' value=', value
                    end if
                    deallocate(seen)
                    return
                end if
                if (seen(value)) then
                    if (verbose) then
                        write(error_unit,'(A,1X,A)') 'Error: duplicated target-site mapping in', trim(target_label)
                        write(error_unit,'(A,I0,A,I0,A,I0)') '  operator=', op, ' position=', pos, ' repeated value=', value
                    end if
                    deallocate(seen)
                    return
                end if
                seen(value) = .true.
            end do
        end do
        deallocate(seen)
        ok = .true.
    end function symmetry_validate_matrix

    subroutine symmetry_restrict_supercell_operators(mgroup1, vgroup1, nop1, na, nb, nc)
        real(dp), intent(inout) :: mgroup1(:,:,:)
        real(dp), intent(inout) :: vgroup1(:,:)
        integer, intent(inout) :: nop1
        integer, intent(in) :: na, nb, nc
        integer :: op_in, op_out

        op_out = 0
        do op_in = 1, nop1
            if (supercell_operator_allowed(mgroup1(op_in,1:3,1:3), na, nb, nc)) then
                op_out = op_out + 1
                if (op_out /= op_in) then
                    mgroup1(op_out,1:3,1:3) = mgroup1(op_in,1:3,1:3)
                    vgroup1(op_out,1:3) = vgroup1(op_in,1:3)
                end if
            end if
        end do

        nop1 = op_out
        if (nop1 <= 0) then
            write(error_unit,'(A)') 'Error: no symmetry operators remain after applying the INSOD supercell constraints.'
            stop 1
        end if
    end subroutine symmetry_restrict_supercell_operators

    pure subroutine symmetry_wrap_fractional_vector(values)
        real(dp), intent(inout) :: values(:)
        integer :: idx

        do idx = 1, size(values)
            values(idx) = wrap_fractional_scalar(values(idx))
        end do
    end subroutine symmetry_wrap_fractional_vector

    subroutine cache_matrix(tmp, nop, npos)
        integer, allocatable, intent(inout) :: tmp(:,:)
        integer, intent(in) :: nop, npos

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
    end subroutine cache_matrix

    pure real(dp) function wrap_fractional_scalar(value) result(wrapped)
        real(dp), intent(in) :: value

        wrapped = value - floor(value)
        if (wrapped > (1.0_dp - wrap_tol)) wrapped = 0.0_dp
        if (wrapped < 0.0_dp .and. abs(wrapped) <= wrap_tol) wrapped = 0.0_dp
    end function wrap_fractional_scalar

    logical function supercell_operator_allowed(rotation, na, nb, nc) result(allowed)
        real(dp), intent(in) :: rotation(3,3)
        integer, intent(in) :: na, nb, nc

        allowed = .true.
        if (na == nb .and. nb == nc) return

        if (na == nb) then
            allowed = abs(rotation(1,3)) <= symop_tol .and. abs(rotation(3,1)) <= symop_tol .and. &
                      abs(rotation(2,3)) <= symop_tol .and. abs(rotation(3,2)) <= symop_tol
        else if (na == nc) then
            allowed = abs(rotation(1,2)) <= symop_tol .and. abs(rotation(2,1)) <= symop_tol .and. &
                      abs(rotation(3,2)) <= symop_tol .and. abs(rotation(2,3)) <= symop_tol
        else if (nb == nc) then
            allowed = abs(rotation(2,1)) <= symop_tol .and. abs(rotation(1,2)) <= symop_tol .and. &
                      abs(rotation(3,1)) <= symop_tol .and. abs(rotation(1,3)) <= symop_tol
        else
            allowed = abs(rotation(2,1)) <= symop_tol .and. abs(rotation(1,2)) <= symop_tol .and. &
                      abs(rotation(3,1)) <= symop_tol .and. abs(rotation(1,3)) <= symop_tol .and. &
                      abs(rotation(3,2)) <= symop_tol .and. abs(rotation(2,3)) <= symop_tol
        end if
    end function supercell_operator_allowed

end module symmetry
