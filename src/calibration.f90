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
!******************************************************************************

!*******************************************************************************
! Common calibration routines shared between SOD ensemble drivers.
! Rutinas comunes de calibración compartidas entre los controladores ensemble de SOD.
!*******************************************************************************
module calibration
    use consts
    use utils
    use settings, only: set_blend_override
    use cli, only: basename_from_path, to_lower_inplace
    use configurations, only: canonicalize_subset, find_subset_index
    use energy_calculations
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    private
    public :: calibrate_level_with_engine, calibrate_level_with_gulp, canonicalize_subset, find_subset_index
    public :: evaluate_subsets_with_engine, evaluate_subsets_with_gulp
    public :: set_calibration_engine_from_filer, set_external_calibration_enabled, external_calibration_is_enabled, calibration_engine_name
    public :: set_calibration_template_gin, set_calibration_protocol_version
    public :: reset_mc_gulp_energy_snapshot
    public :: update_level_blend_override_from_samples

    integer, parameter :: calibration_engine_none = 0
    integer, parameter :: calibration_engine_gulp = 1
    integer, parameter :: calibration_engine_ase  = 2
    logical, save :: external_calibration_enabled = .false.
    integer, parameter :: template_mode_builtin = 0
    integer, parameter :: template_mode_none    = 1
    integer, parameter :: template_mode_custom  = 2
    integer, save :: calibration_engine = calibration_engine_gulp
    integer, save :: template_mode = template_mode_none
    character(len=512), save :: template_gin_override = ''
    character(len=8), save :: calibration_protocol_version = '2.0'

contains

    pure function resolve_level_prefix(level_prefix, default_prefix) result(prefix)
        character(len=*), intent(in), optional :: level_prefix
        character(len=*), intent(in) :: default_prefix
        character(len=32) :: prefix

        prefix = trim(default_prefix)
        if (present(level_prefix)) then
            if (len_trim(level_prefix) > 0) prefix = trim(level_prefix)
        end if
    end function resolve_level_prefix

    pure function to_lower(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: idx, code

        out = text
        do idx = 1, len(text)
            code = iachar(out(idx:idx))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(idx:idx) = achar(code + 32)
            end if
        end do
    end function to_lower

    subroutine set_calibration_engine_from_filer(filer)
        integer, intent(in) :: filer

        select case (filer)
        case (1)
            calibration_engine = calibration_engine_gulp
        case (14)
            calibration_engine = calibration_engine_ase
        case default
            calibration_engine = calibration_engine_none
        end select
    end subroutine set_calibration_engine_from_filer

    ! Subroutine `set_external_calibration_enabled` stores whether external calibration is allowed.
    ! La subrutina `set_external_calibration_enabled` guarda si la calibración externa está permitida.
    subroutine set_external_calibration_enabled(enabled)
        logical, intent(in) :: enabled

        external_calibration_enabled = enabled
    end subroutine set_external_calibration_enabled

    ! Function `external_calibration_is_enabled` reports whether external calibration is currently enabled.
    ! La función `external_calibration_is_enabled` informa de si la calibración externa está habilitada actualmente.
    logical function external_calibration_is_enabled()
        external_calibration_is_enabled = external_calibration_enabled
    end function external_calibration_is_enabled

    function calibration_engine_name() result(name)
        character(len=16) :: name

        name = calibration_engine_name_from_value(calibration_engine)
    end function calibration_engine_name

    pure function calibration_engine_name_from_value(engine) result(name)
        integer, intent(in) :: engine
        character(len=16) :: name

        select case (engine)
        case (calibration_engine_gulp)
            name = 'GULP'
        case (calibration_engine_ase)
            name = 'ASE'
        case default
            name = 'none'
        end select
    end function calibration_engine_name_from_value

    pure function calibration_engine_stub(engine) result(name)
        integer, intent(in) :: engine
        character(len=8) :: name

        select case (engine)
        case (calibration_engine_gulp)
            name = 'gulp'
        case (calibration_engine_ase)
            name = 'ase'
        case default
            name = 'engine'
        end select
    end function calibration_engine_stub

    integer pure function calibration_engine_filer_value(engine)
        integer, intent(in) :: engine

        select case (engine)
        case (calibration_engine_gulp)
            calibration_engine_filer_value = 1
        case (calibration_engine_ase)
            calibration_engine_filer_value = 14
        case default
            calibration_engine_filer_value = 0
        end select
    end function calibration_engine_filer_value

    pure function expected_energy_output_name(base_name, engine) result(name)
        character(len=*), intent(in) :: base_name
        integer, intent(in) :: engine
        character(len=256) :: name

        select case (engine)
        case (calibration_engine_ase)
            name = trim(base_name)//'.out.ase'
        case default
            name = trim(base_name)//'.vasp.gout'
        end select
    end function expected_energy_output_name

    function build_engine_run_command(run_dir, engine) result(command)
        character(len=*), intent(in) :: run_dir
        integer, intent(in) :: engine
        character(len=1024) :: command

        select case (engine)
        case (calibration_engine_gulp)
            command = 'cd ' // trim(run_dir) // ' && SOD_GULP_PROTOCOL_VERSION=' // trim(calibration_protocol_version) // ' bash run_jobs.sh'
        case default
            command = 'cd ' // trim(run_dir) // ' && bash run_jobs.sh'
        end select
    end function build_engine_run_command

    ! Removes the per-level partial GULP energy snapshot generated during MC runs.
! Elimina la instantánea parcial de energías GULP por nivel generada durante las ejecuciones MC.
    subroutine reset_mc_gulp_energy_snapshot(level)
        integer, intent(in) :: level
        character(len=512) :: snapshot_path
        character(len=32) :: level_dir
        logical :: exists
        integer :: unit_file, ios

        level_dir = format_level_directory('mc', level)
        snapshot_path = join_path(trim(level_dir), 'ENERGIES')
        inquire(file=trim(snapshot_path), exist=exists)
        if (.not. exists) return

        open(newunit=unit_file, file=trim(snapshot_path), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        close(unit_file, status='delete')
    end subroutine reset_mc_gulp_energy_snapshot

    ! Appends the ENERGIES file produced in a temporary GULP folder to the MC snapshot.
! Añade al volcado MC el archivo ENERGIES producido en una carpeta temporal de GULP.
    subroutine append_mc_gulp_energy_snapshot(level, source_dir, level_prefix)
        integer, intent(in) :: level
        character(len=*), intent(in) :: source_dir
        character(len=*), intent(in), optional :: level_prefix
        character(len=512) :: snapshot_path, source_path
        character(len=32) :: level_dir
        character(len=32) :: prefix
        character(len=2048) :: line
        integer :: unit_in, unit_out, ios
        logical :: exists

        prefix = resolve_level_prefix(level_prefix, 'mc')
        level_dir = format_level_directory(trim(prefix), level)
        call ensure_directory_exists(trim(level_dir))
        snapshot_path = join_path(trim(level_dir), 'ENERGIES')
        source_path = trim(source_dir)//'/ENERGIES'
        inquire(file=trim(source_path), exist=exists)
        if (.not. exists) return

        open(newunit=unit_in, file=trim(source_path), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        open(newunit=unit_out, file=trim(snapshot_path), status='unknown', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            close(unit_in)
            return
        end if

        do
            read(unit_in,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            write(unit_out,'(A)') trim(line)//' # SOURCE '//trim(source_dir)
        end do

        close(unit_out)
        close(unit_in)
    end subroutine append_mc_gulp_energy_snapshot

    subroutine set_calibration_template_gin(path)
        character(len=*), intent(in) :: path
        character(len=512) :: cleaned
        character(len=512) :: lowered
        logical :: exists

        cleaned = adjustl(path)
        if (len_trim(cleaned) == 0) then
            template_mode = template_mode_none
            template_gin_override = ''
            return
        end if

        lowered = to_lower(trim(cleaned))
        if (trim(lowered) == 'default' .or. trim(lowered) == 'builtin') then
            template_mode = template_mode_builtin
            template_gin_override = ''
        else if (trim(lowered) == 'none' .or. trim(lowered) == 'off' .or. trim(lowered) == 'skip') then
            template_mode = template_mode_none
            template_gin_override = ''
        else
            template_mode = template_mode_custom
            template_gin_override = trim(cleaned)
            inquire(file=trim(template_gin_override), exist=exists)
            if (.not. exists) then
                write(error_unit,'(A)') 'Error: template file not found: '//trim(template_gin_override)
                stop 1
            end if
        end if
    end subroutine set_calibration_template_gin

    subroutine set_calibration_protocol_version(raw_value)
        character(len=*), intent(in) :: raw_value
        character(len=64) :: lowered

        lowered = adjustl(trim(raw_value))
        call to_lower_inplace(lowered)

        select case (trim(lowered))
        case ('', '2', '2.0')
            calibration_protocol_version = '2.0'
        case ('1', '1.0')
            calibration_protocol_version = '1.0'
        case default
            write(error_unit,'(A)') 'Error: invalid protocol version. Use 1.0 or 2.0.'
            stop 1
        end select
    end subroutine set_calibration_protocol_version

    subroutine update_level_blend_override_from_samples(level, total_sites, low_values, high_values, target_values, sample_count, source_label, &
        fitted_lambda_high, fitted_rms_residual, fitted_point_count, level_prefix)
        integer, intent(in) :: level, total_sites, sample_count
        real(dp), intent(in) :: low_values(:), high_values(:), target_values(:)
        character(len=*), intent(in) :: source_label
        character(len=*), intent(in), optional :: level_prefix
        real(dp), intent(out), optional :: fitted_lambda_high, fitted_rms_residual
        integer, intent(out), optional :: fitted_point_count
        real(dp), parameter :: limit = huge(1.0_dp) * 0.5_dp
        real(dp), parameter :: eps = 1.0e-12_dp
        real(dp) :: sum_num, sum_den, lambda_high, delta_val
        real(dp) :: residual, rms_residual, weight_sum, predicted
        character(len=32) :: output_prefix
        integer :: i, valid_count
        logical :: valid_low, valid_high, valid_target

        if (present(fitted_lambda_high)) fitted_lambda_high = -1.0_dp
        if (present(fitted_rms_residual)) fitted_rms_residual = -1.0_dp
        if (present(fitted_point_count)) fitted_point_count = 0
        if (sample_count <= 0) return

        sum_num = 0.0_dp
        sum_den = 0.0_dp
        valid_count = 0

        do i = 1, sample_count
            valid_low = abs(low_values(i)) < limit .and. low_values(i) == low_values(i)
            valid_high = abs(high_values(i)) < limit .and. high_values(i) == high_values(i)
            valid_target = abs(target_values(i)) < limit .and. target_values(i) == target_values(i)
            if (.not. (valid_low .and. valid_high .and. valid_target)) cycle

            delta_val = high_values(i) - low_values(i)
            sum_num = sum_num + delta_val * (target_values(i) - low_values(i))
            sum_den = sum_den + delta_val * delta_val
            valid_count = valid_count + 1
        end do

        if (valid_count < 2 .or. sum_den <= eps) return

        lambda_high = sum_num / sum_den
        lambda_high = max(0.0_dp, min(1.0_dp, lambda_high))

        rms_residual = 0.0_dp
        weight_sum = 0.0_dp
        do i = 1, sample_count
            valid_low = abs(low_values(i)) < limit .and. low_values(i) == low_values(i)
            valid_high = abs(high_values(i)) < limit .and. high_values(i) == high_values(i)
            valid_target = abs(target_values(i)) < limit .and. target_values(i) == target_values(i)
            if (.not. (valid_low .and. valid_high .and. valid_target)) cycle

            predicted = (1.0_dp - lambda_high) * low_values(i) + lambda_high * high_values(i)
            residual = predicted - target_values(i)
            rms_residual = rms_residual + residual * residual
            weight_sum = weight_sum + 1.0_dp
        end do
        if (weight_sum > 0.0_dp) rms_residual = sqrt(rms_residual / weight_sum)

        call set_blend_override(level, total_sites, lambda_high)
        output_prefix = resolve_level_prefix(level_prefix, 'mc')
        call append_blend_override_summary(level, total_sites, output_prefix, trim(source_label), lambda_high, rms_residual, valid_count)
        write(*,'(A,I0,A,F8.5,A,F12.6,A,I0,A)') 'Level ', level, ': fitted Y-side blend weight = ', &
            lambda_high, ' (RMS residual ', rms_residual, ' eV) from ', valid_count, ' '//trim(source_label)//' points'
        flush(output_unit)
        if (present(fitted_lambda_high)) fitted_lambda_high = lambda_high
        if (present(fitted_rms_residual)) fitted_rms_residual = rms_residual
        if (present(fitted_point_count)) fitted_point_count = valid_count
    end subroutine update_level_blend_override_from_samples

    ! Subroutine `calibrate_level_with_engine` dispatches calibration to the active external backend.
    ! La subrutina `calibrate_level_with_engine` delega la calibración al backend externo activo.
    subroutine calibrate_level_with_engine(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, unique_low, &
        unique_high, config, do_low_calibration, do_high_calibration, level_prefix)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_subsets(:,:)
        real(dp), intent(in) :: unique_low_contrib(:,:), unique_high_contrib(:,:)
        real(dp), intent(inout) :: unique_low(:)
        real(dp), intent(inout) :: unique_high(:)
        integer, intent(inout) :: config(:)
        logical, intent(in) :: do_low_calibration, do_high_calibration
        character(len=*), intent(in), optional :: level_prefix

        if (.not. external_calibration_enabled) then
            return
        end if

        if (calibration_engine == calibration_engine_none) then
            write(*,'(A,I0,A)') 'Level ', level, ': no external calibration engine selected; keeping internal Hamiltonian.'
            flush(output_unit)
            return
        end if

        call calibrate_level_with_backend(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, &
            unique_low, unique_high, config, do_low_calibration, do_high_calibration, calibration_engine, level_prefix)
    end subroutine calibrate_level_with_engine

    ! Subroutine `calibrate_level_with_gulp` keeps the legacy GULP-only calibration entry point.
    ! La subrutina `calibrate_level_with_gulp` mantiene el punto de entrada heredado de calibración solo con GULP.
    subroutine calibrate_level_with_gulp(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, unique_low, &
        unique_high, config, do_low_calibration, do_high_calibration, level_prefix)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_subsets(:,:)
        real(dp), intent(in) :: unique_low_contrib(:,:), unique_high_contrib(:,:)
        real(dp), intent(inout) :: unique_low(:)
        real(dp), intent(inout) :: unique_high(:)
        integer, intent(inout) :: config(:)
        logical, intent(in) :: do_low_calibration, do_high_calibration
        character(len=*), intent(in), optional :: level_prefix

        call calibrate_level_with_backend(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, &
            unique_low, unique_high, config, do_low_calibration, do_high_calibration, calibration_engine_gulp, level_prefix)
    end subroutine calibrate_level_with_gulp

    ! Subroutine `calibrate_level_with_backend` implements the shared external calibration workflow.
    subroutine calibrate_level_with_backend(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, unique_low, unique_high, config, do_low_calibration, do_high_calibration, engine, level_prefix)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count, engine
        integer, intent(in) :: unique_subsets(:,:)
        real(dp), intent(in) :: unique_low_contrib(:,:), unique_high_contrib(:,:)
        real(dp), intent(inout) :: unique_low(:)
        real(dp), intent(inout) :: unique_high(:)
        integer, intent(inout) :: config(:)
        logical, intent(in) :: do_low_calibration, do_high_calibration
        character(len=*), intent(in), optional :: level_prefix

        if (engine == calibration_engine_none) then
            write(*,'(A,I0,A)') 'Level ', level, ': no external calibration engine selected; keeping internal Hamiltonian.'
            flush(output_unit)
            return
        end if

        call calibrate_level_with_backend_impl(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, &
            unique_low, unique_high, config, do_low_calibration, do_high_calibration, engine, level_prefix)
    end subroutine calibrate_level_with_backend

    ! Subroutine `calibrate_level_with_backend_impl` performs the actual backend-specific external calibration.
    ! La subrutina `calibrate_level_with_backend_impl` realiza la calibración externa específica del backend.
    subroutine calibrate_level_with_backend_impl(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, unique_low, unique_high, config, do_low_calibration, do_high_calibration, engine, level_prefix)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count, engine
        integer, intent(in) :: unique_subsets(:,:)
        real(dp), intent(in) :: unique_low_contrib(:,:), unique_high_contrib(:,:)
        real(dp), intent(inout) :: unique_low(:)
        real(dp), intent(inout) :: unique_high(:)
        integer, intent(inout) :: config(:)
        logical, intent(in) :: do_low_calibration, do_high_calibration
        character(len=*), intent(in), optional :: level_prefix

        integer, parameter :: max_orders = 4
        integer, parameter :: max_targets = max_orders + 1
        integer, parameter :: max_union = 2 * max_targets
        real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
        character(len=512) :: calib_dir, script_dir
        character(len=64) :: base_name
        character(len=16) :: level_tag
        character(len=32) :: output_prefix
        character(len=128) :: gout_names(max_union)
        integer :: selected_low_idx(max_targets)
        integer :: selected_high_idx(max_targets)
        integer :: union_idx(max_union)
        real(dp) :: measured_energy(max_union)
        logical :: filled_union(max_union)
        integer :: i, j, ios, unit_file, exit_code
        integer :: union_count
        integer :: n_valid_high
        integer :: idx_global
        integer, allocatable :: valid_high_idx(:), tmp_selected(:)
        real(dp), allocatable :: tmp_values(:)
        logical :: ok_scripts, ok_low_selection, ok_high_selection, calibrate_high, ok_system, already
        integer :: low_order_count, high_order_count
        integer :: low_target_count, high_target_count
        real(dp) :: base_energy_low, base_energy_high
        real(dp) :: design_low(max_union, max_orders), target_low(max_union), coeff_low(max_orders)
        real(dp) :: design_high(max_union, max_orders), target_high(max_union), coeff_high(max_orders)
        real(dp) :: energy_val, corrected
        real(dp) :: fitted_energy, residual, rms_low, rms_high, max_res_low, max_res_high
        integer :: fit_count_low, fit_count_high
        character(len=1024) :: command, line
        character(len=256) :: gout_label
        character(len=256) :: engine_label
        integer :: pos_hash, match_idx
        real(dp), allocatable :: blend_low(:), blend_high(:), blend_target(:)
        real(dp) :: blend_lambda_high, blend_rms
        integer :: blend_points

        if (.not. (do_low_calibration .and. do_high_calibration)) return

        low_order_count = max(1, min(max_orders, get_max_low_order()))
        high_order_count = max(1, min(max_orders, get_max_high_order()))
        low_target_count = low_order_count + 1
        high_target_count = high_order_count + 1
        if (unique_count < low_target_count) return

        selected_low_idx = 0
        selected_high_idx = 0

        engine_label = calibration_engine_name_from_value(engine)

        ok_scripts = find_scripts_directory(script_dir, engine)
        if (.not. ok_scripts) then
            write(*,'(A,A,A,I0)') 'Warning: ', trim(engine_label), ' scripts directory not found for level ', level
            return
        end if

        ok_low_selection = .false.
        if (do_low_calibration) then
            call select_calibration_indices(unique_low, unique_count, selected_low_idx(1:low_target_count), ok_low_selection)
            if (.not. ok_low_selection) return
        end if

        allocate(valid_high_idx(unique_count))
        n_valid_high = 0
        do i = 1, unique_count
            if (abs(unique_high(i)) < huge_marker) then
                n_valid_high = n_valid_high + 1
                valid_high_idx(n_valid_high) = i
            end if
        end do

        calibrate_high = .false.
        ok_high_selection = .false.
        if (do_high_calibration .and. n_valid_high >= high_target_count) then
            allocate(tmp_values(n_valid_high))
            allocate(tmp_selected(high_target_count))
            do i = 1, n_valid_high
                tmp_values(i) = unique_high(valid_high_idx(i))
            end do
            call select_calibration_indices(tmp_values, n_valid_high, tmp_selected(1:high_target_count), ok_high_selection)
            if (ok_high_selection) then
                calibrate_high = .true.
                do i = 1, high_target_count
                    selected_high_idx(i) = valid_high_idx(tmp_selected(i))
                end do
            end if
            deallocate(tmp_values)
            deallocate(tmp_selected)
        end if
        if (allocated(valid_high_idx)) deallocate(valid_high_idx)

        output_prefix = resolve_level_prefix(level_prefix, 'mc')
        write(level_tag,'(I4.4)') level
        calib_dir = join_path(trim(format_level_directory(trim(output_prefix), level)), trim(calibration_engine_stub(engine))//'_calib_N'//trim(adjustl(level_tag)))
        if (engine == calibration_engine_ase) then
            call ensure_ase_level_resources(format_level_directory(trim(output_prefix), level))
        end if

        union_idx = 0
        union_count = 0
        if (do_low_calibration) then
            do i = 1, low_target_count
                idx_global = selected_low_idx(i)
                if (idx_global <= 0) cycle
                already = .false.
                do j = 1, union_count
                    if (union_idx(j) == idx_global) then
                        already = .true.
                        exit
                    end if
                end do
                if (.not. already) then
                    if (union_count >= max_union) then
                        write(*,'(A,I0)') 'Warning: insufficient capacity for calibration union at level ', level
                        return
                    end if
                    union_count = union_count + 1
                    union_idx(union_count) = idx_global
                end if
            end do
        end if
        if (calibrate_high) then
            do i = 1, high_target_count
                idx_global = selected_high_idx(i)
                if (idx_global <= 0) cycle
                already = .false.
                do j = 1, union_count
                    if (union_idx(j) == idx_global) then
                        already = .true.
                        exit
                    end if
                end do
                if (.not. already) then
                    if (union_count >= max_union) then
                        write(*,'(A,I0)') 'Warning: insufficient capacity for calibration union at level ', level
                        return
                    end if
                    union_count = union_count + 1
                    union_idx(union_count) = idx_global
                end if
            end do
        end if

        if (union_count == 0) return

        call execute_command_line('rm -rf ' // trim(calib_dir), exitstat=exit_code)
        call execute_command_line('mkdir -p ' // trim(calib_dir), exitstat=exit_code)

        gout_names = ''
        measured_energy = 0.0_dp
        filled_union = .false.

        do i = 1, union_count
            config = 1
            if (level > 0) config(unique_subsets(1:level, union_idx(i))) = 2
            write(base_name,'("calib_N",I4.4,"_c",I2.2)') level, i
            call write_vasp_file(config, total_sites, trim(calib_dir)//'/'//trim(base_name)//'.vasp')
            gout_names(i) = trim(expected_energy_output_name(base_name, engine))
        end do

        call copy_calibration_scripts(script_dir, calib_dir, engine)
        call write_engine_filer_marker(calib_dir, engine)

        command = build_engine_run_command(calib_dir, engine)
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A,A,A,I0,A)') 'Warning: run_jobs.sh reported issues during ', trim(engine_label), ' calibration for level ', &
                level, '; attempting to extract partial energies.'
        end if

        command = 'cd ' // trim(calib_dir) // ' && bash extract.sh'
        call execute_command_line(trim(command), exitstat=exit_code)
        call append_mc_gulp_energy_snapshot(level, calib_dir, trim(output_prefix))
        if (engine == calibration_engine_gulp) call remove_protocol_template_file(calib_dir)
        if (exit_code /= 0) then
            write(*,'(A,A,A,I0)') 'Warning: extract.sh failed during ', trim(engine_label), ' calibration for level ', level
        end if

        open(newunit=unit_file, file=trim(calib_dir)//'/ENERGIES', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,I0)') 'Warning: failed to read ENERGIES for level ', level
            return
        end if

        do
            read(unit_file,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            pos_hash = index(line, '#')
            if (pos_hash <= 0) cycle
            read(line(1:pos_hash-1),*,iostat=ios) energy_val
            if (ios /= 0) cycle
            gout_label = adjustl(trim(line(pos_hash+1:)))
            gout_label = trim(gout_label)
            pos_hash = index(gout_label, '#')
            if (pos_hash > 0) gout_label = adjustl(trim(gout_label(:pos_hash-1)))
            match_idx = 0
            do i = 1, union_count
                if (trim(gout_label) == trim(gout_names(i))) then
                    match_idx = i
                    exit
                end if
            end do
            if (match_idx > 0) then
                measured_energy(match_idx) = energy_val
                filled_union(match_idx) = .true.
            end if
        end do
        close(unit_file)

        do i = 1, union_count
            if (.not. filled_union(i)) then
                write(*,'(A,I0)') 'Warning: incomplete energies in calibration output for level ', level
                return
            end if
        end do

        if (do_low_calibration) then
            base_energy_low = get_base_energy()
            design_low = 0.0_dp
            target_low = 0.0_dp
            fit_count_low = 0
            do i = 1, union_count
                idx_global = union_idx(i)
                if (idx_global <= 0) cycle
                energy_val = measured_energy(i)
                fit_count_low = fit_count_low + 1
                do j = 1, low_order_count
                    design_low(fit_count_low, j) = unique_low_contrib(j, idx_global)
                end do
                target_low(fit_count_low) = energy_val - base_energy_low
            end do

            coeff_low = 0.0_dp
            call solve_least_squares_qr(design_low(1:fit_count_low,1:low_order_count), target_low(1:fit_count_low), &
                coeff_low(1:low_order_count), ok_system)
            if (.not. ok_system) then
                write(*,'(A,I0)') 'Warning: QR least-squares failed in calibration (X side) for level ', level
                return
            end if

            rms_low = 0.0_dp
            max_res_low = 0.0_dp
            if (fit_count_low > 0) then
                do i = 1, union_count
                    idx_global = union_idx(i)
                    if (idx_global <= 0) cycle
                    fitted_energy = base_energy_low + sum(coeff_low(1:low_order_count) * unique_low_contrib(1:low_order_count, idx_global))
                    residual = fitted_energy - measured_energy(i)
                    rms_low = rms_low + residual * residual
                    max_res_low = max(max_res_low, abs(residual))
                end do
                rms_low = sqrt(rms_low / real(fit_count_low, dp))
                    write(*,'(A,I0,A,I0,A,F12.6,A,F12.6)') 'Level ', level, ': X-side calibration uses ', fit_count_low, &
                    ' '//trim(engine_label)//' points; RMS residual = ', rms_low, ' eV, max residual = ', max_res_low
                call print_fit_coefficients(level, 'X', low_order_count, coeff_low)
            end if
        end if

        coeff_high = 0.0_dp
        blend_lambda_high = -1.0_dp
        blend_rms = -1.0_dp
        blend_points = 0
        if (calibrate_high) then
            base_energy_high = get_high_base_energy()
            design_high = 0.0_dp
            target_high = 0.0_dp
            fit_count_high = 0
            do i = 1, union_count
                idx_global = union_idx(i)
                if (idx_global <= 0) cycle
                energy_val = measured_energy(i)
                fit_count_high = fit_count_high + 1
                do j = 1, high_order_count
                    design_high(fit_count_high, j) = unique_high_contrib(j, idx_global)
                end do
                target_high(fit_count_high) = energy_val - base_energy_high
            end do

            call solve_least_squares_qr(design_high(1:fit_count_high,1:high_order_count), target_high(1:fit_count_high), &
                coeff_high(1:high_order_count), ok_system)
            if (.not. ok_system) then
                write(*,'(A,I0)') 'Warning: QR least-squares failed in calibration (Y side) for level ', level
                calibrate_high = .false.
            else
                rms_high = 0.0_dp
                max_res_high = 0.0_dp
                if (fit_count_high > 0) then
                    do i = 1, union_count
                        idx_global = union_idx(i)
                        if (idx_global <= 0) cycle
                        fitted_energy = base_energy_high + sum(coeff_high(1:high_order_count) * unique_high_contrib(1:high_order_count, idx_global))
                        residual = fitted_energy - measured_energy(i)
                        rms_high = rms_high + residual * residual
                        max_res_high = max(max_res_high, abs(residual))
                    end do
                    rms_high = sqrt(rms_high / real(fit_count_high, dp))
                    write(*,'(A,I0,A,I0,A,F12.6,A,F12.6)') 'Level ', level, ': Y-side calibration uses ', fit_count_high, &
                        ' '//trim(engine_label)//' points; RMS residual = ', rms_high, ' eV, max residual = ', max_res_high
                    call print_fit_coefficients(level, 'Y', high_order_count, coeff_high)
                end if
            end if
        end if

        if (do_low_calibration) then
            do i = 1, unique_count
                corrected = base_energy_low + sum(coeff_low(1:low_order_count) * unique_low_contrib(1:low_order_count, i))
                unique_low(i) = corrected
            end do
        end if

        if (calibrate_high) then
            do i = 1, unique_count
                corrected = base_energy_high + sum(coeff_high(1:high_order_count) * unique_high_contrib(1:high_order_count, i))
                unique_high(i) = corrected
            end do
        end if

        if (do_low_calibration .and. calibrate_high) then
            allocate(blend_low(union_count))
            allocate(blend_high(union_count))
            allocate(blend_target(union_count))
            do i = 1, union_count
                idx_global = union_idx(i)
                blend_low(i) = unique_low(idx_global)
                blend_high(i) = unique_high(idx_global)
                blend_target(i) = measured_energy(i)
            end do
            call update_level_blend_override_from_samples(level, total_sites, blend_low, blend_high, blend_target, &
                union_count, 'calibration', blend_lambda_high, blend_rms, blend_points, trim(output_prefix))
            deallocate(blend_low)
            deallocate(blend_high)
            deallocate(blend_target)
        end if

        call write_calibration_coefficients_report(level, total_sites, trim(output_prefix), low_order_count, coeff_low, fit_count_low, rms_low, &
            max_res_low, calibrate_high, high_order_count, coeff_high, fit_count_high, rms_high, max_res_high, &
            blend_lambda_high, blend_rms, blend_points, calib_dir)

        do i = 1, union_count
            idx_global = union_idx(i)
            energy_val = measured_energy(i)
            unique_low(idx_global) = energy_val
            unique_high(idx_global) = energy_val
        end do
    end subroutine calibrate_level_with_backend_impl

    subroutine print_fit_coefficients(level, side_label, order_count, coeff)
        integer, intent(in) :: level, order_count
        character(len=*), intent(in) :: side_label
        real(dp), intent(in) :: coeff(:)
        character(len=512) :: line
        character(len=48) :: value
        integer :: i

        line = 'Level '
        write(value,'(I0)') level
        line = trim(line)//trim(adjustl(value))//': '//trim(side_label)//'-side QR coefficients:'
        do i = 1, order_count
            write(value,'(" c",I0,"=",F14.8)') i, coeff(i)
            line = trim(line)//trim(value)
        end do
        write(*,'(A)') trim(line)
        flush(output_unit)
    end subroutine print_fit_coefficients

    subroutine write_calibration_coefficients_report(level, total_sites, level_prefix, low_order_count, coeff_low, fit_count_low, rms_low, max_res_low, &
        have_high, high_order_count, coeff_high, fit_count_high, rms_high, max_res_high, blend_lambda_high, blend_rms, blend_points, calib_dir)
        integer, intent(in) :: level, total_sites, low_order_count, fit_count_low, high_order_count, fit_count_high, blend_points
        character(len=*), intent(in) :: level_prefix, calib_dir
        logical, intent(in) :: have_high
        real(dp), intent(in) :: coeff_low(:), coeff_high(:)
        real(dp), intent(in) :: rms_low, max_res_low, rms_high, max_res_high
        real(dp), intent(in) :: blend_lambda_high, blend_rms
        character(len=32) :: level_dir
        character(len=512) :: report_path
        integer :: unit_report, ios, i

        level_dir = format_level_directory(trim(level_prefix), level)
        call ensure_directory_exists(trim(level_dir))
        report_path = join_path(trim(level_dir), 'CALIBRATION_COEFFICIENTS.txt')
        open(newunit=unit_report, file=trim(report_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,I0)') 'Warning: failed to write calibration coefficient report for level ', level
            flush(output_unit)
            return
        end if

        write(unit_report,'(A,I0)') 'Level: ', level
        write(unit_report,'(A,A)') 'Level directory: ', trim(level_dir)
        write(unit_report,'(A,A)') 'Calibration folder: ', trim(calib_dir)
        write(unit_report,'(A)') ''
        write(unit_report,'(A)') '[X-side QR fit]'
        write(unit_report,'(A,I0)') 'order_count = ', low_order_count
        write(unit_report,'(A,I0)') 'fit_points = ', fit_count_low
        write(unit_report,'(A,F16.8)') 'rms_residual_eV = ', rms_low
        write(unit_report,'(A,F16.8)') 'max_residual_eV = ', max_res_low
        do i = 1, low_order_count
            write(unit_report,'(A,I0,A,F20.10)') 'c', i, ' = ', coeff_low(i)
        end do

        write(unit_report,'(A)') ''
        write(unit_report,'(A)') '[Y-side QR fit]'
        if (have_high) then
            write(unit_report,'(A,I0)') 'order_count = ', high_order_count
            write(unit_report,'(A,I0)') 'fit_points = ', fit_count_high
            write(unit_report,'(A,F16.8)') 'rms_residual_eV = ', rms_high
            write(unit_report,'(A,F16.8)') 'max_residual_eV = ', max_res_high
            do i = 1, high_order_count
                write(unit_report,'(A,I0,A,F20.10)') 'c', i, ' = ', coeff_high(i)
            end do
        else
            write(unit_report,'(A)') 'not available'
        end if

        write(unit_report,'(A)') ''
        write(unit_report,'(A)') '[Low/High blend]'
        if (blend_points > 0 .and. blend_lambda_high >= 0.0_dp) then
            write(unit_report,'(A,I0)') 'fit_points = ', blend_points
            write(unit_report,'(A,F16.8)') 'lambda_high = ', blend_lambda_high
            write(unit_report,'(A,F16.8)') 'lambda_low = ', 1.0_dp - blend_lambda_high
            write(unit_report,'(A,F16.8)') 'rms_residual_eV = ', blend_rms
        else
            write(unit_report,'(A)') 'not fitted'
        end if

        close(unit_report)
        call append_calibration_coefficients_summary(level, total_sites, level_prefix, low_order_count, coeff_low, fit_count_low, rms_low, &
            max_res_low, have_high, high_order_count, coeff_high, fit_count_high, rms_high, max_res_high, blend_lambda_high, blend_rms, &
            blend_points, calib_dir, level_dir, report_path)
        write(*,'(A,I0,A,A)') 'Level ', level, ': calibration coefficients saved to ', trim(report_path)
        flush(output_unit)
    end subroutine write_calibration_coefficients_report

    subroutine append_calibration_coefficients_summary(level, total_sites, level_prefix, low_order_count, coeff_low, fit_count_low, rms_low, &
        max_res_low, have_high, high_order_count, coeff_high, fit_count_high, rms_high, max_res_high, blend_lambda_high, blend_rms, &
        blend_points, calib_dir, level_dir, report_path)
        integer, intent(in) :: level, total_sites, low_order_count, fit_count_low, high_order_count, fit_count_high, blend_points
        character(len=*), intent(in) :: level_prefix, calib_dir, level_dir, report_path
        logical, intent(in) :: have_high
        real(dp), intent(in) :: coeff_low(:), coeff_high(:)
        real(dp), intent(in) :: rms_low, max_res_low, rms_high, max_res_high, blend_lambda_high, blend_rms
        character(len=*), parameter :: summary_filename = 'CALIBRATION_COEFFICIENTS_SUMMARY.csv'
        character(len=1024) :: case_dir
        character(len=256) :: case_label, template_label, template_payload
        character(len=8192) :: line
        logical :: summary_exists
        integer :: unit_summary, ios

        call get_environment_variable('PWD', case_dir, status=ios)
        if (ios /= 0 .or. len_trim(case_dir) == 0) case_dir = '.'
        case_label = basename_from_path(trim(case_dir))
        template_label = calibration_template_mode_label()
        template_payload = calibration_template_payload_label()

        inquire(file=summary_filename, exist=summary_exists)
        open(newunit=unit_summary, file=summary_filename, status='unknown', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,I0)') 'Warning: failed to append calibration summary CSV for level ', level
            flush(output_unit)
            return
        end if

        if (.not. summary_exists) then
            write(unit_summary,'(A)') 'case_dir,case_label,mode_prefix,level,total_sites,template_mode,template_payload,'// &
                'x_order_count,x_fit_points,x_rms_residual_eV,x_max_residual_eV,x_c1,x_c2,x_c3,x_c4,'// &
                'y_available,y_order_count,y_fit_points,y_rms_residual_eV,y_max_residual_eV,y_c1,y_c2,y_c3,y_c4,'// &
                'blend_fit_points,lambda_high,lambda_low,blend_rms_residual_eV,level_directory,calibration_folder,report_file'
        end if

        line = trim(csv_quote(case_dir))//','//trim(csv_quote(case_label))//','//trim(csv_quote(level_prefix))//','// &
            trim(int_text(level))//','//trim(int_text(total_sites))//','//trim(csv_quote(template_label))//','// &
            trim(csv_quote(template_payload))//','//trim(int_text(low_order_count))//','//trim(int_text(fit_count_low))//','// &
            trim(real_or_blank(rms_low, fit_count_low > 0))//','//trim(real_or_blank(max_res_low, fit_count_low > 0))//','// &
            trim(coeff_or_blank(coeff_low, 1, low_order_count))//','//trim(coeff_or_blank(coeff_low, 2, low_order_count))//','// &
            trim(coeff_or_blank(coeff_low, 3, low_order_count))//','//trim(coeff_or_blank(coeff_low, 4, low_order_count))//','// &
            trim(logical_text(have_high))//','//trim(int_text(high_order_count))//','//trim(int_text(fit_count_high))//','// &
            trim(real_or_blank(rms_high, have_high .and. fit_count_high > 0))//','//trim(real_or_blank(max_res_high, have_high .and. fit_count_high > 0))//','// &
            trim(coeff_or_blank(coeff_high, 1, high_order_count, have_high))//','//trim(coeff_or_blank(coeff_high, 2, high_order_count, have_high))//','// &
            trim(coeff_or_blank(coeff_high, 3, high_order_count, have_high))//','//trim(coeff_or_blank(coeff_high, 4, high_order_count, have_high))//','// &
            trim(int_text(blend_points))//','//trim(real_or_blank(blend_lambda_high, blend_points > 0 .and. blend_lambda_high >= 0.0_dp))//','// &
            trim(real_or_blank(1.0_dp - blend_lambda_high, blend_points > 0 .and. blend_lambda_high >= 0.0_dp))//','// &
            trim(real_or_blank(blend_rms, blend_points > 0 .and. blend_lambda_high >= 0.0_dp))//','//trim(csv_quote(level_dir))//','// &
            trim(csv_quote(calib_dir))//','//trim(csv_quote(report_path))
        write(unit_summary,'(A)') trim(line)
        close(unit_summary)
        write(*,'(A,I0,A,A)') 'Level ', level, ': calibration summary appended to ', summary_filename
        flush(output_unit)
    end subroutine append_calibration_coefficients_summary

    subroutine append_blend_override_summary(level, total_sites, level_prefix, source_label, lambda_high, rms_residual, fit_points)
        integer, intent(in) :: level, total_sites, fit_points
        character(len=*), intent(in) :: level_prefix, source_label
        real(dp), intent(in) :: lambda_high, rms_residual
        character(len=*), parameter :: summary_filename = 'CALIBRATION_BLEND_OVERRIDES_SUMMARY.csv'
        character(len=1024) :: case_dir
        character(len=256) :: case_label, template_label, template_payload
        logical :: summary_exists
        integer :: unit_summary, ios

        call get_environment_variable('PWD', case_dir, status=ios)
        if (ios /= 0 .or. len_trim(case_dir) == 0) case_dir = '.'
        case_label = basename_from_path(trim(case_dir))
        template_label = calibration_template_mode_label()
        template_payload = calibration_template_payload_label()

        inquire(file=summary_filename, exist=summary_exists)
        open(newunit=unit_summary, file=summary_filename, status='unknown', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,I0)') 'Warning: failed to append blend override summary CSV for level ', level
            flush(output_unit)
            return
        end if

        if (.not. summary_exists) then
            write(unit_summary,'(A)') 'case_dir,case_label,mode_prefix,level,total_sites,template_mode,template_payload,'// &
                'source_label,fit_points,lambda_high,lambda_low,rms_residual_eV,level_directory'
        end if

        write(unit_summary,'(A)') trim(csv_quote(case_dir))//','//trim(csv_quote(case_label))//','//trim(csv_quote(level_prefix))//','// &
            trim(int_text(level))//','//trim(int_text(total_sites))//','//trim(csv_quote(template_label))//','//trim(csv_quote(template_payload))//','// &
            trim(csv_quote(source_label))//','//trim(int_text(fit_points))//','//trim(real_text(lambda_high))//','// &
            trim(real_text(1.0_dp - lambda_high))//','//trim(real_text(rms_residual))//','// &
            trim(csv_quote(format_level_directory(trim(level_prefix), level)))
        close(unit_summary)
        write(*,'(A,I0,A,A)') 'Level ', level, ': blend override summary appended to ', summary_filename
        flush(output_unit)
    end subroutine append_blend_override_summary

    function calibration_template_mode_label() result(label)
        character(len=32) :: label

        select case (template_mode)
        case (template_mode_builtin)
            label = 'builtin'
        case (template_mode_custom)
            label = 'custom'
        case default
            label = 'none'
        end select
    end function calibration_template_mode_label

    function calibration_template_payload_label() result(label)
        character(len=256) :: label

        label = ''
        if (template_mode == template_mode_custom) then
            label = basename_from_path(trim(template_gin_override))
        end if
    end function calibration_template_payload_label

    function csv_quote(text) result(out)
        character(len=*), intent(in) :: text
        character(len=2048) :: out
        integer :: i, pos

        out = ' '
        pos = 1
        out(pos:pos) = '"'
        pos = pos + 1
        do i = 1, len_trim(text)
            if (pos > len(out) - 2) exit
            if (text(i:i) == '"') then
                out(pos:pos+1) = '""'
                pos = pos + 2
            else
                out(pos:pos) = text(i:i)
                pos = pos + 1
            end if
        end do
        if (pos <= len(out)) out(pos:pos) = '"'
    end function csv_quote

    function int_text(value) result(out)
        integer, intent(in) :: value
        character(len=32) :: out

        write(out,'(I0)') value
    end function int_text

    function logical_text(flag) result(out)
        logical, intent(in) :: flag
        character(len=5) :: out

        if (flag) then
            out = 'true '
        else
            out = 'false'
        end if
    end function logical_text

    function real_text(value) result(out)
        real(dp), intent(in) :: value
        character(len=48) :: out

        write(out,'(ES24.16E3)') value
    end function real_text

    function real_or_blank(value, enabled) result(out)
        real(dp), intent(in) :: value
        logical, intent(in) :: enabled
        character(len=48) :: out

        out = ''
        if (enabled) out = trim(real_text(value))
    end function real_or_blank

    function coeff_or_blank(coeff, idx, order_count, enabled) result(out)
        real(dp), intent(in) :: coeff(:)
        integer, intent(in) :: idx, order_count
        logical, intent(in), optional :: enabled
        character(len=48) :: out
        logical :: use_coeff

        use_coeff = .true.
        if (present(enabled)) use_coeff = enabled
        out = ''
        if (.not. use_coeff) return
        if (idx < 1 .or. idx > order_count) return
        if (idx > size(coeff)) return
        out = trim(real_text(coeff(idx)))
    end function coeff_or_blank

    ! Subroutine `evaluate_subsets_with_engine` dispatches subset evaluation to the active external backend.
    ! La subrutina `evaluate_subsets_with_engine` delega la evaluación de subconjuntos al backend externo activo.
    subroutine evaluate_subsets_with_engine(level, total_sites, subset_count, subsets, config, energies, success)
        implicit none
        integer, intent(in) :: level, total_sites, subset_count
        integer, intent(in) :: subsets(:,:)
        integer, intent(inout) :: config(:)
        real(dp), intent(out) :: energies(:)
        logical, intent(out) :: success

        if (.not. external_calibration_enabled) then
            success = .false.
            return
        end if

        if (calibration_engine == calibration_engine_none) then
            success = .false.
            return
        end if

        call evaluate_subsets_with_backend(level, total_sites, subset_count, subsets, config, energies, success, calibration_engine)
    end subroutine evaluate_subsets_with_engine

    ! Subroutine `evaluate_subsets_with_gulp` keeps the legacy GULP-only evaluation entry point.
    ! La subrutina `evaluate_subsets_with_gulp` mantiene el punto de entrada heredado de evaluación solo con GULP.
    subroutine evaluate_subsets_with_gulp(level, total_sites, subset_count, subsets, config, energies, success)
        implicit none
        integer, intent(in) :: level, total_sites, subset_count
        integer, intent(in) :: subsets(:,:)
        integer, intent(inout) :: config(:)
        real(dp), intent(out) :: energies(:)
        logical, intent(out) :: success

        call evaluate_subsets_with_backend(level, total_sites, subset_count, subsets, config, energies, success, calibration_engine_gulp)
    end subroutine evaluate_subsets_with_gulp

    ! Subroutine `evaluate_subsets_with_backend` performs external subset evaluation with a selected backend.
    ! La subrutina `evaluate_subsets_with_backend` realiza la evaluación externa de subconjuntos con un backend seleccionado.
    subroutine evaluate_subsets_with_backend(level, total_sites, subset_count, subsets, config, energies, success, engine)
        implicit none
        integer, intent(in) :: level, total_sites, subset_count, engine
        integer, intent(in) :: subsets(:,:)
        integer, intent(inout) :: config(:)
        real(dp), intent(out) :: energies(:)
        logical, intent(out) :: success

        character(len=512) :: script_dir, sample_dir
        character(len=16) :: level_tag
        character(len=64) :: base_name
        character(len=256) :: gout_label
        character(len=256) :: engine_label
        character(len=1024) :: command, line
        character(len=128), allocatable :: gout_files(:)
        logical, allocatable :: filled(:)
        integer :: exit_code, ios, unit_energy
        integer :: i, pos_hash, match_idx
        real(dp) :: energy_val

        success = .false.

        if (subset_count <= 0) return
        if (size(energies) < subset_count) then
            write(*,'(A)') 'Warning: insufficient energy buffer in evaluate_subsets_with_backend.'
            flush(output_unit)
            return
        end if

        engine_label = calibration_engine_name_from_value(engine)

        if (.not. find_scripts_directory(script_dir, engine)) then
            write(*,'(A,A,A)') 'Warning: ', trim(engine_label), ' scripts directory not found.'
            flush(output_unit)
            return
        end if

        write(level_tag,'(I4.4)') level
        sample_dir = join_path(trim(format_level_directory('mc', level)), 'mc_'//trim(calibration_engine_stub(engine))//'_samples_N'//trim(adjustl(level_tag)))
        if (engine == calibration_engine_ase) then
            call ensure_ase_level_resources(format_level_directory('mc', level))
        end if
        call execute_command_line('rm -rf ' // trim(sample_dir), exitstat=exit_code)
        call execute_command_line('mkdir -p ' // trim(sample_dir), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A,A,A)') 'Warning: failed to create the ', trim(engine_label), ' evaluation directory.'
            flush(output_unit)
            return
        end if

        allocate(gout_files(subset_count))
        allocate(filled(subset_count))
        gout_files = ''
        filled = .false.

        do i = 1, subset_count
            config = 1
            if (level > 0 .and. size(subsets,1) >= level) then
                config(subsets(1:level, i)) = 2
            end if
            write(base_name,'("sample_N",I4.4,"_c",I5.5)') level, i
            call write_vasp_file(config, total_sites, trim(sample_dir)//'/'//trim(base_name)//'.vasp')
            gout_files(i) = trim(expected_energy_output_name(base_name, engine))
        end do
        config = 1

        call copy_calibration_scripts(script_dir, sample_dir, engine)
        call write_engine_filer_marker(sample_dir, engine)
        write(*,'(A,I0,A,A,A,A)') 'Level ', level, ': ', trim(engine_label), ' configurations in ', trim(sample_dir)
        flush(output_unit)

        command = build_engine_run_command(sample_dir, engine)
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A,A,A)') 'Warning: run_jobs.sh reported issues while evaluating the ', trim(engine_label), ' samples; attempting to extract partial energies.'
            flush(output_unit)
        end if

        command = 'cd ' // trim(sample_dir) // ' && bash extract.sh'
        call execute_command_line(trim(command), exitstat=exit_code)
        call append_mc_gulp_energy_snapshot(level, sample_dir)
        if (engine == calibration_engine_gulp) call remove_protocol_template_file(sample_dir)
        if (exit_code /= 0) then
            write(*,'(A,A,A)') 'Warning: extract.sh failed while evaluating the ', trim(engine_label), ' samples.'
            flush(output_unit)
        end if

        open(newunit=unit_energy, file=trim(sample_dir)//'/ENERGIES', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,A,A)') 'Warning: failed to open the ENERGIES file generated by ', trim(engine_label), '.'
            flush(output_unit)
            call cleanup()
            return
        end if

        do
            read(unit_energy,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            pos_hash = index(line, '#')
            if (pos_hash <= 0) cycle
            read(line(1:pos_hash-1), *, iostat=ios) energy_val
            if (ios /= 0) cycle
            gout_label = adjustl(trim(line(pos_hash+1:)))
            gout_label = trim(gout_label)
            pos_hash = index(gout_label, '#')
            if (pos_hash > 0) gout_label = adjustl(trim(gout_label(:pos_hash-1)))
            match_idx = 0
            do i = 1, subset_count
                if (gout_label == trim(gout_files(i))) then
                    match_idx = i
                    exit
                end if
            end do
            if (match_idx > 0) then
                energies(match_idx) = energy_val
                filled(match_idx) = .true.
            end if
        end do
        close(unit_energy)

        if (.not. all(filled)) then
            write(*,'(A,A,A)') 'Warning: missing energies in the ', trim(engine_label), ' output.'
            flush(output_unit)
            call cleanup()
            return
        end if

        success = .true.
        call cleanup()
        return

    contains
        subroutine cleanup()
            if (allocated(gout_files)) deallocate(gout_files)
            if (allocated(filled)) deallocate(filled)
        end subroutine cleanup
    end subroutine evaluate_subsets_with_backend

    subroutine remove_protocol_template_file(dir)
        implicit none
        character(len=*), intent(in) :: dir
        character(len=1024) :: path
        integer :: unit_file, ios
        logical :: exists

        path = trim(dir)//'/protocol_template.gin'
        inquire(file=trim(path), exist=exists)
        if (.not. exists) return

        open(newunit=unit_file, file=trim(path), status='old', action='readwrite', iostat=ios)
        if (ios /= 0) return
        close(unit_file, status='delete', iostat=ios)
    end subroutine remove_protocol_template_file

    logical function find_scripts_directory(out_dir, engine)
        implicit none
        character(len=*), intent(out) :: out_dir
        integer, intent(in) :: engine
        character(len=512) :: env_dir
        integer :: lenv, i
        character(len=512), dimension(6) :: candidates

        out_dir = ''
        env_dir = ''
        lenv = 0
        call get_environment_variable('SOD_SCRIPTS', env_dir, length=lenv)
        if (lenv > 0) then
            env_dir = env_dir(1:lenv)
            if (directory_has_scripts(trim(env_dir), engine)) then
                out_dir = trim(env_dir)
                find_scripts_directory = .true.
                return
            end if
        end if

        candidates = ''
        candidates(1) = 'scripts'
        candidates(2) = '../scripts'
        candidates(3) = '../../scripts'
        candidates(4) = '../../../scripts'
        candidates(5) = '../../../../scripts'
        candidates(6) = '/home/salvador/coding/sod/sod/scripts'
        do i = 1, size(candidates)
            if (directory_has_scripts(trim(candidates(i)), engine)) then
                out_dir = trim(candidates(i))
                find_scripts_directory = .true.
                return
            end if
        end do
        find_scripts_directory = .false.
    end function find_scripts_directory

    logical function directory_has_scripts(dir, engine)
        implicit none
        character(len=*), intent(in) :: dir
        integer, intent(in) :: engine
        logical :: exist_vasp, exist_run, exist_extract

        inquire(file=trim(dir)//'/run_jobs.sh', exist=exist_run)
        inquire(file=trim(dir)//'/extract.sh', exist=exist_extract)
        select case (engine)
        case (calibration_engine_ase)
            inquire(file=trim(dir)//'/vasp2ase.py', exist=exist_vasp)
        case default
            inquire(file=trim(dir)//'/vasp2gin.sh', exist=exist_vasp)
        end select
        directory_has_scripts = exist_vasp .and. exist_run .and. exist_extract
    end function directory_has_scripts

    subroutine select_calibration_indices(values, n, indices, ok)
        implicit none
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: n
        integer, intent(out) :: indices(:)
        logical, intent(out) :: ok
        integer, allocatable :: sorted_idx(:)
        integer :: s, m, pick
        real(dp) :: position

        if (n < size(indices)) then
            ok = .false.
            return
        end if

        allocate(sorted_idx(n))
        call sort_indices_by_values(values, sorted_idx, n)

        m = size(indices)
        do s = 1, m
            position = 1.0_dp + floor(real(n - 1, dp) * real(s - 1, dp) / real(m - 1, dp))
            pick = int(position)
            pick = max(1, min(n, pick))
            indices(s) = sorted_idx(pick)
        end do

        deallocate(sorted_idx)
        ok = .true.
    end subroutine select_calibration_indices

    subroutine sort_indices_by_values(values, indices, n)
        implicit none
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: indices(:)
        integer, intent(in) :: n
        integer :: i, j, tmp

        do i = 1, n
            indices(i) = i
        end do

        do i = 1, n - 1
            do j = i + 1, n
                if (values(indices(j)) < values(indices(i))) then
                    tmp = indices(i)
                    indices(i) = indices(j)
                    indices(j) = tmp
                end if
            end do
        end do
    end subroutine sort_indices_by_values

    subroutine copy_calibration_scripts(script_dir, calib_dir, engine)
        implicit none
        character(len=*), intent(in) :: script_dir
        character(len=*), intent(in) :: calib_dir
        integer, intent(in) :: engine
        integer :: exit_code
        character(len=1024) :: command
        logical :: exists

        select case (engine)
        case (calibration_engine_ase)
            inquire(file='template_ase.py', exist=exists)
            if (.not. exists) then
                write(*,'(A)') 'Warning: template_ase.py was not found in the calculation directory; ASE calibration is unavailable.'
                return
            end if
            command = 'cp "' // trim(script_dir)//'/run_jobs.sh" "' // trim(script_dir)//'/extract.sh" "' // &
                trim(script_dir)//'/vasp2ase.py" "template_ase.py" "' // trim(calib_dir) // '"'
            call execute_command_line(trim(command), exitstat=exit_code)
            command = 'chmod +x "' // trim(calib_dir)//'/run_jobs.sh" "' // trim(calib_dir)//'/extract.sh" "' // trim(calib_dir)//'/vasp2ase.py"'
            call execute_command_line(trim(command), exitstat=exit_code)
        case default
            command = 'cp "' // trim(script_dir)//'/vasp2gin.sh" "' // trim(calib_dir)//'/vasp2gin.sh"'
            call execute_command_line(trim(command), exitstat=exit_code)
            command = 'cp "' // trim(script_dir)//'/run_jobs.sh" "' // trim(calib_dir)//'/run_jobs.sh"'
            call execute_command_line(trim(command), exitstat=exit_code)
            command = 'cp "' // trim(script_dir)//'/extract.sh" "' // trim(calib_dir)//'/extract.sh"'
            call execute_command_line(trim(command), exitstat=exit_code)
            command = 'cp "' // trim(script_dir)//'/framework_template.lib" "' // trim(calib_dir)//'/framework_template.lib"'
            call execute_command_line(trim(command), exitstat=exit_code)
            inquire(file=trim(script_dir)//'/reference/protocol_template.gin', exist=exists)
            if (exists) then
                command = 'cp "' // trim(script_dir)//'/reference/protocol_template.gin" "' // trim(calib_dir)//'/protocol_template.gin"'
                call execute_command_line(trim(command), exitstat=exit_code)
            else
                write(*,'(A)') 'Warning: scripts/reference/protocol_template.gin was not found; protocol mode may be unavailable.'
            end if

            select case (template_mode)
            case (template_mode_builtin)
                inquire(file=trim(script_dir)//'/default_template.gin', exist=exists)
                if (exists) then
                    command = 'cp "' // trim(script_dir)//'/default_template.gin" "' // trim(calib_dir)//'/template_payload.gin"'
                    call execute_command_line(trim(command), exitstat=exit_code)
                else
                    inquire(file=trim(script_dir)//'/default_template.include', exist=exists)
                    if (exists) then
                        command = 'cp "' // trim(script_dir)//'/default_template.include" "' // trim(calib_dir)//'/default_template.include"'
                        call execute_command_line(trim(command), exitstat=exit_code)
                    else
                        write(*,'(A)') 'Warning: default_template.include was not found in '//trim(script_dir)//'; skipping template payload.'
                    end if
                end if
            case (template_mode_custom)
                command = 'cp "' // trim(template_gin_override) // '" "' // trim(calib_dir)//'/template_payload.gin"'
                call execute_command_line(trim(command), exitstat=exit_code)
            case (template_mode_none)
                continue
            end select

            command = 'chmod +x "' // trim(calib_dir)//'/vasp2gin.sh" "' // trim(calib_dir)//'/run_jobs.sh" "' // trim(calib_dir)//'/extract.sh"'
            call execute_command_line(trim(command), exitstat=exit_code)
        end select
    end subroutine copy_calibration_scripts

    ! Subroutine `write_engine_filer_marker` writes a helper FILER file for the selected backend.
    ! La subrutina `write_engine_filer_marker` escribe un archivo FILER auxiliar para el backend seleccionado.
    subroutine write_engine_filer_marker(dir, engine)
        implicit none
        character(len=*), intent(in) :: dir
        integer, intent(in) :: engine
        integer :: unit_file, ios

        if (calibration_engine_filer_value(engine) <= 0) return

        open(newunit=unit_file, file=trim(dir)//'/filer', status='replace', action='write', iostat=ios)
        if (ios /= 0) return
        write(unit_file,'(I0)') calibration_engine_filer_value(engine)
        close(unit_file)
    end subroutine write_engine_filer_marker

    ! Subroutine `ensure_ase_level_resources` recreates template-adjacent resources for nested ASE work directories.
    ! La subrutina `ensure_ase_level_resources` recrea los recursos adyacentes al template para directorios ASE anidados.
    subroutine ensure_ase_level_resources(level_dir)
        implicit none
        character(len=*), intent(in) :: level_dir
        logical :: has_tools, level_has_tools
        character(len=1024) :: command
        integer :: exit_code

        inquire(file='tools', exist=has_tools)
        if (.not. has_tools) return

        inquire(file=trim(level_dir)//'/tools', exist=level_has_tools)
        if (level_has_tools) return

        command = 'bash -c "ln -s ../tools ' // trim(level_dir)//'/tools"'
        call execute_command_line(trim(command), exitstat=exit_code)
    end subroutine ensure_ase_level_resources

    subroutine solve_least_squares_qr(design, rhs, coeff, ok)
        implicit none
        real(dp), intent(in) :: design(:,:)
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: coeff(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: q(:,:), r(:,:), v(:), y(:)
        real(dp) :: norm_v, proj
        integer :: m, n, i, j, k
        real(dp), parameter :: tol = 1.0e-10_dp

        m = size(design, 1)
        n = size(design, 2)
        if (m <= 0 .or. n <= 0) then
            ok = .false.
            return
        end if
        if (size(rhs) /= m .or. size(coeff) /= n) then
            ok = .false.
            return
        end if
        if (m < n) then
            ok = .false.
            return
        end if

        allocate(q(m, n))
        allocate(r(n, n))
        allocate(v(m))
        allocate(y(n))
        q = 0.0_dp
        r = 0.0_dp
        y = 0.0_dp
        coeff = 0.0_dp

        do k = 1, n
            v = design(:, k)
            do j = 1, k - 1
                proj = dot_product(q(:, j), v)
                r(j, k) = proj
                v = v - proj * q(:, j)
            end do

            norm_v = sqrt(max(0.0_dp, dot_product(v, v)))
            if (norm_v < tol) then
                deallocate(q, r, v, y)
                ok = .false.
                return
            end if

            r(k, k) = norm_v
            q(:, k) = v / norm_v
        end do

        do k = 1, n
            y(k) = dot_product(q(:, k), rhs)
        end do

        do i = n, 1, -1
            if (abs(r(i, i)) < tol) then
                deallocate(q, r, v, y)
                ok = .false.
                return
            end if
            coeff(i) = y(i)
            do j = i + 1, n
                coeff(i) = coeff(i) - r(i, j) * coeff(j)
            end do
            coeff(i) = coeff(i) / r(i, i)
        end do

        deallocate(q, r, v, y)
        ok = .true.
    end subroutine solve_least_squares_qr

end module calibration
