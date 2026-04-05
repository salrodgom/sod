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
! two folders and emits gnuplot-ready data plus a fitting script.
! Comparison helper that combines MC and configurational-entropy summaries from
! Utilidad de comparación que combina los resúmenes de MC y de entropía configuracional de
! dos carpetas y genera datos listos para gnuplot junto con un script de ajuste.
!*******************************************************************************
module compare
    use consts
    use cli, only: basename_from_path, compose_mode_command, is_help_token, to_lower_inplace
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    private

    real(dp), parameter :: ev_to_kjmol = 96.48533212331002_dp
    real(dp), parameter :: invalid_projection_marker = 9.0e29_dp

    type :: phase_curve_t
        integer :: count = 0
        character(len=512) :: folder = ''
        character(len=256) :: label = ''
        integer, allocatable :: levels(:)
        real(dp), allocatable :: xy_fraction(:)
        real(dp), allocatable :: delta_h(:)
        real(dp), allocatable :: delta_h_low(:)
        real(dp), allocatable :: delta_h_high(:)
        real(dp), allocatable :: s_conf_site(:)
        real(dp), allocatable :: tdelta_s(:)
    end type phase_curve_t

    public :: run_sod_ensemble_compare

contains

    subroutine run_sod_ensemble_compare(arg_offset)
        integer, intent(in), optional :: arg_offset
        character(len=512) :: system_folder, reference_folder
        character(len=256) :: system_label, reference_label
        character(len=512) :: output_prefix
        real(dp) :: temperature
        logical :: has_temperature
        type(phase_curve_t) :: system_curve, reference_curve
        character(len=512) :: system_data_file, reference_data_file, script_file

        call parse_arguments_compare(system_folder, reference_folder, system_label, reference_label, output_prefix, &
            temperature, has_temperature, arg_offset)

        if (.not. has_temperature) then
            write(error_unit,'(A)') 'Error: the comparison mode requires --temperature <K>.'
            call print_usage_compare()
            stop 1
        end if

        call normalize_folder_path(system_folder)
        call normalize_folder_path(reference_folder)

        if (len_trim(system_label) == 0) system_label = basename_from_path(system_folder)
        if (len_trim(reference_label) == 0) reference_label = basename_from_path(reference_folder)
        if (len_trim(output_prefix) == 0) then
            output_prefix = trim(basename_from_path(system_folder))//'_vs_'//trim(basename_from_path(reference_folder))
        end if

        call load_phase_curve(system_folder, system_label, temperature, system_curve, require_mc=.true.)
        call load_phase_curve(reference_folder, reference_label, temperature, reference_curve, require_mc=.false.)

        system_data_file = trim(output_prefix)//'_system.dat'
        reference_data_file = trim(output_prefix)//'_reference.dat'
        script_file = trim(output_prefix)//'.gnuplot'

        call write_phase_data(system_data_file, system_curve, temperature)
        call write_phase_data(reference_data_file, reference_curve, temperature)
        call write_gnuplot_script(script_file, trim(system_data_file), trim(reference_data_file), trim(output_prefix), &
            trim(system_curve%label), trim(reference_curve%label), temperature)

        write(*,'(A)') '--- Configurational DeltaG comparison ---'
        write(*,'(A,F10.3,A)') 'Temperature: ', temperature, ' K'
        write(*,'(A)') 'System folder:    '//trim(system_curve%folder)
        write(*,'(A)') 'Reference folder: '//trim(reference_curve%folder)
        write(*,'(A,I0)') 'System points loaded: ', system_curve%count
        write(*,'(A,I0)') 'Reference points loaded: ', reference_curve%count
        write(*,'(A)') 'System data file:    '//trim(system_data_file)
        write(*,'(A)') 'Reference data file: '//trim(reference_data_file)
        write(*,'(A)') 'Gnuplot script:      '//trim(script_file)
        write(*,'(A)') 'Run "gnuplot '//trim(script_file)//'" to generate the comparison plots.'
        flush(output_unit)
    end subroutine run_sod_ensemble_compare

    subroutine parse_arguments_compare(system_folder, reference_folder, system_label, reference_label, output_prefix, &
        temperature, has_temperature, arg_offset)
        character(len=*), intent(out) :: system_folder, reference_folder
        character(len=*), intent(out) :: system_label, reference_label
        character(len=*), intent(out) :: output_prefix
        real(dp), intent(out) :: temperature
        logical, intent(out) :: has_temperature
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset, eq_pos, ios
        character(len=512) :: arg, lowered, value

        system_folder = ''
        reference_folder = ''
        system_label = ''
        reference_label = ''
        output_prefix = ''
        temperature = 0.0_dp
        has_temperature = .false.

        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)
        if (argc <= offset) then
            call print_usage_compare()
            stop 1
        end if

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            if (is_help_token(arg)) then
                call print_usage_compare()
                stop
            end if

            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--system=') == 1) then
                eq_pos = index(arg, '=')
                system_folder = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--reference=') == 1) then
                eq_pos = index(arg, '=')
                reference_folder = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--temperature=') == 1) then
                eq_pos = index(arg, '=')
                read(arg(eq_pos+1:), *, iostat=ios) temperature
                if (ios /= 0 .or. temperature <= 0.0_dp) then
                    write(error_unit,'(A)') 'Error: invalid value provided to --temperature.'
                    stop 1
                end if
                has_temperature = .true.
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--output-prefix=') == 1) then
                eq_pos = index(arg, '=')
                output_prefix = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--system-label=') == 1) then
                eq_pos = index(arg, '=')
                system_label = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--reference-label=') == 1) then
                eq_pos = index(arg, '=')
                reference_label = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('--system')
                call require_following_value(iarg, argc, '--system', value)
                system_folder = trim(value)
                iarg = iarg + 2
            case ('--reference')
                call require_following_value(iarg, argc, '--reference', value)
                reference_folder = trim(value)
                iarg = iarg + 2
            case ('-t','--temperature')
                call require_following_value(iarg, argc, '--temperature', value)
                read(value, *, iostat=ios) temperature
                if (ios /= 0 .or. temperature <= 0.0_dp) then
                    write(error_unit,'(A)') 'Error: invalid value provided to --temperature.'
                    stop 1
                end if
                has_temperature = .true.
                iarg = iarg + 2
            case ('-o','--output-prefix','--output')
                call require_following_value(iarg, argc, '--output-prefix', value)
                output_prefix = trim(value)
                iarg = iarg + 2
            case ('--system-label')
                call require_following_value(iarg, argc, '--system-label', value)
                system_label = trim(value)
                iarg = iarg + 2
            case ('--reference-label')
                call require_following_value(iarg, argc, '--reference-label', value)
                reference_label = trim(value)
                iarg = iarg + 2
            case default
                write(error_unit,'(A)') 'Error: unrecognized argument in compare mode.'
                call print_usage_compare()
                stop 1
            end select
        end do

        if (len_trim(system_folder) == 0) then
            write(error_unit,'(A)') 'Error: --system <folder> is required.'
            call print_usage_compare()
            stop 1
        end if
        if (len_trim(reference_folder) == 0) then
            write(error_unit,'(A)') 'Error: --reference <folder> is required.'
            call print_usage_compare()
            stop 1
        end if
    end subroutine parse_arguments_compare

    subroutine require_following_value(iarg, argc, opt_name, value)
        integer, intent(in) :: iarg, argc
        character(len=*), intent(in) :: opt_name
        character(len=*), intent(out) :: value

        if (iarg + 1 > argc) then
            write(error_unit,'(A)') 'Error: missing value after '//trim(opt_name)//'.'
            stop 1
        end if
        call get_command_argument(iarg + 1, value)
        value = adjustl(value)
        if (len_trim(value) == 0) then
            write(error_unit,'(A)') 'Error: missing value after '//trim(opt_name)//'.'
            stop 1
        end if
    end subroutine require_following_value

    subroutine print_usage_compare()
        character(len=256) :: command_name

        command_name = compose_mode_command('compare')
        write(*,'(A)') 'Usage: '//trim(command_name)//' --system <folder> --reference <folder> --temperature <K> [options]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: compare, diff, reference'
        write(*,'(A)') ''
        write(*,'(A)') 'Required arguments:'
        write(*,'(A)') '  --system <folder>        Folder containing sod_ensemble_summary.csv and sod_entropy_summary.csv.'
        write(*,'(A)') '  --reference <folder>     Reference folder containing sod_entropy_summary.csv.'
        write(*,'(A)') '  -t, --temperature <K>    Temperature used to convert S_conf into TDeltaS.'
        write(*,'(A)') ''
        write(*,'(A)') 'Optional arguments:'
        write(*,'(A)') '  -o, --output-prefix <p>  Prefix for the generated .dat and .gnuplot files.'
        write(*,'(A)') '  --system-label <text>    Plot label for the system folder.'
        write(*,'(A)') '  --reference-label <text> Plot label for the reference folder.'
        write(*,'(A)') ''
        write(*,'(A)') 'The generated gnuplot script fits the reference TDeltaS curve with'
        write(*,'(A)') 'order-2/3/4 x*(1-x) polynomials and evaluates total, low-branch, and'
        write(*,'(A)') 'high-branch DeltaG estimates for the study system.'
    end subroutine print_usage_compare

    subroutine normalize_folder_path(path)
        character(len=*), intent(inout) :: path

        path = adjustl(path)
        do while (len_trim(path) > 1 .and. path(len_trim(path):len_trim(path)) == '/')
            path(len_trim(path):len_trim(path)) = ' '
        end do
    end subroutine normalize_folder_path

    subroutine load_phase_curve(folder, label, temperature, curve, require_mc)
        character(len=*), intent(in) :: folder, label
        real(dp), intent(in) :: temperature
        type(phase_curve_t), intent(out) :: curve
        logical, intent(in) :: require_mc
        integer, allocatable :: mc_levels(:), entropy_levels(:)
        real(dp), allocatable :: mc_x(:), mc_delta_h(:), mc_delta_h_low(:), mc_delta_h_high(:)
        real(dp), allocatable :: entropy_x(:), entropy_sconf_site(:)
        integer :: mc_count, entropy_count
        integer :: idx, entropy_idx

        curve%folder = trim(folder)
        curve%label = trim(label)

        call read_entropy_summary(trim(folder), entropy_levels, entropy_x, entropy_sconf_site, entropy_count)
        if (require_mc) then
            call read_mc_summary(trim(folder), mc_levels, mc_x, mc_delta_h, mc_delta_h_low, mc_delta_h_high, mc_count)

            do idx = 1, mc_count
                entropy_idx = find_level_index(entropy_levels, mc_levels(idx))
                if (entropy_idx <= 0) then
                    write(error_unit,'(A,I0,A)') 'Warning: level ', mc_levels(idx), ' is missing in sod_entropy_summary.csv; skipping.'
                    cycle
                end if

                if (abs(mc_x(idx) - entropy_x(entropy_idx)) > 1.0e-8_dp) then
                    write(error_unit,'(A,I0,A,ES12.5,A,ES12.5,A)') 'Warning: inconsistent xY for level ', mc_levels(idx), &
                        ' (MC=', mc_x(idx), ', entropy=', entropy_x(entropy_idx), '); using the MC value.'
                end if

                call append_curve_point(curve, mc_levels(idx), mc_x(idx), mc_delta_h(idx), mc_delta_h_low(idx), mc_delta_h_high(idx), &
                    entropy_sconf_site(entropy_idx), temperature)
            end do
        else
            do idx = 1, entropy_count
                call append_curve_point(curve, entropy_levels(idx), entropy_x(idx), 0.0_dp, 0.0_dp, 0.0_dp, &
                    entropy_sconf_site(idx), temperature)
            end do
        end if

        if (curve%count == 0) then
            if (require_mc) then
                write(error_unit,'(A)') 'Error: no common levels were found between the MC and entropy summaries in '//trim(folder)//'.'
            else
                write(error_unit,'(A)') 'Error: no valid entropy levels were found in '//trim(folder)//'.'
            end if
            stop 1
        end if

    end subroutine load_phase_curve

    subroutine read_mc_summary(folder, levels, xy_fraction, delta_h, delta_h_low, delta_h_high, count)
        character(len=*), intent(in) :: folder
        integer, allocatable, intent(out) :: levels(:)
        real(dp), allocatable, intent(out) :: xy_fraction(:), delta_h(:), delta_h_low(:), delta_h_high(:)
        integer, intent(out) :: count
        character(len=1024) :: filename, line
        integer :: unit_mc, ios, level
        real(dp) :: frac_y, delta_total, delta_low_val, delta_high_val
        logical :: ok

        count = 0
        allocate(levels(0))
        allocate(xy_fraction(0))
        allocate(delta_h(0))
        allocate(delta_h_low(0))
        allocate(delta_h_high(0))

        filename = trim(folder)//'/sod_ensemble_summary.csv'
        open(newunit=unit_mc, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to open '//trim(filename)
            stop 1
        end if

        do
            read(unit_mc,'(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            call parse_mc_summary_row(line, level, frac_y, delta_total, delta_low_val, delta_high_val, ok)
            if (.not. ok) then
                write(error_unit,'(A)') 'Warning: skipping malformed MC summary line in '//trim(filename)
                cycle
            end if

            call append_integer(levels, count, level)
            call append_real(xy_fraction, count, frac_y)
            call append_real(delta_h, count, delta_total)
            call append_real(delta_h_low, count, delta_low_val)
            call append_real(delta_h_high, count, delta_high_val)
        end do
        close(unit_mc)

        if (count == 0) then
            write(error_unit,'(A)') 'Error: no valid data rows were found in '//trim(filename)
            stop 1
        end if
    end subroutine read_mc_summary

    subroutine read_entropy_summary(folder, levels, xy_fraction, s_conf_site, count)
        character(len=*), intent(in) :: folder
        integer, allocatable, intent(out) :: levels(:)
        real(dp), allocatable, intent(out) :: xy_fraction(:), s_conf_site(:)
        integer, intent(out) :: count
        character(len=1024) :: filename, line, work
        integer :: unit_entropy, ios, level, unique_configs
        real(dp) :: frac_y, total_weight, log10_weight
        real(dp) :: s_conf_total, s_conf_site_val, s_conf_mol, s_max, s_ideal, reduction

        count = 0
        allocate(levels(0))
        allocate(xy_fraction(0))
        allocate(s_conf_site(0))

        filename = trim(folder)//'/sod_entropy_summary.csv'
        open(newunit=unit_entropy, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to open '//trim(filename)
            stop 1
        end if

        do
            read(unit_entropy,'(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            work = line
            call replace_char_inplace(work, ',', ' ')
            read(work, *, iostat=ios) level, frac_y, unique_configs, total_weight, log10_weight, s_conf_total, &
                s_conf_site_val, s_conf_mol, s_max, s_ideal, reduction
            if (ios /= 0) then
                write(error_unit,'(A)') 'Warning: skipping malformed entropy summary line in '//trim(filename)
                cycle
            end if

            call append_integer(levels, count, level)
            call append_real(xy_fraction, count, frac_y)
            call append_real(s_conf_site, count, s_conf_site_val)
        end do
        close(unit_entropy)

        if (count == 0) then
            write(error_unit,'(A)') 'Error: no valid data rows were found in '//trim(filename)
            stop 1
        end if
    end subroutine read_entropy_summary

    integer function find_level_index(levels, level) result(idx)
        integer, intent(in) :: levels(:)
        integer, intent(in) :: level
        integer :: i

        idx = 0
        do i = 1, size(levels)
            if (levels(i) == level) then
                idx = i
                return
            end if
        end do
    end function find_level_index

    subroutine append_curve_point(curve, level, xy_fraction, delta_h, delta_h_low, delta_h_high, s_conf_site, temperature)
        type(phase_curve_t), intent(inout) :: curve
        integer, intent(in) :: level
        real(dp), intent(in) :: xy_fraction, delta_h, delta_h_low, delta_h_high, s_conf_site, temperature
        integer :: new_count
        integer, allocatable :: tmp_levels(:)
        real(dp), allocatable :: tmp_xy_fraction(:), tmp_delta_h(:), tmp_delta_h_low(:), tmp_delta_h_high(:), tmp_s_conf(:), tmp_tdelta_s(:)
        real(dp) :: tdelta_s_val

        new_count = curve%count + 1
        tdelta_s_val = temperature * s_conf_site * ev_to_kjmol

        allocate(tmp_levels(new_count))
        allocate(tmp_xy_fraction(new_count))
        allocate(tmp_delta_h(new_count))
        allocate(tmp_delta_h_low(new_count))
        allocate(tmp_delta_h_high(new_count))
        allocate(tmp_s_conf(new_count))
        allocate(tmp_tdelta_s(new_count))

        if (curve%count > 0) then
            tmp_levels(1:curve%count) = curve%levels
            tmp_xy_fraction(1:curve%count) = curve%xy_fraction
            tmp_delta_h(1:curve%count) = curve%delta_h
            tmp_delta_h_low(1:curve%count) = curve%delta_h_low
            tmp_delta_h_high(1:curve%count) = curve%delta_h_high
            tmp_s_conf(1:curve%count) = curve%s_conf_site
            tmp_tdelta_s(1:curve%count) = curve%tdelta_s
        end if

        tmp_levels(new_count) = level
        tmp_xy_fraction(new_count) = xy_fraction
        tmp_delta_h(new_count) = delta_h
        tmp_delta_h_low(new_count) = delta_h_low
        tmp_delta_h_high(new_count) = delta_h_high
        tmp_s_conf(new_count) = s_conf_site
        tmp_tdelta_s(new_count) = tdelta_s_val

        if (allocated(curve%levels)) deallocate(curve%levels)
        if (allocated(curve%xy_fraction)) deallocate(curve%xy_fraction)
        if (allocated(curve%delta_h)) deallocate(curve%delta_h)
        if (allocated(curve%delta_h_low)) deallocate(curve%delta_h_low)
        if (allocated(curve%delta_h_high)) deallocate(curve%delta_h_high)
        if (allocated(curve%s_conf_site)) deallocate(curve%s_conf_site)
        if (allocated(curve%tdelta_s)) deallocate(curve%tdelta_s)

        call move_alloc(tmp_levels, curve%levels)
        call move_alloc(tmp_xy_fraction, curve%xy_fraction)
        call move_alloc(tmp_delta_h, curve%delta_h)
        call move_alloc(tmp_delta_h_low, curve%delta_h_low)
        call move_alloc(tmp_delta_h_high, curve%delta_h_high)
        call move_alloc(tmp_s_conf, curve%s_conf_site)
        call move_alloc(tmp_tdelta_s, curve%tdelta_s)

        curve%count = new_count
    end subroutine append_curve_point

    subroutine append_integer(vec, count, value)
        integer, allocatable, intent(inout) :: vec(:)
        integer, intent(inout) :: count
        integer, intent(in) :: value
        integer, allocatable :: tmp(:)

        allocate(tmp(count + 1))
        if (count > 0) tmp(1:count) = vec(1:count)
        tmp(count + 1) = value
        if (allocated(vec)) deallocate(vec)
        call move_alloc(tmp, vec)
        count = count + 1
    end subroutine append_integer

    subroutine append_real(vec, count, value)
        real(dp), allocatable, intent(inout) :: vec(:)
        integer, intent(in) :: count
        real(dp), intent(in) :: value
        real(dp), allocatable :: tmp(:)

        allocate(tmp(count))
        if (count > 1) tmp(1:count-1) = vec(1:count-1)
        tmp(count) = value
        if (allocated(vec)) deallocate(vec)
        call move_alloc(tmp, vec)
    end subroutine append_real

    subroutine replace_char_inplace(text, old_char, new_char)
        character(len=*), intent(inout) :: text
        character(len=1), intent(in) :: old_char, new_char
        integer :: idx

        do idx = 1, len_trim(text)
            if (text(idx:idx) == old_char) text(idx:idx) = new_char
        end do
    end subroutine replace_char_inplace

    subroutine write_phase_data(filename, curve, temperature)
        character(len=*), intent(in) :: filename
        type(phase_curve_t), intent(in) :: curve
        real(dp), intent(in) :: temperature
        integer :: unit_data, ios, idx

        open(newunit=unit_data, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to create '//trim(filename)
            stop 1
        end if

        write(unit_data,'(A)') '# folder='//trim(curve%folder)
        write(unit_data,'(A)') '# label='//trim(curve%label)
        write(unit_data,'(A,F10.4)') '# temperature_K=', temperature
        write(unit_data,'(A)') '# N xY DeltaH_kJmol DeltaH_low_kJmol DeltaH_high_kJmol Sconf_site_eV_per_K TDeltaS_kJmol'

        do idx = 1, curve%count
            write(unit_data,'(I0,1X,F10.6,5(1X,ES20.10))') curve%levels(idx), curve%xy_fraction(idx), curve%delta_h(idx), &
                curve%delta_h_low(idx), curve%delta_h_high(idx), curve%s_conf_site(idx), curve%tdelta_s(idx)
        end do

        close(unit_data)
    end subroutine write_phase_data

    subroutine write_gnuplot_script(filename, system_data_file, reference_data_file, output_prefix, system_label, reference_label, temperature)
        character(len=*), intent(in) :: filename, system_data_file, reference_data_file, output_prefix
        character(len=*), intent(in) :: system_label, reference_label
        real(dp), intent(in) :: temperature
        integer :: unit_script, ios
        character(len=512) :: safe_system_label, safe_reference_label
        character(len=512) :: plot_file, fit_log_file

        safe_system_label = gnuplot_safe_text(system_label)
        safe_reference_label = gnuplot_safe_text(reference_label)
        plot_file = trim(output_prefix)//'_deltaG_compare.png'
        fit_log_file = trim(output_prefix)//'_fit.log'

        open(newunit=unit_script, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to create '//trim(filename)
            stop 1
        end if

        write(unit_script,'(A)') '#!/usr/bin/env gnuplot'
        write(unit_script,'(A)') '# Generated by sod compare.'
        write(unit_script,'(A)') '# Column 7 stores TDeltaS_kJmol. Columns 3/4/5 store total/low/high DeltaH_kJmol for the system.'
        write(unit_script,'(A)') ''
        write(unit_script,'(A)') 'system_file = "'//trim(system_data_file)//'"'
        write(unit_script,'(A)') 'reference_file = "'//trim(reference_data_file)//'"'
        write(unit_script,'(A)') 'set fit logfile "'//trim(fit_log_file)//'"'
        write(unit_script,'(A)') 'set fit quiet'
        write(unit_script,'(A,F10.4)') 'temperature = ', temperature
        write(unit_script,'(A)') ''
        write(unit_script,'(A)') 'g2(x) = x*(1.0-x)*(a0)'
        write(unit_script,'(A)') 'g3(x) = x*(1.0-x)*(b0 + b1*x)'
        write(unit_script,'(A)') 'g4(x) = x*(1.0-x)*(c0 + c1*x + c2*x*x)'
        write(unit_script,'(A)') 'a0 = 1.0'
        write(unit_script,'(A)') 'b0 = 1.0'
        write(unit_script,'(A)') 'b1 = 1.0'
        write(unit_script,'(A)') 'c0 = 1.0'
        write(unit_script,'(A)') 'c1 = 1.0'
        write(unit_script,'(A)') 'c2 = 1.0'
        write(unit_script,'(A)') 'fit g2(x) reference_file using 2:7 via a0'
        write(unit_script,'(A)') 'std2 = FIT_STDFIT'
        write(unit_script,'(A)') 'fit g3(x) reference_file using 2:7 via b0,b1'
        write(unit_script,'(A)') 'std3 = FIT_STDFIT'
        write(unit_script,'(A)') 'fit g4(x) reference_file using 2:7 via c0,c1,c2'
        write(unit_script,'(A)') 'std4 = FIT_STDFIT'
        write(unit_script,'(A)') ''
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_reference_tds_fit_order2.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(g2($2)) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_reference_tds_fit_order3.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(g3($2)) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_reference_tds_fit_order4.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(g4($2)) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') ''
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_order2.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:($3-$7+g2($2)) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_order3.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:($3-$7+g3($2)) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_order4.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:($3-$7+g4($2)) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') ''
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_low_order2.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($4) < 1.0e20 ? $4-$7+g2($2) : 1/0) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_low_order3.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($4) < 1.0e20 ? $4-$7+g3($2) : 1/0) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_low_order4.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($4) < 1.0e20 ? $4-$7+g4($2) : 1/0) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_high_order2.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($5) < 1.0e20 ? $5-$7+g2($2) : 1/0) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_high_order3.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($5) < 1.0e20 ? $5-$7+g3($2) : 1/0) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') 'set table "'//trim(output_prefix)//'_deltaG_high_order4.dat"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($5) < 1.0e20 ? $5-$7+g4($2) : 1/0) with table'
        write(unit_script,'(A)') 'unset table'
        write(unit_script,'(A)') ''
        write(unit_script,'(A)') 'set terminal pngcairo size 1600,1000 enhanced'
        write(unit_script,'(A)') 'set output "'//trim(plot_file)//'"'
        write(unit_script,'(A)') 'set multiplot layout 3,1 title sprintf("Configurational DeltaG comparison at T = %.1f K", temperature)'
        write(unit_script,'(A)') 'set key outside'
        write(unit_script,'(A)') 'set xlabel "xY"'
        write(unit_script,'(A)') 'set ylabel "TDeltaS_reference (kJ/mol)"'
        write(unit_script,'(A)') 'plot reference_file using 2:7 with points pt 7 ps 1.2 title "'//trim(safe_reference_label)//' TDeltaS data", \'
        write(unit_script,'(A)') '     g2(x) with lines lw 2 title sprintf("'//trim(safe_reference_label)//' TDeltaS fit O2 (std=%.4g)", std2), \'
        write(unit_script,'(A)') '     g3(x) with lines lw 2 title sprintf("'//trim(safe_reference_label)//' TDeltaS fit O3 (std=%.4g)", std3), \'
        write(unit_script,'(A)') '     g4(x) with lines lw 2 title sprintf("'//trim(safe_reference_label)//' TDeltaS fit O4 (std=%.4g)", std4)'
        write(unit_script,'(A)') 'set xlabel "xY"'
        write(unit_script,'(A)') 'set ylabel "DeltaG_system (kJ/mol)"'
        write(unit_script,'(A)') 'plot system_file using 2:($3-$7+g2($2)) with linespoints lw 2 pt 5 title "'//trim(safe_system_label)//' DeltaG O2", \'
        write(unit_script,'(A)') '     system_file using 2:($3-$7+g3($2)) with linespoints lw 2 pt 7 title "'//trim(safe_system_label)//' DeltaG O3", \'
        write(unit_script,'(A)') '     system_file using 2:($3-$7+g4($2)) with linespoints lw 2 pt 9 title "'//trim(safe_system_label)//' DeltaG O4"'
        write(unit_script,'(A)') 'set xlabel "xY"'
        write(unit_script,'(A)') 'set ylabel "Projected DeltaG (kJ/mol)"'
        write(unit_script,'(A)') 'plot system_file using 2:(abs($4) < 1.0e20 ? $4-$7+g2($2) : 1/0) with linespoints lw 2 dt 2 pt 5 title "'//trim(safe_system_label)//' low-branch O2", \'
        write(unit_script,'(A)') '     system_file using 2:(abs($5) < 1.0e20 ? $5-$7+g2($2) : 1/0) with linespoints lw 2 dt 1 pt 5 title "'//trim(safe_system_label)//' high-branch O2", \'
        write(unit_script,'(A)') '     system_file using 2:(abs($4) < 1.0e20 ? $4-$7+g3($2) : 1/0) with linespoints lw 2 dt 2 pt 7 title "'//trim(safe_system_label)//' low-branch O3", \'
        write(unit_script,'(A)') '     system_file using 2:(abs($5) < 1.0e20 ? $5-$7+g3($2) : 1/0) with linespoints lw 2 dt 1 pt 7 title "'//trim(safe_system_label)//' high-branch O3", \'
        write(unit_script,'(A)') '     system_file using 2:(abs($4) < 1.0e20 ? $4-$7+g4($2) : 1/0) with linespoints lw 2 dt 2 pt 9 title "'//trim(safe_system_label)//' low-branch O4", \'
        write(unit_script,'(A)') '     system_file using 2:(abs($5) < 1.0e20 ? $5-$7+g4($2) : 1/0) with linespoints lw 2 dt 1 pt 9 title "'//trim(safe_system_label)//' high-branch O4"'
        write(unit_script,'(A)') 'unset multiplot'
        write(unit_script,'(A)') 'unset output'

        close(unit_script)
    end subroutine write_gnuplot_script

    subroutine parse_mc_summary_row(line, level, frac_y, delta_total, delta_low, delta_high, ok)
        character(len=*), intent(in) :: line
        integer, intent(out) :: level
        real(dp), intent(out) :: frac_y, delta_total, delta_low, delta_high
        logical, intent(out) :: ok
        character(len=128) :: token
        logical :: ok_level, ok_frac, ok_total, ok_low, ok_high

        call extract_delimited_field(line, ';', 1, token)
        call parse_integer_token(token, level, ok_level)
        call extract_delimited_field(line, ';', 2, token)
        call parse_real_token(token, frac_y, ok_frac)
        call extract_delimited_field(line, ';', 11, token)
        call parse_real_token(token, delta_total, ok_total)
        call extract_delimited_field(line, ';', 13, token)
        call parse_optional_real_token(token, delta_low, ok_low)
        call extract_delimited_field(line, ';', 15, token)
        call parse_optional_real_token(token, delta_high, ok_high)

        if (.not. ok_low) delta_low = invalid_projection_marker
        if (.not. ok_high) delta_high = invalid_projection_marker
        ok = ok_level .and. ok_frac .and. ok_total
    end subroutine parse_mc_summary_row

    subroutine extract_delimited_field(line, delimiter, field_index, token)
        character(len=*), intent(in) :: line
        character(len=1), intent(in) :: delimiter
        integer, intent(in) :: field_index
        character(len=*), intent(out) :: token
        integer :: idx, start_pos, end_pos, current_field, line_len

        token = ''
        line_len = len_trim(line)
        if (field_index <= 0 .or. line_len <= 0) return

        start_pos = 1
        current_field = 1
        do idx = 1, line_len + 1
            if (idx > line_len .or. line(idx:idx) == delimiter) then
                if (current_field == field_index) then
                    end_pos = idx - 1
                    if (end_pos >= start_pos) then
                        token = adjustl(line(start_pos:end_pos))
                    else
                        token = ''
                    end if
                    return
                end if
                current_field = current_field + 1
                start_pos = idx + 1
            end if
        end do
    end subroutine extract_delimited_field

    subroutine parse_integer_token(token, value, ok)
        character(len=*), intent(in) :: token
        integer, intent(out) :: value
        logical, intent(out) :: ok
        integer :: ios

        read(token, *, iostat=ios) value
        ok = (ios == 0)
        if (.not. ok) value = 0
    end subroutine parse_integer_token

    subroutine parse_real_token(token, value, ok)
        character(len=*), intent(in) :: token
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        integer :: ios

        read(token, *, iostat=ios) value
        ok = (ios == 0)
        if (.not. ok) value = 0.0_dp
    end subroutine parse_real_token

    subroutine parse_optional_real_token(token, value, ok)
        character(len=*), intent(in) :: token
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        character(len=128) :: lowered

        lowered = adjustl(token)
        call to_lower_inplace(lowered)
        if (len_trim(lowered) == 0 .or. trim(lowered) == '--' .or. trim(lowered) == '-') then
            value = invalid_projection_marker
            ok = .false.
            return
        end if

        call parse_real_token(token, value, ok)
    end subroutine parse_optional_real_token

    function gnuplot_safe_text(text) result(out)
        character(len=*), intent(in) :: text
        character(len=512) :: out
        integer :: idx

        out = ' '
        out(1:min(len_trim(text), len(out))) = text(1:min(len_trim(text), len(out)))
        do idx = 1, len_trim(out)
            if (out(idx:idx) == '"') out(idx:idx) = ''''
        end do
        out = trim(adjustl(out))
    end function gnuplot_safe_text

end module compare
