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
! Exhaustive enumeration driver that prints all minimum-energy configurations
! Controlador de enumeración exhaustiva que imprime todas las configuraciones de energía mínima
! for each substitution level found via the SOD energy calculator.
! para cada nivel de sustitución encontrado mediante el calculador de energía de SOD.
!*******************************************************************************
! Module `exact` performs exhaustive enumeration using the modern effective Hamiltonian.
! El módulo `exact` realiza la enumeración exhaustiva usando el Hamiltoniano efectivo moderno.
module exact
    use consts
    use cli, only: build_level_sequence, classify_template_file, compose_mode_command, &
        engine_name_to_filer, is_help_token, parse_level_spec, to_lower_inplace
    use configurations, only: enumerate_unique_subsets, write_outsod_unit
    use inputs, only: insod_file_data, read_insod_file
    use utils
    use settings, only: blend_low_high_energy_level, reference_relative, reset_blend_overrides
    use symmetry
    use energy_calculations
    use calibration
    use structure_io, only: motif_data_type, read_motif_file
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    private
    public :: run_sod_ensemble_exact

    character(len=512), save :: exact_motif_file = ''
    type(motif_data_type), save :: exact_motif
    logical, save :: exact_have_motif = .false.
    integer, save :: exact_filer = 0

    contains
    subroutine run_sod_ensemble_exact(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: level_min, level_max
        integer, allocatable :: requested_levels(:)
        integer, allocatable :: levels(:)
        real(dp) :: energy_tolerance
        logical :: just_outsod
        logical :: has_level_list
        character(len=512) :: template_gin_option
        character(len=16) :: protocol_option
        character(len=64) :: engine_option
        logical :: external_calibration_enabled
        integer :: nop, total_sites
        integer, pointer :: eqmatrix(:,:)
        integer, allocatable :: scratch_config(:)
        integer :: level_index
        type(insod_file_data) :: insod
        if (present(arg_offset)) then
            call parse_arguments_exact(level_min, level_max, requested_levels, has_level_list, energy_tolerance, &
                just_outsod, template_gin_option, protocol_option, external_calibration_enabled, exact_motif_file, engine_option, arg_offset)
        else
            call parse_arguments_exact(level_min, level_max, requested_levels, has_level_list, energy_tolerance, &
                just_outsod, template_gin_option, protocol_option, external_calibration_enabled, exact_motif_file, engine_option)
        end if
        call set_calibration_template_gin(trim(template_gin_option))
        call set_calibration_protocol_version(trim(protocol_option))
        call set_external_calibration_enabled(external_calibration_enabled)
        call read_insod_file(insod)

        ! Apply engine override
        if (len_trim(engine_option) > 0) then
            exact_filer = engine_name_to_filer(engine_option)
            if (exact_filer < 0) then
                write(error_unit,'(A)') 'Error: unrecognized engine: '//trim(engine_option)
                write(error_unit,'(A)') 'Accepted engines: gulp, lammps, vasp, castep, qe, ase'
                stop 1
            end if
            write(*,'(A,A,A,I0,A)') 'Engine override: --engine ', trim(engine_option), ' (FILER=', exact_filer, ')'
        else
            exact_filer = insod%filer
        end if
        call set_calibration_engine_from_filer(exact_filer)

        exact_have_motif = .false.
        if (len_trim(exact_motif_file) > 0) then
            call read_motif_file(trim(exact_motif_file), exact_motif)
            exact_have_motif = (exact_motif%natoms > 0)
        end if

        call init_energy_calc(skip_energy_files=just_outsod)
        call symmetry_initialize()
        call symmetry_get_matrix(eqmatrix, nop, total_sites)
        if (.not. associated(eqmatrix) .or. total_sites <= 0) then
            write(error_unit,'(A)') 'Error: unable to obtain EQMATRIX or no substitutable sites are available.'
            stop 1
        end if
        nullify(eqmatrix)
        call reset_blend_overrides(total_sites)

        call build_level_sequence(level_min, level_max, requested_levels, has_level_list, total_sites, levels)
        if (.not. allocated(levels) .or. size(levels) == 0) then
            write(error_unit,'(A)') 'Error: the level selection did not produce any valid levels.'
            stop 1
        end if

        allocate(scratch_config(total_sites))
        write(*,'(A)') '--- Exhaustive configuration enumeration ---'
        write(*,'(A,I0)') 'Substitutable sites (npos): ', total_sites
        if (has_level_list) then
            write(*,'(A)', advance='no') 'Requested levels: '
            do level_index = 1, size(levels)
                if (level_index > 1) write(*,'(A)', advance='no') ', '
                write(*,'(I0)', advance='no') levels(level_index)
            end do
            write(*,*)
        else
            write(*,'(A,I0,A,I0)') 'Levels evaluated: ', levels(1), ' .. ', levels(size(levels))
        end if
        write(*,'(A,ES12.5)') 'Energy tolerance for grouping minima (eV): ', energy_tolerance
        if (just_outsod) then
            write(*,'(A)') 'Mode --just-outsod enabled: only xNN/OUTSOD files will be generated.'
            write(*,'(A)') 'ENERGIES, POSCAR_* and energy summaries will be skipped.'
        end if
        if (external_calibration_enabled) then
            write(*,'(A,A)') 'External calibration: enabled with backend ', trim(calibration_engine_name())
        else
            write(*,'(A)') 'External calibration: disabled (default)'
        end if
        write(*,*)
        flush(output_unit)

        do level_index = 1, size(levels)
            call process_level_exact(levels(level_index), total_sites, scratch_config, energy_tolerance, just_outsod)
        end do
        deallocate(scratch_config)
        call symmetry_finalize()
        call cleanup_energy_calc()
    end subroutine run_sod_ensemble_exact

    subroutine parse_arguments_exact(level_min, level_max, level_list, has_level_list, tol, just_outsod, template_gin_option, protocol_option, external_calibration_enabled, motif_file, engine_option, arg_offset)
        integer, intent(out) :: level_min, level_max
        integer, allocatable, intent(out) :: level_list(:)
        logical, intent(out) :: has_level_list
        real(dp), intent(out) :: tol
        logical, intent(out) :: just_outsod
        character(len=*), intent(out) :: template_gin_option
        character(len=*), intent(out) :: protocol_option
        logical, intent(out) :: external_calibration_enabled
        character(len=*), intent(out) :: motif_file
        character(len=*), intent(out) :: engine_option
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, ios
        character(len=256) :: arg, spec, lowered
        integer :: level_candidate
        real(dp) :: tol_candidate
        logical :: level_specified
        integer :: eq_pos
        integer :: offset
        character(len=512) :: tpl_path
        character(len=8) :: tpl_cat
        level_min = 0
        level_max = -1
        has_level_list = .false.
        tol = 1.0e-6_dp
        just_outsod = .false.
        external_calibration_enabled = .false.
        level_specified = .false.
        template_gin_option = 'none'
        protocol_option = '2.0'
        motif_file = ''
        engine_option = ''
        allocate(level_list(0))

        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)
        if (argc <= offset) return

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            if (is_help_token(arg)) then
                call print_usage_exact()
                stop
            end if

            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--engine=') == 1) then
                eq_pos = index(arg, '=')
                engine_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--template=') == 1) then
                eq_pos = index(arg, '=')
                tpl_path = adjustl(arg(eq_pos+1:))
                tpl_cat = classify_template_file(tpl_path)
                if (trim(tpl_cat) == 'lib' .or. trim(tpl_cat) == 'gin') then
                    template_gin_option = tpl_path
                else
                    write(error_unit,'(A)') 'Warning: --template '//trim(tpl_path)//' ignored in exact mode (only .lib/.gin/.include are used).'
                end if
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--template-gin=') == 1 .or. index(trim(lowered), '--template_gin=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid argument for --template-gin.'
                    call print_usage_exact()
                    stop 1
                end if
                template_gin_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--protocol=') == 1 .or. index(trim(lowered), '--protocole=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid argument for --protocol.'
                    call print_usage_exact()
                    stop 1
                end if
                protocol_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (trim(lowered) == '--template-gin' .or. trim(lowered) == '--template_gin') then
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --template-gin.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, template_gin_option)
                template_gin_option = adjustl(template_gin_option)
                iarg = iarg + 2
                cycle
            else if (trim(lowered) == '--protocol' .or. trim(lowered) == '--protocole') then
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --protocol.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, protocol_option)
                protocol_option = adjustl(protocol_option)
                iarg = iarg + 2
                cycle
            else if (trim(lowered) == '--no-template-gin' .or. trim(lowered) == '--skip-template' .or. trim(lowered) == '--skip_template') then
                template_gin_option = 'none'
                iarg = iarg + 1
                cycle
            else if (trim(lowered) == '--calibration' .or. trim(lowered) == '--external-calibration') then
                external_calibration_enabled = .true.
                iarg = iarg + 1
                cycle
            else if (trim(lowered) == '--no-calibration' .or. trim(lowered) == '--no-external-calibration') then
                external_calibration_enabled = .false.
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--motif=') == 1) then
                eq_pos = index(arg, '=')
                motif_file = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (trim(lowered) == '--motif') then
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --motif.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, motif_file)
                motif_file = adjustl(motif_file)
                iarg = iarg + 2
                cycle
            else if (trim(lowered) == '--engine') then
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing engine name after --engine.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, engine_option)
                engine_option = adjustl(engine_option)
                iarg = iarg + 2
                cycle
            else if (trim(lowered) == '--template') then
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --template.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, tpl_path)
                tpl_path = adjustl(tpl_path)
                tpl_cat = classify_template_file(tpl_path)
                if (trim(tpl_cat) == 'lib' .or. trim(tpl_cat) == 'gin') then
                    template_gin_option = tpl_path
                else
                    write(error_unit,'(A)') 'Warning: --template '//trim(tpl_path)//' ignored in exact mode (only .lib/.gin/.include are used).'
                end if
                iarg = iarg + 2
                cycle
            else if (trim(lowered) == '-t' .or. trim(lowered) == '--tolerance') then
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --tolerance.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, spec)
                read(spec, *, iostat=ios) tol_candidate
                if (ios /= 0 .or. tol_candidate <= 0.0_dp) then
                    write(error_unit,'(A)') 'Error: invalid tolerance provided to --tolerance.'
                    call print_usage_exact()
                    stop 1
                end if
                tol = tol_candidate
                iarg = iarg + 2
                cycle
            else if (index(trim(lowered), '--tolerance=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid argument for --tolerance.'
                    call print_usage_exact()
                    stop 1
                end if
                read(arg(eq_pos+1:), *, iostat=ios) tol_candidate
                if (ios /= 0 .or. tol_candidate <= 0.0_dp) then
                    write(error_unit,'(A)') 'Error: invalid tolerance provided to --tolerance.'
                    call print_usage_exact()
                    stop 1
                end if
                tol = tol_candidate
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-N','-n')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing specification after -N.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, spec)
                call parse_level_spec(spec, level_min, level_max, level_list, has_level_list)
                level_specified = .true.
                iarg = iarg + 2
            case ('--just-outsod','--solo-outsod','--only-outsod')
                just_outsod = .true.
                iarg = iarg + 1
            case default
                read(arg, *, iostat=ios) level_candidate
                if (ios == 0 .and. .not. level_specified) then
                    if (level_candidate < 0) then
                        level_min = 0
                        level_max = -1
                    else
                        level_min = 0
                        level_max = level_candidate
                    end if
                    level_specified = .true.
                    iarg = iarg + 1
                else
                    read(arg, *, iostat=ios) tol_candidate
                    if (ios /= 0) then
                        write(error_unit,'(A)') 'Error: unrecognized argument.'
                        call print_usage_exact()
                        stop 1
                    end if
                    if (tol_candidate <= 0.0_dp) then
                        write(error_unit,'(A)') 'Warning: invalid tolerance; using 1e-6 instead.'
                        tol = 1.0e-6_dp
                    else
                        tol = tol_candidate
                    end if
                    iarg = iarg + 1
                end if
            end select
        end do
    end subroutine parse_arguments_exact

    subroutine print_usage_exact()
        character(len=256) :: command_name

        command_name = compose_mode_command('exact')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [-N <spec>] [-t <tol_eV>] [--engine <name>] [--template <file>...] [--just-outsod]'
        write(*,'(A)') '       '//trim(command_name)//' [Nmax] [tol_eV]  (legacy compatibility mode)'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: exact, exhaustive, enumerate, full'
        write(*,'(A)') ''
        write(*,'(A)') 'Optional arguments (legacy positional syntax is preserved for compatibility):'
        write(*,'(A)') '  -N <spec>             Selects the substitution levels to evaluate:'
        write(*,'(A)') '                        -N -1       -> all levels (0..npos)'
        write(*,'(A)') '                        -N 12       -> only level 12'
        write(*,'(A)') '                        -N 1:12     -> levels 1 through 12'
        write(*,'(A)') '                        -N 3,6,9    -> explicit list (order preserved)'
        write(*,'(A)') '  -t, --tolerance <tol> Energy tolerance (eV) used to group minima [1e-6].'
        write(*,'(A)') '  --engine <name>       Select the energy engine (overrides INSOD FILER).'
        write(*,'(A)') '                         gulp, lammps, vasp, castep, qe, ase'
        write(*,'(A)') '  --template <file>     Provide a template/model file (repeatable). See setup --help.'
        write(*,'(A)') '  --just-outsod         Writes only xNN/OUTSOD (no ENERGIES or POSCAR files).'
        write(*,'(A)') '  --template-gin <file> Alias for --template with a .gin/.include file.'
        write(*,'(A)') '  --no-template-gin     Skip template fragments when creating .gin files (default).'
        write(*,'(A)') '  --calibration          Enable external calibration/evaluation with the backend selected by FILER.'
        write(*,'(A)') '  --external-calibration Alias for --calibration. Disabled by default.'
        write(*,'(A)') '  --protocol <ver>      Select the GULP workflow protocol: 1.0 or 2.0 [2.0].'
        write(*,'(A)') '  --motif <file>        Include a molecular motif (.include file) in GULP output (FILER=1).'
        write(*,'(A)') ''
        write(*,'(A)') 'Examples:'
        write(*,'(A)') '  '//trim(command_name)//' -N 5:10 -t 1e-5'
        write(*,'(A)') '  '//trim(command_name)//' -N 8'
        write(*,'(A)') '  '//trim(command_name)//' -N 3,6,9 --just-outsod'
        write(*,'(A)') '  '//trim(command_name)//' 12 1e-6'
        write(*,'(A)') ''
        write(*,'(A)') 'The program exhaustively enumerates every configuration for each requested level'
        write(*,'(A)') 'and stores its auxiliary files inside xNN folders.'
    end subroutine print_usage_exact

    subroutine process_level_exact(level, total_sites, config, tol, just_outsod)
        implicit none
        integer, intent(in) :: level, total_sites
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: tol
        logical, intent(in) :: just_outsod

        integer(ip) :: total_comb
        real(dp) :: energy, low_estimate, high_estimate
        integer :: idx
        integer :: i, nop_local
        integer :: unique_count
        integer :: total_degeneracy_weight
    integer, allocatable :: best_subsets_si(:,:), best_subsets_ge(:,:)
    integer, allocatable :: unique_subsets(:,:)
    integer, allocatable :: unique_deg(:)
    integer, allocatable :: best_deg_si(:), best_deg_ge(:)
    real(dp), allocatable :: best_values_si(:), best_values_ge(:)
    real(dp), allocatable :: unique_low(:), unique_high(:)
    real(dp), allocatable :: unique_low_contrib(:,:), unique_high_contrib(:,:)
    integer, allocatable :: tmp_subsets(:,:)
    integer, allocatable :: tmp_deg(:)
    real(dp), allocatable :: tmp_low(:), tmp_high(:)
    logical :: allow_low_estimate, allow_high_estimate
    logical :: has_low_data, has_high_data
    integer, pointer :: eqmatrix_local(:,:)
    character(len=80), parameter :: separator = '------------------------------------------------------------------------'
    real(dp), parameter :: sentinel_energy = huge(1.0_dp)
    real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
    real(dp), parameter :: temp_targets(3) = (/300.0_dp, 800.0_dp, 1200.0_dp/)
    real(dp), parameter :: boltzmann_temperature = 300.0_dp
        real(dp) :: entropy_total
        real(dp) :: max_entropy, ideal_entropy, x
    real(dp) :: mean_low_all, mean_high_all
    real(dp) :: variance_low_all, variance_high_all
    real(dp) :: weighted_low_sum, weighted_high_sum
    real(dp) :: weighted_low_sq_sum, weighted_high_sq_sum
    real(dp) :: weighted_total_sum, weighted_total_sq_sum
        integer :: total_low_weight, total_high_weight
        integer :: total_total_weight
        real(dp) :: min_low_energy, min_high_energy
        real(dp) :: min_total_energy
        real(dp) :: degeneracy_weight
        real(dp) :: boltzmann_low_weight_sum, boltzmann_high_weight_sum
        real(dp) :: boltzmann_total_weight_sum
        real(dp) :: boltzmann_low_energy_sum, boltzmann_high_energy_sum
        real(dp) :: boltzmann_total_energy_sum
        real(dp) :: boltzmann_low_energy_sq_sum, boltzmann_high_energy_sq_sum
        real(dp) :: boltzmann_total_energy_sq_sum
        real(dp) :: boltzmann_reference_low, boltzmann_reference_high
        real(dp) :: boltzmann_reference_total
        real(dp) :: boltzmann_factor
        real(dp) :: boltzmann_mean_low, boltzmann_variance_low
        real(dp) :: boltzmann_mean_high, boltzmann_variance_high
        real(dp) :: boltzmann_mean_total, boltzmann_variance_total
        real(dp) :: mean_total_all, variance_total_all
        real(dp) :: entropy_low, entropy_high
        real(dp) :: ts_low(3), ts_high(3)
        real(dp) :: free_energy_low(3), free_energy_high(3)
        real(dp) :: deltaF_reference_x, deltaF_reference_y
        integer :: capacity_si, capacity_ge
    integer :: best_count_si, best_count_ge
        real(dp) :: best_energy_si, best_energy_ge
        integer :: hole_count
        logical :: need_low_calibration, need_high_calibration
    real(dp) :: y_fraction
    real(dp) :: expected_mix
    real(dp) :: delta_exp_total, delta_min_total, delta_exp_low, delta_min_low, delta_exp_high, delta_min_high, delta_exp_mix
    real(dp) :: accept_ratio
    real(dp) :: combined_energy
    logical :: valid_low, valid_high, valid_total

        ! Get EQMATRIX for symmetry checking
! Obtener EQMATRIX para comprobar la simetría
        call symmetry_get_matrix(eqmatrix_local, nop_local, i)
        if (.not. associated(eqmatrix_local)) then
            write(*,'(A)') 'Error: EQMATRIX is not available for exact enumeration.'
            return
        end if

        total_comb = binomial_int64(total_sites, level)

        write(*,'(A)') separator
        write(*,'(A,I0)') 'Level: ', level
        write(*,'(A,I0)') 'Total combinations: ', total_comb

        if (total_comb > 0_int64) then
            ! kB_eVk: Boltzmann constant in eV/K; dp: double precision kind
! kB_eVk: constante de Boltzmann en eV/K; dp: tipo de doble precisión
            max_entropy = kB_eVk * log(real(total_comb, dp))
        else
            max_entropy = 0.0_dp
        end if
        write(*,'(A,F16.6)') 'Maximum entropy (k_B ln W) [eV/K]: ', max_entropy

        if (total_sites > 0) then
            x = real(level, dp) / real(total_sites, dp)
            if (x <= 0.0_dp .or. x >= 1.0_dp) then
                ideal_entropy = 0.0_dp
            else
                ! Binary mixing entropy: S = -N k_B [x ln x + (1-x) ln(1-x)]
! Entropía de mezcla binaria: S = -N k_B [x ln x + (1-x) ln(1-x)]
                ! Reference: Kittel, "Introduction to Solid State Physics", mixing entropy equation
! Referencia: Kittel, "Introduction to Solid State Physics", ecuación de la entropía de mezcla
                ideal_entropy = -real(total_sites, dp) * kB_eVk * (x * log(x) + (1.0_dp - x) * log(1.0_dp - x))
            end if
        else
            ideal_entropy = 0.0_dp
        end if
        write(*,'(A,F16.6)') 'Ideal binary-mixing entropy [eV/K]: ', ideal_entropy
        flush(output_unit)

        hole_count = total_sites - level
        allow_low_estimate = .true.
        allow_high_estimate = .true.
        if (hole_count <= get_max_high_order()) allow_low_estimate = .false.
        if (level <= get_max_low_order()) allow_high_estimate = .false.

        if (level == 0 .and. just_outsod) then
            allocate(tmp_subsets(0,1))
            allocate(tmp_deg(1))
            allocate(tmp_low(1))
            allocate(tmp_high(1))
            tmp_deg(1) = 1
            tmp_low(1) = sentinel_energy
            tmp_high(1) = sentinel_energy
            call write_level_outputs(level, total_sites, 1, tmp_subsets, tmp_deg, tmp_low, tmp_high, .false.)
            deallocate(tmp_subsets)
            deallocate(tmp_deg)
            deallocate(tmp_low)
            deallocate(tmp_high)
            nullify(eqmatrix_local)
            return
        end if

        if (level == 0) then
            config = 1
            call calculate_structure_energy(config, total_sites, energy, low_estimate, high_estimate)
            if (.not. allow_low_estimate) low_estimate = sentinel_energy
            if (.not. allow_high_estimate) high_estimate = sentinel_energy

            total_degeneracy_weight = 1

            total_low_weight = 0
            total_high_weight = 0
            mean_low_all = 0.0_dp
            mean_high_all = 0.0_dp
            variance_low_all = 0.0_dp
            variance_high_all = 0.0_dp
            min_low_energy = sentinel_energy
            min_high_energy = sentinel_energy
            boltzmann_mean_low = 0.0_dp
            boltzmann_mean_high = 0.0_dp
            boltzmann_variance_low = 0.0_dp
            boltzmann_variance_high = 0.0_dp
            entropy_total = 0.0_dp
            entropy_low = 0.0_dp
            entropy_high = 0.0_dp
            ts_low = 0.0_dp
            ts_high = 0.0_dp
            free_energy_low = 0.0_dp
            free_energy_high = 0.0_dp
            deltaF_reference_x = 0.0_dp
            deltaF_reference_y = 0.0_dp

            y_fraction = 0.0_dp
            if (total_sites > 0) y_fraction = real(level, dp) / real(total_sites, dp)
            valid_low = allow_low_estimate
            valid_high = allow_high_estimate
            valid_total = valid_low .or. valid_high

            if (valid_low .and. valid_high) then
                combined_energy = blend_low_high_energy_level(level, total_sites, low_estimate, high_estimate)
            else if (valid_low) then
                combined_energy = low_estimate
            else if (valid_high) then
                combined_energy = high_estimate
            else
                combined_energy = 0.0_dp
            end if

            min_total_energy = combined_energy
            mean_total_all = combined_energy
            variance_total_all = 0.0_dp
            boltzmann_mean_total = combined_energy
            boltzmann_variance_total = 0.0_dp
            expected_mix = combined_energy

            if (.not. valid_low) then
                min_low_energy = 0.0_dp
                boltzmann_mean_low = 0.0_dp
            else
                min_low_energy = low_estimate
            end if

            if (.not. valid_high) then
                min_high_energy = 0.0_dp
                boltzmann_mean_high = 0.0_dp
            else
                min_high_energy = high_estimate
            end if

            if (.not. valid_total) then
                min_total_energy = 0.0_dp
                boltzmann_mean_total = 0.0_dp
                expected_mix = 0.0_dp
            end if

            delta_exp_total = merge(reference_relative(level, total_sites, boltzmann_mean_total), 0.0_dp, valid_total)
            delta_min_total = merge(reference_relative(level, total_sites, min_total_energy), 0.0_dp, valid_total)
            delta_exp_low = merge(reference_relative(level, total_sites, boltzmann_mean_low), 0.0_dp, valid_low)
            delta_min_low = merge(reference_relative(level, total_sites, min_low_energy), 0.0_dp, valid_low)
            delta_exp_high = merge(reference_relative(level, total_sites, boltzmann_mean_high), 0.0_dp, valid_high)
            delta_min_high = merge(reference_relative(level, total_sites, min_high_energy), 0.0_dp, valid_high)
            delta_exp_mix = merge(reference_relative(level, total_sites, expected_mix), 0.0_dp, valid_total)
            accept_ratio = 1.0_dp


            if (allow_high_estimate) then
                total_high_weight = 1
                mean_high_all = high_estimate
                min_high_energy = high_estimate
                boltzmann_mean_high = high_estimate
                free_energy_high = (/high_estimate, high_estimate, high_estimate/)
                deltaF_reference_y = reference_relative(level, total_sites, free_energy_high(1))
            end if

            has_low_data = allow_low_estimate
            has_high_data = allow_high_estimate

            allocate(tmp_subsets(0,1))
            allocate(tmp_deg(1))
            allocate(tmp_low(1))
            allocate(tmp_high(1))
            tmp_deg(1) = 1
            if (allow_low_estimate) then
                tmp_low(1) = low_estimate
            else
                tmp_low(1) = sentinel_energy
            end if
            if (allow_high_estimate) then
                tmp_high(1) = high_estimate
            else
                tmp_high(1) = sentinel_energy
            end if
            call write_level_outputs(level, total_sites, 1, tmp_subsets, tmp_deg, tmp_low, tmp_high, .not. just_outsod)
            deallocate(tmp_subsets)
            deallocate(tmp_deg)
            deallocate(tmp_low)
            deallocate(tmp_high)

            if (just_outsod) then
                nullify(eqmatrix_local)
                return
            end if

            allocate(best_subsets_si(0,1))
            allocate(best_subsets_ge(0,1))
            allocate(best_values_si(1))
            allocate(best_values_ge(1))
            allocate(best_deg_si(1))
            allocate(best_deg_ge(1))

            if (allow_low_estimate) then
                best_values_si(1) = low_estimate
                best_deg_si(1) = 1
                best_energy_si = low_estimate
                best_count_si = 1
            else
                best_values_si(1) = 0.0_dp
                best_deg_si(1) = 0
                best_energy_si = sentinel_energy
                best_count_si = 0
            end if

            if (allow_high_estimate) then
                best_values_ge(1) = high_estimate
                best_deg_ge(1) = 1
                best_energy_ge = high_estimate
                best_count_ge = 1
            else
                best_values_ge(1) = 0.0_dp
                best_deg_ge(1) = 0
                best_energy_ge = sentinel_energy
                best_count_ge = 0
            end if

            call emit_side_statistics('X side', 'X', level, best_count_si, best_energy_si, best_subsets_si, best_values_si, best_deg_si, total_sites, config, min_low_energy, mean_low_all, variance_low_all, entropy_low, has_low_data)
            call emit_side_statistics('Y side', 'Y', level, best_count_ge, best_energy_ge, best_subsets_ge, best_values_ge, best_deg_ge, total_sites, config, min_high_energy, mean_high_all, variance_high_all, entropy_high, has_high_data)

            call append_normalized_summary('sod_ensemble_exact.txt', level, y_fraction, boltzmann_mean_total, min_total_energy, &
                 boltzmann_variance_total, boltzmann_mean_low, min_low_energy, boltzmann_mean_high, min_high_energy, expected_mix, &
                 delta_exp_total, delta_min_total, delta_exp_low, delta_min_low, delta_exp_high, delta_min_high, delta_exp_mix, &
                 accept_ratio)

            deallocate(best_subsets_si)
            deallocate(best_subsets_ge)
            deallocate(best_values_si)
            deallocate(best_values_ge)
            deallocate(best_deg_si)
            deallocate(best_deg_ge)
            nullify(eqmatrix_local)
            return
        end if

        if (.not. just_outsod) then
            capacity_si = max(8, level)
            capacity_ge = max(8, level)
            allocate(best_subsets_si(level, capacity_si))
            allocate(best_subsets_ge(level, capacity_ge))
            allocate(best_values_si(capacity_si))
            allocate(best_values_ge(capacity_ge))
            allocate(best_deg_si(capacity_si))
            allocate(best_deg_ge(capacity_ge))
            best_values_si = 0.0_dp
            best_values_ge = 0.0_dp
            best_deg_si = 0
            best_deg_ge = 0
            best_count_si = 0
            best_count_ge = 0
            best_energy_si = huge(1.0_dp)
            best_energy_ge = huge(1.0_dp)
        end if

        call enumerate_unique_subsets(total_sites, level, eqmatrix_local, nop_local, unique_subsets, unique_deg, unique_count)

        allocate(unique_low(unique_count))
        allocate(unique_high(unique_count))
        allocate(unique_low_contrib(4, unique_count))
        allocate(unique_high_contrib(4, unique_count))
        unique_low = sentinel_energy
        unique_high = sentinel_energy
        unique_low_contrib = 0.0_dp
        unique_high_contrib = 0.0_dp
        config = 1

        write(*,'(A,I0)') 'Unique configurations evaluated: ', unique_count
        flush(output_unit)

        if (just_outsod) then
            call write_level_outputs(level, total_sites, unique_count, unique_subsets, unique_deg, unique_low, unique_high, .false.)
        else
            do idx = 1, unique_count
                if (unique_low(idx) == sentinel_energy) then
                    config = 1
                    config(unique_subsets(1:level, idx)) = 2
                    call calculate_structure_energy(config, total_sites, energy, energy_low_side=low_estimate, &
                         energy_high_side=high_estimate, low_contrib=unique_low_contrib(:, idx), &
                         high_contrib=unique_high_contrib(:, idx))
                    if (.not. allow_low_estimate) then
                        low_estimate = sentinel_energy
                        unique_low_contrib(:, idx) = 0.0_dp
                    end if
                    if (.not. allow_high_estimate) then
                        high_estimate = sentinel_energy
                        unique_high_contrib(:, idx) = 0.0_dp
                    end if
                    unique_low(idx) = low_estimate
                    unique_high(idx) = high_estimate
                end if
            end do

            hole_count = total_sites - level
            need_low_calibration = allow_low_estimate .and. (level > get_max_low_order())
            need_high_calibration = allow_high_estimate .and. (hole_count > get_max_high_order())

            if (need_low_calibration .and. need_high_calibration) then
                call calibrate_level_with_engine(level, total_sites, unique_count, unique_subsets, unique_low_contrib, &
                    unique_high_contrib, unique_low, unique_high, config, need_low_calibration, need_high_calibration, 'x')
            end if

            y_fraction = 0.0_dp
            if (total_sites > 0) y_fraction = real(level, dp) / real(total_sites, dp)
            weighted_low_sum = 0.0_dp
            weighted_high_sum = 0.0_dp
            weighted_low_sq_sum = 0.0_dp
            weighted_high_sq_sum = 0.0_dp
            weighted_total_sum = 0.0_dp
            weighted_total_sq_sum = 0.0_dp
            boltzmann_low_weight_sum = 0.0_dp
            boltzmann_high_weight_sum = 0.0_dp
            boltzmann_total_weight_sum = 0.0_dp
            boltzmann_low_energy_sum = 0.0_dp
            boltzmann_high_energy_sum = 0.0_dp
            boltzmann_total_energy_sum = 0.0_dp
            boltzmann_low_energy_sq_sum = 0.0_dp
            boltzmann_high_energy_sq_sum = 0.0_dp
            boltzmann_total_energy_sq_sum = 0.0_dp
            total_degeneracy_weight = 0
            total_low_weight = 0
            total_high_weight = 0
            total_total_weight = 0
            min_low_energy = huge_marker
            min_high_energy = huge_marker
            min_total_energy = huge_marker

            do idx = 1, unique_count
                degeneracy_weight = real(unique_deg(idx), dp)
                total_degeneracy_weight = total_degeneracy_weight + unique_deg(idx)

                if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                    total_low_weight = total_low_weight + unique_deg(idx)
                    weighted_low_sum = weighted_low_sum + degeneracy_weight * unique_low(idx)
                    weighted_low_sq_sum = weighted_low_sq_sum + degeneracy_weight * unique_low(idx) * unique_low(idx)
                    if (unique_low(idx) < min_low_energy) min_low_energy = unique_low(idx)
                end if

                if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                    total_high_weight = total_high_weight + unique_deg(idx)
                    weighted_high_sum = weighted_high_sum + degeneracy_weight * unique_high(idx)
                    weighted_high_sq_sum = weighted_high_sq_sum + degeneracy_weight * unique_high(idx) * unique_high(idx)
                    if (unique_high(idx) < min_high_energy) min_high_energy = unique_high(idx)
                end if

                combined_energy = blend_low_high_energy_level(level, total_sites, unique_low(idx), unique_high(idx))
                if (combined_energy < huge_marker) then
                    total_total_weight = total_total_weight + unique_deg(idx)
                    weighted_total_sum = weighted_total_sum + degeneracy_weight * combined_energy
                    weighted_total_sq_sum = weighted_total_sq_sum + degeneracy_weight * combined_energy * combined_energy
                    if (combined_energy < min_total_energy) min_total_energy = combined_energy
                end if
            end do

            if (total_low_weight > 0) then
                mean_low_all = weighted_low_sum / real(total_low_weight, dp)
                variance_low_all = max(0.0_dp, weighted_low_sq_sum / real(total_low_weight, dp) - mean_low_all * mean_low_all)
            else
                mean_low_all = 0.0_dp
                variance_low_all = 0.0_dp
                min_low_energy = huge_marker
            end if

            if (total_high_weight > 0) then
                mean_high_all = weighted_high_sum / real(total_high_weight, dp)
                variance_high_all = max(0.0_dp, weighted_high_sq_sum / real(total_high_weight, dp) - mean_high_all * mean_high_all)
            else
                mean_high_all = 0.0_dp
                variance_high_all = 0.0_dp
                min_high_energy = huge_marker
            end if

            if (total_total_weight > 0) then
                mean_total_all = weighted_total_sum / real(total_total_weight, dp)
                variance_total_all = max(0.0_dp, weighted_total_sq_sum / real(total_total_weight, dp) - mean_total_all * mean_total_all)
            else
                mean_total_all = 0.0_dp
                variance_total_all = 0.0_dp
                min_total_energy = huge_marker
            end if

            call write_level_outputs(level, total_sites, unique_count, unique_subsets, unique_deg, unique_low, unique_high, .true.)

            boltzmann_mean_low = mean_low_all
            boltzmann_variance_low = variance_low_all
            boltzmann_mean_high = mean_high_all
            boltzmann_variance_high = variance_high_all
            boltzmann_mean_total = mean_total_all
            boltzmann_variance_total = variance_total_all

            if (total_low_weight > 0) then
                boltzmann_reference_low = min_low_energy
                if (boltzmann_reference_low >= huge_marker) boltzmann_reference_low = mean_low_all
                boltzmann_low_weight_sum = 0.0_dp
                boltzmann_low_energy_sum = 0.0_dp
                boltzmann_low_energy_sq_sum = 0.0_dp
                do idx = 1, unique_count
                    if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                        degeneracy_weight = real(unique_deg(idx), dp)
                        boltzmann_factor = exp(-(unique_low(idx) - boltzmann_reference_low) / (kB_eVk * boltzmann_temperature))
                        boltzmann_low_weight_sum = boltzmann_low_weight_sum + degeneracy_weight * boltzmann_factor
                        boltzmann_low_energy_sum = boltzmann_low_energy_sum + degeneracy_weight * boltzmann_factor * unique_low(idx)
                        boltzmann_low_energy_sq_sum = boltzmann_low_energy_sq_sum + degeneracy_weight * boltzmann_factor * unique_low(idx) * unique_low(idx)
                    end if
                end do
                if (boltzmann_low_weight_sum > 0.0_dp) then
                    boltzmann_mean_low = boltzmann_low_energy_sum / boltzmann_low_weight_sum
                    boltzmann_variance_low = max(0.0_dp, boltzmann_low_energy_sq_sum / boltzmann_low_weight_sum - boltzmann_mean_low * boltzmann_mean_low)
                end if
            end if

            if (total_high_weight > 0) then
                boltzmann_reference_high = min_high_energy
                if (boltzmann_reference_high >= huge_marker) boltzmann_reference_high = mean_high_all
                boltzmann_high_weight_sum = 0.0_dp
                boltzmann_high_energy_sum = 0.0_dp
                boltzmann_high_energy_sq_sum = 0.0_dp
                do idx = 1, unique_count
                    if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                        degeneracy_weight = real(unique_deg(idx), dp)
                        boltzmann_factor = exp(-(unique_high(idx) - boltzmann_reference_high) / (kB_eVk * boltzmann_temperature))
                        boltzmann_high_weight_sum = boltzmann_high_weight_sum + degeneracy_weight * boltzmann_factor
                        boltzmann_high_energy_sum = boltzmann_high_energy_sum + degeneracy_weight * boltzmann_factor * unique_high(idx)
                        boltzmann_high_energy_sq_sum = boltzmann_high_energy_sq_sum + degeneracy_weight * boltzmann_factor * unique_high(idx) * unique_high(idx)
                    end if
                end do
                if (boltzmann_high_weight_sum > 0.0_dp) then
                    boltzmann_mean_high = boltzmann_high_energy_sum / boltzmann_high_weight_sum
                    boltzmann_variance_high = max(0.0_dp, boltzmann_high_energy_sq_sum / boltzmann_high_weight_sum - boltzmann_mean_high * boltzmann_mean_high)
                end if
            end if

            if (total_total_weight > 0) then
                boltzmann_reference_total = min_total_energy
                if (boltzmann_reference_total >= huge_marker) boltzmann_reference_total = mean_total_all
                boltzmann_total_weight_sum = 0.0_dp
                boltzmann_total_energy_sum = 0.0_dp
                boltzmann_total_energy_sq_sum = 0.0_dp
                do idx = 1, unique_count
                    degeneracy_weight = real(unique_deg(idx), dp)
                    combined_energy = blend_low_high_energy_level(level, total_sites, unique_low(idx), unique_high(idx))
                    if (combined_energy < huge_marker) then
                        boltzmann_factor = exp(-(combined_energy - boltzmann_reference_total) / (kB_eVk * boltzmann_temperature))
                        boltzmann_total_weight_sum = boltzmann_total_weight_sum + degeneracy_weight * boltzmann_factor
                        boltzmann_total_energy_sum = boltzmann_total_energy_sum + degeneracy_weight * boltzmann_factor * combined_energy
                        boltzmann_total_energy_sq_sum = boltzmann_total_energy_sq_sum + degeneracy_weight * boltzmann_factor * combined_energy * combined_energy
                    end if
                end do
                if (boltzmann_total_weight_sum > 0.0_dp) then
                    boltzmann_mean_total = boltzmann_total_energy_sum / boltzmann_total_weight_sum
                    boltzmann_variance_total = max(0.0_dp, boltzmann_total_energy_sq_sum / boltzmann_total_weight_sum - boltzmann_mean_total * boltzmann_mean_total)
                end if
            end if

            if (total_low_weight > 0) then
                entropy_low = kB_eVk * log(real(total_low_weight, dp))
                ts_low = entropy_low * temp_targets
                free_energy_low = boltzmann_mean_low - ts_low
                deltaF_reference_x = reference_relative(level, total_sites, free_energy_low(1))
            else
                entropy_low = 0.0_dp
                ts_low = 0.0_dp
                free_energy_low = 0.0_dp
                deltaF_reference_x = 0.0_dp
            end if

            if (total_high_weight > 0) then
                entropy_high = kB_eVk * log(real(total_high_weight, dp))
                ts_high = entropy_high * temp_targets
                free_energy_high = boltzmann_mean_high - ts_high
                deltaF_reference_y = reference_relative(level, total_sites, free_energy_high(1))
            else
                entropy_high = 0.0_dp
                ts_high = 0.0_dp
                free_energy_high = 0.0_dp
                deltaF_reference_y = 0.0_dp
            end if

            if (total_degeneracy_weight > 0) then
                entropy_total = kB_eVk * log(real(total_degeneracy_weight, dp))
            else
                entropy_total = 0.0_dp
            end if

            if (min_low_energy >= huge_marker) min_low_energy = mean_low_all
            if (min_high_energy >= huge_marker) min_high_energy = mean_high_all
            if (min_total_energy >= huge_marker) min_total_energy = mean_total_all

            if (total_low_weight > 0 .and. total_high_weight > 0) then
                expected_mix = blend_low_high_energy_level(level, total_sites, boltzmann_mean_low, boltzmann_mean_high)
            else if (total_low_weight > 0) then
                expected_mix = boltzmann_mean_low
            else if (total_high_weight > 0) then
                expected_mix = boltzmann_mean_high
            else
                expected_mix = boltzmann_mean_total
            end if

            valid_low = (total_low_weight > 0)
            valid_high = (total_high_weight > 0)
            valid_total = (total_total_weight > 0)
            delta_exp_total = merge(reference_relative(level, total_sites, boltzmann_mean_total), 0.0_dp, valid_total)
            delta_min_total = merge(reference_relative(level, total_sites, min_total_energy), 0.0_dp, valid_total)
            delta_exp_low = merge(reference_relative(level, total_sites, boltzmann_mean_low), 0.0_dp, valid_low)
            delta_min_low = merge(reference_relative(level, total_sites, min_low_energy), 0.0_dp, valid_low)
            delta_exp_high = merge(reference_relative(level, total_sites, boltzmann_mean_high), 0.0_dp, valid_high)
            delta_min_high = merge(reference_relative(level, total_sites, min_high_energy), 0.0_dp, valid_high)
            delta_exp_mix = reference_relative(level, total_sites, expected_mix)
            accept_ratio = 1.0_dp

            do idx = 1, unique_count
                if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                    call update_best_structures_side(level, unique_subsets(1:level, idx), unique_low(idx), unique_deg(idx), tol, best_energy_si, best_subsets_si, best_values_si, best_deg_si, best_count_si, capacity_si, eqmatrix_local, nop_local)
                end if
                if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                    call update_best_structures_side(level, unique_subsets(1:level, idx), unique_high(idx), unique_deg(idx), tol, best_energy_ge, best_subsets_ge, best_values_ge, best_deg_ge, best_count_ge, capacity_ge, eqmatrix_local, nop_local)
                end if
            end do

            has_low_data = (total_low_weight > 0)
            has_high_data = (total_high_weight > 0)

            call emit_side_statistics('X side', 'X', level, best_count_si, best_energy_si, best_subsets_si, best_values_si, best_deg_si, total_sites, config, min_low_energy, mean_low_all, variance_low_all, entropy_low, has_low_data)
            call emit_side_statistics('Y side', 'Y', level, best_count_ge, best_energy_ge, best_subsets_ge, best_values_ge, best_deg_ge, total_sites, config, min_high_energy, mean_high_all, variance_high_all, entropy_high, has_high_data)

            call append_normalized_summary('sod_ensemble_exact.txt', level, y_fraction, boltzmann_mean_total, min_total_energy, &
                 boltzmann_variance_total, boltzmann_mean_low, min_low_energy, boltzmann_mean_high, min_high_energy, expected_mix, &
                 delta_exp_total, delta_min_total, delta_exp_low, delta_min_low, delta_exp_high, delta_min_high, delta_exp_mix, &
                 accept_ratio)
        end if

    if (allocated(unique_subsets)) deallocate(unique_subsets)
    if (allocated(unique_low)) deallocate(unique_low)
    if (allocated(unique_high)) deallocate(unique_high)
    if (allocated(unique_deg)) deallocate(unique_deg)
    if (allocated(unique_low_contrib)) deallocate(unique_low_contrib)
    if (allocated(unique_high_contrib)) deallocate(unique_high_contrib)
    if (allocated(best_subsets_si)) deallocate(best_subsets_si)
    if (allocated(best_subsets_ge)) deallocate(best_subsets_ge)
    if (allocated(best_values_si)) deallocate(best_values_si)
    if (allocated(best_values_ge)) deallocate(best_values_ge)
    if (allocated(best_deg_si)) deallocate(best_deg_si)
    if (allocated(best_deg_ge)) deallocate(best_deg_ge)
    nullify(eqmatrix_local)
    end subroutine process_level_exact

    subroutine emit_side_statistics(label, tag, level, count, best_energy, best_subsets, best_values, best_deg, total_sites, config, min_energy, mean_energy, variance_energy, entropy_side, has_data)
        implicit none
        character(len=*), intent(in) :: label
        character(len=*), intent(in) :: tag
        integer, intent(in) :: level, count, total_sites
        real(dp), intent(in) :: best_energy, min_energy, mean_energy, variance_energy, entropy_side
        integer, intent(inout) :: config(:)
        integer, intent(in) :: best_subsets(:,:)
        real(dp), intent(in) :: best_values(:)
        integer, intent(in) :: best_deg(:)
        logical, intent(in) :: has_data
        integer :: idx
        integer :: total_deg
        real(dp) :: entropy_min

        write(*,'(A)') trim(label)//' - statistics'
        if (.not. has_data) then
            write(*,'(A)') '  Minimum energy (eV): -'
            write(*,'(A)') '  Mean energy (eV): -'
            write(*,'(A)') '  Energy variance [eV^2]: -'
            write(*,'(A)') '  Reference minimum energy (eV): -'
            write(*,'(A)') '  Minimum-energy configurations: 0'
            write(*,'(A)') '  Total degeneracy of the minima: 0'
            write(*,'(A)') '  Entropy of the minima (eV/K): -'
            write(*,'(A)') '  Total side entropy ('//trim(tag)//') (eV/K): -'
            flush(output_unit)
            return
        end if

        if (count <= 0) then
            write(*,'(A)') '  No measurable configurations are available for this side.'
            flush(output_unit)
            return
        end if

        write(*,'(A,F16.6)') '  Minimum energy (eV): ', best_energy
        write(*,'(A,F16.6)') '  Mean energy (eV): ', mean_energy
        write(*,'(A,F16.6)') '  Energy variance [eV^2]: ', variance_energy
        write(*,'(A,F16.6)') '  Reference minimum energy (eV): ', min_energy

        total_deg = 0
        do idx = 1, count
            total_deg = total_deg + best_deg(idx)
        end do
        write(*,'(A,I0)') '  Minimum-energy configurations: ', count
        write(*,'(A,I0)') '  Total degeneracy of the minima: ', total_deg
        if (total_deg > 0) then
            entropy_min = kB_eVk * log(real(total_deg, dp))
        else
            entropy_min = 0.0_dp
        end if
        write(*,'(A,F16.6)') '  Entropy of the minima (eV/K): ', entropy_min
        write(*,'(A,F16.6)') '  Total side entropy ('//trim(tag)//') (eV/K): ', entropy_side
        flush(output_unit)

        do idx = 1, count
            config = 1
            if (level > 0) config(best_subsets(1:level, idx)) = 2
            call emit_configuration_info_side(trim(label), idx, best_values(idx), best_subsets(1:level, idx), best_deg(idx))
            call write_exact_poscar_with_tag(level, idx, config, total_sites, trim(tag))
        end do
    end subroutine emit_side_statistics

    subroutine emit_configuration_info_side(label, index, energy, subset, degeneracy)
        implicit none
        character(len=*), intent(in) :: label
        integer, intent(in) :: index
        real(dp), intent(in) :: energy
        integer, intent(in) :: subset(:)
        integer, intent(in) :: degeneracy

        write(*,'(A,I0)') '  '//trim(label)//' - configuration ', index
        if (size(subset) > 0) then
            write(*,'(A)', advance='no') '    Y sites: '
            if (size(subset) > 0) write(*,'( *(1X,I0) )') subset
        else
            write(*,'(A)') '    Base configuration without substitutions.'
        end if
        write(*,'(A,I0)') '    Degeneracy (g): ', degeneracy
        write(*,'(A,F16.6)') '    Corrected energy (eV): ', energy
        flush(output_unit)
    end subroutine emit_configuration_info_side

    subroutine write_exact_poscar_with_tag(level, index, config, total_sites, tag)
        implicit none
        integer, intent(in) :: level, index, total_sites
        integer, intent(in) :: config(:)
        character(len=*), intent(in) :: tag
        character(len=512) :: filename
        character(len=32) :: level_dir
        character(len=16) :: config_tag

        level_dir = format_level_directory('x', level)
        call ensure_directory_exists(trim(level_dir))
        write(filename,'("POSCAR_",A,"_N",I4.4,"_",I4.4,".vasp")') trim(tag), level, index
        filename = join_path(trim(level_dir), trim(filename))
        if (exact_have_motif) then
            call write_vasp_file(config, total_sites, trim(filename), exact_motif%atoms, exact_motif%natoms)
        else
            call write_vasp_file(config, total_sites, trim(filename))
        end if
        write(*,'(A)') '    POSCAR written to '//trim(filename)

        if (exact_filer == 1) then
            write(filename,'(A,"_N",I4.4,"_",I4.4,".gin")') trim(tag), level, index
            filename = join_path(trim(level_dir), trim(filename))
            write(config_tag,'("c",I5.5)') index
            if (exact_have_motif) then
                call write_gulp_output_file(config, total_sites, trim(filename), 'template_input.gin', &
                    trim(config_tag), motif=exact_motif)
            else
                call write_gulp_output_file(config, total_sites, trim(filename), 'template_input.gin', &
                    trim(config_tag))
            end if
            write(*,'(A)') '    GULP file written to '//trim(filename)
        end if

        if (exact_filer == 14) then
            write(filename,'("POSCAR_",A,"_N",I4.4,"_",I4.4,".vasp")') trim(tag), level, index
            filename = join_path(trim(level_dir), trim(filename))
            write(*,'(A)') '    ASE VASP file written to '//trim(filename)
        end if
        flush(output_unit)
    end subroutine write_exact_poscar_with_tag

    subroutine write_level_outputs(level, total_sites, unique_count, unique_subsets, unique_deg, unique_low, unique_high, write_energy_file)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_subsets(:,:), unique_deg(:)
        real(dp), intent(in) :: unique_low(:), unique_high(:)
        logical, intent(in) :: write_energy_file
        character(len=512) :: outsod_name, energies_name
        character(len=32) :: level_dir
        character(len=32) :: low_field, high_field, mix_field
        integer :: unit_outsod, unit_energy
        integer :: idx
        real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
        real(dp) :: mixed_energy
        logical :: valid_low, valid_high
        logical :: energy_open

        level_dir = format_level_directory('x', level)
        call ensure_directory_exists(trim(level_dir))
        outsod_name = join_path(trim(level_dir), 'OUTSOD')
        energies_name = join_path(trim(level_dir), 'ENERGIES')

        open(newunit=unit_outsod, file=trim(outsod_name), status='replace', action='write')
        energy_open = .false.
        if (write_energy_file) then
            open(newunit=unit_energy, file=trim(energies_name), status='replace', action='write')
            energy_open = .true.
        end if

        call write_outsod_unit(unit_outsod, level, total_sites, unique_count, unique_deg, unique_subsets)
        if (energy_open) write(unit_energy,'(A)') '# idx degeneracy energy_Si_eV energy_Ge_eV energy_mix_eV'

        if (unique_count <= 0) then
            close(unit_outsod)
            if (energy_open) close(unit_energy)
            return
        end if

        do idx = 1, unique_count
            if (energy_open) then
                valid_low = abs(unique_low(idx)) < huge_marker
                valid_high = abs(unique_high(idx)) < huge_marker

                if (valid_low) then
                    write(low_field,'(F18.8)') unique_low(idx)
                else
                    low_field = 'NaN'
                end if
                if (valid_high) then
                    write(high_field,'(F18.8)') unique_high(idx)
                else
                    high_field = 'NaN'
                end if

                if (valid_low .and. valid_high) then
                    mixed_energy = blend_low_high_energy_level(level, total_sites, unique_low(idx), unique_high(idx))
                    write(mix_field,'(F18.8)') mixed_energy
                else if (valid_low) then
                    write(mix_field,'(F18.8)') unique_low(idx)
                else if (valid_high) then
                    write(mix_field,'(F18.8)') unique_high(idx)
                else
                    mix_field = 'NaN'
                end if

                write(unit_energy,'(I6,1X,I12,3(1X,A))') idx, unique_deg(idx), adjustl(low_field), adjustl(high_field), adjustl(mix_field)
            end if
        end do

        close(unit_outsod)
        if (energy_open) close(unit_energy)
    end subroutine write_level_outputs

    subroutine append_normalized_summary(filename, level, frac_ge, e_exp_total, e_min_total, var_total, e_exp_low, e_min_low, &
            e_exp_high, e_min_high, e_exp_mix, delta_exp_total, delta_min_total, delta_exp_low, delta_min_low, delta_exp_high, &
            delta_min_high, delta_exp_mix, ratio)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: level
        real(dp), intent(in) :: frac_ge, e_exp_total, e_min_total, var_total
        real(dp), intent(in) :: e_exp_low, e_min_low, e_exp_high, e_min_high, e_exp_mix
        real(dp), intent(in) :: delta_exp_total, delta_min_total, delta_exp_low, delta_min_low
        real(dp), intent(in) :: delta_exp_high, delta_min_high, delta_exp_mix, ratio
        logical :: summary_exists
        integer :: unit_summary, ios
        real(dp) :: frac_val, exp_total_val, min_total_val, var_total_val
        real(dp) :: exp_low_val, min_low_val, exp_high_val, min_high_val
        real(dp) :: exp_mix_val
        real(dp) :: delta_exp_total_val, delta_min_total_val
        real(dp) :: delta_exp_low_val, delta_min_low_val
        real(dp) :: delta_exp_high_val, delta_min_high_val
        real(dp) :: delta_exp_mix_val, ratio_val

        frac_val = max(0.0_dp, min(1.0_dp, sanitize_value(frac_ge)))
        exp_total_val = sanitize_value(e_exp_total)
        min_total_val = sanitize_value(e_min_total)
        var_total_val = max(0.0_dp, sanitize_value(var_total))
        exp_low_val = sanitize_value(e_exp_low)
        min_low_val = sanitize_value(e_min_low)
        exp_high_val = sanitize_value(e_exp_high)
        min_high_val = sanitize_value(e_min_high)
        exp_mix_val = sanitize_value(e_exp_mix)
        delta_exp_total_val = sanitize_value(delta_exp_total)
        delta_min_total_val = sanitize_value(delta_min_total)
        delta_exp_low_val = sanitize_value(delta_exp_low)
        delta_min_low_val = sanitize_value(delta_min_low)
        delta_exp_high_val = sanitize_value(delta_exp_high)
        delta_min_high_val = sanitize_value(delta_min_high)
        delta_exp_mix_val = sanitize_value(delta_exp_mix)
        ratio_val = max(0.0_dp, min(1.0_dp, sanitize_value(ratio)))

        inquire(file=filename, exist=summary_exists)
        open(newunit=unit_summary, file=filename, status='unknown', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') 'Warning: failed to write the normalized summary.'
            flush(output_unit)
            return
        end if

        if (.not. summary_exists) then
            write(unit_summary,'(A)') '#N FracY E_exp_total E_min_total Var_total E_exp_X_side E_min_X_side ' // &
                'E_exp_Y_side E_min_Y_side E_exp_combined Delta_exp_total Delta_min_total ' // &
                'Delta_exp_X_side Delta_min_X_side Delta_exp_Y_side Delta_min_Y_side ' // &
                'Delta_exp_combined Acceptance_ratio'
        end if

        write(unit_summary,'(I0,1X,F7.4,16(1X,F12.6))') level, frac_val, exp_total_val, min_total_val, var_total_val, &
            exp_low_val, min_low_val, exp_high_val, min_high_val, exp_mix_val, delta_exp_total_val, delta_min_total_val, &
            delta_exp_low_val, delta_min_low_val, delta_exp_high_val, delta_min_high_val, delta_exp_mix_val, ratio_val

        close(unit_summary)
    end subroutine append_normalized_summary

    pure real(dp) function sanitize_value(value) result(clean_val)
        implicit none
        real(dp), intent(in) :: value
        real(dp), parameter :: limit = huge(1.0_dp) * 0.5_dp

        if (abs(value) >= limit .or. value /= value) then
            clean_val = 0.0_dp
        else
            clean_val = value
        end if
    end function sanitize_value

    subroutine update_best_structures_side(level, subset, energy, degeneracy, tol, best_energy, best_subsets, best_values, best_deg, best_count, capacity, eqmatrix, nop)
        implicit none
        integer, intent(in) :: level, degeneracy, nop
        integer, intent(in) :: subset(:)
        real(dp), intent(in) :: energy, tol
        real(dp), intent(inout) :: best_energy
        integer, allocatable, intent(inout) :: best_subsets(:,:)
        real(dp), allocatable, intent(inout) :: best_values(:)
        integer, allocatable, intent(inout) :: best_deg(:)
        integer, intent(inout) :: best_count, capacity
        integer, intent(in) :: eqmatrix(:,:)
        real(dp) :: diff
        logical :: equivalent

        diff = energy - best_energy

        if (best_count == 0 .or. diff < -tol) then
            best_energy = energy
            best_count = 1
            call ensure_side_capacity(level, best_subsets, best_values, best_deg, capacity, best_count)
            if (level > 0) best_subsets(1:level,1) = subset(1:level)
            best_values(1) = energy
            best_deg(1) = degeneracy
            return
        end if

        if (abs(diff) <= tol) then
            if (best_count > 0) then
                equivalent = is_symmetry_equivalent(subset, level, best_subsets, best_count, eqmatrix, nop)
                if (equivalent) return
            end if

            best_count = best_count + 1
            call ensure_side_capacity(level, best_subsets, best_values, best_deg, capacity, best_count)
            if (level > 0) best_subsets(1:level, best_count) = subset(1:level)
            best_values(best_count) = energy
            best_deg(best_count) = degeneracy
        end if
    end subroutine update_best_structures_side

    subroutine ensure_side_capacity(level, best_subsets, best_values, best_deg, capacity, required)
        implicit none
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: best_subsets(:,:)
        real(dp), allocatable, intent(inout) :: best_values(:)
        integer, allocatable, intent(inout) :: best_deg(:)
        integer :: old_capacity, new_capacity
        integer, allocatable :: tmp_subsets(:,:)
        real(dp), allocatable :: tmp_values(:)
        integer, allocatable :: tmp_int(:)

        if (required <= capacity) return

        old_capacity = capacity
        new_capacity = max(required, max(2*max(1, old_capacity), 8))

        if (level > 0) then
            allocate(tmp_subsets(level, new_capacity))
            if (old_capacity > 0) tmp_subsets(:,1:old_capacity) = best_subsets(:,1:old_capacity)
            if (allocated(best_subsets)) deallocate(best_subsets)
            call move_alloc(tmp_subsets, best_subsets)
        end if

        allocate(tmp_values(new_capacity))
        if (old_capacity > 0) tmp_values(1:old_capacity) = best_values(1:old_capacity)
        if (allocated(best_values)) deallocate(best_values)
        call move_alloc(tmp_values, best_values)

        allocate(tmp_int(new_capacity))
        if (old_capacity > 0) tmp_int(1:old_capacity) = best_deg(1:old_capacity)
        if (allocated(best_deg)) deallocate(best_deg)
        call move_alloc(tmp_int, best_deg)

        capacity = new_capacity
    end subroutine ensure_side_capacity

    logical function is_symmetry_equivalent(subset, level, best_subsets, best_count, eqmatrix, nop)
        implicit none
        integer, intent(in) :: level, best_count, nop
        integer, intent(in) :: subset(:)
        integer, intent(in) :: best_subsets(:,:)
        integer, intent(in) :: eqmatrix(:,:)
        integer :: idx, op, i, j
        integer :: mapped_subset(level), sorted_mapped(level), sorted_stored(level)
        logical :: match

        is_symmetry_equivalent = .false.

        ! Check against each stored configuration
! Comprobar frente a cada configuración almacenada
        do idx = 1, best_count
            ! Try all symmetry operations
! Probar todas las operaciones de simetría
            do op = 1, nop
                ! Apply symmetry operation to current subset
! Aplicar la operación de simetría al subconjunto actual
                do i = 1, level
                    mapped_subset(i) = eqmatrix(op, subset(i))
                end do

                ! Sort both arrays for comparison
! Ordenar ambos vectores para compararlos
                sorted_mapped = mapped_subset
                call sort_int_ascending(sorted_mapped, level)
                sorted_stored(1:level) = best_subsets(1:level, idx)
                call sort_int_ascending(sorted_stored, level)

                ! Compare
! Comparar
                match = .true.
                do j = 1, level
                    if (sorted_mapped(j) /= sorted_stored(j)) then
                        match = .false.
                        exit
                    end if
                end do

                if (match) then
                    is_symmetry_equivalent = .true.
                    return
                end if
            end do
        end do
    end function is_symmetry_equivalent

end module exact
