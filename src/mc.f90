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

! Module `mc` implements the Monte Carlo sampling workflows.
! Monte Carlo driver that samples substitution levels, calls the SOD energy
! calculator, and aggregates Boltzmann-weighted statistics per level.
! Controlador Monte Carlo que muestrea niveles de sustitución y llama a la energía de SOD
! y agrega estadísticas ponderadas de Boltzmann por nivel.
! El módulo `mc` implementa los flujos de muestreo Monte Carlo.
module mc
    use consts
    use cli, only: classify_template_file, compose_mode_command, engine_name_to_filer, &
        is_help_token, parse_level_spec, to_lower_inplace
    use inputs, only: insod_file_data, read_insod_file
    use utils
    use settings, only: default_max_exact_combos, uniform_unique_cap, uniform_unique_min_cap, &
         uniform_cap_shrink, default_summary_filename, default_summary_txt_filename, blend_low_high_energy_level, &
         reference_relative, reset_blend_overrides
    use symmetry
    use workspace
    use calibration
    use energy_calculations
    use structure_io, only: motif_data_type, read_motif_file
    use omp_lib, only: omp_in_parallel
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none

    character(len=512), save :: summary_filename = ''
    character(len=512), save :: summary_txt_filename = ''
    character(len=512), save :: mc_motif_file = ''
    type(motif_data_type), save :: mc_motif
    logical, save :: mc_have_motif = .false.
    integer, save :: mc_filer = 0

    contains
    ! Drives Monte Carlo mode: parses options, initializes state, and iterates levels.
! Ejecuta el modo Monte Carlo: analiza opciones, inicializa el estado y recorre los niveles.
    subroutine run_sod_ensemble_mc(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: max_exact_combos
        real(dp) :: temperature
        integer :: max_substitutions, samples_per_level
        integer :: seed_value
        character(len=16) :: sampling_mode
        integer :: summary_unit, summary_txt_unit
        integer :: tempering_replicas, tempering_swap_every
        logical :: use_parallel
        logical :: omp_available
        logical :: force_restart_accept
        logical :: force_mc_sampling
        integer :: level_min, level_max
        integer :: level_start, level_end
        integer, pointer :: eqmatrix(:,:)
        integer, allocatable :: config(:)
        integer, allocatable :: local_config(:)
        integer, allocatable :: level_overrides(:)
        integer, allocatable :: level_targets(:)
        integer :: nop, total_sites
        integer :: level, effective_max
        integer :: level_idx
        logical :: has_level_overrides
        logical :: allow_parallel_levels
        logical :: effective_use_parallel
        logical :: force_parallel_lists
        character(len=512) :: template_gin_option
        character(len=64) :: engine_option
        character(len=16) :: protocol_option
        real(dp) :: tempering_max_temp_factor
        logical :: external_calibration_enabled
        integer :: offset
        type(insod_file_data) :: insod

        use_parallel = .false.
        omp_available = .false.
        ! $  use_parallel = .true.
        ! $  omp_available = .true.

        max_exact_combos = default_max_exact_combos
        summary_filename = default_summary_filename
        summary_txt_filename = default_summary_txt_filename

        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)

        template_gin_option = 'none'
        protocol_option = '2.0'
        mc_motif_file = ''
        engine_option = ''
        call parse_arguments(temperature, level_min, level_max, max_substitutions, samples_per_level, seed_value, sampling_mode, &
            use_parallel, omp_available, force_mc_sampling, level_overrides, has_level_overrides, force_parallel_lists, &
            template_gin_option, protocol_option, external_calibration_enabled, tempering_replicas, tempering_swap_every, &
            tempering_max_temp_factor, mc_motif_file, engine_option, offset)
        call set_calibration_template_gin(trim(template_gin_option))
        if (len_trim(mc_motif_file) > 0) then
            call read_motif_file(trim(mc_motif_file), mc_motif)
            mc_have_motif = (mc_motif%natoms > 0)
            if (mc_have_motif) then
                write(*,'(A,I0,A,1X,A)') 'Loaded ', mc_motif%natoms, ' motif atoms from', trim(mc_motif_file)
            end if
        end if
        call set_calibration_protocol_version(trim(protocol_option))
        call set_external_calibration_enabled(external_calibration_enabled)
        call read_insod_file(insod)

        ! Apply engine override
        if (len_trim(engine_option) > 0) then
            mc_filer = engine_name_to_filer(engine_option)
            if (mc_filer < 0) then
                write(error_unit,'(A)') 'Error: unrecognized engine: '//trim(engine_option)
                write(error_unit,'(A)') 'Accepted engines: gulp, lammps, vasp, castep, qe, ase'
                stop 1
            end if
            write(*,'(A,A,A,I0,A)') 'Engine override: --engine ', trim(engine_option), ' (FILER=', mc_filer, ')'
        else
            mc_filer = insod%filer
        end if
        call set_calibration_engine_from_filer(mc_filer)
        call configure_random_seed(seed_value)
        call configure_restart_mode(force_restart_accept)

        if (trim(sampling_mode) == 'tempering') then
            if (tempering_replicas < 2) then
                write(*,'(A)') 'Warning: parallel tempering requires at least 2 replicas; using 2.'
                tempering_replicas = 2
            end if
            if (tempering_swap_every < 1) then
                write(*,'(A)') 'Warning: invalid swap interval for tempering; using 25.'
                tempering_swap_every = 25
            end if
            if (tempering_max_temp_factor <= 1.0_dp) then
                write(*,'(A)') 'Warning: invalid maximum temperature factor for tempering; using 2.0.'
                tempering_max_temp_factor = 2.0_dp
            end if
        end if

        call init_energy_calc()
        call symmetry_initialize()
        call symmetry_get_matrix(eqmatrix, nop, total_sites)
        if (.not. associated(eqmatrix) .or. nop <= 0 .or. total_sites <= 0) then
            write(*,'(A)') 'Error: unable to obtain EQMATRIX or the number of positions/operators is invalid.'
            stop 1
        end if
        nullify(eqmatrix)
        call reset_blend_overrides(total_sites)

        call init_summary_files(summary_unit, summary_txt_unit)
        allow_parallel_levels = use_parallel
        effective_use_parallel = use_parallel

        if (allow_parallel_levels .and. .not. force_parallel_lists) then
            allow_parallel_levels = .false.
            if (use_parallel) then
                write(*,'(A)') 'Warning: disabling outer per-level parallelism; inner parallelism per level remains active.'
            end if
        end if

        if (max_substitutions < 0) then
            effective_max = total_sites
        else
            effective_max = min(max_substitutions, total_sites)
        end if

        if (has_level_overrides) then
            call finalize_level_overrides(level_overrides, level_targets, total_sites)
            if (.not. allocated(level_targets)) then
                write(*,'(A)') 'Error: the level list provided via -N contains no valid values.'
                stop 1
            end if
            if (allow_parallel_levels .and. .not. force_parallel_lists) then
                write(*,'(A)') 'Warning: disabling outer parallelism to preserve the order specified by -N.'
                allow_parallel_levels = .false.
                if (effective_use_parallel) then
                    write(*,'(A)') 'Warning: per-level sampling will execute sequentially for stability.'
                    effective_use_parallel = .false.
                end if
            else if (force_parallel_lists .and. allow_parallel_levels) then
                write(*,'(A)') 'Warning: retaining parallel execution with explicit list (--parallel-lists); results may appear out of order.'
            end if
        else
            if (level_max < 0) then
                level_start = max(0, level_min)
                level_end = effective_max
            else
                level_start = max(0, min(level_min, total_sites))
                level_end = min(level_max, total_sites)
            end if
            if (level_end > effective_max) level_end = effective_max
            if (level_start > level_end) level_start = level_end
        end if

        write(*,'(A)') '--- Calculation parameters ---'
        write(*,'(A,F10.2)') 'Temperature (K): ', temperature
        write(*,'(A,I6)') 'Substitutable sites (npos): ', total_sites
        write(*,'(A,I6)') 'Max evaluated substitutions: ', effective_max
        if (has_level_overrides) then
            call print_level_overrides(level_targets)
        else
            write(*,'(A,I0,A,I0)') 'Levels evaluated: ', level_start, ' .. ', level_end
        end if
        write(*,'(A,I8)') 'Exact enumeration threshold: ', max_exact_combos
        write(*,'(A,I8)') 'Random samples (if threshold exceeded): ', samples_per_level
        write(*,'(A)') 'Per-level results stored in: '//trim(summary_filename)
        write(*,'(A)') '                          and: '//trim(summary_txt_filename)
        write(*,'(A)') 'MC sampling method (if applicable): '//trim(sampling_mode)
        if (trim(sampling_mode) == 'tempering') then
            write(*,'(A,I6)') 'Tempering replicas: ', tempering_replicas
            write(*,'(A,I6)') 'Tempering swap interval: ', tempering_swap_every
            write(*,'(A,F8.3)') 'Tempering Tmax/T target: ', tempering_max_temp_factor
        end if
        if (external_calibration_enabled) then
            write(*,'(A,A)') 'External calibration: enabled with backend ', trim(calibration_engine_name())
        else
            write(*,'(A)') 'External calibration: disabled (default)'
        end if
        if (effective_use_parallel) then
            write(*,'(A)') 'OpenMP parallel: Yes'
        else
            write(*,'(A)') 'OpenMP parallel: No'
        end if
        if (force_mc_sampling) then
            write(*,'(A)') 'Force MC sampling: Yes'
        else
            write(*,'(A)') 'Force MC sampling: No'
        end if
        if (use_parallel) then
            write(*,'(A)') 'Note: per-level outputs may print non-sequentially during parallel execution.'
            write(*,'(A)') '      Summary files are reordered at the end to restore ascending levels.'
        end if
        write(*,*)
        flush(output_unit)

        call mc_workspace_prepare(use_parallel)

        if (allow_parallel_levels) then
            ! $omp parallel default(shared) private(local_config, level, level_idx)
            allocate(local_config(total_sites))
            if (has_level_overrides) then
                ! $omp do schedule(dynamic)
                do level_idx = 1, size(level_targets)
                    level = level_targets(level_idx)
                    call process_level(level, total_sites, local_config, temperature, samples_per_level, &
                        max_exact_combos, sampling_mode, force_mc_sampling, force_restart_accept, effective_use_parallel, summary_unit, summary_txt_unit, &
                        tempering_replicas, tempering_swap_every, tempering_max_temp_factor)
                end do
                ! $omp end do
            else
                ! $omp do schedule(dynamic)
                do level = level_start, level_end
                    call process_level(level, total_sites, local_config, temperature, samples_per_level, &
                        max_exact_combos, sampling_mode, force_mc_sampling, force_restart_accept, effective_use_parallel, summary_unit, summary_txt_unit, &
                        tempering_replicas, tempering_swap_every, tempering_max_temp_factor)
                end do
                ! $omp end do
            end if
            deallocate(local_config)
            ! $omp end parallel
        else
            allocate(config(total_sites))
            if (has_level_overrides) then
                do level_idx = 1, size(level_targets)
                    level = level_targets(level_idx)
                    call process_level(level, total_sites, config, temperature, samples_per_level, &
                        max_exact_combos, sampling_mode, force_mc_sampling, force_restart_accept, effective_use_parallel, summary_unit, summary_txt_unit, &
                        tempering_replicas, tempering_swap_every, tempering_max_temp_factor)
                end do
            else
                do level = level_start, level_end
                    call process_level(level, total_sites, config, temperature, samples_per_level, &
                        max_exact_combos, sampling_mode, force_mc_sampling, force_restart_accept, effective_use_parallel, summary_unit, summary_txt_unit, &
                        tempering_replicas, tempering_swap_every, tempering_max_temp_factor)
                end do
            end if
            deallocate(config)
        end if

        call close_summary_files(summary_unit, summary_txt_unit)
        call reorder_summary_outputs()
        call mc_workspace_finalize()
        call symmetry_finalize()
        call cleanup_energy_calc()
    end subroutine run_sod_ensemble_mc

    ! Parses command-line arguments for Monte Carlo mode and populates runtime parameters.
! Analiza los argumentos de línea de comandos del modo Monte Carlo y rellena los parámetros de ejecución.
    subroutine parse_arguments(temp, level_min, level_max, max_subs, samples_level, seed, sampler, use_parallel, omp_available, force_mc, level_list, has_level_list, force_parallel_lists, template_gin_option, protocol_option, external_calibration_enabled, tempering_replicas, tempering_swap_every, tempering_max_temp_factor, motif_file, engine_option, arg_offset)
        real(dp), intent(out) :: temp
        integer, intent(out) :: level_min, level_max, max_subs, samples_level, seed
        character(len=*), intent(out) :: sampler
        logical, intent(inout) :: use_parallel
        logical, intent(in) :: omp_available
        logical, intent(out) :: force_mc
        integer, allocatable, intent(out) :: level_list(:)
        logical, intent(out) :: has_level_list
        logical, intent(out) :: force_parallel_lists
        character(len=*), intent(out) :: template_gin_option
        character(len=*), intent(out) :: protocol_option
        logical, intent(out) :: external_calibration_enabled
        integer, intent(out) :: tempering_replicas, tempering_swap_every
        real(dp), intent(out) :: tempering_max_temp_factor
        character(len=*), intent(out) :: motif_file
        character(len=*), intent(out) :: engine_option
        integer, intent(in), optional :: arg_offset
        integer :: argc, ios, i
        character(len=256) :: carg, lowered, spec
        character(len=256), allocatable :: args(:)
        logical, allocatable :: skip(:)
        logical :: found_range, handled
        logical :: temp_set, max_subs_set, samples_set, seed_set, sampler_set, omp_set, template_set, protocol_set
        character(len=512) :: spec_trim
        integer :: eq_pos
        integer :: offset
        character(len=512) :: tpl_path
        character(len=16)  :: tpl_cat
        
        temp = 1000.0_dp
        level_min = 0
        level_max = -1
        max_subs = -1
        samples_level = 5000
        seed = -1
        sampler = 'uniform'
        force_mc = .false.
        has_level_list = .false.
        force_parallel_lists = .false.
        template_gin_option = 'none'
        protocol_option = '2.0'
        external_calibration_enabled = .false.
        tempering_replicas = 8
        tempering_swap_every = 25
        tempering_max_temp_factor = 2.0_dp
        motif_file = ''
        engine_option = ''
        template_set = .false.
        protocol_set = .false.
        temp_set = .false.
        max_subs_set = .false.
        samples_set = .false.
        seed_set = .false.
        sampler_set = .false.
        omp_set = .false.
        
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)

        argc = command_argument_count()
        if (argc <= offset) return
        
        allocate(args(argc))
        allocate(skip(argc))
        skip = .false.
        if (offset > 0 .and. argc > 0) then
            skip(1:min(offset, argc)) = .true.
        end if
        
        do i = 1, argc
            call get_command_argument(i, args(i))
            if (is_help_token(args(i))) then
                call print_usage(omp_available)
                stop
            end if
        end do
        
        do i = 1, argc
            if (skip(i)) cycle
            lowered = adjustl(args(i))
            call to_lower_inplace(lowered)
            if (trim(lowered) == '--force-mc' .or. trim(lowered) == '--force_mc' .or. &
            trim(lowered) == 'force-mc' .or. trim(lowered) == 'forcemc') then
                force_mc = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--parallel-lists' .or. trim(lowered) == '--parallel_lists' .or. &
            trim(lowered) == 'parallel-lists' .or. trim(lowered) == 'parallellists') then
                force_parallel_lists = .true.
                skip(i) = .true.
            else if (index(trim(lowered), '--template-gin=') == 1 .or. index(trim(lowered), '--template_gin=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --template-gin.'
                    call print_usage(omp_available)
                    stop 1
                end if
                template_gin_option = adjustl(args(i)(eq_pos+1:))
                template_set = .true.
                skip(i) = .true.
            else if (index(trim(lowered), '--protocol=') == 1 .or. index(trim(lowered), '--protocole=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --protocol.'
                    call print_usage(omp_available)
                    stop 1
                end if
                protocol_option = adjustl(args(i)(eq_pos+1:))
                protocol_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--template-gin' .or. trim(lowered) == '--template_gin') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing path after --template-gin.'
                    call print_usage(omp_available)
                    stop 1
                end if
                template_gin_option = adjustl(args(i+1))
                template_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (trim(lowered) == '--protocol' .or. trim(lowered) == '--protocole') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --protocol.'
                    call print_usage(omp_available)
                    stop 1
                end if
                protocol_option = adjustl(args(i+1))
                protocol_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (trim(lowered) == '--no-template-gin' .or. trim(lowered) == '--skip-template' .or. trim(lowered) == '--skip_template') then
                template_gin_option = 'none'
                template_set = .true.
                skip(i) = .true.
            else if (index(trim(lowered), '--motif=') == 1) then
                eq_pos = index(args(i), '=')
                motif_file = adjustl(args(i)(eq_pos+1:))
                skip(i) = .true.
            else if (trim(lowered) == '--motif') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing path after --motif.'
                    call print_usage(omp_available)
                    stop 1
                end if
                motif_file = adjustl(args(i+1))
                skip(i) = .true.
                skip(i+1) = .true.
            else if (trim(lowered) == '--calibration' .or. trim(lowered) == '--external-calibration') then
                external_calibration_enabled = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--no-calibration' .or. trim(lowered) == '--no-external-calibration') then
                external_calibration_enabled = .false.
                skip(i) = .true.
            else if (index(trim(lowered), '--engine=') == 1) then
                eq_pos = index(args(i), '=')
                engine_option = adjustl(args(i)(eq_pos+1:))
                skip(i) = .true.
            else if (trim(lowered) == '--engine') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --engine.'
                    call print_usage(omp_available)
                    stop 1
                end if
                engine_option = adjustl(args(i+1))
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--template=') == 1) then
                eq_pos = index(args(i), '=')
                tpl_path = adjustl(args(i)(eq_pos+1:))
                tpl_cat = classify_template_file(tpl_path)
                select case (trim(tpl_cat))
                case ('gin')
                    template_gin_option = trim(tpl_path)
                    template_set = .true.
                case default
                    write(error_unit,'(A)') 'Warning: --template file "'//trim(tpl_path)// &
                        '" (type: '//trim(tpl_cat)//') is not used in mc mode; only .gin/.include are supported.'
                end select
                skip(i) = .true.
            else if (trim(lowered) == '--template') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing path after --template.'
                    call print_usage(omp_available)
                    stop 1
                end if
                tpl_path = adjustl(args(i+1))
                tpl_cat = classify_template_file(tpl_path)
                select case (trim(tpl_cat))
                case ('gin')
                    template_gin_option = trim(tpl_path)
                    template_set = .true.
                case default
                    write(error_unit,'(A)') 'Warning: --template file "'//trim(tpl_path)// &
                        '" (type: '//trim(tpl_cat)//') is not used in mc mode; only .gin/.include are supported.'
                end select
                skip(i) = .true.
                skip(i+1) = .true.
            else if (trim(lowered) == '-t' .or. trim(lowered) == '--temperature') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --temperature.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_real_option(args(i+1), temp, 'temperature (--temperature/-T)')
                temp_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--temperature=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --temperature.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_real_option(adjustl(args(i)(eq_pos+1:)), temp, 'temperature (--temperature/-T)')
                temp_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '-m' .or. trim(lowered) == '--max-substitutions') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --max-substitutions.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_int_option(args(i+1), max_subs, 'maximum substitutions (--max-substitutions/-M)')
                max_subs_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--max-substitutions=') == 1 .or. index(trim(lowered), '--max_substitutions=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --max-substitutions.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_int_option(adjustl(args(i)(eq_pos+1:)), max_subs, 'maximum substitutions (--max-substitutions/-M)')
                max_subs_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '-s' .or. trim(lowered) == '--seed') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --seed.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_int_option(args(i+1), seed, 'seed (--seed/-s)')
                seed_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--seed=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --seed.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_int_option(adjustl(args(i)(eq_pos+1:)), seed, 'seed (--seed/-s)')
                seed_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '-c' .or. trim(lowered) == '--samples') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --samples.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_int_option(args(i+1), samples_level, 'samples per level (--samples/-S)')
                samples_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--samples=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --samples.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_int_option(adjustl(args(i)(eq_pos+1:)), samples_level, 'samples per level (--samples/-S)')
                samples_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--replicas' .or. trim(lowered) == '--tempering-replicas' .or. trim(lowered) == '--tempering_replicas') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --replicas.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_int_option(args(i+1), tempering_replicas, 'tempering replicas (--replicas)')
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--replicas=') == 1 .or. index(trim(lowered), '--tempering-replicas=') == 1 .or. &
                index(trim(lowered), '--tempering_replicas=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --replicas.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_int_option(adjustl(args(i)(eq_pos+1:)), tempering_replicas, 'tempering replicas (--replicas)')
                skip(i) = .true.
            else if (trim(lowered) == '--swap-every' .or. trim(lowered) == '--swap_every' .or. &
                trim(lowered) == '--tempering-swap-every' .or. trim(lowered) == '--tempering_swap_every') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --swap-every.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_int_option(args(i+1), tempering_swap_every, 'tempering swap interval (--swap-every)')
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--swap-every=') == 1 .or. index(trim(lowered), '--swap_every=') == 1 .or. &
                index(trim(lowered), '--tempering-swap-every=') == 1 .or. index(trim(lowered), '--tempering_swap_every=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --swap-every.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_int_option(adjustl(args(i)(eq_pos+1:)), tempering_swap_every, 'tempering swap interval (--swap-every)')
                skip(i) = .true.
            else if (trim(lowered) == '--max-temperature-factor' .or. trim(lowered) == '--max_temperature_factor') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --max-temperature-factor.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_real_option(args(i+1), tempering_max_temp_factor, 'tempering maximum temperature factor (--max-temperature-factor)')
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--max-temperature-factor=') == 1 .or. index(trim(lowered), '--max_temperature_factor=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --max-temperature-factor.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_positive_real_option(adjustl(args(i)(eq_pos+1:)), tempering_max_temp_factor, 'tempering maximum temperature factor (--max-temperature-factor)')
                skip(i) = .true.
            else if (trim(lowered) == '-a' .or. trim(lowered) == '--sampler') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing value after --sampler.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_sampler_option(args(i+1), sampler)
                sampler_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (index(trim(lowered), '--sampler=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: invalid argument for --sampler.'
                    call print_usage(omp_available)
                    stop 1
                end if
                call parse_sampler_option(adjustl(args(i)(eq_pos+1:)), sampler)
                sampler_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--omp' .or. trim(lowered) == '--enable-omp') then
                use_parallel = .true.
                omp_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--no-omp' .or. trim(lowered) == '--disable-omp') then
                use_parallel = .false.
                omp_set = .true.
                skip(i) = .true.
            end if
        end do
        
        found_range = .false.
        do i = 1, argc
            if (skip(i)) cycle
            lowered = adjustl(args(i))
            call to_lower_inplace(lowered)
            if (trim(lowered) == '-n') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: missing specification after -N.'
                    call print_usage(omp_available)
                    stop 1
                end if
                spec = adjustl(args(i+1))
                spec_trim = adjustl(spec)
                call parse_level_spec(spec_trim, level_min, level_max, level_list, has_level_list)
                skip(i) = .true.
                if (i < argc) skip(i+1) = .true.
                found_range = .true.
            end if
        end do
        if (.not. template_set) template_gin_option = 'none'
        if (.not. protocol_set) protocol_option = '2.0'
        do i = 1, argc
            if (skip(i)) cycle
            carg = args(i)
            handled = .false.
            
            lowered = adjustl(carg)
            call to_lower_inplace(lowered)
            
            if (.not. sampler_set) then
                if (trim(lowered) == 'uniform' .or. trim(lowered) == 'metropolis' .or. trim(lowered) == 'tempering') then
                    sampler = trim(lowered)
                    sampler_set = .true.
                    handled = .true.
                end if
            end if
            
            if (.not. handled .and. .not. omp_set) then
                if (trim(lowered) == 'omp' .or. trim(lowered) == 'noomp' .or. trim(lowered) == 'o' .or. trim(lowered) == 'no' .or. &
                trim(lowered) == 'true' .or. trim(lowered) == 'false' .or. trim(lowered) == '1' .or. trim(lowered) == '0') then
                    call parse_omp_flag(carg, use_parallel, omp_available)
                    omp_set = .true.
                    handled = .true.
                end if
            end if
            
            if (.not. handled .and. .not. temp_set) then
                read(carg,*,iostat=ios) temp
                if (ios == 0) then
                    temp_set = .true.
                    handled = .true.
                end if
            end if
            
            if (.not. handled .and. .not. max_subs_set) then
                read(carg,*,iostat=ios) max_subs
                if (ios == 0) then
                    max_subs_set = .true.
                    handled = .true.
                end if
            end if
            
            if (.not. handled .and. .not. samples_set) then
                read(carg,*,iostat=ios) samples_level
                if (ios == 0) then
                    if (samples_level <= 0) then
                        write(error_unit,'(A)') 'Warning: Nsamples must be positive; using 5000 samples instead.'
                        samples_level = 5000
                    end if
                    samples_set = .true.
                    handled = .true.
                end if
            end if
            
            if (.not. handled .and. .not. seed_set) then
                read(carg,*,iostat=ios) seed
                if (ios == 0) then
                    seed_set = .true.
                    handled = .true.
                end if
            end if
            
            if (.not. handled) then
                write(error_unit,'(A)') 'Warning: ignored argument -> '//trim(carg)
            end if
        end do
        
        if (.not. found_range .and. level_max >= 0) then
            if (level_min < 0) level_min = 0
        end if
        
        if (allocated(args)) deallocate(args)
        if (allocated(skip)) deallocate(skip)
        return
    end subroutine parse_arguments

    ! Parses a floating-point option value and raises an error when invalid.
! Analiza un valor de opción de punto flotante y lanza un error si es inválido.
    subroutine parse_real_option(raw_value, target, label)
        character(len=*), intent(in) :: raw_value
        real(dp), intent(out) :: target
        character(len=*), intent(in) :: label
        integer :: ios

        character(len=256) :: text
        text = adjustl(raw_value)
        read(text, *, iostat=ios) target
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: invalid value for '//trim(label)//'.'
            stop 1
        end if
    end subroutine parse_real_option

    ! Parses an integer option value.
! Analiza un valor de opción entera.
    subroutine parse_int_option(raw_value, target, label)
        character(len=*), intent(in) :: raw_value
        integer, intent(out) :: target
        character(len=*), intent(in) :: label
        integer :: ios

        character(len=256) :: text
        text = adjustl(raw_value)
        read(text, *, iostat=ios) target
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: invalid value for '//trim(label)//'.'
            stop 1
        end if
    end subroutine parse_int_option

    ! Parses a strictly positive integer option value.
! Analiza un valor de opción entera estrictamente positivo.
    subroutine parse_positive_int_option(raw_value, target, label)
        character(len=*), intent(in) :: raw_value
        integer, intent(out) :: target
        character(len=*), intent(in) :: label
        integer :: ios

        character(len=256) :: text
        text = adjustl(raw_value)
        read(text, *, iostat=ios) target
        if (ios /= 0 .or. target <= 0) then
            write(error_unit,'(A)') 'Error: invalid value for '//trim(label)//'.'
            stop 1
        end if
    end subroutine parse_positive_int_option

    ! Parses a strictly positive floating-point option value.
! Analiza un valor de opción real estrictamente positivo.
    subroutine parse_positive_real_option(raw_value, target, label)
        character(len=*), intent(in) :: raw_value
        real(dp), intent(out) :: target
        character(len=*), intent(in) :: label
        integer :: ios

        character(len=256) :: text
        text = adjustl(raw_value)
        read(text, *, iostat=ios) target
        if (ios /= 0 .or. target <= 0.0_dp) then
            write(error_unit,'(A)') 'Error: invalid value for '//trim(label)//'.'
            stop 1
        end if
    end subroutine parse_positive_real_option

    ! Normalizes the sampler argument to a supported keyword.
! Normaliza el argumento del muestreador a una palabra clave admitida.
    subroutine parse_sampler_option(raw_value, sampler)
        character(len=*), intent(in) :: raw_value
        character(len=*), intent(inout) :: sampler
        character(len=256) :: lowered
        integer :: idx

        lowered = adjustl(raw_value)
        do idx = 1, len_trim(lowered)
            if (lowered(idx:idx) >= 'A' .and. lowered(idx:idx) <= 'Z') lowered(idx:idx) = achar(iachar(lowered(idx:idx)) + 32)
        end do
        lowered = trim(lowered)
        select case (lowered)
        case ('uniform', 'metropolis', 'tempering')
            sampler = lowered
        case default
            write(error_unit,'(A)') 'Error: invalid sampler. Use "uniform", "metropolis", or "tempering".'
            stop 1
        end select
    end subroutine parse_sampler_option

    ! Interprets an OpenMP flag token and updates the parallel execution mode.
! Interpreta un token de bandera OpenMP y actualiza el modo de ejecución en paralelo.
subroutine parse_omp_flag(raw, use_parallel, omp_available)
    character(len=*), intent(in) :: raw
    logical, intent(inout) :: use_parallel
    logical, intent(in) :: omp_available
    character(len=256) :: token
    integer :: i, lt
    
    token = adjustl(raw)
    lt = len_trim(token)
    do i = 1, lt
        if (token(i:i) >= 'A' .and. token(i:i) <= 'Z') token(i:i) = achar(iachar(token(i:i)) + 32)
    end do
    
    select case (trim(token))
    case ('omp', 'o', '1', 'true')
        if (omp_available) then
            use_parallel = .true.
        else
            write(error_unit,'(A)') 'Warning: OpenMP not available in this build; ignoring argument.'
            use_parallel = .false.
        end if
    case ('noomp', 'no', '0', 'false')
        use_parallel = .false.
    case default
        write(error_unit,'(A)') 'Warning: unknown argument for parallel mode; keeping previous configuration.'
    end select
end subroutine parse_omp_flag

    ! Filters and deduplicates explicit level overrides provided via -N lists.
! Filtra y elimina duplicados de los niveles explícitos proporcionados mediante listas con -N.
    subroutine finalize_level_overrides(raw_levels, targets, total_sites)
        integer, allocatable, intent(inout) :: raw_levels(:)
        integer, allocatable, intent(out) :: targets(:)
        integer, intent(in) :: total_sites
        integer, allocatable :: temp(:)
        integer :: idx, valid_count, val
        logical :: duplicate
        
        if (allocated(targets)) deallocate(targets)
        if (.not. allocated(raw_levels)) return
        allocate(temp(size(raw_levels)))
        temp = 0
        valid_count = 0
        do idx = 1, size(raw_levels)
            val = raw_levels(idx)
            if (val < 0 .or. val > total_sites) then
                write(*,'(A,I0,A)') 'Warning: level ', val, ' is out of range; skipping.'
                flush(output_unit)
                cycle
            end if
            duplicate = .false.
            if (valid_count > 0) then
                if (any(temp(1:valid_count) == val)) duplicate = .true.
            end if
            if (duplicate) cycle
            valid_count = valid_count + 1
            temp(valid_count) = val
        end do
        if (valid_count > 0) then
            allocate(targets(valid_count))
            targets(1:valid_count) = temp(1:valid_count)
        end if
        deallocate(temp)
        deallocate(raw_levels)
    end subroutine finalize_level_overrides

    ! Prints the explicit level list selected via -N.
! Imprime la lista explícita de niveles seleccionada mediante -N.
    subroutine print_level_overrides(levels)
        integer, intent(in) :: levels(:)
        integer :: idx
        write(*,'(A)') 'Levels evaluated: explicit list'
        write(*,'(A)', advance='no') 'Values: '
        do idx = 1, size(levels)
            if (idx > 1) write(*,'(A)', advance='no') ', '
            write(*,'(I0)', advance='no') levels(idx)
        end do
        write(*,*)
    end subroutine print_level_overrides

! Emits program usage information including optional OpenMP note.
! Emite la información de uso del programa, incluida una nota opcional sobre OpenMP.
subroutine print_usage(omp_available)
    logical, intent(in) :: omp_available
    
    character(len=32) :: cap_str
    character(len=256) :: command_name
    character(len=32) :: shrink_str
    real(dp) :: shrink_percent
    
    command_name = compose_mode_command('mc')
    write(cap_str,'(I0)') uniform_unique_min_cap
    shrink_percent = uniform_cap_shrink * 100.0_dp
    write(shrink_str,'(F5.1)') shrink_percent
    
    write(*,'(A)') 'Usage: '//trim(command_name)//' [-T <K>] [-M <Nmax>] [-C <Nsamples>] [-s <seed>] [-a <sampler>] [--omp|--no-omp] [--force-mc] [-N range] [--calibration|--external-calibration]'
    write(*,'(A)') '       '//trim(command_name)//' --help'
    write(*,'(A)') ''
    write(*,'(A)') 'Aliases: mc, sampling, sample, spbemc'
    write(*,'(A)') ''
    write(*,'(A)') 'Optional arguments (defaults in brackets; legacy positional syntax remains available):'
    write(*,'(A)') '  -T, --temperature <K>     Temperature in Kelvin for Boltzmann weights [1000].'
    write(*,'(A)') '  -M, --max-substitutions N Maximum substitutions evaluated when -N is absent [-1 -> all].'
    write(*,'(A)') '  -C, --samples N           MC samples per level when C(N,npos) exceeds the threshold [5000].'
    write(*,'(A)') '  -s, --seed value          Random-number seed [-1 -> derived from system_clock].'
    write(*,'(A)') '  -a, --sampler mode        "uniform", "metropolis", or "tempering" [uniform].'
    write(*,'(A)') '      --replicas N          Replica count for parallel tempering [8].'
    write(*,'(A)') '      --swap-every N        Local steps per replica between swap attempts [25].'
    write(*,'(A)') '      --max-temperature-factor X  Maximum replica temperature as X*T [2.0].'
    write(*,'(A)') '      --omp / --no-omp      Explicitly enable or disable OpenMP.'
    write(*,'(A)') '  -N spec                   Range or list: -N 5 (only 5), -N 3:8 (3 to 8), -N 12,30,45 (explicit list).'
    write(*,'(A)') '  --parallel-lists          Keep OpenMP enabled even with -N lists (results may appear out of order).'
    write(*,'(A)') '  --force-mc                Force Monte Carlo sampling even if C(N,npos) <= exact threshold.'
    write(*,'(A)') '  --template-gin <file>     Append the specified template fragment when generating .gin files.'
    write(*,'(A)') '                             Use "default" to copy the bundled template include.'
    write(*,'(A)') '  --no-template-gin         Skip template fragments when creating .gin files (default).'
    write(*,'(A)') '  --engine <name>           Override the INSOD FILER engine: gulp, lammps, vasp, castep, qe, ase.'
    write(*,'(A)') '  --template <file>         Provide a template file (auto-classified by extension).'
    write(*,'(A)') '                             .gin / .include -> GULP optimization template.'
    write(*,'(A)') '  --calibration             Enable external calibration/evaluation with the backend selected by FILER.'
    write(*,'(A)') '  --external-calibration    Alias for --calibration. Disabled by default.'
    write(*,'(A)') '  --no-calibration          Explicitly disable external calibration/evaluation.'
    write(*,'(A)') '  --protocol <ver>          Select the GULP workflow protocol: 1.0 or 2.0 [2.0].'
    write(*,'(A)') ''
    write(*,'(A)') 'Additional details:'
    write(*,'(A)') '  - By default every level from N=0 up to Nmax is evaluated unless -N restricts the set.'
    write(*,'(A)') '  - If C(N,npos) <= 200000 all configurations are enumerated; otherwise Monte Carlo sampling is used.'
    write(*,'(A)') '  - Uniform sampling applies adaptive control to the unique configuration target ' // &
    '(minimum '//trim(cap_str)//', factor '//trim(adjustl(shrink_str))//'%).'
    write(*,'(A)') '  - Aggregated results for each level are written to '//trim(summary_filename)//'.'
    write(*,'(A)') '  - The plain-text summary is written to '//trim(summary_txt_filename)//'.'
    if (omp_available) then
        write(*,'(A)') '  - This binary was compiled with OpenMP support.'
    else
        write(*,'(A)') '  - This binary was compiled without OpenMP support.'
    end if
    write(*,'(A)') ''
    write(*,'(A)') 'Examples:'
    write(*,'(A)') '  '//trim(command_name)//'                     # Run with defaults.'
    write(*,'(A)') '  '//trim(command_name)//' 800 6 2000 1234 metropolis omp'
    write(*,'(A)') '                                        # 800 K, up to 6 substitutions, Metropolis and OpenMP.'
    write(*,'(A)') '  '//trim(command_name)//' -N 12,30,45         # Evaluate only levels 12, 30, and 45.'
    write(*,'(A)') '  '//trim(command_name)//' -a tempering --replicas 8 --swap-every 25 -T 1500'
    write(*,'(A)') '                                        # Parallel tempering at 1500 K.'
end subroutine print_usage

! Initializes the intrinsic random generator with a deterministic or clock-based seed.
! Inicializa el generador aleatorio intrínseco con una semilla determinista o basada en el reloj.
subroutine configure_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: seed_array(:)
    integer :: count, count_rate
    
    call random_seed(size=n)
    allocate(seed_array(n))
    
    if (seed >= 0) then
        seed_array = seed + 37 * [(i-1, i=1,n)]
    else
        call system_clock(count, count_rate)
        if (count == 0 .and. count_rate == 0) then
            seed_array = 123456 + 37 * [(i-1, i=1,n)]
        else
            seed_array = count + 37 * [(i-1, i=1,n)]
        end if
    end if
    
    call random_seed(put=seed_array)
    deallocate(seed_array)
end subroutine configure_random_seed

! Reads the restart environment toggle and reports whether restart moves are forced.
! Lee el interruptor de entorno de reinicio e informa si los movimientos de reinicio están forzados.
subroutine configure_restart_mode(force_restart)
    logical, intent(out) :: force_restart
    character(len=32) :: env_value
    integer :: status, len_env, i
    character(len=32) :: token
    
    force_restart = .true.
    call get_environment_variable('SOD_FORCE_RESTART_ACCEPT', env_value, length=len_env, status=status)
    if (status /= 0 .or. len_env <= 0) then
        write(*,'(A)') 'Forced restart acceptance: enabled (default)'
        return
    end if
    
    token = adjustl(env_value(1:len_env))
    do i = 1, len_trim(token)
        if (token(i:i) >= 'A' .and. token(i:i) <= 'Z') then
            token(i:i) = achar(iachar(token(i:i)) + 32)
        end if
    end do
    token = token(1:len_trim(token))
    
    select case (trim(token))
    case ('0', 'no', 'false', 'off', 'disable', 'disabled')
        force_restart = .false.
        write(*,'(A)') 'Forced restart acceptance: disabled (SOD_FORCE_RESTART_ACCEPT='//trim(token)//')'
    case ('1', 'yes', 'true', 'on', 'enable', 'enabled')
        force_restart = .true.
        write(*,'(A)') 'Forced restart acceptance: enabled (SOD_FORCE_RESTART_ACCEPT='//trim(token)//')'
    case default
        force_restart = .true.
        write(*,'(A)') 'Forced restart acceptance: enabled (unrecognized value, using default)'
    end select
end subroutine configure_restart_mode

! Opens the summary CSV and text files and writes their headers.
! Abre los archivos de resumen CSV y de texto y escribe sus cabeceras.
subroutine init_summary_files(unit_csv, unit_txt)
    integer, intent(out) :: unit_csv, unit_txt
    integer :: ios
    
    open(newunit=unit_csv, file=summary_filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A)') 'Error: failed to create summary file '//trim(summary_filename)
        stop 1
    end if
    open(newunit=unit_txt, file=summary_txt_filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A)') 'Error: failed to create summary file '//trim(summary_txt_filename)
        stop 1
    end if
    write(unit_csv,'(A)') '#N;FracY;E_exp_total;E_min_total;Var_total;E_exp_X_side;E_min_X_side;'// &
    'E_exp_Y_side;E_min_Y_side;E_exp_mixed;Delta_exp_total;Delta_min_total;'// &
    'Delta_exp_X_side;Delta_min_X_side;Delta_exp_Y_side;Delta_min_Y_side;Delta_exp_mixed;'// &
    'Acceptance_ratio'
    write(unit_txt,'(A)') '#N FracY E_exp_total E_min_total Var_total E_exp_X_side E_min_X_side '// &
    'E_exp_Y_side E_min_Y_side E_exp_mixed Delta_exp_total Delta_min_total '// &
    'Delta_exp_X_side Delta_min_X_side Delta_exp_Y_side Delta_min_Y_side Delta_exp_mixed '// &
    'Acceptance_ratio'
end subroutine init_summary_files

! Closes the summary file units if they are valid.
! Cierra las unidades de los archivos de resumen si son válidas.
subroutine close_summary_files(unit_csv, unit_txt)
    integer, intent(in) :: unit_csv, unit_txt
    
    if (unit_csv /= 0) close(unit_csv)
    if (unit_txt /= 0) close(unit_txt)
end subroutine close_summary_files

! Chooses exhaustive or stochastic evaluation for a substitution level and dispatches it.
! Elige una evaluación exhaustiva o estocástica para un nivel de sustitución y la despacha.
subroutine process_level(level, total_sites, config, temperature, samples_level, max_exact, sampler, &
    force_sampling, force_restart_accept, use_parallel, summary_unit, summary_txt_unit, tempering_replicas, &
    tempering_swap_every, tempering_max_temp_factor)
    integer, intent(in) :: level, total_sites, samples_level, max_exact
    integer, intent(inout) :: config(:)
    real(dp), intent(in) :: temperature
    character(len=*), intent(in) :: sampler
    logical, intent(in) :: force_sampling
    logical, intent(in) :: force_restart_accept
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit
    integer, intent(in) :: tempering_replicas, tempering_swap_every
    real(dp), intent(in) :: tempering_max_temp_factor
    
    integer(ip) :: total_comb
    integer :: ncomb_int
    logical :: use_sampling
    logical :: force_this_level
    
    if (level > total_sites) return
    
    total_comb = binomial_int64(total_sites, level)
    if (total_comb == 0_ip) then
        write(*,'(A,I3)') 'Level ', level, ': no valid combinations.'
        return
    end if

    ! Reset the per-level partial GULP snapshot before collecting new MC points.
! Reinicia la instantánea parcial GULP por nivel antes de recoger nuevos puntos MC.
    call reset_mc_gulp_energy_snapshot(level)
    
    force_this_level = force_sampling .and. (level > 0 .and. level < total_sites)
    
    if (level == 0 .or. level == total_sites) then
        use_sampling = .false.
    else
        use_sampling = force_this_level .or. (total_comb > max_exact)
    end if
    
    if (use_sampling) then
        if (force_this_level .and. total_comb <= max_exact) then
            write(*,'(A,I0,A,I0,A,I0,A)') 'Level ', level, ': forced Monte Carlo sampling (combinations=', int(total_comb, kind=kind(max_exact)), ', exact threshold=', max_exact, ').'
        else if (sampler == 'metropolis') then
            write(*,'(A,I0,A)') 'Level ', level, ': using Metropolis-Hastings Monte Carlo sampling.'
        else if (sampler == 'tempering') then
            write(*,'(A,I0,A,I0,A,I0,A,F6.2,A)') 'Level ', level, ': using parallel tempering with ', &
                tempering_replicas, ' replicas, swap interval ', tempering_swap_every, ', Tmax/T=', &
                tempering_max_temp_factor, '.'
        else
            write(*,'(A,I0,A)') 'Level ', level, ': using uniform Monte Carlo sampling.'
        end if
        flush(output_unit)
        if (sampler == 'metropolis') then
            call metropolis_level(level, total_sites, config, temperature, samples_level, total_comb, &
            force_restart_accept, use_parallel, summary_unit, summary_txt_unit)
        else if (sampler == 'tempering') then
            call tempering_level(level, total_sites, config, temperature, samples_level, total_comb, &
            force_restart_accept, use_parallel, summary_unit, summary_txt_unit, tempering_replicas, &
            tempering_swap_every, tempering_max_temp_factor)
        else
            call monte_carlo_level(level, total_sites, config, temperature, samples_level, total_comb, &
            use_parallel, summary_unit, summary_txt_unit)
        end if
    else
        ncomb_int = int(total_comb, kind=kind(1))
        call exhaustive_level(level, total_sites, config, temperature, total_comb, ncomb_int, &
        use_parallel, summary_unit, summary_txt_unit)
    end if
end subroutine process_level

! Enumerates every configuration for a level and accumulates energy statistics.
! Enumera cada configuración de un nivel y acumula estadísticas de energía.
subroutine exhaustive_level(level, total_sites, config, temperature, total_comb, ncomb_int, use_parallel, &
    summary_unit, summary_txt_unit)
    integer, intent(in) :: level, total_sites, ncomb_int
    integer(ip), intent(in) :: total_comb
    integer, intent(inout) :: config(:)
    real(dp), intent(in) :: temperature
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit
    
    real(dp), pointer :: energies(:)
    real(dp), pointer :: energies_low(:)
    real(dp), pointer :: energies_high(:)
    integer :: idx
    integer :: best_positions(max(1, level))
    integer :: best_count
    real(dp) :: best_energy
    integer :: comb(max(1, level))
    integer :: j
    
    allocate(energies(ncomb_int))
    allocate(energies_low(ncomb_int))
    allocate(energies_high(ncomb_int))
    idx = 0
    best_energy = huge(1.0_dp)
    best_positions = 0
    best_count = 0

    if (level == 0) then
        config = 1
        idx = 1
        call calculate_structure_energy(config, total_sites, energies(idx), energies_low(idx), energies_high(idx))
        best_energy = energies(idx)
    else
        do j = 1, level
            comb(j) = j
        end do
        do
            config = 1
            config(comb(1:level)) = 2
            idx = idx + 1
            call calculate_structure_energy(config, total_sites, energies(idx), energies_low(idx), energies_high(idx))
            if (energies(idx) < best_energy) then
                best_energy = energies(idx)
                best_positions(1:level) = comb(1:level)
                best_count = level
            end if
            if (.not. next_combination(comb, total_sites)) exit
        end do
    end if
    
    if (idx == 0) then
        write(*,'(A,I3)') 'Level ', level, ': could not generate valid configurations.'
        flush(output_unit)
        deallocate(energies, energies_low, energies_high)
        return
    end if
    
    call summarize_level(level, total_sites, temperature, total_comb, real(idx, dp), energies(1:idx), best_energy, &
    best_positions, best_count, 0, &
    energies_low(1:idx), energies_high(1:idx), &
    use_parallel, summary_unit, summary_txt_unit)
    
    deallocate(energies, energies_low, energies_high)
end subroutine exhaustive_level

! Samples random configurations uniformly when exhaustive enumeration is infeasible.
! Muestrea configuraciones aleatorias de forma uniforme cuando la enumeración exhaustiva no es viable.
subroutine monte_carlo_level(level, total_sites, config, temperature, samples_level, total_comb, use_parallel, &
    summary_unit, summary_txt_unit)
    integer, intent(in) :: level, total_sites, samples_level
    integer(ip), intent(in) :: total_comb
    integer, intent(inout) :: config(:)
    real(dp), intent(in) :: temperature
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit
    
    integer :: subset(max(1, level))
    integer :: unique_count
    integer :: unique_cap, initial_unique_cap
    integer :: samples_target
    integer :: attempt_count, max_attempts
    integer :: adaptive_window, adaptive_threshold, stall_counter, prev_unique_snapshot
    integer :: new_cap
    real(dp), pointer :: energies(:) => null()
    real(dp), pointer :: energies_low(:) => null()
    real(dp), pointer :: energies_high(:) => null()
    real(dp) :: best_energy
    integer :: best_subset(max(1, level))
    integer :: best_count
    integer :: best_idx
    integer :: trace_unit
    integer :: calib_best_idx
    integer, pointer :: unique_subsets(:,:) => null()
    integer, pointer :: accept_attempt(:) => null()
    integer, allocatable :: config_local(:)
    real(dp), pointer :: low_contribs(:,:)
    real(dp), pointer :: high_contribs(:,:)
    real(dp) :: low_contrib_tmp(4), high_contrib_tmp(4)
    integer, pointer :: eqmatrix(:,:)
    integer, pointer :: canonical_subset(:)
    integer :: existing
    integer :: nop, npos
    integer(ip) :: cap_ip
    real(dp), pointer :: gulp_energies(:)
    logical :: gulp_success
    real(dp) :: sum_energy, sumsq_energy, mean_energy, variance_energy, std_energy
    logical :: cap_reduced, allow_inner_parallel
    integer :: idx
    
    if (level == 0) then
        call exhaustive_level(level, total_sites, config, temperature, total_comb, 1, use_parallel, &
        summary_unit, summary_txt_unit)
        return
    end if
    
    samples_target = max(1, max(samples_level, uniform_unique_cap))
    if (total_comb > 0_ip) then
        cap_ip = min(total_comb, int(uniform_unique_cap, kind=ip))
    else
        cap_ip = int(uniform_unique_cap, kind=ip)
    end if
    unique_cap = max(1, int(cap_ip))
    initial_unique_cap = unique_cap
    adaptive_window = max(100, samples_target / 4)
    adaptive_threshold = max(1, adaptive_window / 8)
    stall_counter = 0
    prev_unique_snapshot = 0
    cap_reduced = .false.
    call mc_workspace_checkout(level, unique_cap, total_sites, unique_subsets, energies, energies_low, energies_high, &
         low_contribs, high_contribs, accept_attempt, gulp_energies, canonical_subset)
    call symmetry_get_matrix(eqmatrix, nop, npos)
    if (.not. associated(eqmatrix)) then
        write(*,'(A)') 'Warning: could not obtain EQMATRIX for sample deduplication.'
        flush(output_unit)
        nullify(unique_subsets)
        nullify(energies)
        nullify(energies_low)
        nullify(energies_high)
        nullify(low_contribs)
        nullify(high_contribs)
        nullify(accept_attempt)
        nullify(gulp_energies)
        nullify(canonical_subset)
        nullify(eqmatrix)
        return
    end if
    calib_best_idx = 0
    unique_count = 0
    unique_subsets = 0
    energies = 0.0_dp
    energies_low = huge(1.0_dp)
    energies_high = huge(1.0_dp)
    low_contribs = 0.0_dp
    high_contribs = 0.0_dp
    accept_attempt = 0
    best_energy = huge(1.0_dp)
    best_subset = 0
    best_count = 0
    config = 1
    call open_mc_trace_file(level, trace_unit)
    attempt_count = 0
    max_attempts = max(50, samples_target * 50)
    
    ! Deduplicate uniform samples until enough symmetry-unique configurations are collected.
! Elimina duplicados de las muestras uniformes hasta reunir suficientes configuraciones únicas por simetría.
    do while (unique_count < unique_cap .and. attempt_count < max_attempts)
        call random_subset(total_sites, level, subset)
        call canonicalize_subset(subset, level, eqmatrix, nop, canonical_subset)
        existing = find_subset_index(canonical_subset, level, unique_subsets, unique_count)
        attempt_count = attempt_count + 1
        if (existing /= 0) cycle
        
        unique_count = unique_count + 1
        unique_subsets(:, unique_count) = canonical_subset(1:level)
        accept_attempt(unique_count) = attempt_count
        
        if (adaptive_window > 0) then
            if (mod(attempt_count, adaptive_window) == 0) then
                if (unique_count - prev_unique_snapshot < adaptive_threshold) then
                    stall_counter = stall_counter + 1
                    if (stall_counter >= 2) then
                        new_cap = max(unique_count, int(real(unique_cap, dp) * uniform_cap_shrink))
                        if (initial_unique_cap > uniform_unique_min_cap) then
                            new_cap = max(new_cap, uniform_unique_min_cap)
                        end if
                        if (new_cap < unique_cap) then
                            unique_cap = new_cap
                            cap_reduced = .true.
                            write(*,'(A,I0,A,I0,A)') 'Level ', level, ': adaptive control shrinks target to ', unique_cap, ' unique configurations.'
                            flush(output_unit)
                        end if
                    end if
                else
                    stall_counter = 0
                end if
                prev_unique_snapshot = unique_count
            end if
        end if
    end do
    
    if (attempt_count >= max_attempts .and. unique_count < unique_cap) then
        write(*,'(A,I0,A,I0,A)') 'Level ', level, ': random sampling reached ', attempt_count, ' attempts without meeting the unique quota.'
        flush(output_unit)
    end if
    
    if (unique_count == 0) then
        write(*,'(A,I0)') 'Level ', level, ': no symmetry-unique configuration could be found.'
        flush(output_unit)
        if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
        nullify(eqmatrix, unique_subsets, energies, energies_low, energies_high, low_contribs, high_contribs, &
               gulp_energies, canonical_subset, accept_attempt)
        return
    else
        write(*,'(A,I0,A,I0)') 'Level ', level, ': symmetry-unique configurations collected: ', unique_count
        flush(output_unit)
        if (cap_reduced) then
            write(*,'(A,I0,A,I0)') 'Level ', level, ': final adaptive target: ', unique_cap
            flush(output_unit)
        else if (unique_count < initial_unique_cap) then
            write(*,'(A,I0,A,I0,A)') 'Level ', level, ': target of ', initial_unique_cap, ' unique samples not reached.'
            flush(output_unit)
        end if
    end if
    
    if (unique_count > 0) then
        allow_inner_parallel = use_parallel
        if (allow_inner_parallel) then
            if (omp_in_parallel()) allow_inner_parallel = .false.
        end if
        
        if (allow_inner_parallel) then
            ! $omp parallel default(shared) private(idx, low_contrib_tmp, high_contrib_tmp, config_local)
            allocate(config_local(total_sites))
            config_local = 1
            ! $omp do schedule(dynamic)
            do idx = 1, unique_count
                config_local = 1
                config_local(unique_subsets(1:level, idx)) = 2
                call calculate_structure_energy(config_local, total_sites, energies(idx), &
                energy_low_side=energies_low(idx), energy_high_side=energies_high(idx), &
                low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)
                low_contribs(:, idx) = low_contrib_tmp
                high_contribs(:, idx) = high_contrib_tmp
            end do
            ! $omp end do
            deallocate(config_local)
            ! $omp end parallel
        else
            do idx = 1, unique_count
                config = 1
                config(unique_subsets(1:level, idx)) = 2
                call calculate_structure_energy(config, total_sites, energies(idx), &
                energy_low_side=energies_low(idx), energy_high_side=energies_high(idx), &
                low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)
                low_contribs(:, idx) = low_contrib_tmp
                high_contribs(:, idx) = high_contrib_tmp
            end do
        end if
        
        if (trace_unit /= 0) then
            do idx = 1, unique_count
                call write_mc_trace_step(trace_unit, idx, accept_attempt(idx), energies(idx), energies_low(idx), energies_high(idx))
            end do
            flush(trace_unit)
        end if
        
        best_idx = 1
        best_energy = energies(1)
        do idx = 2, unique_count
            if (energies(idx) < best_energy) then
                best_energy = energies(idx)
                best_idx = idx
            end if
        end do
        best_subset(1:level) = unique_subsets(:, best_idx)
        best_count = level
    end if
    config = 1
    
    if (external_calibration_is_enabled()) then
        call attempt_calibration_from_samples(level, total_sites, 1, unique_count, unique_subsets, &
        energies, energies_low, energies_high, low_contribs, high_contribs, calib_best_idx)
        if (calib_best_idx > 0) then
            best_energy = energies(calib_best_idx)
            best_subset(1:level) = unique_subsets(:, calib_best_idx)
            best_count = level
        end if
        
        gulp_success = .false.
        gulp_energies = 0.0_dp
        ! Evaluate precise energies with the selected external engine for the stored unique configurations.
! Evalúa energías precisas con el motor externo seleccionado para las configuraciones únicas almacenadas.
        call evaluate_subsets_with_engine(level, total_sites, unique_count, unique_subsets(:,1:unique_count), config, gulp_energies, gulp_success)
        if (gulp_success) then
            call update_level_blend_override_from_samples(level, total_sites, energies_low(1:unique_count), &
                energies_high(1:unique_count), gulp_energies(1:unique_count), unique_count, 'MC external sample', level_prefix='mc')
            energies(1:unique_count) = gulp_energies(1:unique_count)
            best_energy = minval(energies(1:unique_count))
            existing = 1
            do while (existing <= unique_count)
                if (energies(existing) == best_energy) then
                    best_subset(1:level) = unique_subsets(:, existing)
                    exit
                end if
                existing = existing + 1
            end do
            sum_energy = 0.0_dp
            sumsq_energy = 0.0_dp
            do existing = 1, unique_count
                sum_energy = sum_energy + energies(existing)
                sumsq_energy = sumsq_energy + energies(existing)**2
            end do
            mean_energy = sum_energy / real(unique_count, dp)
            variance_energy = max(0.0_dp, (sumsq_energy / real(unique_count, dp)) - mean_energy**2)
            std_energy = sqrt(variance_energy)
            write(*,'(A,I0,A,A,A,F16.6,A,F16.6)') 'Level ', level, ': ', trim(calibration_engine_name()), ' mean = ', mean_energy, ' eV, deviation = ', std_energy
            flush(output_unit)
        else
            write(*,'(A,I0,A,A,A)') 'Level ', level, ': some samples could not be evaluated with ', trim(calibration_engine_name()), '; keeping internal energies.'
            flush(output_unit)
        end if
    end if
    call write_monitored_outsod(level, total_sites, unique_count, unique_subsets(:,1:unique_count))
    config = 1
    
    call summarize_level(level, total_sites, temperature, total_comb, real(unique_count, dp), energies(1:unique_count), &
    best_energy, best_subset, best_count, 0, &
    energies_low(1:unique_count), energies_high(1:unique_count), &
    use_parallel, summary_unit, summary_txt_unit)
    
    if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
    nullify(eqmatrix)
    nullify(unique_subsets)
    nullify(energies)
    nullify(energies_low)
    nullify(energies_high)
    nullify(low_contribs)
    nullify(high_contribs)
    nullify(gulp_energies)
    nullify(accept_attempt)
    nullify(canonical_subset)
end subroutine monte_carlo_level

! Runs a Metropolis-Hastings walk with optional restart moves for a substitution level.
! Ejecuta una caminata de Metropolis-Hastings con movimientos de reinicio opcionales para un nivel de sustitución.
subroutine metropolis_level(level, total_sites, config, temperature, samples_level, total_comb, force_restart_accept, use_parallel, &
    summary_unit, summary_txt_unit)
    integer, intent(in) :: level, total_sites, samples_level
    integer(ip), intent(in) :: total_comb
    integer, intent(inout) :: config(:)
    real(dp), intent(in) :: temperature
    logical, intent(in) :: force_restart_accept
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit
    
    integer :: current_subset(max(1, level))
    integer :: trial_subset(max(1, level))
    integer :: best_subset(max(1, level))
    integer :: canonical_buf(max(1, level))
    integer :: canonical_subset_buf(max(1, level))
    integer :: best_count
    integer :: accept_count, max_trials, remove_idx, add_site
    integer :: i_attempt
    integer :: skipped
    integer :: tries
    real(dp) :: beta, delta_e, rand_num
    real(dp) :: current_energy, current_low, current_high
    real(dp) :: trial_energy, trial_low, trial_high
    real(dp) :: best_energy
    real(dp), pointer :: energies(:)
    real(dp), pointer :: energies_low(:)
    real(dp), pointer :: energies_high(:)
    logical :: accepted, valid_trial
    integer :: trace_unit
    integer :: best_step
    integer :: burn_keep, burn_start
    logical :: best_in_subset
    real(dp), parameter :: restart_move_prob = 0.01_dp
    logical :: use_restart_move
    integer :: restart_attempts, restart_accepts
    integer :: flip_attempts, flip_accepts
    integer :: calib_best_idx
    integer :: burn_sample_count, unique_count
    integer :: sample_idx, existing
    integer :: nop, npos
    integer, pointer :: sampled_subsets(:,:)
    real(dp), pointer :: low_contribs(:,:)
    real(dp), pointer :: high_contribs(:,:)
    integer, pointer :: accept_attempt(:)
    real(dp), pointer :: gulp_energies(:)
    integer, pointer :: canonical_subset(:)
    integer, pointer :: eqmatrix(:,:)
    integer, allocatable :: unique_subsets(:,:), unique_deg(:)
    real(dp), allocatable :: unique_low_eval(:), unique_high_eval(:)
    real(dp) :: low_contrib_tmp(4), high_contrib_tmp(4)
    logical :: gulp_success
    real(dp) :: sum_energy, sumsq_energy, mean_energy, std_energy
    
    if (level == 0) then
        call exhaustive_level(level, total_sites, config, temperature, total_comb, 1, use_parallel, &
        summary_unit, summary_txt_unit)
        return
    end if
    
    call mc_workspace_checkout(level, samples_level, total_sites, sampled_subsets, energies, energies_low, energies_high, &
         low_contribs, high_contribs, accept_attempt, gulp_energies, canonical_subset)
    energies = 0.0_dp
    energies_low = 0.0_dp
    energies_high = 0.0_dp
    if (level > 0) then
        sampled_subsets = 0
        low_contribs = 0.0_dp
        high_contribs = 0.0_dp
    end if
    accept_attempt = 0
    gulp_energies = 0.0_dp
    canonical_subset = 0
    calib_best_idx = 0
    
    beta = 1.0_dp / (kB_eVk * temperature)
    skipped = 0
    
    call random_subset(total_sites, level, current_subset)
    config = 1
    config(current_subset(1:level)) = 2
    call calculate_structure_energy(config, total_sites, current_energy, energy_low_side=current_low, &
    energy_high_side=current_high, low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)
    best_energy = current_energy
    best_subset = current_subset
    best_count = level
    accept_count = 1
    energies(accept_count) = current_energy
    energies_low(accept_count) = current_low
    energies_high(accept_count) = current_high
    if (level > 0) then
        sampled_subsets(:, accept_count) = current_subset(1:level)
        low_contribs(:, accept_count) = low_contrib_tmp
        high_contribs(:, accept_count) = high_contrib_tmp
    end if
    best_step = accept_count
    call open_mc_trace_file(level, trace_unit)
    if (trace_unit /= 0) then
        call write_mc_trace_step(trace_unit, accept_count, 0, current_energy, current_low, current_high)
        flush(trace_unit)
    end if
    
    max_trials = samples_level * 100
    i_attempt = 0
    restart_attempts = 0
    restart_accepts = 0
    flip_attempts = 0
    flip_accepts = 0
    
    do while (accept_count < samples_level .and. i_attempt < max_trials)
        i_attempt = i_attempt + 1
        trial_subset = current_subset
        
        call random_number(rand_num)
        use_restart_move = (rand_num < restart_move_prob)
        
        if (use_restart_move) then
            call random_subset(total_sites, level, trial_subset)
            if (all(trial_subset(1:level) == current_subset(1:level))) cycle
            restart_attempts = restart_attempts + 1
        else
            call random_number(rand_num)
            remove_idx = int(rand_num * real(level, dp)) + 1
            
            tries = 0
            valid_trial = .false.
            do while (.not. valid_trial .and. tries < total_sites * 2)
                tries = tries + 1
                call random_number(rand_num)
                add_site = int(rand_num * real(total_sites, dp)) + 1
                if (.not. any(trial_subset == add_site)) then
                    trial_subset(remove_idx) = add_site
                    call sort_int_ascending(trial_subset, level)
                    valid_trial = .true.
                end if
            end do
            if (.not. valid_trial) cycle
            flip_attempts = flip_attempts + 1
        end if
        
        config = 1
        config(trial_subset(1:level)) = 2
        call calculate_structure_energy(config, total_sites, trial_energy, energy_low_side=trial_low, &
        energy_high_side=trial_high, low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)
        
        if (use_restart_move .and. force_restart_accept) then
            accepted = .true.
        else
            delta_e = trial_energy - current_energy
            if (delta_e <= 0.0_dp) then
                accepted = .true.
            else
                call random_number(rand_num)
                accepted = (rand_num < exp(-beta * delta_e))
            end if
        end if
        
        if (accepted) then
            accept_count = accept_count + 1
            current_subset = trial_subset
            current_energy = trial_energy
            current_low = trial_low
            current_high = trial_high
            if (use_restart_move) then
                restart_accepts = restart_accepts + 1
            else
                flip_accepts = flip_accepts + 1
            end if
            energies(accept_count) = current_energy
            energies_low(accept_count) = current_low
            energies_high(accept_count) = current_high
            if (level > 0) then
                sampled_subsets(:, accept_count) = current_subset(1:level)
                low_contribs(:, accept_count) = low_contrib_tmp
                high_contribs(:, accept_count) = high_contrib_tmp
            end if
            if (current_energy < best_energy) then
                best_energy = current_energy
                best_subset = current_subset
                best_count = level
                best_step = accept_count
            end if
            if (trace_unit /= 0) then
                call write_mc_trace_step(trace_unit, accept_count, i_attempt, current_energy, current_low, current_high)
                flush(trace_unit)
            end if
        else
            skipped = skipped + 1
        end if
    end do
    
    burn_keep = max(1, ceiling(0.8_dp * real(accept_count, dp)))
    burn_start = max(1, accept_count - burn_keep + 1)
    best_in_subset = (best_step >= burn_start)
    if (burn_start > 1) then
        write(*,'(A,I0)') 'Metropolis: configurations discarded during burn-in = ', burn_start - 1
        flush(output_unit)
    end if
    
    if (level > 0 .and. accept_count >= burn_start) then
        if (external_calibration_is_enabled()) then
            call attempt_calibration_from_samples(level, total_sites, burn_start, accept_count, sampled_subsets, &
            energies, energies_low, energies_high, low_contribs, high_contribs, calib_best_idx)
            if (calib_best_idx > 0) then
                best_energy = energies(calib_best_idx)
                if (level > 0) best_subset(1:level) = sampled_subsets(:, calib_best_idx)
                best_count = level
                best_step = calib_best_idx
                best_in_subset = (best_step >= burn_start)
            end if
        end if
    end if

    gulp_success = .false.
    if (external_calibration_is_enabled() .and. level > 0 .and. accept_count >= burn_start) then
        burn_sample_count = accept_count - burn_start + 1
        if (burn_sample_count > 0) then
            call symmetry_get_matrix(eqmatrix, nop, npos)
            if (associated(eqmatrix)) then
                allocate(unique_subsets(level, burn_sample_count))
                allocate(unique_deg(burn_sample_count))
                allocate(unique_low_eval(burn_sample_count))
                allocate(unique_high_eval(burn_sample_count))
                unique_subsets = 0
                unique_deg = 0
                unique_low_eval = huge(1.0_dp)
                unique_high_eval = huge(1.0_dp)
                unique_count = 0
                do sample_idx = burn_start, accept_count
                    canonical_buf(1:level) = sampled_subsets(:, sample_idx)
                    call canonicalize_subset(canonical_buf, level, eqmatrix, nop, canonical_subset_buf)
                    existing = find_subset_index(canonical_subset_buf, level, unique_subsets, unique_count)
                    if (existing == 0) then
                        unique_count = unique_count + 1
                        unique_subsets(:, unique_count) = canonical_subset_buf(1:level)
                        unique_deg(unique_count) = 1
                    else
                        unique_deg(existing) = unique_deg(existing) + 1
                    end if
                    if (abs(energies_low(sample_idx)) < huge(1.0_dp) * 0.5_dp) then
                        unique_low_eval(existing) = min(unique_low_eval(existing), energies_low(sample_idx))
                    end if
                    if (abs(energies_high(sample_idx)) < huge(1.0_dp) * 0.5_dp) then
                        unique_high_eval(existing) = min(unique_high_eval(existing), energies_high(sample_idx))
                    end if
                end do
                call write_monitored_outsod(level, total_sites, unique_count, unique_subsets(:,1:unique_count), &
                    unique_deg(1:unique_count))
                if (unique_count > 0) then
                    gulp_energies(1:unique_count) = 0.0_dp
                    call evaluate_subsets_with_engine(level, total_sites, unique_count, unique_subsets(:,1:unique_count), config, &
                        gulp_energies, gulp_success)
                    if (gulp_success) then
                        call update_level_blend_override_from_samples(level, total_sites, unique_low_eval(1:unique_count), &
                            unique_high_eval(1:unique_count), gulp_energies(1:unique_count), unique_count, 'Metropolis external sample', &
                            level_prefix='mc')
                        do sample_idx = burn_start, accept_count
                            canonical_buf(1:level) = sampled_subsets(:, sample_idx)
                            call canonicalize_subset(canonical_buf, level, eqmatrix, nop, canonical_subset_buf)
                            existing = find_subset_index(canonical_subset_buf, level, unique_subsets, unique_count)
                            if (existing > 0) then
                                energies(sample_idx) = gulp_energies(existing)
                            end if
                        end do
                        best_energy = minval(energies(burn_start:accept_count))
                        do sample_idx = burn_start, accept_count
                            if (energies(sample_idx) == best_energy) then
                                best_subset(1:level) = sampled_subsets(:, sample_idx)
                                best_step = sample_idx
                                best_in_subset = (best_step >= burn_start)
                                exit
                            end if
                        end do
                        best_count = level
                        sum_energy = sum(energies(burn_start:accept_count))
                        sumsq_energy = sum(energies(burn_start:accept_count)**2)
                        mean_energy = sum_energy / real(burn_sample_count, dp)
                        std_energy = sqrt(max(0.0_dp, (sumsq_energy / real(burn_sample_count, dp)) - mean_energy**2))
                        write(*,'(A,I0,A,A,A,F16.6,A,F16.6)') 'Level ', level, ': ', trim(calibration_engine_name()), ' mean (Metropolis) = ', mean_energy, ' eV, deviation = ', std_energy
                        flush(output_unit)
                    else
                        write(*,'(A,I0,A,A,A)') 'Level ', level, ': Metropolis samples could not be evaluated with ', trim(calibration_engine_name()), '.'
                        flush(output_unit)
                    end if
                end if
                if (allocated(unique_subsets)) deallocate(unique_subsets)
                if (allocated(unique_deg)) deallocate(unique_deg)
                if (allocated(unique_low_eval)) deallocate(unique_low_eval)
                if (allocated(unique_high_eval)) deallocate(unique_high_eval)
            else
                write(*,'(A,A,A)') 'Warning: EQMATRIX not available to evaluate Metropolis results with ', trim(calibration_engine_name()), '.'
                flush(output_unit)
            end if
            nullify(eqmatrix)
        end if
    end if
    
    call summarize_level(level, total_sites, temperature, total_comb, real(accept_count, dp), &
    energies(burn_start:accept_count), best_energy, best_subset, best_count, skipped, &
    energies_low(burn_start:accept_count), energies_high(burn_start:accept_count), &
    use_parallel, summary_unit, summary_txt_unit, best_in_subset, &
    restart_attempts, restart_accepts, flip_attempts, flip_accepts)
    
    if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
    nullify(sampled_subsets)
    nullify(low_contribs)
    nullify(high_contribs)
    nullify(accept_attempt)
    nullify(gulp_energies)
    nullify(canonical_subset)
    nullify(energies)
    nullify(energies_low)
    nullify(energies_high)
end subroutine metropolis_level

! Runs parallel tempering (replica exchange) at fixed substitution level.
! Ejecuta parallel tempering (intercambio de réplicas) a nivel de sustitución fijo.
subroutine tempering_level(level, total_sites, config, temperature, samples_level, total_comb, force_restart_accept, use_parallel, &
    summary_unit, summary_txt_unit, replica_count, swap_every, max_temp_factor)
    integer, intent(in) :: level, total_sites, samples_level, replica_count, swap_every
    integer(ip), intent(in) :: total_comb
    integer, intent(inout) :: config(:)
    real(dp), intent(in) :: temperature, max_temp_factor
    logical, intent(in) :: force_restart_accept
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit

    integer :: best_subset(max(1, level))
    integer :: trial_subset(max(1, level))
    integer :: canonical_buf(max(1, level))
    integer :: canonical_subset_buf(max(1, level))
    integer :: best_count
    integer :: trace_unit, swap_unit
    integer :: burn_keep, burn_start, burn_sample_count, unique_count
    integer :: record_count, sweep_idx, rep_idx, local_step
    integer :: remove_idx, add_site, tries
    integer :: sample_idx, existing, nop, npos
    integer :: best_step, calib_best_idx
    integer :: swap_phase_start, left_idx, right_idx
    integer :: swap_attempts, swap_accepts
    integer, pointer :: sampled_subsets(:,:)
    real(dp), pointer :: energies(:)
    real(dp), pointer :: energies_low(:)
    real(dp), pointer :: energies_high(:)
    real(dp), pointer :: low_contribs(:,:)
    real(dp), pointer :: high_contribs(:,:)
    integer, pointer :: accept_attempt(:)
    real(dp), pointer :: gulp_energies(:)
    integer, pointer :: canonical_subset(:)
    integer, pointer :: eqmatrix(:,:)
    integer, allocatable :: unique_subsets(:,:), unique_deg(:)
    integer, allocatable :: per_temp_unique_subsets(:,:), per_temp_unique_deg(:)
    integer, allocatable :: all_unique_subsets(:,:), all_unique_deg(:)
    integer, allocatable :: replica_subsets(:,:)
    integer, allocatable :: trace_units(:)
    integer, allocatable :: tempering_subsets(:,:,:)
    integer, allocatable :: aggregate_subsets(:,:)
    integer, allocatable :: local_attempts(:), local_accepts(:)
    integer, allocatable :: restart_attempts(:), restart_accepts(:)
    integer, allocatable :: flip_attempts(:), flip_accepts(:)
    real(dp), allocatable :: betas(:), temps(:)
    real(dp), allocatable :: replica_energies(:), replica_low(:), replica_high(:)
    real(dp), allocatable :: replica_low_contrib(:,:), replica_high_contrib(:,:)
    real(dp), allocatable :: unique_low_eval(:), unique_high_eval(:)
    real(dp) :: current_energy, current_low, current_high
    real(dp) :: trial_energy, trial_low, trial_high
    real(dp) :: best_energy, delta_e, rand_num, delta_swap
    real(dp) :: swap_ratio
    real(dp) :: low_contrib_tmp(4), high_contrib_tmp(4)
    logical :: accepted, valid_trial, use_restart_move
    logical :: best_in_subset, gulp_success

    if (level == 0) then
        call exhaustive_level(level, total_sites, config, temperature, total_comb, 1, use_parallel, &
        summary_unit, summary_txt_unit)
        return
    end if

    if (replica_count < 2) then
        call metropolis_level(level, total_sites, config, temperature, samples_level, total_comb, force_restart_accept, &
            use_parallel, summary_unit, summary_txt_unit)
        return
    end if

    call mc_workspace_checkout(level, samples_level, total_sites, sampled_subsets, energies, energies_low, energies_high, &
         low_contribs, high_contribs, accept_attempt, gulp_energies, canonical_subset)
    energies = 0.0_dp
    energies_low = 0.0_dp
    energies_high = 0.0_dp
    sampled_subsets = 0
    low_contribs = 0.0_dp
    high_contribs = 0.0_dp
    accept_attempt = 0
    gulp_energies = 0.0_dp
    canonical_subset = 0
    calib_best_idx = 0

    allocate(temps(replica_count), betas(replica_count))
    allocate(replica_subsets(level, replica_count))
    allocate(tempering_subsets(level, samples_level, replica_count))
    allocate(replica_energies(replica_count), replica_low(replica_count), replica_high(replica_count))
    allocate(replica_low_contrib(4, replica_count), replica_high_contrib(4, replica_count))
    allocate(local_attempts(replica_count), local_accepts(replica_count))
    allocate(restart_attempts(replica_count), restart_accepts(replica_count))
    allocate(flip_attempts(replica_count), flip_accepts(replica_count))
    allocate(trace_units(replica_count))

    call build_tempering_schedule(temperature, replica_count, max_temp_factor, temps)
    betas = 1.0_dp / (kB_eVk * temps)
    local_attempts = 0
    local_accepts = 0
    restart_attempts = 0
    restart_accepts = 0
    flip_attempts = 0
    flip_accepts = 0
    swap_attempts = 0
    swap_accepts = 0
    best_energy = huge(1.0_dp)
    best_subset = 0
    best_count = level
    best_step = 0
    trace_units = 0
    tempering_subsets = 0

    call write_tempering_temperature_file(level, temps)
    call open_mc_trace_file(level, trace_unit)
    call open_tempering_trace_files(level, temps, trace_units)
    call open_tempering_swap_file(level, swap_unit)

    do rep_idx = 1, replica_count
        call random_subset(total_sites, level, replica_subsets(:, rep_idx))
        config = 1
        config(replica_subsets(1:level, rep_idx)) = 2
        call calculate_structure_energy(config, total_sites, current_energy, energy_low_side=current_low, &
            energy_high_side=current_high, low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)
        replica_energies(rep_idx) = current_energy
        replica_low(rep_idx) = current_low
        replica_high(rep_idx) = current_high
        replica_low_contrib(:, rep_idx) = low_contrib_tmp
        replica_high_contrib(:, rep_idx) = high_contrib_tmp
        tempering_subsets(:, 1, rep_idx) = replica_subsets(:, rep_idx)
    end do

    record_count = 1
    energies(record_count) = replica_energies(1)
    energies_low(record_count) = replica_low(1)
    energies_high(record_count) = replica_high(1)
    sampled_subsets(:, record_count) = replica_subsets(:, 1)
    low_contribs(:, record_count) = replica_low_contrib(:, 1)
    high_contribs(:, record_count) = replica_high_contrib(:, 1)
    accept_attempt(record_count) = 0
    best_energy = energies(record_count)
    best_subset(1:level) = sampled_subsets(:, record_count)
    best_step = record_count
    if (trace_unit /= 0) then
        call write_mc_trace_step(trace_unit, record_count, 0, energies(record_count), energies_low(record_count), energies_high(record_count))
        flush(trace_unit)
    end if
    call write_tempering_trace_step_all(trace_units, replica_count, record_count, 0, replica_energies, replica_low, replica_high)

    sweep_idx = 0
    swap_phase_start = 1
    do while (record_count < samples_level)
        sweep_idx = sweep_idx + 1

        do rep_idx = 1, replica_count
            do local_step = 1, swap_every
                trial_subset = replica_subsets(:, rep_idx)

                call random_number(rand_num)
                use_restart_move = (rand_num < 0.01_dp)
                valid_trial = .false.

                if (use_restart_move) then
                    call random_subset(total_sites, level, trial_subset)
                    if (.not. all(trial_subset(1:level) == replica_subsets(1:level, rep_idx))) then
                        valid_trial = .true.
                        restart_attempts(rep_idx) = restart_attempts(rep_idx) + 1
                    end if
                else
                    call random_number(rand_num)
                    remove_idx = int(rand_num * real(level, dp)) + 1
                    tries = 0
                    do while (.not. valid_trial .and. tries < total_sites * 2)
                        tries = tries + 1
                        call random_number(rand_num)
                        add_site = int(rand_num * real(total_sites, dp)) + 1
                        if (.not. any(trial_subset == add_site)) then
                            trial_subset(remove_idx) = add_site
                            call sort_int_ascending(trial_subset, level)
                            valid_trial = .true.
                        end if
                    end do
                    if (valid_trial) flip_attempts(rep_idx) = flip_attempts(rep_idx) + 1
                end if

                if (.not. valid_trial) cycle
                local_attempts(rep_idx) = local_attempts(rep_idx) + 1

                config = 1
                config(trial_subset(1:level)) = 2
                call calculate_structure_energy(config, total_sites, trial_energy, energy_low_side=trial_low, &
                    energy_high_side=trial_high, low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)

                if (use_restart_move .and. force_restart_accept) then
                    accepted = .true.
                else
                    delta_e = trial_energy - replica_energies(rep_idx)
                    if (delta_e <= 0.0_dp) then
                        accepted = .true.
                    else
                        call random_number(rand_num)
                        accepted = (rand_num < exp(-betas(rep_idx) * delta_e))
                    end if
                end if

                if (accepted) then
                    local_accepts(rep_idx) = local_accepts(rep_idx) + 1
                    if (use_restart_move) then
                        restart_accepts(rep_idx) = restart_accepts(rep_idx) + 1
                    else
                        flip_accepts(rep_idx) = flip_accepts(rep_idx) + 1
                    end if
                    replica_subsets(:, rep_idx) = trial_subset(1:level)
                    replica_energies(rep_idx) = trial_energy
                    replica_low(rep_idx) = trial_low
                    replica_high(rep_idx) = trial_high
                    replica_low_contrib(:, rep_idx) = low_contrib_tmp
                    replica_high_contrib(:, rep_idx) = high_contrib_tmp
                end if
            end do
        end do

        do left_idx = swap_phase_start, replica_count - 1, 2
            right_idx = left_idx + 1
            delta_swap = (betas(left_idx) - betas(right_idx)) * (replica_energies(right_idx) - replica_energies(left_idx))
            swap_attempts = swap_attempts + 1
            accepted = .false.
            if (delta_swap >= 0.0_dp) then
                accepted = .true.
            else
                call random_number(rand_num)
                accepted = (rand_num < exp(delta_swap))
            end if
            if (accepted) then
                swap_accepts = swap_accepts + 1
                call swap_replica_states(level, replica_subsets(:, left_idx), replica_subsets(:, right_idx), &
                    replica_energies(left_idx), replica_energies(right_idx), replica_low(left_idx), replica_low(right_idx), &
                    replica_high(left_idx), replica_high(right_idx), replica_low_contrib(:, left_idx), &
                    replica_low_contrib(:, right_idx), replica_high_contrib(:, left_idx), replica_high_contrib(:, right_idx))
            end if
            call write_tempering_swap_step(swap_unit, sweep_idx, left_idx, right_idx, temps(left_idx), temps(right_idx), &
                delta_swap, accepted)
        end do
        if (swap_phase_start == 1) then
            swap_phase_start = 2
        else
            swap_phase_start = 1
        end if

        record_count = record_count + 1
        energies(record_count) = replica_energies(1)
        energies_low(record_count) = replica_low(1)
        energies_high(record_count) = replica_high(1)
        sampled_subsets(:, record_count) = replica_subsets(:, 1)
        low_contribs(:, record_count) = replica_low_contrib(:, 1)
        high_contribs(:, record_count) = replica_high_contrib(:, 1)
        accept_attempt(record_count) = sweep_idx
        if (energies(record_count) < best_energy) then
            best_energy = energies(record_count)
            best_subset(1:level) = sampled_subsets(:, record_count)
            best_step = record_count
        end if
        if (trace_unit /= 0) then
            call write_mc_trace_step(trace_unit, record_count, sweep_idx, energies(record_count), energies_low(record_count), energies_high(record_count))
            flush(trace_unit)
        end if
        call write_tempering_trace_step_all(trace_units, replica_count, record_count, sweep_idx, replica_energies, replica_low, replica_high)
        do rep_idx = 1, replica_count
            tempering_subsets(:, record_count, rep_idx) = replica_subsets(:, rep_idx)
        end do
    end do

    burn_keep = max(1, ceiling(0.8_dp * real(record_count, dp)))
    burn_start = max(1, record_count - burn_keep + 1)
    best_in_subset = (best_step >= burn_start)
    if (burn_start > 1) then
        write(*,'(A,I0)') 'Tempering: configurations discarded during burn-in = ', burn_start - 1
        flush(output_unit)
    end if

    if (external_calibration_is_enabled()) then
        call attempt_calibration_from_samples(level, total_sites, burn_start, record_count, sampled_subsets, &
            energies, energies_low, energies_high, low_contribs, high_contribs, calib_best_idx)
        if (calib_best_idx > 0) then
            best_energy = energies(calib_best_idx)
            best_subset(1:level) = sampled_subsets(:, calib_best_idx)
            best_count = level
            best_step = calib_best_idx
            best_in_subset = (best_step >= burn_start)
        end if
    end if

    gulp_success = .false.
    burn_sample_count = record_count - burn_start + 1
    if (external_calibration_is_enabled() .and. burn_sample_count > 0) then
        call symmetry_get_matrix(eqmatrix, nop, npos)
        if (associated(eqmatrix)) then
            allocate(per_temp_unique_subsets(level, burn_sample_count))
            allocate(per_temp_unique_deg(burn_sample_count))
            do rep_idx = 1, replica_count
                call collect_unique_canonical_subsets(level, burn_sample_count, tempering_subsets(:, burn_start:record_count, rep_idx), &
                    eqmatrix, nop, per_temp_unique_subsets, per_temp_unique_deg, unique_count)
                call write_monitored_outsod(level, total_sites, unique_count, per_temp_unique_subsets(:,1:unique_count), &
                    per_temp_unique_deg(1:unique_count), output_name=tempering_monitored_name(rep_idx))
            end do

            allocate(aggregate_subsets(level, burn_sample_count * replica_count))
            do rep_idx = 1, replica_count
                do sample_idx = 1, burn_sample_count
                    aggregate_subsets(:, (rep_idx - 1) * burn_sample_count + sample_idx) = &
                        tempering_subsets(:, burn_start + sample_idx - 1, rep_idx)
                end do
            end do
            allocate(all_unique_subsets(level, burn_sample_count * replica_count))
            allocate(all_unique_deg(burn_sample_count * replica_count))
            call collect_unique_canonical_subsets(level, burn_sample_count * replica_count, aggregate_subsets, eqmatrix, nop, &
                all_unique_subsets, all_unique_deg, unique_count)
            call write_monitored_outsod(level, total_sites, unique_count, all_unique_subsets(:,1:unique_count), &
                all_unique_deg(1:unique_count), output_name='monitored_OUTSOD_ALL_TEMPERATURES')

            allocate(unique_subsets(level, burn_sample_count))
            allocate(unique_deg(burn_sample_count))
            allocate(unique_low_eval(burn_sample_count))
            allocate(unique_high_eval(burn_sample_count))
            call collect_unique_canonical_subsets(level, burn_sample_count, sampled_subsets(:, burn_start:record_count), eqmatrix, nop, &
                unique_subsets, unique_deg, unique_count)
            unique_low_eval = huge(1.0_dp)
            unique_high_eval = huge(1.0_dp)
            do sample_idx = burn_start, record_count
                canonical_buf(1:level) = sampled_subsets(:, sample_idx)
                call canonicalize_subset(canonical_buf, level, eqmatrix, nop, canonical_subset_buf)
                existing = find_subset_index(canonical_subset_buf, level, unique_subsets, unique_count)
                if (existing <= 0) cycle
                if (abs(energies_low(sample_idx)) < huge(1.0_dp) * 0.5_dp) then
                    unique_low_eval(existing) = min(unique_low_eval(existing), energies_low(sample_idx))
                end if
                if (abs(energies_high(sample_idx)) < huge(1.0_dp) * 0.5_dp) then
                    unique_high_eval(existing) = min(unique_high_eval(existing), energies_high(sample_idx))
                end if
            end do
            call write_monitored_outsod(level, total_sites, unique_count, unique_subsets(:,1:unique_count), unique_deg(1:unique_count))
            if (unique_count > 0) then
                gulp_energies(1:unique_count) = 0.0_dp
                call evaluate_subsets_with_engine(level, total_sites, unique_count, unique_subsets(:,1:unique_count), config, &
                    gulp_energies, gulp_success)
                if (gulp_success) then
                    call update_level_blend_override_from_samples(level, total_sites, unique_low_eval(1:unique_count), &
                        unique_high_eval(1:unique_count), gulp_energies(1:unique_count), unique_count, 'Tempering external sample', &
                        level_prefix='mc')
                    do sample_idx = burn_start, record_count
                        canonical_buf(1:level) = sampled_subsets(:, sample_idx)
                        call canonicalize_subset(canonical_buf, level, eqmatrix, nop, canonical_subset_buf)
                        existing = find_subset_index(canonical_subset_buf, level, unique_subsets, unique_count)
                        if (existing > 0) energies(sample_idx) = gulp_energies(existing)
                    end do
                    best_energy = minval(energies(burn_start:record_count))
                    do sample_idx = burn_start, record_count
                        if (energies(sample_idx) == best_energy) then
                            best_subset(1:level) = sampled_subsets(:, sample_idx)
                            best_step = sample_idx
                            best_in_subset = (best_step >= burn_start)
                            exit
                        end if
                    end do
                    best_count = level
                else
                    write(*,'(A,I0,A,A,A)') 'Level ', level, ': tempering samples could not be evaluated with ', trim(calibration_engine_name()), '.'
                    flush(output_unit)
                end if
            end if
            if (allocated(per_temp_unique_subsets)) deallocate(per_temp_unique_subsets)
            if (allocated(per_temp_unique_deg)) deallocate(per_temp_unique_deg)
            if (allocated(aggregate_subsets)) deallocate(aggregate_subsets)
            if (allocated(all_unique_subsets)) deallocate(all_unique_subsets)
            if (allocated(all_unique_deg)) deallocate(all_unique_deg)
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            if (allocated(unique_deg)) deallocate(unique_deg)
            if (allocated(unique_low_eval)) deallocate(unique_low_eval)
            if (allocated(unique_high_eval)) deallocate(unique_high_eval)
        else
            write(*,'(A,A,A)') 'Warning: EQMATRIX not available to evaluate tempering results with ', trim(calibration_engine_name()), '.'
            flush(output_unit)
        end if
        nullify(eqmatrix)
    end if

    call summarize_level(level, total_sites, temperature, total_comb, real(record_count, dp), energies(burn_start:record_count), &
        best_energy, best_subset, best_count, max(0, local_attempts(1) - local_accepts(1)), &
        energies_low(burn_start:record_count), energies_high(burn_start:record_count), &
        use_parallel, summary_unit, summary_txt_unit, best_in_subset, &
        restart_attempts(1), restart_accepts(1), flip_attempts(1), flip_accepts(1))

    if (swap_attempts > 0) then
        swap_ratio = real(swap_accepts, dp) / real(swap_attempts, dp)
    else
        swap_ratio = 0.0_dp
    end if
    write(*,'(A,I8,A,I8,A,F8.4)') 'Tempering swaps accepted: ', swap_accepts, ' / ', swap_attempts, '  ratio: ', swap_ratio
    flush(output_unit)

    if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
    call close_tempering_trace_files(trace_units)
    call close_tempering_swap_file(swap_unit)
    deallocate(temps, betas, replica_subsets, replica_energies, replica_low, replica_high)
    deallocate(replica_low_contrib, replica_high_contrib)
    deallocate(local_attempts, local_accepts, restart_attempts, restart_accepts, flip_attempts, flip_accepts)
    deallocate(trace_units)
    deallocate(tempering_subsets)
    nullify(sampled_subsets, energies, energies_low, energies_high, low_contribs, high_contribs, accept_attempt, gulp_energies, canonical_subset)
end subroutine tempering_level

! Attempts to calibrate side-specific energy estimates from sampled configurations.
! Intenta calibrar estimaciones de energía específicas de cada lado a partir de configuraciones muestreadas.
subroutine attempt_calibration_from_samples(level, total_sites, sample_start, sample_end, subsets, energies, energies_low, energies_high, low_contribs, high_contribs, best_index)
    integer, intent(in) :: level, total_sites, sample_start, sample_end
    integer, intent(in) :: subsets(:,:)
    real(dp), intent(inout) :: energies(:), energies_low(:), energies_high(:)
    real(dp), intent(in) :: low_contribs(:,:), high_contribs(:,:)
    integer, intent(out), optional :: best_index
    
    integer :: sample_count, idx, actual_idx, best_global_idx
    integer :: hole_count
    integer :: nop, npos
    integer, pointer :: eqmatrix(:,:)
    integer, allocatable :: unique_subsets(:,:), unique_deg(:)
    integer, allocatable :: sample_map(:), subset_buf(:), canonical(:), config(:)
    real(dp), allocatable :: unique_low(:), unique_high(:)
    real(dp), allocatable :: unique_low_contrib(:,:), unique_high_contrib(:,:)
    real(dp), allocatable :: unique_low_max(:), unique_high_max(:)
    real(dp) :: val_low, val_high
    real(dp) :: huge_marker
    integer :: existing, unique_count
    logical :: need_low_calib, need_high_calib, do_calibrate
    logical :: have_low, have_high
    real(dp) :: min_low, max_low, min_high, max_high
    integer :: best_loc(1)
    
    if (present(best_index)) best_index = 0
    
    if (.not. external_calibration_is_enabled()) return
    if (level <= 0) return
    if (sample_end < sample_start) return
    
    sample_count = sample_end - sample_start + 1
    if (sample_count <= 0) return
    
    hole_count = total_sites - level
    need_low_calib = (level > get_max_low_order())
    need_high_calib = (hole_count > get_max_high_order())
    if (.not. (need_low_calib .or. need_high_calib)) return
    
    call symmetry_get_matrix(eqmatrix, nop, npos)
    if (.not. associated(eqmatrix)) return
    if (level > size(subsets,1)) then
        nullify(eqmatrix)
        return
    end if
    
    huge_marker = huge(1.0_dp) * 0.5_dp
    allocate(unique_subsets(level, sample_count))
    allocate(unique_low(sample_count))
    allocate(unique_high(sample_count))
    allocate(unique_low_contrib(4, sample_count))
    allocate(unique_high_contrib(4, sample_count))
    allocate(unique_low_max(sample_count))
    allocate(unique_high_max(sample_count))
    allocate(unique_deg(sample_count))
    allocate(sample_map(sample_count))
    allocate(subset_buf(level))
    allocate(canonical(level))
    
    unique_low = huge_marker
    unique_high = huge_marker
    unique_low_max = -huge_marker
    unique_high_max = -huge_marker
    unique_low_contrib = 0.0_dp
    unique_high_contrib = 0.0_dp
    unique_deg = 0
    unique_count = 0
    sample_map = 0
    
    do idx = 1, sample_count
        actual_idx = sample_start + idx - 1
        subset_buf = subsets(1:level, actual_idx)
        call canonicalize_subset(subset_buf, level, eqmatrix, nop, canonical)
        existing = find_subset_index(canonical, level, unique_subsets, unique_count)
        if (existing == 0) then
            unique_count = unique_count + 1
            unique_subsets(:, unique_count) = canonical
            unique_deg(unique_count) = 1
            existing = unique_count
        else
            unique_deg(existing) = unique_deg(existing) + 1
        end if
        sample_map(idx) = existing
        
        val_low = energies_low(actual_idx)
        if (abs(val_low) < huge_marker) then
            if (val_low < unique_low(existing)) then
                unique_low(existing) = val_low
                unique_low_contrib(:, existing) = low_contribs(:, actual_idx)
            end if
            if (val_low > unique_low_max(existing)) unique_low_max(existing) = val_low
        end if
        
        val_high = energies_high(actual_idx)
        if (abs(val_high) < huge_marker) then
            if (val_high < unique_high(existing)) then
                unique_high(existing) = val_high
                unique_high_contrib(:, existing) = high_contribs(:, actual_idx)
            end if
            if (val_high > unique_high_max(existing)) unique_high_max(existing) = val_high
        end if
    end do
    
    have_low = any(abs(unique_low(1:unique_count)) < huge_marker)
    have_high = any(abs(unique_high(1:unique_count)) < huge_marker)
    need_low_calib = need_low_calib .and. have_low
    need_high_calib = need_high_calib .and. have_high
    do_calibrate = need_low_calib .and. need_high_calib
    
    if (have_low) then
        min_low = huge_marker
        max_low = -huge_marker
        do idx = 1, unique_count
            if (abs(unique_low(idx)) < huge_marker) then
                if (unique_low(idx) < min_low) min_low = unique_low(idx)
            end if
            if (unique_low_max(idx) > max_low) max_low = unique_low_max(idx)
        end do
        write(*,'(A,I0,A,F16.6,A,F16.6)') 'Level ', level, ': X-side energy min=', min_low, ', max=', max_low
        flush(output_unit)
    end if
    if (have_high) then
        min_high = huge_marker
        max_high = -huge_marker
        do idx = 1, unique_count
            if (abs(unique_high(idx)) < huge_marker) then
                if (unique_high(idx) < min_high) min_high = unique_high(idx)
            end if
            if (unique_high_max(idx) > max_high) max_high = unique_high_max(idx)
        end do
        write(*,'(A,I0,A,F16.6,A,F16.6)') 'Level ', level, ': Y-side energy min=', min_high, ', max=', max_high
        flush(output_unit)
    end if
    
    if (do_calibrate .and. unique_count >= 5) then
        allocate(config(total_sites))
        config = 1
        call calibrate_level_with_engine(level, total_sites, unique_count, unique_subsets(:,1:unique_count), &
            unique_low_contrib(:,1:unique_count), unique_high_contrib(:,1:unique_count), unique_low(1:unique_count), &
            unique_high(1:unique_count), config, need_low_calib, need_high_calib)
        deallocate(config)
        
        do idx = 1, sample_count
            existing = sample_map(idx)
            if (existing <= 0) cycle
            actual_idx = sample_start + idx - 1
            if (need_low_calib .and. abs(unique_low(existing)) < huge_marker) then
                energies_low(actual_idx) = unique_low(existing)
            end if
            if (need_high_calib .and. abs(unique_high(existing)) < huge_marker) then
                energies_high(actual_idx) = unique_high(existing)
            end if
        end do
        
        ! Recompose total energies using the same shared low/high blending rule.
! Recompone las energías totales usando la misma regla compartida de mezcla low/high.
        do idx = 1, sample_count
            actual_idx = sample_start + idx - 1
            energies(actual_idx) = blend_low_high_energy_level(level, total_sites, energies_low(actual_idx), &
                energies_high(actual_idx))
        end do
        
        ! Propagate the updated best-sample index to the caller if needed.
! Propaga al llamador el índice actualizado de la mejor muestra si hace falta.
        if (present(best_index)) then
            best_loc = minloc(energies(sample_start:sample_end))
            best_global_idx = sample_start + best_loc(1) - 1
            best_index = best_global_idx
        end if
        
        write(*,'(A,I0,A,I0)') 'MC calibration completed at level ', level, ' with ', unique_count
        flush(output_unit)
    else
        if (do_calibrate) then
            write(*,'(A,I0)') 'MC calibration skipped due to insufficient unique configurations at level ', level
        else if (need_low_calib .neqv. need_high_calib) then
            write(*,'(A,I0)') 'MC calibration skipped because only one side requires adjustment at level ', level
        end if
        flush(output_unit)
    end if
    
    nullify(eqmatrix)
    deallocate(unique_subsets, unique_low, unique_high, unique_low_contrib, unique_high_contrib)
    deallocate(unique_low_max, unique_high_max)
    deallocate(unique_deg, sample_map, subset_buf, canonical)
    
end subroutine attempt_calibration_from_samples

! Computes Boltzmann statistics for a level and records them to screen and summaries.
! Calcula estadísticas de Boltzmann para un nivel y las registra en pantalla y en los resúmenes.
subroutine summarize_level(level, total_sites, temperature, total_comb, processed, energies, best_energy, best_positions, &
    best_count, skipped, energies_low, energies_high, use_parallel, summary_unit, summary_txt_unit, &
    best_included, restart_attempts, restart_accepts, flip_attempts, flip_accepts)
    integer, intent(in) :: level, total_sites
    real(dp), intent(in) :: temperature, processed
    integer(ip), intent(in) :: total_comb
    real(dp), intent(in) :: energies(:)
    real(dp), intent(in) :: best_energy
    integer, intent(in) :: best_positions(:)
    integer, intent(in) :: best_count
    integer, intent(in) :: skipped
    real(dp), intent(in), optional :: energies_low(:), energies_high(:)
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit
    logical, intent(in), optional :: best_included
    integer, intent(in), optional :: restart_attempts, restart_accepts
    integer, intent(in), optional :: flip_attempts, flip_accepts
    
    real(dp) :: beta, emin, expected, wsum, variance, prob_best
    real(dp) :: expected_low, expected_high, wsum_low, wsum_high
    real(dp) :: low_ref, high_ref, expected_mix
    real(dp) :: ge_fraction
    real(dp) :: expected_sum, variance_sum
    real(dp) :: expected_low_sum, expected_high_sum
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: weights_low(:), weights_high(:)
    integer :: i, ncomb
    character(len=80), parameter :: separator = '------------------------------------------------------------------------'
    character(len=32) :: total_str
    character(len=32) :: level_str
    logical :: have_low, have_high, valid_low, valid_high
    real(dp) :: low_min, high_min
    character(len=32) :: label_low, label_high
    character(len=32) :: exp_low_str, exp_high_str, mix_exp_str, frac_str
    character(len=32) :: exp_total_str, min_total_str
    character(len=32) :: variance_str
    character(len=32) :: delta_exp_total_str, delta_min_total_str
    character(len=32) :: delta_exp_low_str, delta_min_low_str
    character(len=32) :: delta_exp_high_str, delta_min_high_str
    character(len=32) :: delta_mix_str
    character(len=768) :: csv_line, txt_line
    character(len=32) :: accept_ratio_str
    character(len=32) :: processed_str, rejected_str, stddev_str
    real(dp) :: huge_marker
    real(dp) :: rel_exp_total, rel_min_total
    real(dp) :: rel_exp_low, rel_min_low
    real(dp) :: rel_exp_high, rel_min_high
    real(dp) :: rel_exp_mix
    real(dp) :: accept_ratio, total_trials
    real(dp) :: emin_ref
    logical :: best_sampled
    logical :: have_move_stats
    integer :: restart_attempts_val, restart_accepts_val
    integer :: flip_attempts_val, flip_accepts_val
    real(dp) :: restart_ratio, flip_ratio
    
    if (.not. use_parallel) then
        continue
    end if
    
    ncomb = size(energies)
    beta = 1.0_dp / (kB_eVk * temperature)
    emin = minval(energies)
    best_sampled = .true.
    if (present(best_included)) best_sampled = best_included
    emin_ref = min(emin, best_energy)
    allocate(weights(ncomb))
    wsum = 0.0_dp
    expected_sum = 0.0_dp
    ! $omp parallel do default(shared) private(i) reduction(+:wsum, expected_sum) if(use_parallel)
    do i = 1, ncomb
        weights(i) = exp(-beta * (energies(i) - emin_ref))
        wsum = wsum + weights(i)
        expected_sum = expected_sum + weights(i) * energies(i)
    end do
    ! $omp end parallel do
    if (wsum > 0.0_dp) then
        expected = expected_sum / wsum
    else
        expected = emin
    end if
    variance_sum = 0.0_dp
    if (wsum > 0.0_dp) then
        ! $omp parallel do default(shared) private(i) reduction(+:variance_sum) if(use_parallel)
        do i = 1, ncomb
            variance_sum = variance_sum + weights(i) * (energies(i) - expected)**2
        end do
        ! $omp end parallel do
        variance = variance_sum / wsum
    else
        variance = 0.0_dp
    end if
    if (variance < 0.0_dp) variance = 0.0_dp
    if (wsum > 0.0_dp .and. best_sampled) then
        prob_best = exp(-beta * (best_energy - emin_ref)) / wsum
        if (prob_best > 1.0_dp) prob_best = 1.0_dp
    else
        prob_best = 0.0_dp
    end if
    
    huge_marker = huge(1.0_dp)
    have_low = present(energies_low)
    have_high = present(energies_high)
    valid_low = .false.
    valid_high = .false.
    label_low = ' --'
    label_high = ' --'
    exp_low_str = ' --'
    exp_high_str = ' --'
    mix_exp_str = ' --'
    frac_str = ' --'
    delta_exp_total_str = ' --'
    delta_min_total_str = ' --'
    delta_exp_low_str = ' --'
    delta_min_low_str = ' --'
    delta_exp_high_str = ' --'
    delta_min_high_str = ' --'
    delta_mix_str = ' --'
    expected_low = huge_marker
    expected_high = huge_marker
    low_min = 0.0_dp
    high_min = 0.0_dp
    rel_exp_total = huge_marker
    rel_min_total = huge_marker
    rel_exp_low = huge_marker
    rel_min_low = huge_marker
    rel_exp_high = huge_marker
    rel_min_high = huge_marker
    rel_exp_mix = huge_marker
    
    if (have_low) then
        low_min = minval(energies_low)
        call format_summary_real(low_min, label_low)
        valid_low = .true.
        low_ref = low_min
        allocate(weights_low(ncomb))
        wsum_low = 0.0_dp
        expected_low_sum = 0.0_dp
        ! $omp parallel do default(shared) private(i) reduction(+:wsum_low, expected_low_sum) if(use_parallel)
        do i = 1, ncomb
            weights_low(i) = exp(-beta * (energies_low(i) - low_ref))
            wsum_low = wsum_low + weights_low(i)
            expected_low_sum = expected_low_sum + weights_low(i) * energies_low(i)
        end do
        ! $omp end parallel do
        if (wsum_low > 0.0_dp) then
            expected_low = expected_low_sum / wsum_low
        else
            expected_low = low_ref
        end if
        call format_summary_real(expected_low, exp_low_str)
        deallocate(weights_low)
    end if
    
    if (have_high) then
        if (.not. all(energies_high > 0.5_dp * huge_marker)) then
            high_min = minval(energies_high)
            call format_summary_real(high_min, label_high)
            valid_high = .true.
            high_ref = high_min
            allocate(weights_high(ncomb))
            wsum_high = 0.0_dp
            expected_high_sum = 0.0_dp
            ! $omp parallel do default(shared) private(i) reduction(+:wsum_high, expected_high_sum) if(use_parallel)
            do i = 1, ncomb
                weights_high(i) = exp(-beta * (energies_high(i) - high_ref))
                wsum_high = wsum_high + weights_high(i)
                expected_high_sum = expected_high_sum + weights_high(i) * energies_high(i)
            end do
            ! $omp end parallel do
            if (wsum_high > 0.0_dp) then
                expected_high = expected_high_sum / wsum_high
            else
                expected_high = high_ref
            end if
            call format_summary_real(expected_high, exp_high_str)
            deallocate(weights_high)
        end if
    end if
    
    ge_fraction = 0.0_dp
    if (total_sites > 0) ge_fraction = real(level, dp) / real(total_sites, dp)
    write(frac_str,'(F7.4)') ge_fraction
    
    if (valid_low .and. valid_high) then
        expected_mix = blend_low_high_energy_level(level, total_sites, expected_low, expected_high)
    else if (valid_low) then
        expected_mix = expected_low
    else if (valid_high) then
        expected_mix = expected_high
    else
        expected_mix = expected
    end if
    if (valid_low .or. valid_high) then
        call format_summary_real(expected_mix, mix_exp_str)
    end if
    
    call format_summary_real(expected, exp_total_str)
    call format_summary_real(best_energy, min_total_str)
    call format_summary_real(variance, variance_str)
    rel_exp_total = reference_relative(level, total_sites, expected)
    rel_min_total = reference_relative(level, total_sites, best_energy)
    call format_summary_real(rel_exp_total, delta_exp_total_str)
    call format_summary_real(rel_min_total, delta_min_total_str)
    if (valid_low) then
        rel_exp_low = reference_relative(level, total_sites, expected_low)
        rel_min_low = reference_relative(level, total_sites, low_min)
        call format_summary_real(rel_exp_low, delta_exp_low_str)
        call format_summary_real(rel_min_low, delta_min_low_str)
    end if
    if (valid_high) then
        rel_exp_high = reference_relative(level, total_sites, expected_high)
        rel_min_high = reference_relative(level, total_sites, high_min)
        call format_summary_real(rel_exp_high, delta_exp_high_str)
        call format_summary_real(rel_min_high, delta_min_high_str)
    end if
    if (valid_low .or. valid_high) then
        rel_exp_mix = reference_relative(level, total_sites, expected_mix)
        call format_summary_real(rel_exp_mix, delta_mix_str)
    end if
    write(total_str,'(I0)') total_comb
    write(level_str,'(I0)') level
    total_trials = processed + real(skipped, dp)
    if (total_trials > 0.0_dp) then
        accept_ratio = processed / total_trials
    else
        accept_ratio = 1.0_dp
    end if
    write(accept_ratio_str,'(F12.6)') accept_ratio
    call format_summary_count(processed, processed_str)
    write(rejected_str,'(I0)') max(skipped, 0)
    call format_summary_real(sqrt(variance), stddev_str)
    ! $omp critical(summary_io)
    write(*,'(A)') separator
    write(*,'(A,I0)') 'Substitutions (N): ', level
    write(*,'(A)') 'Total combinations: '//trim(total_str)
    write(*,'(A,A)') 'Configurations processed: ', trim(processed_str)
    write(*,'(A,A)') 'Rejected attempts: ', trim(rejected_str)
    write(*,'(A,A)') 'Minimum energy (eV): ', trim(min_total_str)
    write(*,'(A,A)') 'Expected energy <E> (eV): ', trim(exp_total_str)
    write(*,'(A,A)') 'Variance around <E> (eV^2): ', trim(variance_str)
    write(*,'(A,A)') 'Standard deviation (eV): ', trim(stddev_str)
    write(*,'(A,A)') 'Minimum X-side energy (eV): ', trim(label_low)
    write(*,'(A,A)') 'Expected X-side energy (eV): ', trim(exp_low_str)
    write(*,'(A,A)') 'Minimum Y-side energy (eV): ', trim(label_high)
    write(*,'(A,A)') 'Expected Y-side energy (eV): ', trim(exp_high_str)
    write(*,'(A,A)') 'Expected mixed energy (eV): ', trim(mix_exp_str)
    write(*,'(A,F10.6)') 'Boltzmann probability of minimum: ', prob_best
    write(*,'(A,F8.4)') 'Acceptance ratio: ', accept_ratio
    have_move_stats = present(restart_attempts) .and. present(restart_accepts) .and. &
    present(flip_attempts) .and. present(flip_accepts)
    if (have_move_stats) then
        restart_attempts_val = restart_attempts
        restart_accepts_val = restart_accepts
        flip_attempts_val = flip_attempts
        flip_accepts_val = flip_accepts
        if (restart_attempts_val > 0) then
            restart_ratio = real(restart_accepts_val, dp) / real(restart_attempts_val, dp)
        else
            restart_ratio = 0.0_dp
        end if
        if (flip_attempts_val > 0) then
            flip_ratio = real(flip_accepts_val, dp) / real(flip_attempts_val, dp)
        else
            flip_ratio = 0.0_dp
        end if
        write(*,'(A,I8,A,I8,A,F8.4)') 'Restart accepted: ', restart_accepts_val, ' / ', restart_attempts_val, &
        '  ratio: ', restart_ratio
        write(*,'(A,I8,A,I8,A,F8.4)') 'Flip accepted:    ', flip_accepts_val, ' / ', flip_attempts_val, &
        '  ratio: ', flip_ratio
    end if
    call print_best_positions(best_positions, best_count)
    call save_best_structure_poscar(level, total_sites, best_positions, best_count)
    if (processed < real(total_comb, dp)) then
        write(*,'(A)') 'Warning: not all combinations were covered (sampling applied).'
    end if
    if (level == 0) then
        write(*,'(A)') 'FracY; E_exp(X_side); E_min(X_side); E_exp(Y_side); E_min(Y_side)'
    end if
    write(*,'(A)') trim(frac_str)//'; '//trim(exp_low_str)//'; '//trim(label_low)//'; ' &
    //trim(exp_high_str)//'; '//trim(label_high)
    flush(output_unit)
    if (summary_unit /= 0) then
        csv_line = trim(adjustl(level_str))//';'//trim(adjustl(frac_str))//';'//trim(adjustl(exp_total_str))//';'// &
        trim(adjustl(min_total_str))//';'//trim(adjustl(variance_str))//';'//trim(adjustl(exp_low_str))//';'//trim(adjustl(label_low))//';'// &
        trim(adjustl(exp_high_str))//';'//trim(adjustl(label_high))//';'//trim(adjustl(mix_exp_str))//';'// &
        trim(adjustl(delta_exp_total_str))//';'//trim(adjustl(delta_min_total_str))//';'// &
        trim(adjustl(delta_exp_low_str))//';'//trim(adjustl(delta_min_low_str))//';'// &
        trim(adjustl(delta_exp_high_str))//';'//trim(adjustl(delta_min_high_str))//';'// &
        trim(adjustl(delta_mix_str))//';'//trim(adjustl(accept_ratio_str))
        write(summary_unit,'(A)') trim(csv_line)
        flush(summary_unit)
    end if
    if (summary_txt_unit /= 0) then
        txt_line = trim(adjustl(level_str))//' '//trim(adjustl(frac_str))//' '// &
        trim(adjustl(exp_total_str))//' '//trim(adjustl(min_total_str))//' '//trim(adjustl(variance_str))//' '// &
        trim(adjustl(exp_low_str))//' '//trim(adjustl(label_low))//' '// &
        trim(adjustl(exp_high_str))//' '//trim(adjustl(label_high))//' '// &
        trim(adjustl(mix_exp_str))//' '//trim(adjustl(delta_exp_total_str))//' '// &
        trim(adjustl(delta_min_total_str))//' '//trim(adjustl(delta_exp_low_str))//' '// &
        trim(adjustl(delta_min_low_str))//' '//trim(adjustl(delta_exp_high_str))//' '// &
        trim(adjustl(delta_min_high_str))//' '//trim(adjustl(delta_mix_str))//' '// &
        trim(adjustl(accept_ratio_str))
        write(summary_txt_unit,'(A)') trim(txt_line)
        flush(summary_txt_unit)
    end if
    ! $omp end critical(summary_io)
    deallocate(weights)
end subroutine summarize_level

! Sorts the summary files by substitution level to keep outputs monotonic.
! Ordena los archivos de resumen por nivel de sustitución para mantener las salidas monótonas.
subroutine reorder_summary_outputs()
    call reorder_single_summary(summary_filename, ';')
    call reorder_single_summary(summary_txt_filename, ' ')
end subroutine reorder_summary_outputs

! Loads a summary file, sorts its lines by level, and writes it back.
! Carga un archivo de resumen, ordena sus líneas por nivel y lo vuelve a escribir.
subroutine reorder_single_summary(filename, delimiter)
    character(len=*), intent(in) :: filename
    character(len=1), intent(in) :: delimiter
    integer :: unit, ios, raw_count, idx, level_val
    character(len=512) :: header_line
    character(len=512) :: line
    character(len=512), allocatable :: lines(:)
    integer, allocatable :: levels(:)
    integer :: actual_count
    
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        return
    end if
    
    read(unit,'(A)',iostat=ios) header_line
    if (ios /= 0) then
        close(unit)
        return
    end if
    
    raw_count = 0
    do
        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        raw_count = raw_count + 1
    end do
    
    if (raw_count <= 1) then
        close(unit)
        return
    end if
    
    rewind(unit)
    read(unit,'(A)',iostat=ios) header_line
    if (ios /= 0) then
        close(unit)
        return
    end if
    
    allocate(lines(raw_count))
    allocate(levels(raw_count))
    actual_count = 0
    do
        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        level_val = extract_level_from_line(line, delimiter)
        if (level_val == huge(1)) cycle
        actual_count = actual_count + 1
        lines(actual_count) = trim(line)
        levels(actual_count) = level_val
    end do
    close(unit)
    
    if (actual_count <= 1) then
        deallocate(lines, levels)
        return
    end if
    
    call sort_lines_by_level(levels, lines, actual_count)
    
    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        deallocate(lines, levels)
        return
    end if
    
    write(unit,'(A)') trim(header_line)
    do idx = 1, actual_count
        write(unit,'(A)') trim(lines(idx))
    end do
    close(unit)
    deallocate(lines, levels)
end subroutine reorder_single_summary

! Parses the leading level integer from a delimited summary line.
! Analiza el entero de nivel inicial de una línea de resumen delimitada.
integer function extract_level_from_line(line, delimiter) result(level_val)
character(len=*), intent(in) :: line
character(len=1), intent(in) :: delimiter
character(len=512) :: buffer
character(len=64) :: token
integer :: pos, ios, lt, i

buffer = adjustl(line)
lt = len_trim(buffer)
if (lt == 0) then
    level_val = huge(1)
    return
end if

if (delimiter == ' ') then
    pos = 0
    do i = 1, lt
        if (buffer(i:i) == ' ' .or. buffer(i:i) == char(9)) then
            pos = i
            exit
        end if
    end do
    if (pos == 0) then
        pos = lt + 1
    end if
else
    pos = index(buffer(1:lt), delimiter)
    if (pos == 0) then
        pos = lt + 1
    end if
end if

if (pos <= 1) then
    token = buffer(1:lt)
else
    token = buffer(1:pos-1)
end if

read(token,*,iostat=ios) level_val
if (ios /= 0) then
    level_val = huge(1)
end if
end function extract_level_from_line

! Performs insertion sort on parallel arrays of levels and lines.
! Realiza una ordenación por inserción sobre vectores paralelos de niveles y líneas.
subroutine sort_lines_by_level(levels, lines, n)
    integer, intent(inout) :: levels(:)
    character(len=*), intent(inout) :: lines(:)
    integer, intent(in) :: n
    integer :: i, j, key_level
    character(len=512) :: key_line
    
    do i = 2, n
        key_level = levels(i)
        key_line = lines(i)
        j = i - 1
        do while (j >= 1 .and. levels(j) > key_level)
            levels(j+1) = levels(j)
            lines(j+1) = lines(j)
            j = j - 1
        end do
        levels(j+1) = key_level
        lines(j+1) = key_line
    end do
end subroutine sort_lines_by_level

! Writes the minimum-energy configuration for a level as a VASP POSCAR file.
! Escribe la configuración de energía mínima de un nivel como un archivo POSCAR de VASP.
subroutine save_best_structure_poscar(level, total_sites, best_positions, best_count)
    integer, intent(in) :: level, total_sites, best_count
    integer, intent(in) :: best_positions(:)
    
    integer, allocatable :: best_config(:)
    character(len=512) :: filename
    character(len=32) :: level_dir
    character(len=16) :: config_tag
    logical :: valid_positions
    integer :: i
    
    if (total_sites <= 0) then
        return
    end if
    
    allocate(best_config(total_sites))
    best_config = 1
    valid_positions = .true.
    
    if (level > 0 .and. best_count < level) then
        valid_positions = .false.
    end if
    
    if (valid_positions .and. best_count > 0) then
        do i = 1, best_count
            if (best_positions(i) < 1 .or. best_positions(i) > total_sites) then
                valid_positions = .false.
                exit
            end if
            best_config(best_positions(i)) = 2
        end do
    end if
    
    if (valid_positions) then
        level_dir = format_level_directory('mc', level)
        call ensure_directory_exists(trim(level_dir))
        filename = join_path(trim(level_dir), 'POSCAR.vasp')
        if (mc_have_motif) then
            call write_vasp_file(best_config, total_sites, trim(filename), mc_motif%atoms, mc_motif%natoms)
        else
            call write_vasp_file(best_config, total_sites, trim(filename))
        end if
        write(*,'(A)') 'Minimum POSCAR saved to '//trim(filename)

        if (mc_filer == 1) then
            filename = join_path(trim(level_dir), 'best.gin')
            write(config_tag,'("c",I5.5)') level
            if (mc_have_motif) then
                call write_gulp_output_file(best_config, total_sites, trim(filename), 'template_input.gin', &
                    trim(config_tag), motif=mc_motif)
            else
                call write_gulp_output_file(best_config, total_sites, trim(filename), 'template_input.gin', &
                    trim(config_tag))
            end if
            write(*,'(A)') 'GULP file saved to '//trim(filename)
        end if

        if (mc_filer == 14) then
            write(*,'(A)') 'ASE VASP file saved to '//trim(join_path(trim(level_dir), 'POSCAR.vasp'))
        end if
    else
        write(*,'(A)') 'Warning: could not generate minimum POSCAR due to invalid positions.'
    end if
    flush(output_unit)
    
    deallocate(best_config)
end subroutine save_best_structure_poscar

! Builds a geometric temperature ladder for replica exchange runs.
! Construye una escalera geométrica de temperaturas para replica exchange.
subroutine build_tempering_schedule(base_temperature, replica_count, max_temp_factor, temperatures)
    real(dp), intent(in) :: base_temperature, max_temp_factor
    integer, intent(in) :: replica_count
    real(dp), intent(out) :: temperatures(:)
    integer :: idx
    real(dp) :: exponent

    if (replica_count <= 1) then
        temperatures(1) = base_temperature
        return
    end if

    do idx = 1, replica_count
        exponent = real(idx - 1, dp) / real(replica_count - 1, dp)
        temperatures(idx) = base_temperature * max_temp_factor**exponent
    end do
end subroutine build_tempering_schedule

! Swaps two replica states in-place during replica exchange.
! Intercambia in situ dos estados de réplica durante replica exchange.
subroutine swap_replica_states(level, subset_a, subset_b, energy_a, energy_b, low_a, low_b, high_a, high_b, &
    low_contrib_a, low_contrib_b, high_contrib_a, high_contrib_b)
    integer, intent(in) :: level
    integer, intent(inout) :: subset_a(:), subset_b(:)
    real(dp), intent(inout) :: energy_a, energy_b, low_a, low_b, high_a, high_b
    real(dp), intent(inout) :: low_contrib_a(:), low_contrib_b(:), high_contrib_a(:), high_contrib_b(:)
    integer :: tmp_subset(max(1, level))
    real(dp) :: tmp_energy, tmp_low, tmp_high
    real(dp) :: tmp_low_contrib(4), tmp_high_contrib(4)

    tmp_subset = subset_a(1:level)
    subset_a(1:level) = subset_b(1:level)
    subset_b(1:level) = tmp_subset(1:level)

    tmp_energy = energy_a
    energy_a = energy_b
    energy_b = tmp_energy

    tmp_low = low_a
    low_a = low_b
    low_b = tmp_low

    tmp_high = high_a
    high_a = high_b
    high_b = tmp_high

    tmp_low_contrib = low_contrib_a(1:4)
    low_contrib_a(1:4) = low_contrib_b(1:4)
    low_contrib_b(1:4) = tmp_low_contrib

    tmp_high_contrib = high_contrib_a(1:4)
    high_contrib_a(1:4) = high_contrib_b(1:4)
    high_contrib_b(1:4) = tmp_high_contrib
end subroutine swap_replica_states

! Writes the temperature ladder used by parallel tempering.
! Escribe la escalera de temperaturas usada por parallel tempering.
subroutine write_tempering_temperature_file(level, temperatures)
    integer, intent(in) :: level
    real(dp), intent(in) :: temperatures(:)
    character(len=32) :: level_dir
    character(len=512) :: filename
    integer :: unit_id, ios, idx

    level_dir = format_level_directory('mc', level)
    call ensure_directory_exists(trim(level_dir))
    filename = join_path(trim(level_dir), 'TEMPERING_TEMPERATURES.csv')
    open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(unit_id,'(A)') '#replica;temperature_K'
    do idx = 1, size(temperatures)
        write(unit_id,'(I0,";",F18.8)') idx, temperatures(idx)
    end do
    close(unit_id)
end subroutine write_tempering_temperature_file

! Opens the swap log for a parallel tempering level.
! Abre el registro de intercambios para un nivel con parallel tempering.
subroutine open_tempering_swap_file(level, unit_id)
    integer, intent(in) :: level
    integer, intent(out) :: unit_id
    character(len=32) :: level_dir
    character(len=512) :: filename
    integer :: ios

    level_dir = format_level_directory('mc', level)
    call ensure_directory_exists(trim(level_dir))
    filename = join_path(trim(level_dir), 'TEMPERING_SWAPS.csv')
    open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        unit_id = 0
        return
    end if
    write(unit_id,'(A)') '#sweep;left_replica;right_replica;temp_left_K;temp_right_K;delta;accepted'
end subroutine open_tempering_swap_file

! Appends one swap attempt to the parallel tempering swap log.
! Añade un intento de intercambio al registro de swaps de parallel tempering.
subroutine write_tempering_swap_step(unit_id, sweep_idx, left_idx, right_idx, temp_left, temp_right, delta_val, accepted)
    integer, intent(in) :: unit_id, sweep_idx, left_idx, right_idx
    real(dp), intent(in) :: temp_left, temp_right, delta_val
    logical, intent(in) :: accepted

    if (unit_id == 0) return
    write(unit_id,'(I0,";",I0,";",I0,";",F18.8,";",F18.8,";",ES16.8,";",L1)') &
        sweep_idx, left_idx, right_idx, temp_left, temp_right, delta_val, accepted
end subroutine write_tempering_swap_step

! Closes the swap log for a parallel tempering run.
! Cierra el registro de intercambios de una ejecución con parallel tempering.
subroutine close_tempering_swap_file(unit_id)
    integer, intent(in) :: unit_id

    if (unit_id /= 0) close(unit_id)
end subroutine close_tempering_swap_file

! Opens one trace CSV per temperature slot for a tempering run.
! Abre un CSV de traza por slot de temperatura para una ejecución con tempering.
subroutine open_tempering_trace_files(level, temperatures, unit_ids)
    integer, intent(in) :: level
    real(dp), intent(in) :: temperatures(:)
    integer, intent(inout) :: unit_ids(:)
    character(len=32) :: level_dir
    character(len=512) :: filename
    integer :: idx, ios

    level_dir = format_level_directory('mc', level)
    call ensure_directory_exists(trim(level_dir))
    do idx = 1, min(size(temperatures), size(unit_ids))
        write(filename,'("MC_TRACE_T",I4.4,".csv")') idx
        filename = join_path(trim(level_dir), trim(filename))
        open(newunit=unit_ids(idx), file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,I0,A,I0)') 'Warning: failed to open tempering trace file for N=', level, ', slot=', idx
            unit_ids(idx) = 0
        else
            write(unit_ids(idx),'(A,F18.8)') '#temperature_K=', temperatures(idx)
            write(unit_ids(idx),'(A)') '#step;attempt;energy_eV;energy_low_eV;energy_high_eV'
        end if
    end do
end subroutine open_tempering_trace_files

! Appends one trace entry for every temperature slot in a tempering run.
! Añade una entrada de traza para cada slot de temperatura en una ejecución con tempering.
subroutine write_tempering_trace_step_all(unit_ids, replica_count, step_idx, attempt_idx, energies, energies_low, energies_high)
    integer, intent(in) :: unit_ids(:)
    integer, intent(in) :: replica_count, step_idx, attempt_idx
    real(dp), intent(in) :: energies(:), energies_low(:), energies_high(:)
    integer :: idx, nslots

    nslots = min(replica_count, min(size(unit_ids), min(size(energies), min(size(energies_low), size(energies_high)))))
    do idx = 1, nslots
        call write_mc_trace_step(unit_ids(idx), step_idx, attempt_idx, energies(idx), energies_low(idx), energies_high(idx))
        if (unit_ids(idx) /= 0) flush(unit_ids(idx))
    end do
end subroutine write_tempering_trace_step_all

! Closes all per-temperature trace CSVs for a tempering run.
! Cierra todos los CSV de traza por temperatura de una ejecución con tempering.
subroutine close_tempering_trace_files(unit_ids)
    integer, intent(in) :: unit_ids(:)
    integer :: idx

    do idx = 1, size(unit_ids)
        if (unit_ids(idx) /= 0) close(unit_ids(idx))
    end do
end subroutine close_tempering_trace_files

! Builds the monitored_OUTSOD file name for one tempering temperature slot.
! Construye el nombre del archivo monitored_OUTSOD para un slot térmico de tempering.
pure function tempering_monitored_name(slot_index) result(filename)
    integer, intent(in) :: slot_index
    character(len=64) :: filename

    write(filename,'("monitored_OUTSOD_T",I4.4)') slot_index
end function tempering_monitored_name

! Canonicalizes and deduplicates sampled subsets, counting their occurrences.
! Canonicaliza y elimina duplicados de subconjuntos muestreados, contando sus ocurrencias.
subroutine collect_unique_canonical_subsets(level, sample_count, subsets, eqmatrix, nop, unique_subsets, unique_deg, unique_count)
    integer, intent(in) :: level, sample_count, nop
    integer, intent(in) :: subsets(:,:), eqmatrix(:,:)
    integer, intent(out) :: unique_subsets(:,:), unique_deg(:), unique_count
    integer :: idx, existing
    integer :: subset_buf(max(1, level))
    integer :: canonical_buf(max(1, level))

    unique_subsets = 0
    unique_deg = 0
    unique_count = 0
    if (sample_count <= 0) return

    do idx = 1, sample_count
        if (level > 0) subset_buf(1:level) = subsets(1:level, idx)
        call canonicalize_subset(subset_buf, level, eqmatrix, nop, canonical_buf)
        existing = find_subset_index(canonical_buf, level, unique_subsets, unique_count)
        if (existing == 0) then
            unique_count = unique_count + 1
            unique_subsets(:, unique_count) = canonical_buf(1:level)
            unique_deg(unique_count) = 1
        else
            unique_deg(existing) = unique_deg(existing) + 1
        end if
    end do
end subroutine collect_unique_canonical_subsets

! Creates a trace CSV for Monte Carlo steps at the specified level.
! Crea un CSV de traza para los pasos Monte Carlo en el nivel especificado.
subroutine open_mc_trace_file(level, unit_id)
    integer, intent(in) :: level
    integer, intent(out) :: unit_id
    character(len=512) :: filename
    character(len=32) :: level_dir
    integer :: ios
    
    level_dir = format_level_directory('mc', level)
    call ensure_directory_exists(trim(level_dir))
    filename = join_path(trim(level_dir), 'MC_TRACE.csv')
    open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,I0)') 'Warning: failed to open MC trace file for N=', level
        unit_id = 0
    else
        write(unit_id,'(A)') '#step;attempt;energy_eV;energy_low_eV;energy_high_eV'
    end if
end subroutine open_mc_trace_file

! Appends one Monte Carlo sample entry to the trace CSV if the unit is valid.
! Añade una entrada de muestra Monte Carlo al CSV de traza si la unidad es válida.
subroutine write_mc_trace_step(unit_id, step_idx, attempt_idx, energy, energy_low, energy_high)
    integer, intent(in) :: unit_id
    integer, intent(in) :: step_idx, attempt_idx
    real(dp), intent(in) :: energy, energy_low, energy_high
    
    if (unit_id == 0) return
    write(unit_id,'(I0,";",I0,";",F18.10,";",F18.10,";",F18.10)') step_idx, attempt_idx, energy, energy_low, energy_high
end subroutine write_mc_trace_step

! Closes the trace CSV unit when tracing has finished.
! Cierra la unidad del CSV de traza cuando la traza ha finalizado.
subroutine close_mc_trace_file(unit_id)
    integer, intent(in) :: unit_id
    
    if (unit_id /= 0) close(unit_id)
end subroutine close_mc_trace_file

! Writes the monitored symmetry-unique subsets sampled by MC for a level.
! Escribe los subconjuntos únicos por simetría monitorizados por MC para un nivel.
subroutine write_monitored_outsod(level, total_sites, unique_count, unique_subsets, counts, output_name)
    integer, intent(in) :: level, total_sites, unique_count
    integer, intent(in) :: unique_subsets(:,:)
    integer, intent(in), optional :: counts(:)
    character(len=*), intent(in), optional :: output_name

    character(len=32) :: level_dir
    character(len=512) :: filename
    integer :: unit_outsod
    integer :: idx, site
    integer :: count_value

    level_dir = format_level_directory('mc', level)
    call ensure_directory_exists(trim(level_dir))
    if (present(output_name)) then
        filename = join_path(trim(level_dir), trim(output_name))
    else
        filename = join_path(trim(level_dir), 'monitored_OUTSOD')
    end if

    open(newunit=unit_outsod, file=trim(filename), status='replace', action='write')
    write(unit_outsod,'(I12,"  substitutions in",I12," sites")') level, total_sites
    write(unit_outsod,'(I12,"  monitored configurations")') unique_count
    write(unit_outsod,'(A)') '# second column = monitored occurrences in the MC sample'
    do idx = 1, unique_count
        count_value = 1
        if (present(counts)) count_value = counts(idx)
        write(unit_outsod,'(I6,1X,I12)', advance='no') idx, count_value
        if (level > 0) then
            do site = 1, level
                write(unit_outsod,'(1X,I6)', advance='no') unique_subsets(site, idx)
            end do
        end if
        write(unit_outsod,*)
    end do
    close(unit_outsod)
end subroutine write_monitored_outsod

! Prints the indices of the lowest-energy configuration sampled for a level.
! Imprime los índices de la configuración de menor energía muestreada para un nivel.
subroutine print_best_positions(best_positions, best_count)
    integer, intent(in) :: best_positions(:)
    integer, intent(in) :: best_count
    
    if (best_count <= 0) then
        write(*,'(A)') 'Minimum-energy configuration: base structure'
    else
        write(*,'(A)', advance='no') 'Minimum-energy configuration (positions 1..n): '
        write(*,'(*(1X,I3))') best_positions(1:best_count)
    end if
    flush(output_unit)
end subroutine print_best_positions

! Formats a real value for summaries without overflowing fixed-width fields.
! Formatea un valor real para resúmenes sin desbordar campos de anchura fija.
subroutine format_summary_real(value, text)
    real(dp), intent(in) :: value
    character(len=*), intent(out) :: text
    real(dp) :: abs_value

    abs_value = abs(value)
    if ((abs_value >= 1.0e11_dp) .or. (abs_value > 0.0_dp .and. abs_value < 1.0e-4_dp)) then
        write(text,'(ES22.12E3)') value
    else
        write(text,'(F20.6)') value
    end if
    text = adjustl(text)
end subroutine format_summary_real

! Formats counters as integers when they are effectively whole numbers.
! Formatea contadores como enteros cuando son efectivamente números enteros.
subroutine format_summary_count(value, text)
    real(dp), intent(in) :: value
    character(len=*), intent(out) :: text

    if (abs(value - anint(value)) <= 1.0e-9_dp * max(1.0_dp, abs(value))) then
        write(text,'(I0)') nint(value, kind=ip)
    else
        call format_summary_real(value, text)
    end if
    text = adjustl(text)
end subroutine format_summary_count
end module mc
