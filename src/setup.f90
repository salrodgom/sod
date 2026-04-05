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
! Setup mode: enumerates symmetry-inequivalent configurations, prepares POSCAR
! El modo setup enumera configuraciones no equivalentes por simetría y prepara los POSCAR
! files per level, runs the external job pipeline, and collects energies.
! por nivel, ejecuta el flujo externo de trabajos y recopila energías.
!*******************************************************************************
! Module `setup` prepares per-level folders, structures, and external workflow inputs.
! El módulo `setup` prepara carpetas por nivel, estructuras y entradas para el flujo externo.
module setup
    use consts
    use cli, only: compose_mode_command, classify_template_file, engine_name_to_filer, &
        is_help_token, parse_level_spec, build_level_sequence, to_lower_inplace
    use clean, only: ensure_gulp_job_sender, prune_level_runtime_files
    use comb, only: run_combination_workflow
    use utils, only: print_sod_banner
    use core, only: binomial_int64, enumerate_unique_subsets, symmetry_finalize, symmetry_get_matrix, &
        symmetry_initialize, write_outsod_file
    use calibration, only: set_calibration_template_gin, set_calibration_protocol_version
    use energy_calculations, only: init_energy_calc, cleanup_energy_calc, write_vasp_file
    use structure_io, only: motif_data_type, motif_atom_type, read_motif_file
    use inputs, only: insod_file_data, read_insod_file
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    character(len=512), save :: scripts_directory = ''
    logical, save :: scripts_ready = .false.
    contains

    ! Orchestrates setup mode: parses CLI options, enumerates levels, and launches outputs per level.
! Orquesta el modo setup: analiza las opciones de la CLI, enumera niveles y lanza las salidas por nivel.
    subroutine run_sod_ensemble_setup(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: level_min, level_max
        logical :: has_list
        integer, allocatable :: level_list(:)
        integer, allocatable :: levels(:)
        character(len=512) :: template_gin_option
        character(len=512) :: template_lib_option
        character(len=512) :: template_lammps_option
        character(len=512) :: template_ase_option
        character(len=64) :: engine_option
        character(len=16) :: protocol_option
        character(len=512) :: motif_option
        integer :: nop, total_sites, active_filer
        integer, pointer :: eqmatrix(:,:)
        integer :: i
        logical :: success
        type(insod_file_data) :: insod

        level_min = 0
        level_max = -1
        has_list = .false.
        template_gin_option = 'none'
        template_lib_option = ''
        template_lammps_option = ''
        template_ase_option = ''
        engine_option = ''
        protocol_option = '2.0'
        motif_option = ''
        allocate(level_list(0))

        call print_sod_banner()
        call parse_setup_arguments(level_min, level_max, level_list, has_list, template_gin_option, protocol_option, &
            motif_option, engine_option, template_lib_option, template_lammps_option, template_ase_option, arg_offset)
        call read_insod_file(insod)

        ! Determine the active FILER: --engine overrides INSOD
        if (len_trim(engine_option) > 0) then
            active_filer = engine_name_to_filer(engine_option)
            if (active_filer < 0) then
                write(error_unit,'(A)') 'Error: unrecognized engine: '//trim(engine_option)
                write(error_unit,'(A)') 'Accepted engines: gulp, lammps, vasp, castep, qe, ase'
                stop 1
            end if
            write(*,'(A,A,A,I0,A)') 'Engine override: --engine ', trim(engine_option), ' (FILER=', active_filer, ')'
        else
            active_filer = insod%filer
        end if

        call init_energy_calc(skip_energy_files=.true.)
        call symmetry_initialize()
        call symmetry_get_matrix(eqmatrix, nop, total_sites)
        if (.not. associated(eqmatrix) .or. total_sites <= 0) then
            write(error_unit,'(A)') 'Error: unable to obtain EQMATRIX or no substitutable sites are available.'
            call symmetry_finalize()
            call cleanup_energy_calc()
            stop 1
        end if

        call build_level_sequence(level_min, level_max, level_list, has_list, total_sites, levels)
        if (.not. allocated(levels) .or. size(levels) == 0) then
            write(error_unit,'(A)') 'Error: the -N specification did not generate valid levels.'
            nullify(eqmatrix)
            call symmetry_finalize()
            call cleanup_energy_calc()
            stop 1
        end if

        if (active_filer == 2) then
            if (trim(template_gin_option) /= 'none') then
                write(*,'(A)') 'Warning: --template-gin is ignored when engine=lammps.'
            end if
            call verify_support_scripts()
            write(*,'(A)') 'Setup mode: preparing and launching the LAMMPS workflow.'
            flush(output_unit)

            nullify(eqmatrix)
            call symmetry_finalize()
            call cleanup_energy_calc()

            success = .true.
            do i = 1, size(levels)
                if (.not. process_level_lammps_setup(levels(i), motif_option, template_lammps_option)) then
                    success = .false.
                    exit
                end if
            end do
            if (allocated(level_list)) deallocate(level_list)
            deallocate(levels)
            if (.not. success) stop 1
            return
        end if

        if (active_filer == 14) then
            if (trim(template_gin_option) /= 'none') then
                write(*,'(A)') 'Warning: --template-gin is ignored when engine=ase.'
            end if
            call verify_support_scripts()
            write(*,'(A)') 'Setup mode: preparing and launching the ASE workflow.'
            flush(output_unit)

            nullify(eqmatrix)
            call symmetry_finalize()
            call cleanup_energy_calc()

            success = .true.
            do i = 1, size(levels)
                if (.not. process_level_ase_setup(levels(i), motif_option, template_ase_option)) then
                    success = .false.
                    exit
                end if
            end do
            if (allocated(level_list)) deallocate(level_list)
            deallocate(levels)
            if (.not. success) stop 1
            return
        end if

        if (active_filer /= 1) then
            if (trim(template_gin_option) /= 'none') then
                write(*,'(A)') 'Warning: --template-gin is ignored when engine is not GULP.'
            end if
            write(*,'(A,I0,A)') 'Setup mode detected FILER=', active_filer, '; using combinatorial file generation only.'
            write(*,'(A)') 'No external energy engine will be launched in this setup run.'
            flush(output_unit)

            nullify(eqmatrix)
            call symmetry_finalize()
            call cleanup_energy_calc()

            success = .true.
            do i = 1, size(levels)
                if (len_trim(motif_option) > 0) then
                    call run_combination_workflow(levels(i), motif_file=trim(motif_option))
                else
                    call run_combination_workflow(levels(i))
                end if
            end do
            if (allocated(level_list)) deallocate(level_list)
            deallocate(levels)
            return
        end if

        call set_calibration_template_gin(trim(template_gin_option))
        call set_calibration_protocol_version(trim(protocol_option))
        call verify_support_scripts()

        write(*,'(A)') '--- Exact configuration setup ---'
        write(*,'(A,I0)') 'Substitutable sites (npos): ', total_sites
        write(*,'(A)', advance='no') 'Levels to process: '
        do i = 1, size(levels)
            write(*,'(I0)', advance='no') levels(i)
            if (i < size(levels)) write(*,'(A)', advance='no') ', '
        end do
        write(*,*)
        flush(output_unit)

        success = .true.
        do i = 1, size(levels)
            if (.not. process_level_setup(levels(i), total_sites, eqmatrix, nop, template_gin_option, protocol_option, &
                motif_option, template_lib_option)) then
                success = .false.
                exit
            end if
        end do

        if (allocated(level_list)) deallocate(level_list)
        deallocate(levels)
        nullify(eqmatrix)
        call symmetry_finalize()
        call cleanup_energy_calc()

        if (.not. success) then
            stop 1
        end if
    end subroutine run_sod_ensemble_setup

    ! Parses CLI arguments specific to setup mode, capturing level selectors.
! Analiza los argumentos de la CLI específicos del modo setup y recoge los selectores de nivel.
    subroutine parse_setup_arguments(level_min, level_max, level_list, has_list, template_gin_option, protocol_option, &
        motif_option, engine_option, template_lib_option, template_lammps_option, template_ase_option, arg_offset)
        integer, intent(inout) :: level_min, level_max
        integer, allocatable, intent(inout) :: level_list(:)
        logical, intent(inout) :: has_list
        character(len=*), intent(out) :: template_gin_option
        character(len=*), intent(out) :: protocol_option
        character(len=*), intent(out) :: motif_option
        character(len=*), intent(out) :: engine_option
        character(len=*), intent(out) :: template_lib_option
        character(len=*), intent(out) :: template_lammps_option
        character(len=*), intent(out) :: template_ase_option
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset, eq_pos
        character(len=256) :: arg, spec
        character(len=256) :: lowered
        character(len=512) :: tpl_path
        character(len=8) :: tpl_cat

        if (allocated(level_list)) deallocate(level_list)
        allocate(level_list(0))
        has_list = .false.
        level_min = 0
        level_max = -1
        template_gin_option = 'none'
        protocol_option = '2.0'
        motif_option = ''
        engine_option = ''
        template_lib_option = ''
        template_lammps_option = ''
        template_ase_option = ''

        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)
        if (argc <= offset) return

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--engine=') == 1) then
                eq_pos = index(arg, '=')
                engine_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            end if

            if (index(trim(lowered), '--template=') == 1) then
                eq_pos = index(arg, '=')
                tpl_path = adjustl(arg(eq_pos+1:))
                tpl_cat = classify_template_file(tpl_path)
                select case (trim(tpl_cat))
                case ('lib')
                    template_lib_option = tpl_path
                case ('gin')
                    template_gin_option = tpl_path
                case ('lammps')
                    template_lammps_option = tpl_path
                case ('py')
                    template_ase_option = tpl_path
                case default
                    write(error_unit,'(A)') 'Error: unrecognized template file extension: '//trim(tpl_path)
                    write(error_unit,'(A)') 'Accepted extensions: .lib, .gin, .include, .lammps, .py'
                    call print_setup_usage()
                    stop 1
                end select
                iarg = iarg + 1
                cycle
            end if

            if (index(trim(lowered), '--motif=') == 1) then
                eq_pos = index(arg, '=')
                motif_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            end if

            if (index(trim(lowered), '--template-gin=') == 1 .or. index(trim(lowered), '--template_gin=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid argument for --template-gin.'
                    call print_setup_usage()
                    stop 1
                end if
                template_gin_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            else if (index(trim(lowered), '--protocol=') == 1 .or. index(trim(lowered), '--protocole=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid argument for --protocol.'
                    call print_setup_usage()
                    stop 1
                end if
                protocol_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-n')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing specification after -N.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, spec)
                call parse_level_spec(spec, level_min, level_max, level_list, has_list)
                iarg = iarg + 2
            case ('--engine')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing engine name after --engine.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, engine_option)
                engine_option = adjustl(engine_option)
                iarg = iarg + 2
            case ('--template')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --template.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, tpl_path)
                tpl_path = adjustl(tpl_path)
                tpl_cat = classify_template_file(tpl_path)
                select case (trim(tpl_cat))
                case ('lib')
                    template_lib_option = tpl_path
                case ('gin')
                    template_gin_option = tpl_path
                case ('lammps')
                    template_lammps_option = tpl_path
                case ('py')
                    template_ase_option = tpl_path
                case default
                    write(error_unit,'(A)') 'Error: unrecognized template file extension: '//trim(tpl_path)
                    write(error_unit,'(A)') 'Accepted extensions: .lib, .gin, .include, .lammps, .py'
                    call print_setup_usage()
                    stop 1
                end select
                iarg = iarg + 2
            case ('--motif')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --motif.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, motif_option)
                motif_option = adjustl(motif_option)
                iarg = iarg + 2
            case ('--template-gin','--template_gin')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --template-gin.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, template_gin_option)
                template_gin_option = adjustl(template_gin_option)
                iarg = iarg + 2
            case ('--protocol','--protocole')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --protocol.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, protocol_option)
                protocol_option = adjustl(protocol_option)
                iarg = iarg + 2
            case ('--no-template-gin','--skip-template','--skip_template')
                template_gin_option = 'none'
                iarg = iarg + 1
            case ('--help','-h','--ayuda','-ayuda','ayuda')
                call print_setup_usage()
                stop
            case default
                write(error_unit,'(A)') 'Error: unrecognized argument in setup mode. Use --help for more information.'
                call print_setup_usage()
                stop 1
            end select
        end do
    end subroutine parse_setup_arguments

    ! Displays usage instructions for setup mode.
! Muestra las instrucciones de uso del modo setup.
    subroutine print_setup_usage()
        character(len=256) :: command_name

        command_name = compose_mode_command('setup')
        write(*,'(A)') 'Usage: '//trim(command_name)//' -N <spec> [--engine <name>] [--template <file>...] [--motif <file>]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: setup, prepare, inputs, energy'
        write(*,'(A)') ''
        write(*,'(A)') '  -N, -n    Defines the substitution levels (Y atoms) to prepare. Accepts:'
        write(*,'(A)') '            * a single value (e.g. 4)'
        write(*,'(A)') '            * a colon-separated range (e.g. 2:6)'
        write(*,'(A)') '            * a comma-separated list (e.g. 0,3,5,7)'
        write(*,'(A)') ''
        write(*,'(A)') '  --engine <name>        Select the energy engine (overrides INSOD FILER).'
        write(*,'(A)') '                         gulp, lammps, vasp, castep, qe, ase'
        write(*,'(A)') ''
        write(*,'(A)') '  --template <file>      Provide a template/model file for the engine (repeatable).'
        write(*,'(A)') '                         Auto-classified by extension:'
        write(*,'(A)') '                           .lib      GULP potential library'
        write(*,'(A)') '                           .gin      GULP optimization fragment'
        write(*,'(A)') '                           .include  GULP optimization fragment'
        write(*,'(A)') '                           .lammps   LAMMPS input template'
        write(*,'(A)') '                           .py       ASE calculator script'
        write(*,'(A)') ''
        write(*,'(A)') '  --template-gin <file>  Alias for --template with a .gin/.include file.'
        write(*,'(A)') '                         Use "default" to copy the bundled template.'
        write(*,'(A)') '  --no-template-gin      Skip template fragments when creating .gin files.'
        write(*,'(A)') '  --motif <file>         Include a molecular motif (.include) in every configuration.'
        write(*,'(A)') '  --protocol <ver>       Select the GULP workflow protocol: 1.0 or 2.0 [2.0].'
        write(*,'(A)') ''
        write(*,'(A)') 'Examples:'
        write(*,'(A)') '  sod setup -N 4 --engine gulp --template custom.lib --template optim.gin --motif OSDA.include'
        write(*,'(A)') '  sod setup -N 4 --engine lammps --template template_in.lammps'
        write(*,'(A)') '  sod setup -N 4 --engine ase --template template_ase_mace.py'
        write(*,'(A)') ''
        write(*,'(A)') 'If --engine is not given, the engine is selected from the INSOD FILER value.'
        write(*,'(A)') 'If --template is not given, defaults are used (bundled .lib for GULP,'
        write(*,'(A)') 'template_in.lammps from CWD for LAMMPS, template_ase.py from CWD for ASE).'
    end subroutine print_setup_usage

    ! Processes a single level: enumerates unique subsets, writes artifacts, and runs external scripts.
! Procesa un único nivel: enumera subconjuntos únicos, escribe artefactos y ejecuta scripts externos.
    logical function process_level_setup(level, total_sites, eqmatrix, nop, template_gin_option, protocol_option, &
        motif_option, template_lib_option) result(success)
        integer, intent(in) :: level, total_sites, nop
        integer, intent(in) :: eqmatrix(:,:)
        character(len=*), intent(in) :: template_gin_option
        character(len=*), intent(in) :: protocol_option
        character(len=*), intent(in) :: motif_option
        character(len=*), intent(in) :: template_lib_option
        type(motif_data_type) :: motif
        logical :: have_motif
        integer :: unique_count
        integer, allocatable :: unique_subsets(:,:)
        integer, allocatable :: unique_deg(:)
        integer :: idx
        integer(ip) :: total_comb
        character(len=16) :: level_dir
        logical :: ok
        character(len=512) :: command
        integer, allocatable :: config(:)

        success = .false.

        total_comb = binomial_int64(total_sites, level)
        write(*,'(A,I0)') 'Level: ', level
        write(*,'(A,I0)') 'Total combinations: ', total_comb
        flush(output_unit)

        level_dir = format_level_directory(level)
        command = 'mkdir -p ' // trim(level_dir)
        if (.not. run_shell_command(command, 'create directory ' // trim(level_dir))) return

        command = 'bash -c "rm -f ' // trim(level_dir)//'/c*.vasp ' // trim(level_dir)//'/*.vasp.gin ' // trim(level_dir)//'/*.vasp.gout ' // trim(level_dir)//'/*.vasp.grs ' // trim(level_dir)//'/OUTSOD ' // trim(level_dir)//'/ENERGIES"'
        if (.not. run_shell_command(command, 'clean directory ' // trim(level_dir))) return

        call enumerate_unique_subsets(total_sites, level, eqmatrix, nop, unique_subsets, unique_deg, unique_count)

        if (level > 0 .and. unique_count <= 0) then
            write(error_unit,'(A)') 'Error: no unique configurations were found.'
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            if (allocated(unique_deg)) deallocate(unique_deg)
            return
        end if

        write(*,'(A,I0)') 'Unique configurations: ', max(1, unique_count)
        flush(output_unit)

        ok = write_outsod_file(trim(level_dir)//'/OUTSOD', level, total_sites, unique_count, unique_deg, unique_subsets)
        if (.not. ok) then
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        allocate(config(total_sites))
        config = 1

        have_motif = (len_trim(motif_option) > 0)
        if (have_motif) call read_motif_file(trim(motif_option), motif)

        if (level > 0) then
            do idx = 1, unique_count
                config = 1
                config(unique_subsets(1:level, idx)) = 2
                if (have_motif) then
                    if (.not. write_poscar(level_dir, idx, config, total_sites, motif%atoms, motif%natoms)) then
                        deallocate(config)
                        if (allocated(unique_subsets)) deallocate(unique_subsets)
                        deallocate(unique_deg)
                        return
                    end if
                else
                    if (.not. write_poscar(level_dir, idx, config, total_sites)) then
                        deallocate(config)
                        if (allocated(unique_subsets)) deallocate(unique_subsets)
                        deallocate(unique_deg)
                        return
                    end if
                end if
            end do
        else
            if (have_motif) then
                if (.not. write_poscar(level_dir, 1, config, total_sites, motif%atoms, motif%natoms)) then
                    deallocate(config)
                    deallocate(unique_deg)
                    return
                end if
            else
                if (.not. write_poscar(level_dir, 1, config, total_sites)) then
                    deallocate(config)
                    deallocate(unique_deg)
                    return
                end if
            end if
        end if
        deallocate(config)

        if (.not. copy_support_scripts(level_dir, template_gin_option, template_lib_option)) then
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        command = 'bash -c "cd ' // trim(level_dir) // ' && SOD_GULP_PROTOCOL_VERSION=' // trim(protocol_option) // ' ./run_jobs.sh"'
        if (.not. run_shell_command(command, 'run_jobs.sh in '//trim(level_dir))) then
            call remove_protocol_template_file(level_dir)
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        command = 'bash -c "cd ' // trim(level_dir) // ' && ./extract.sh"'
        if (.not. run_shell_command(command, 'extract.sh in '//trim(level_dir))) then
            call remove_protocol_template_file(level_dir)
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if
        call remove_protocol_template_file(level_dir)

        if (.not. energies_generated(level_dir)) then
            write(error_unit,'(A)') 'Error: the ENERGIES file was not found after extract.sh finished.'
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        call ensure_gulp_job_sender(level_dir)
        call prune_level_runtime_files(level_dir, 1)

        if (allocated(unique_subsets)) deallocate(unique_subsets)
        deallocate(unique_deg)

        success = .true.
    end function process_level_setup

    logical function process_level_lammps_setup(level, motif_option, template_lammps_option) result(ok)
        integer, intent(in) :: level
        character(len=*), intent(in) :: motif_option
        character(len=*), intent(in) :: template_lammps_option
        character(len=16) :: level_dir
        character(len=512) :: command

        level_dir = format_level_directory(level)
        ok = .false.

        if (len_trim(motif_option) > 0) then
            call run_combination_workflow(level, motif_file=trim(motif_option))
        else
            call run_combination_workflow(level)
        end if
        if (.not. copy_lammps_support_files(level_dir, template_lammps_option)) return

        write(*,'(A)') 'Running LAMMPS jobs...'
        flush(output_unit)
        command = 'bash -c "cd ' // trim(level_dir) // ' && ./run_jobs.sh"'
        if (.not. run_shell_command(command, 'run_jobs.sh in '//trim(level_dir))) then
            return
        end if

        write(*,'(A)') 'Extracting energies and relaxed structures...'
        flush(output_unit)
        command = 'bash -c "cd ' // trim(level_dir) // ' && ./extract.sh"'
        if (.not. run_shell_command(command, 'extract.sh in '//trim(level_dir))) then
            return
        end if

        if (.not. energies_generated(level_dir)) then
            write(error_unit,'(A)') 'Error: the ENERGIES file was not found after extract.sh finished.'
            return
        end if

        call prune_level_runtime_files(level_dir, 2)

        ok = .true.
    end function process_level_lammps_setup

    logical function process_level_ase_setup(level, motif_option, template_ase_option) result(ok)
        integer, intent(in) :: level
        character(len=*), intent(in) :: motif_option
        character(len=*), intent(in) :: template_ase_option
        character(len=16) :: level_dir
        character(len=512) :: command

        level_dir = format_level_directory(level)
        ok = .false.

        if (len_trim(motif_option) > 0) then
            call run_combination_workflow(level, motif_file=trim(motif_option))
        else
            call run_combination_workflow(level)
        end if
        if (.not. copy_ase_support_files(level_dir, template_ase_option)) return

        write(*,'(A)') 'Running ASE jobs...'
        flush(output_unit)
        command = 'bash -c "cd ' // trim(level_dir) // ' && ./run_jobs.sh"'
        if (.not. run_shell_command(command, 'run_jobs.sh in '//trim(level_dir))) then
            return
        end if

        write(*,'(A)') 'Extracting energies and relaxed structures...'
        flush(output_unit)
        command = 'bash -c "cd ' // trim(level_dir) // ' && ./extract.sh"'
        if (.not. run_shell_command(command, 'extract.sh in '//trim(level_dir))) then
            return
        end if

        if (.not. energies_generated(level_dir)) then
            write(error_unit,'(A)') 'Error: the ENERGIES file was not found after extract.sh finished.'
            return
        end if

        call prune_level_runtime_files(level_dir, 14)

        ok = .true.
    end function process_level_ase_setup

    ! Serializes a configuration into a POSCAR file inside the level directory.
! Serializa una configuración en un archivo POSCAR dentro del directorio del nivel.
    logical function write_poscar(dir, index, config, total_sites, motif_atoms, motif_count) result(ok)
        character(len=*), intent(in) :: dir
        integer, intent(in) :: index, total_sites
        integer, intent(in) :: config(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        character(len=256) :: filename
        character(len=32) :: index_str
        integer :: pad_len

        ok = .false.
        index_str = ''
        write(index_str,'(I0)') index
        pad_len = max(0, 5 - len_trim(index_str))
        filename = trim(dir)//'/c'//repeat('0', pad_len)//trim(index_str)//'.vasp'
        if (present(motif_atoms) .and. present(motif_count)) then
            call write_vasp_file(config, total_sites, trim(filename), motif_atoms=motif_atoms, motif_count=motif_count)
        else
            call write_vasp_file(config, total_sites, trim(filename))
        end if
        ok = .true.
    end function write_poscar

    ! Copies helper shell scripts into the working level directory and ensures they are executable.
! Copia los scripts auxiliares de shell al directorio de trabajo del nivel y asegura que sean ejecutables.
    logical function copy_support_scripts(dir, template_gin_option, template_lib_option) result(ok)
        character(len=*), intent(in) :: dir
        character(len=*), intent(in) :: template_gin_option
        character(len=*), intent(in) :: template_lib_option
        character(len=512) :: command
        character(len=512) :: lowered_template
        logical :: exists

        if (.not. scripts_ready) call verify_support_scripts()

        ! Copy standard scripts
        command = 'cp "' // trim(scripts_directory)//'/run_jobs.sh" "' // trim(scripts_directory)//'/extract.sh" "' // &
            trim(scripts_directory)//'/vasp2gin.sh" "' // trim(dir) // '"'
        if (.not. run_shell_command(command, 'copy scripts to '//trim(dir))) then
            ok = .false.
            return
        end if

        ! Copy potential library: custom --template .lib or bundled default
        if (len_trim(template_lib_option) > 0) then
            command = 'cp "' // trim(template_lib_option) // '" "' // trim(dir)//'/framework_template.lib"'
            if (.not. run_shell_command(command, 'copy custom .lib to '//trim(dir))) then
                ok = .false.
                return
            end if
        else
            command = 'cp "' // trim(scripts_directory)//'/framework_template.lib" "' // trim(dir) // '"'
            if (.not. run_shell_command(command, 'copy default .lib to '//trim(dir))) then
                ok = .false.
                return
            end if
        end if

        inquire(file=trim(scripts_directory)//'/reference/protocol_template.gin', exist=exists)
        if (exists) then
            command = 'cp "' // trim(scripts_directory)//'/reference/protocol_template.gin" "' // trim(dir)//'/protocol_template.gin"'
            if (.not. run_shell_command(command, 'copy protocol template to '//trim(dir))) then
                ok = .false.
                return
            end if
        else
            write(*,'(A)') 'Warning: scripts/reference/protocol_template.gin was not found; protocol mode may be unavailable.'
        end if

        lowered_template = adjustl(template_gin_option)
        call to_lower_inplace(lowered_template)

        select case (trim(lowered_template))
        case ('', 'none', 'off', 'skip')
            continue
        case ('default', 'builtin')
            inquire(file=trim(scripts_directory)//'/default_template.gin', exist=exists)
            if (exists) then
                command = 'cp "' // trim(scripts_directory)//'/default_template.gin" "' // trim(dir)//'/template_payload.gin"'
                if (.not. run_shell_command(command, 'copy default template payload to '//trim(dir))) then
                    ok = .false.
                    return
                end if
            else
                inquire(file=trim(scripts_directory)//'/default_template.include', exist=exists)
                if (exists) then
                    command = 'cp "' // trim(scripts_directory)//'/default_template.include" "' // trim(dir)//'/default_template.include"'
                    if (.not. run_shell_command(command, 'copy default template include to '//trim(dir))) then
                        ok = .false.
                        return
                    end if
                else
                    write(*,'(A)') 'Warning: scripts/default_template.include was not found; template fragment will be omitted.'
                end if
            end if
        case default
            command = 'cp "' // trim(template_gin_option) // '" "' // trim(dir)//'/template_payload.gin"'
            if (.not. run_shell_command(command, 'copy custom template payload to '//trim(dir))) then
                ok = .false.
                return
            end if
        end select

        command = 'chmod +x "' // trim(dir)//'/run_jobs.sh" "' // trim(dir)//'/extract.sh" "' // trim(dir)//'/vasp2gin.sh"'
        if (.not. run_shell_command(command, 'grant permissions in '//trim(dir))) then
            ok = .false.
            return
        end if

        ok = .true.
    end function copy_support_scripts

    logical function copy_lammps_support_files(dir, template_lammps_option) result(ok)
        character(len=*), intent(in) :: dir
        character(len=*), intent(in) :: template_lammps_option
        character(len=512) :: command, lammps_source
        logical :: exists

        if (.not. scripts_ready) call verify_support_scripts()

        ! Resolve LAMMPS template: --template .lammps, or CWD fallback
        if (len_trim(template_lammps_option) > 0) then
            lammps_source = template_lammps_option
        else
            lammps_source = 'template_in.lammps'
        end if

        inquire(file=trim(lammps_source), exist=exists)
        if (.not. exists) then
            write(error_unit,'(A)') 'Error: LAMMPS template not found: '//trim(lammps_source)
            write(error_unit,'(A)') 'Provide one via --template <file>.lammps or place template_in.lammps in the working directory.'
            ok = .false.
            return
        end if

        command = 'cp "' // trim(scripts_directory)//'/run_jobs.sh" "' // trim(scripts_directory)//'/extract.sh" "' // trim(dir) // '"'
        if (.not. run_shell_command(command, 'copy LAMMPS scripts to '//trim(dir))) then
            ok = .false.
            return
        end if

        command = 'cp "' // trim(lammps_source) // '" "' // trim(dir)//'/template_in.lammps"'
        if (.not. run_shell_command(command, 'copy LAMMPS template to '//trim(dir))) then
            ok = .false.
            return
        end if

        command = 'chmod +x "' // trim(dir)//'/run_jobs.sh" "' // trim(dir)//'/extract.sh"'
        if (.not. run_shell_command(command, 'grant LAMMPS permissions in '//trim(dir))) then
            ok = .false.
            return
        end if

        ok = .true.
    end function copy_lammps_support_files

    logical function copy_ase_support_files(dir, template_ase_option) result(ok)
        character(len=*), intent(in) :: dir
        character(len=*), intent(in) :: template_ase_option
        character(len=512) :: command, ase_source
        logical :: exists

        if (.not. scripts_ready) call verify_support_scripts()

        ! Resolve ASE template: --template .py, or CWD fallback
        if (len_trim(template_ase_option) > 0) then
            ase_source = template_ase_option
        else
            ase_source = 'template_ase.py'
        end if

        inquire(file=trim(ase_source), exist=exists)
        if (.not. exists) then
            write(error_unit,'(A)') 'Error: ASE template not found: '//trim(ase_source)
            write(error_unit,'(A)') 'Provide one via --template <file>.py or place template_ase.py in the working directory.'
            ok = .false.
            return
        end if

        command = 'cp "' // trim(scripts_directory)//'/run_jobs.sh" "' // trim(scripts_directory)//'/extract.sh" "' // &
            trim(scripts_directory)//'/vasp2ase.py" "' // trim(dir) // '"'
        if (.not. run_shell_command(command, 'copy ASE scripts to '//trim(dir))) then
            ok = .false.
            return
        end if

        command = 'cp "' // trim(ase_source) // '" "' // trim(dir)//'/template_ase.py"'
        if (.not. run_shell_command(command, 'copy ASE template to '//trim(dir))) then
            ok = .false.
            return
        end if

        command = 'chmod +x "' // trim(dir)//'/run_jobs.sh" "' // trim(dir)//'/extract.sh" "' // trim(dir)//'/vasp2ase.py"'
        if (.not. run_shell_command(command, 'grant ASE permissions in '//trim(dir))) then
            ok = .false.
            return
        end if

        ok = .true.
    end function copy_ase_support_files

    subroutine remove_protocol_template_file(dir)
        character(len=*), intent(in) :: dir
        character(len=512) :: path
        integer :: unit_file, ios
        logical :: exists

        path = trim(dir)//'/protocol_template.gin'
        inquire(file=trim(path), exist=exists)
        if (.not. exists) return

        open(newunit=unit_file, file=trim(path), status='old', action='readwrite', iostat=ios)
        if (ios /= 0) return
        close(unit_file, status='delete', iostat=ios)
    end subroutine remove_protocol_template_file

    ! Executes a shell command and reports contextualized failures.
! Ejecuta un comando de shell e informa de fallos con contexto.
    logical function run_shell_command(command, context) result(ok)
        character(len=*), intent(in) :: command
        character(len=*), intent(in) :: context
        integer :: cmdstat, exitstat

        call execute_command_line(trim(command), wait=.true., cmdstat=cmdstat, exitstat=exitstat)
        if (cmdstat /= 0) then
            write(error_unit,'(A,": failed to launch command (cmdstat=",I0,")")') trim(context), cmdstat
            ok = .false.
            return
        end if
        if (exitstat /= 0) then
            write(error_unit,'(A,": the command returned exit code ",I0)') trim(context), exitstat
            ok = .false.
            return
        end if
        ok = .true.
    end function run_shell_command

    ! Checks whether extract.sh produced the ENERGIES file in the level directory.
! Comprueba si extract.sh produjo el archivo ENERGIES en el directorio del nivel.
    logical function energies_generated(dir) result(ok)
        character(len=*), intent(in) :: dir
        logical :: exists
        inquire(file=trim(dir)//'/ENERGIES', exist=exists)
        ok = exists
    end function energies_generated

    ! Locates the scripts directory and caches it for later copy operations.
! Localiza el directorio de scripts y lo guarda en caché para operaciones de copia posteriores.
    subroutine verify_support_scripts()
        logical :: ok

        if (scripts_ready) return

        ok = locate_scripts_directory(scripts_directory)
        if (.not. ok) then
            write(error_unit,'(A)') 'Error: scripts directory with run_jobs.sh, extract.sh, vasp2gin.sh, and vasp2ase.py was not found.'
            stop 1
        end if

        scripts_ready = .true.
    end subroutine verify_support_scripts

    ! Searches parent directories for the scripts folder containing required helpers.
! Busca en los directorios padre la carpeta scripts que contiene los auxiliares requeridos.
    logical function locate_scripts_directory(out_dir) result(found)
        character(len=*), intent(out) :: out_dir
        integer :: depth, step
        character(len=512) :: candidate
        logical :: exists_run, exists_extract, exists_vasp2gin, exists_vasp2ase

        out_dir = ''
        found = .false.

        do depth = 0, 6
            if (depth == 0) then
                candidate = 'scripts'
            else
                candidate = ''
                do step = 1, depth
                    candidate = trim(candidate)//'../'
                end do
                candidate = trim(candidate)//'scripts'
            end if

            inquire(file=trim(candidate)//'/run_jobs.sh', exist=exists_run)
            if (.not. exists_run) cycle
            inquire(file=trim(candidate)//'/extract.sh', exist=exists_extract)
            if (.not. exists_extract) cycle
            inquire(file=trim(candidate)//'/vasp2gin.sh', exist=exists_vasp2gin)
            if (.not. exists_vasp2gin) cycle
            inquire(file=trim(candidate)//'/vasp2ase.py', exist=exists_vasp2ase)
            if (.not. exists_vasp2ase) cycle

            out_dir = trim(candidate)
            found = .true.
            exit
        end do
    end function locate_scripts_directory

    ! Formats the directory name used for a given substitution level (e.g., n03).
! Da formato al nombre del directorio usado para un nivel de sustitución dado (por ejemplo, n03).
    pure function format_level_directory(level) result(dir)
        integer, intent(in) :: level
        character(len=16) :: dir
        write(dir,'("n0",I0)') level
    end function format_level_directory

end module setup
