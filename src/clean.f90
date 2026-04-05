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
! Cleanup mode: prunes transient helper files from completed per-level folders.
! El modo clean elimina ficheros auxiliares transitorios de carpetas por nivel ya completadas.
!*******************************************************************************
! Module `clean` centralizes conservative runtime cleanup helpers and a dedicated
! command-line mode that only prunes finished setup directories.
! El módulo `clean` centraliza utilidades conservadoras de limpieza en tiempo de
! ejecución y un modo dedicado de línea de comandos que solo poda directorios de
! setup ya terminados.
module clean
    use cli, only: compose_mode_command, is_help_token, parse_level_spec, build_level_sequence, to_lower_inplace
    use inputs, only: insod_file_data, read_insod_file
    use utils, only: format_level_directory, join_path
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    private

    public :: run_sod_clean
    public :: prune_level_runtime_files
    public :: ensure_gulp_job_sender

contains

    ! Orchestrates the cleanup mode, pruning only finished level directories.
    ! Orquesta el modo clean y poda solo los directorios de nivel terminados.
    subroutine run_sod_clean(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: offset, level_min, level_max, total_sites
        integer :: filer, idx, cleaned_count, skipped_count, missing_count
        integer, allocatable :: level_list(:)
        integer, allocatable :: levels(:)
        logical :: has_list, level_specified
        character(len=32) :: level_dir
        type(insod_file_data) :: insod

        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)

        level_min = 0
        level_max = -1
        has_list = .false.
        level_specified = .false.
        allocate(level_list(0))

        call parse_clean_arguments(level_min, level_max, level_list, has_list, level_specified, offset)
        call read_insod_file(insod)

        if (.not. level_specified) then
            level_min = insod%nsubs_min
            level_max = insod%nsubs_max
        end if

        total_sites = infer_total_sites(insod)
        call build_level_sequence(level_min, level_max, level_list, has_list, total_sites, levels)
        if (.not. allocated(levels) .or. size(levels) == 0) then
            write(error_unit,'(A)') 'Error: clean mode did not resolve any valid levels.'
            if (allocated(level_list)) deallocate(level_list)
            stop 1
        end if

        cleaned_count = 0
        skipped_count = 0
        missing_count = 0

        write(*,'(A)') '------------------------------------------------------------------------'
        write(*,'(A)') 'Cleaning finished setup directories (nNN only)...'

        do idx = 1, size(levels)
            level_dir = format_level_directory('n', levels(idx))
            if (.not. file_exists(join_path(level_dir, 'OUTSOD'))) then
                missing_count = missing_count + 1
                write(*,'(A,1X,A)') 'Skipping missing level directory:', trim(level_dir)
                cycle
            end if

            filer = detect_level_filer(level_dir, insod%filer)
            if (.not. level_cleanup_eligible(level_dir, filer)) then
                skipped_count = skipped_count + 1
                write(*,'(A,1X,A)') 'Skipping unfinished or unsupported level directory:', trim(level_dir)
                cycle
            end if

            if (filer == 1) call ensure_gulp_job_sender(level_dir)
            call prune_level_runtime_files(level_dir, filer)
            cleaned_count = cleaned_count + 1
            write(*,'(A,1X,A,1X,A,I0)') 'Cleaned finished level directory:', trim(level_dir), 'FILER=', filer
        end do

        write(*,'(A)') '------------------------------------------------------------------------'
        write(*,'(A,I0)') 'Finished levels cleaned: ', cleaned_count
        write(*,'(A,I0)') 'Skipped levels: ', skipped_count
        write(*,'(A,I0)') 'Missing levels: ', missing_count

        if (allocated(level_list)) deallocate(level_list)
        if (allocated(levels)) deallocate(levels)
    end subroutine run_sod_clean

    ! Parses CLI arguments accepted by clean mode.
    ! Analiza los argumentos de la CLI aceptados por el modo clean.
    subroutine parse_clean_arguments(level_min, level_max, level_list, has_list, level_specified, arg_offset)
        integer, intent(inout) :: level_min, level_max
        integer, allocatable, intent(inout) :: level_list(:)
        logical, intent(inout) :: has_list
        logical, intent(out) :: level_specified
        integer, intent(in) :: arg_offset
        integer :: argc, iarg, eq_pos
        character(len=256) :: arg, lowered, spec

        level_specified = .false.
        argc = command_argument_count()
        iarg = 1 + arg_offset

        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '-n=') == 1 .or. index(trim(lowered), '--level=') == 1 .or. &
                index(trim(lowered), '--levels=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid level specification in clean mode.'
                    call print_clean_usage()
                    stop 1
                end if
                spec = adjustl(arg(eq_pos+1:))
                call parse_level_spec(trim(spec), level_min, level_max, level_list, has_list)
                level_specified = .true.
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-n','--level','--levels')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after -N/--level in clean mode.'
                    call print_clean_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, spec)
                call parse_level_spec(trim(spec), level_min, level_max, level_list, has_list)
                level_specified = .true.
                iarg = iarg + 2
            case ('--help','-h','--ayuda','-ayuda','ayuda')
                call print_clean_usage()
                stop
            case default
                write(error_unit,'(A,1X,A)') 'Error: unrecognized argument in clean mode:', trim(arg)
                call print_clean_usage()
                stop 1
            end select
        end do
    end subroutine parse_clean_arguments

    ! Prints the dedicated help text for clean mode.
    ! Imprime la ayuda específica del modo clean.
    subroutine print_clean_usage()
        character(len=256) :: command_name

        command_name = compose_mode_command('clean')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [-N <level-spec>]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: clean, prune, tidy'
        write(*,'(A)') ''
        write(*,'(A)') 'Clean mode prunes transient helper files from finished nNN directories only.'
        write(*,'(A)') 'It never deletes finished outputs such as OUTSOD, ENERGIES, relaxed structures,'
        write(*,'(A)') 'or the minimal files required to relaunch a completed calculation.'
        write(*,'(A)') ''
        write(*,'(A)') '  -N, --level, --levels <spec>  Optional level, list, or range to clean.'
        write(*,'(A)') '                                If omitted, INSOD nsubs_min:nsubs_max is used.'
        write(*,'(A)') ''
        write(*,'(A)') 'Examples:'
        write(*,'(A)') '  sod clean'
        write(*,'(A)') '  sod clean -N 4'
        write(*,'(A)') '  sod clean -N 3,5,7'
        write(*,'(A)') '  sod clean -N 2:6'
    end subroutine print_clean_usage

    ! Ensures a rerunnable GULP job_sender exists for a finished level directory.
    ! Garantiza que exista un job_sender relanzable de GULP para un directorio de nivel terminado.
    subroutine ensure_gulp_job_sender(level_dir)
        character(len=*), intent(in) :: level_dir
        logical :: has_job_sender
        integer :: unit_job, idx, unique_count
        character(len=32) :: basename

        has_job_sender = file_exists(join_path(level_dir, 'job_sender'))
        if (has_job_sender) return

        unique_count = read_outsod_unique_count(join_path(level_dir, 'OUTSOD'))
        if (unique_count <= 0) then
            write(error_unit,'(A,1X,A)') 'Warning: unable to reconstruct job_sender for', trim(level_dir)
            return
        end if

        open(newunit=unit_job, file=trim(join_path(level_dir, 'job_sender')), status='replace', action='write')
        do idx = 1, unique_count
            write(basename,'("c",I5.5)') idx
            write(unit_job,'(A)') 'gulp < ' // trim(basename)//'.vasp.gin > ' // trim(basename)//'.vasp.gout'
        end do
        close(unit_job)
    end subroutine ensure_gulp_job_sender

    ! Removes helper scripts and transient support files once the level finished successfully.
    ! Elimina scripts auxiliares y ficheros transitorios de soporte cuando el nivel terminó con éxito.
    subroutine prune_level_runtime_files(level_dir, filer)
        character(len=*), intent(in) :: level_dir
        integer, intent(in) :: filer
        character(len=2048) :: command

        select case (filer)
        case (1)
            command = 'bash -c "rm -f ' // trim(level_dir)//'/run_jobs.sh ' // trim(level_dir)//'/extract.sh ' // &
                trim(level_dir)//'/vasp2gin.sh ' // trim(level_dir)//'/framework_template.lib ' // &
                trim(level_dir)//'/protocol_template.gin ' // trim(level_dir)//'/template_payload.gin ' // &
                trim(level_dir)//'/default_template.gin ' // trim(level_dir)//'/default_template.include ; ' // &
                'rm -rf ' // trim(level_dir)//'/slot_* ' // trim(level_dir)//'/__pycache__"'
        case (2)
            command = 'bash -c "rm -f ' // trim(level_dir)//'/run_jobs.sh ' // trim(level_dir)//'/extract.sh ' // &
                trim(level_dir)//'/template_in.lammps ; rm -rf ' // trim(level_dir)//'/__pycache__"'
        case (14)
            command = 'bash -c "rm -f ' // trim(level_dir)//'/run_jobs.sh ' // trim(level_dir)//'/extract.sh ' // &
                trim(level_dir)//'/*.relaxed.vasp ; ' // &
                'rm -rf ' // trim(level_dir)//'/__pycache__"'
        case default
            return
        end select

        if (.not. run_shell_command(command, 'prune runtime files in '//trim(level_dir))) then
            write(error_unit,'(A,1X,A)') 'Warning: failed to prune runtime helpers in', trim(level_dir)
        end if
    end subroutine prune_level_runtime_files

    ! Returns the number of substitutable sites inferred directly from INSOD.
    ! Devuelve el número de sitios sustituibles inferido directamente desde INSOD.
    integer function infer_total_sites(insod) result(total_sites)
        type(insod_file_data), intent(in) :: insod

        total_sites = 0
        if (allocated(insod%natsp0)) then
            if (insod%sptarget >= 1 .and. insod%sptarget <= size(insod%natsp0)) then
                total_sites = insod%natsp0(insod%sptarget) * max(1, insod%na) * max(1, insod%nb) * max(1, insod%nc)
            end if
        end if
        if (total_sites <= 0) total_sites = max(insod%nsubs_max, insod%nsubs)
    end function infer_total_sites

    ! Detects the engine FILER associated with an existing level directory.
    ! Detecta el engine FILER asociado a un directorio de nivel existente.
    integer function detect_level_filer(level_dir, default_filer) result(filer)
        character(len=*), intent(in) :: level_dir
        integer, intent(in) :: default_filer
        integer :: unit_filer, ios
        logical :: exists

        filer = default_filer

        inquire(file=trim(join_path(level_dir, 'filer')), exist=exists)
        if (exists) then
            open(newunit=unit_filer, file=trim(join_path(level_dir, 'filer')), status='old', action='read', iostat=ios)
            if (ios == 0) then
                read(unit_filer, *, iostat=ios) filer
                close(unit_filer)
                if (ios == 0) return
            end if
        end if

        if (file_exists(join_path(level_dir, 'vasp2ase.py')) .or. file_exists(join_path(level_dir, 'template_ase.py')) .or. &
            file_exists(join_path(level_dir, 'c00001.out.ase'))) then
            filer = 14
        else if (file_exists(join_path(level_dir, 'template_in.lammps')) .or. file_exists(join_path(level_dir, 'c00001.in.lammps')) .or. &
            file_exists(join_path(level_dir, 'c00001.out.lammps'))) then
            filer = 2
        else if (file_exists(join_path(level_dir, 'c00001.vasp.gin')) .or. file_exists(join_path(level_dir, 'c00001.vasp.gout'))) then
            filer = 1
        end if
    end function detect_level_filer

    ! Returns whether a level directory is finished and safe to prune.
    ! Devuelve si un directorio de nivel está terminado y es seguro podarlo.
    logical function level_cleanup_eligible(level_dir, filer)
        character(len=*), intent(in) :: level_dir
        integer, intent(in) :: filer

        if (.not. file_exists(join_path(level_dir, 'OUTSOD'))) then
            level_cleanup_eligible = .false.
            return
        end if

        select case (filer)
        case (1, 2, 14)
            level_cleanup_eligible = file_exists(join_path(level_dir, 'ENERGIES'))
        case default
            level_cleanup_eligible = .false.
        end select
    end function level_cleanup_eligible

    ! Reads the unique configuration count stored on the second line of OUTSOD.
    ! Lee el número de configuraciones únicas almacenado en la segunda línea de OUTSOD.
    integer function read_outsod_unique_count(filename) result(unique_count)
        character(len=*), intent(in) :: filename
        integer :: unit_outsod, ios
        character(len=512) :: line

        unique_count = 0
        if (.not. file_exists(filename)) return

        open(newunit=unit_outsod, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        read(unit_outsod, '(A)', iostat=ios) line
        if (ios /= 0) then
            close(unit_outsod)
            return
        end if
        read(unit_outsod, '(A)', iostat=ios) line
        if (ios /= 0) then
            close(unit_outsod)
            return
        end if
        close(unit_outsod)

        read(line, *, iostat=ios) unique_count
        if (ios /= 0) unique_count = 0
    end function read_outsod_unique_count

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

    ! Checks whether a file exists on disk.
    ! Comprueba si un fichero existe en disco.
    logical function file_exists(path)
        character(len=*), intent(in) :: path
        inquire(file=trim(path), exist=file_exists)
    end function file_exists

end module clean
