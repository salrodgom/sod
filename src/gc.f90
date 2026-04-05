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
! Grand-canonical post-processing over a range of substitution levels.
! This mode reuses OUTSOD and ENERGIES files generated for several levels and
! computes grand-canonical probabilities, cumulative populations by level, and
! finite-temperature thermodynamics.
!
! Postproceso grand-canónico sobre un rango de niveles de sustitución.
! Este modo reutiliza archivos OUTSOD y ENERGIES generados para varios niveles
! y calcula probabilidades grand-canónicas, poblaciones acumuladas por nivel y
! termodinámica a temperatura finita.
!*******************************************************************************
module gc
    use consts
    use cli, only: compose_mode_command, is_help_token, to_lower_inplace
    use utils, only: format_level_directory, join_path, print_sod_banner
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    private

    integer, parameter :: max_temps_default = 2
    integer, parameter :: max_tokens_per_line = 16
    real(dp), parameter :: tolmu = 1.0d-10
    real(dp), parameter :: tolq = 1.0d-12
    real(dp), parameter :: eva3togpa = 160.2176621_dp
    real(dp), parameter :: tiny_q = 1.0d-16

    type :: ingc_data_t
        integer :: nsubsmin = 0
        integer :: nsubsmax = 0
        character(len=2) :: xormu = 'x '
        real(dp) :: xormuvalue = 0.0_dp
        real(dp) :: lambda = 0.0_dp
        real(dp) :: v0 = 0.0_dp
        real(dp) :: v1 = 0.0_dp
        real(dp) :: bv = 0.0_dp
        real(dp) :: bm0 = 0.0_dp
        real(dp) :: bm1 = 0.0_dp
        real(dp) :: bb = 0.0_dp
    end type ingc_data_t

    public :: run_sod_gc

contains

    ! Runs the grand-canonical post-processing workflow from INGC and level data.
    ! Ejecuta el flujo de postproceso grand-canónico a partir de INGC y de los datos por nivel.
    subroutine run_sod_gc(arg_offset)
        integer, intent(in), optional :: arg_offset
        type(ingc_data_t) :: ingc
        character(len=16) :: source_mode
        integer :: ntt, tt, npos, nconfigmax
        integer :: nsubsmin, nsubsmax
        integer, allocatable :: mm(:), omega(:,:)
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: ene(:,:), enemun(:,:), enemunrel(:,:)
        real(dp), allocatable :: p(:,:), pinf(:,:)
        real(dp), allocatable :: pn(:,:), tau(:)
        real(dp), allocatable :: z(:), e(:), f(:), s(:), mu_values(:), xeq_values(:)
        integer :: unit_prob, unit_thermo
        real(dp) :: x, mu, xeq
        real(dp) :: einf, sinf
        character(len=16) :: header_source

        call parse_arguments_gc(source_mode, arg_offset)
        call print_sod_banner()
        call read_ingc_file('INGC', ingc)
        call read_temperatures_file('TEMPERATURES', t, ntt)

        nsubsmin = ingc%nsubsmin
        nsubsmax = ingc%nsubsmax

        call load_gc_level_dataset(nsubsmin, nsubsmax, trim(source_mode), mm, omega, ene, npos, nconfigmax, header_source)

        if (trim(ingc%xormu) == 'x') then
            x = ingc%xormuvalue
            if (x < real(nsubsmin, dp) / real(npos, dp) .or. x > real(nsubsmax, dp) / real(npos, dp)) then
                write(error_unit,'(A)') 'Error: x is out of range for the levels available in INGC.'
                stop 1
            end if
            mu = 0.0_dp
        else
            mu = ingc%xormuvalue
            x = 0.0_dp
        end if

        allocate(enemun(nsubsmin:nsubsmax, 1:nconfigmax))
        allocate(enemunrel(nsubsmin:nsubsmax, 1:nconfigmax))
        allocate(p(nsubsmin:nsubsmax, 1:nconfigmax))
        allocate(pinf(nsubsmin:nsubsmax, 1:nconfigmax))
        allocate(pn(nsubsmin:nsubsmax, 1:ntt))
        allocate(tau(nsubsmin:nsubsmax))
        allocate(z(1:ntt))
        allocate(e(1:ntt))
        allocate(f(1:ntt))
        allocate(s(1:ntt))
        allocate(mu_values(1:ntt))
        allocate(xeq_values(1:ntt))

        call open_gc_outputs(unit_prob, unit_thermo)

        write(*,'(A)') '--- Grand-canonical analysis ---'
        write(*,'(A,I0,A,I0)') 'Levels considered: ', nsubsmin, ' .. ', nsubsmax
        write(*,'(A,I0)') 'Substitutable sites (npos): ', npos
        write(*,'(A,A)') 'Energy source: ', trim(header_source)
        if (trim(ingc%xormu) == 'x') then
            write(*,'(A,F14.6)') 'Target composition x: ', x
        else
            write(*,'(A,F14.6,A)') 'Chemical potential mu: ', mu, ' eV'
        end if
        write(*,'(A,I0)') 'Temperature points: ', ntt
        write(*,*)
        flush(output_unit)

        do tt = 1, ntt
            call compute_gc_for_temperature(t(tt), nsubsmin, nsubsmax, npos, nconfigmax, mm, omega, ene, ingc, &
                x, mu, enemun, enemunrel, p, pn(:,tt), z(tt), e(tt), f(tt), s(tt), mu_values(tt), xeq_values(tt))
            call write_probability_block(unit_prob, t(tt), nsubsmin, nsubsmax, mm, omega, ene, enemun, p, pn(:,tt), &
                mu_values(tt), xeq_values(tt))
            write(unit_thermo,'(F10.1,2X,F14.4,2X,F14.4,2X,E12.6)') t(tt), e(tt), f(tt), s(tt)
        end do

        if (trim(ingc%xormu) == 'mu') x = xeq_values(ntt)

        call compute_infinite_limit(nsubsmin, nsubsmax, npos, nconfigmax, mm, omega, ene, x, pinf, tau, xeq, einf, sinf)
        call write_infinite_probability_block(unit_prob, nsubsmin, nsubsmax, mm, omega, ene, pinf, xeq)
        write(unit_thermo,'(A10,2X,F14.4,8X,A3,7X,E12.6)') 'Infinite', einf, ' - ', sinf

        close(unit_prob)
        close(unit_thermo)

        write(*,'(A)') 'Grand-canonical analysis completed.'
        write(*,'(A)') 'Outputs written to probabilities.dat and thermodynamics.dat.'
        flush(output_unit)
    end subroutine run_sod_gc

    ! Parses the gc-mode command line.
    ! Analiza la línea de comandos del modo gc.
    subroutine parse_arguments_gc(source_mode, arg_offset)
        character(len=*), intent(out) :: source_mode
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset, eq_pos
        character(len=256) :: arg, lowered

        source_mode = 'auto'
        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)
        if (argc <= offset) return

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            if (is_help_token(arg)) then
                call print_usage_gc()
                stop
            end if

            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--source=') == 1) then
                eq_pos = index(arg, '=')
                source_mode = adjustl(arg(eq_pos+1:))
                call normalize_source_mode(source_mode)
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('--source')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --source.'
                    call print_usage_gc()
                    stop 1
                end if
                call get_command_argument(iarg + 1, source_mode)
                call normalize_source_mode(source_mode)
                iarg = iarg + 2
            case default
                write(error_unit,'(A)') 'Error: unrecognized argument in gc mode.'
                call print_usage_gc()
                stop 1
            end select
        end do
    end subroutine parse_arguments_gc

    ! Prints the command-line usage for gc mode.
    ! Imprime la ayuda de línea de comandos para el modo gc.
    subroutine print_usage_gc()
        character(len=256) :: command_name

        command_name = compose_mode_command('gc')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [--source auto|legacy|n|x]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: gc, grand, populations'
        write(*,'(A)') ''
        write(*,'(A)') 'The mode reads INGC in the current directory and combines OUTSOD/ENERGIES'
        write(*,'(A)') 'over the requested range of substitution levels.'
        write(*,'(A)') ''
        write(*,'(A)') 'Optional arguments:'
        write(*,'(A)') '  --source auto    Prefer legacy OUTSOD_XX/ENERGIES_XX, then nNN, then xNN (default).'
        write(*,'(A)') '  --source legacy  Use OUTSOD_XX and ENERGIES_XX in the current directory.'
        write(*,'(A)') '  --source n       Use nNN/OUTSOD and nNN/ENERGIES.'
        write(*,'(A)') '  --source x       Use xNN/OUTSOD and xNN/ENERGIES.'
        write(*,'(A)') ''
        write(*,'(A)') 'Outputs: probabilities.dat and thermodynamics.dat'
    end subroutine print_usage_gc

    ! Normalizes the requested gc data source.
    ! Normaliza la fuente de datos solicitada para gc.
    subroutine normalize_source_mode(source_mode)
        character(len=*), intent(inout) :: source_mode

        source_mode = adjustl(source_mode)
        call to_lower_inplace(source_mode)
        select case (trim(source_mode))
        case ('auto', 'legacy', 'n', 'x')
            continue
        case default
            write(error_unit,'(A)') 'Error: invalid --source value. Use auto, legacy, n, or x.'
            stop 1
        end select
    end subroutine normalize_source_mode

    ! Reads INGC using the legacy gcstatsod field order.
    ! Lee INGC usando el orden de campos heredado de gcstatsod.
    subroutine read_ingc_file(filename, ingc)
        character(len=*), intent(in) :: filename
        type(ingc_data_t), intent(out) :: ingc
        integer :: unit_in, ios
        character(len=1024) :: line

        open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: INGC file does not exist, but it is needed for grand-canonical analysis.'
            stop 1
        end if

        call read_next_data_line(unit_in, line, ios)
        if (ios /= 0) stop 'Error: failed to read nsubs range from INGC.'
        read(line, *, iostat=ios) ingc%nsubsmin, ingc%nsubsmax
        if (ios /= 0) stop 'Error: invalid nsubs range in INGC.'

        call read_next_data_line(unit_in, line, ios)
        if (ios /= 0) stop 'Error: failed to read x/mu selector from INGC.'
        read(line, *, iostat=ios) ingc%xormu, ingc%xormuvalue
        if (ios /= 0) stop 'Error: invalid x/mu selector in INGC.'
        ingc%xormu = adjustl(ingc%xormu)
        call to_lower_inplace(ingc%xormu)
        if (trim(ingc%xormu) == 'x') then
            ingc%xormu = 'x '
        else if (trim(ingc%xormu) == 'mu') then
            ingc%xormu = 'mu'
        else
            stop 'Error: INGC must specify either x or mu.'
        end if

        call read_next_data_line(unit_in, line, ios)
        if (ios /= 0) stop 'Error: failed to read lambda from INGC.'
        read(line, *, iostat=ios) ingc%lambda
        if (ios /= 0) stop 'Error: invalid lambda in INGC.'

        call read_next_data_line(unit_in, line, ios)
        if (ios /= 0) stop 'Error: failed to read volume parameters from INGC.'
        read(line, *, iostat=ios) ingc%v0, ingc%v1, ingc%bv
        if (ios /= 0) stop 'Error: invalid volume parameters in INGC.'

        call read_next_data_line(unit_in, line, ios)
        if (ios /= 0) stop 'Error: failed to read bulk-modulus parameters from INGC.'
        read(line, *, iostat=ios) ingc%bm0, ingc%bm1, ingc%bb
        if (ios /= 0) stop 'Error: invalid bulk-modulus parameters in INGC.'

        close(unit_in)
    end subroutine read_ingc_file

    ! Reads TEMPERATURES or falls back to 300 K and 1000 K.
    ! Lee TEMPERATURES o usa 300 K y 1000 K como valores por defecto.
    subroutine read_temperatures_file(filename, temperatures, count)
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: temperatures(:)
        integer, intent(out) :: count
        integer :: unit_in, ios, ntt
        character(len=1024) :: line
        logical :: exists

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            count = max_temps_default
            allocate(temperatures(count))
            temperatures(1) = 300.0_dp
            temperatures(2) = 1000.0_dp
            return
        end if

        open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to open TEMPERATURES.'
            stop 1
        end if

        ntt = 0
        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = strip_comment(line)
            if (len_trim(line) == 0) cycle
            ntt = ntt + 1
        end do
        rewind(unit_in)

        if (ntt <= 0) then
            close(unit_in)
            count = max_temps_default
            allocate(temperatures(count))
            temperatures(1) = 300.0_dp
            temperatures(2) = 1000.0_dp
            return
        end if

        count = ntt
        allocate(temperatures(count))
        ntt = 0
        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = strip_comment(line)
            if (len_trim(line) == 0) cycle
            ntt = ntt + 1
            read(line, *, iostat=ios) temperatures(ntt)
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid temperature in TEMPERATURES.'
                stop 1
            end if
        end do
        close(unit_in)
    end subroutine read_temperatures_file

    ! Loads OUTSOD and ENERGIES for the requested gc level range.
    ! Carga OUTSOD y ENERGIES para el rango de niveles solicitado por gc.
    subroutine load_gc_level_dataset(nsubsmin, nsubsmax, source_mode, mm, omega, ene, npos, nconfigmax, header_source)
        integer, intent(in) :: nsubsmin, nsubsmax
        character(len=*), intent(in) :: source_mode
        integer, allocatable, intent(out) :: mm(:), omega(:,:)
        real(dp), allocatable, intent(out) :: ene(:,:)
        integer, intent(out) :: npos, nconfigmax
        character(len=*), intent(out) :: header_source
        integer :: nsubs, mm_level, level_npos
        integer, allocatable :: omega_level(:)
        real(dp), allocatable :: ene_level(:)
        character(len=512) :: outsod_path, energies_path
        logical :: found
        character(len=16) :: level_source

        allocate(mm(nsubsmin:nsubsmax))
        mm = 0
        npos = -1
        nconfigmax = 0
        header_source = 'unknown'

        do nsubs = nsubsmin, nsubsmax
            call locate_level_files(nsubs, trim(source_mode), outsod_path, energies_path, level_source, found)
            if (.not. found) then
                write(error_unit,'(A,I0,A)') 'Error: could not locate OUTSOD/ENERGIES for level ', nsubs, '.'
                stop 1
            end if
            call read_outsod_level(outsod_path, level_npos, mm_level, omega_level)
            call read_energy_level(energies_path, mm_level, ene_level)
            if (npos < 0) then
                npos = level_npos
                header_source = trim(level_source)
            else if (level_npos /= npos) then
                write(error_unit,'(A,I0,A)') 'Error: inconsistent number of substitutable sites detected while reading level ', nsubs, '.'
                stop 1
            end if
            mm(nsubs) = mm_level
            nconfigmax = max(nconfigmax, mm_level)
            deallocate(omega_level)
            deallocate(ene_level)
        end do

        allocate(omega(nsubsmin:nsubsmax, 1:nconfigmax))
        allocate(ene(nsubsmin:nsubsmax, 1:nconfigmax))
        omega = 0
        ene = 0.0_dp

        do nsubs = nsubsmin, nsubsmax
            call locate_level_files(nsubs, trim(source_mode), outsod_path, energies_path, level_source, found)
            call read_outsod_level(outsod_path, level_npos, mm_level, omega_level)
            call read_energy_level(energies_path, mm_level, ene_level)
            omega(nsubs, 1:mm_level) = omega_level(1:mm_level)
            ene(nsubs, 1:mm_level) = ene_level(1:mm_level)
            deallocate(omega_level)
            deallocate(ene_level)
        end do
    end subroutine load_gc_level_dataset

    ! Locates the OUTSOD and ENERGIES files for one level using the selected source preference.
    ! Localiza los archivos OUTSOD y ENERGIES de un nivel usando la preferencia de fuente seleccionada.
    subroutine locate_level_files(level, source_mode, outsod_path, energies_path, source_label, found)
        integer, intent(in) :: level
        character(len=*), intent(in) :: source_mode
        character(len=*), intent(out) :: outsod_path, energies_path, source_label
        logical, intent(out) :: found
        character(len=32) :: level_tag
        logical :: outsod_exists, energies_exists

        write(level_tag,'(I2.2)') level
        outsod_path = ''
        energies_path = ''
        source_label = ''
        found = .false.

        select case (trim(source_mode))
        case ('legacy')
            call probe_candidate('OUTSOD_'//trim(level_tag), 'ENERGIES_'//trim(level_tag), 'legacy', found, outsod_path, energies_path, source_label)
        case ('n')
            call probe_candidate(join_path(format_level_directory('n', level), 'OUTSOD'), &
                join_path(format_level_directory('n', level), 'ENERGIES'), 'n', found, outsod_path, energies_path, source_label)
        case ('x')
            call probe_candidate(join_path(format_level_directory('x', level), 'OUTSOD'), &
                join_path(format_level_directory('x', level), 'ENERGIES'), 'x', found, outsod_path, energies_path, source_label)
        case default
            call probe_candidate('OUTSOD_'//trim(level_tag), 'ENERGIES_'//trim(level_tag), 'legacy', found, outsod_path, energies_path, source_label)
            if (.not. found) then
                call probe_candidate(join_path(format_level_directory('n', level), 'OUTSOD'), &
                    join_path(format_level_directory('n', level), 'ENERGIES'), 'n', found, outsod_path, energies_path, source_label)
            end if
            if (.not. found) then
                call probe_candidate(join_path(format_level_directory('x', level), 'OUTSOD'), &
                    join_path(format_level_directory('x', level), 'ENERGIES'), 'x', found, outsod_path, energies_path, source_label)
            end if
        end select

    contains
        subroutine probe_candidate(outsod_candidate, energies_candidate, label, located, outsod_final, energies_final, label_final)
            character(len=*), intent(in) :: outsod_candidate, energies_candidate, label
            logical, intent(inout) :: located
            character(len=*), intent(inout) :: outsod_final, energies_final, label_final

            if (located) return
            inquire(file=trim(outsod_candidate), exist=outsod_exists)
            inquire(file=trim(energies_candidate), exist=energies_exists)
            if (outsod_exists .and. energies_exists) then
                outsod_final = trim(outsod_candidate)
                energies_final = trim(energies_candidate)
                label_final = trim(label)
                located = .true.
            end if
        end subroutine probe_candidate
    end subroutine locate_level_files

    ! Reads one OUTSOD file and extracts npos, the number of unique states, and their degeneracies.
    ! Lee un archivo OUTSOD y extrae npos, el número de estados únicos y sus degeneraciones.
    subroutine read_outsod_level(filename, npos, mm_level, omega_level)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: npos, mm_level
        integer, allocatable, intent(out) :: omega_level(:)
        integer :: unit_in, ios, idx, dummy_nsubs, dummy_index, deg
        character(len=1024) :: line
        character(len=32) :: word

        open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: failed to open OUTSOD file', trim(filename)
            stop 1
        end if

        read(unit_in, '(A)', iostat=ios) line
        if (ios /= 0) stop 'Error: failed to read OUTSOD header.'
        read(line, *, iostat=ios) dummy_nsubs, word, word, npos
        if (ios /= 0) stop 'Error: invalid OUTSOD header.'

        read(unit_in, *, iostat=ios) mm_level
        if (ios /= 0 .or. mm_level <= 0) stop 'Error: invalid unique-count line in OUTSOD.'

        allocate(omega_level(mm_level))
        do idx = 1, mm_level
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) stop 'Error: failed to read a configuration line from OUTSOD.'
            read(line, *, iostat=ios) dummy_index, deg
            if (ios /= 0) stop 'Error: invalid configuration line in OUTSOD.'
            omega_level(idx) = deg
        end do
        close(unit_in)
    end subroutine read_outsod_level

    ! Reads one ENERGIES file and returns a single energy per inequivalent configuration.
    ! Lee un archivo ENERGIES y devuelve una sola energía por configuración inequivalente.
    subroutine read_energy_level(filename, mm_level, ene_level)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: mm_level
        real(dp), allocatable, intent(out) :: ene_level(:)
        integer :: unit_in, ios, idx, numeric_count
        character(len=1024) :: line
        real(dp) :: numeric_tokens(max_tokens_per_line)

        allocate(ene_level(mm_level))
        open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: failed to open ENERGIES file', trim(filename)
            stop 1
        end if

        idx = 0
        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = strip_comment(line)
            if (len_trim(line) == 0) cycle
            call extract_numeric_tokens(line, numeric_tokens, numeric_count)
            if (numeric_count <= 0) cycle
            idx = idx + 1
            if (idx > mm_level) exit
            ene_level(idx) = numeric_tokens(numeric_count)
        end do
        close(unit_in)

        if (idx /= mm_level) then
            write(error_unit,'(A,1X,A)') 'Error: ENERGIES count does not match OUTSOD count for', trim(filename)
            stop 1
        end if
    end subroutine read_energy_level

    ! Computes finite-temperature grand-canonical probabilities and thermodynamics for one temperature.
    ! Calcula probabilidades grand-canónicas y termodinámica a temperatura finita para una temperatura dada.
    subroutine compute_gc_for_temperature(temp, nsubsmin, nsubsmax, npos, nconfigmax, mm, omega, ene, ingc, x_target, mu_input, &
        enemun, enemunrel, p, pn, z, e, f, s, mu_out, xeq_out)
        real(dp), intent(in) :: temp
        integer, intent(in) :: nsubsmin, nsubsmax, npos, nconfigmax
        integer, intent(in) :: mm(nsubsmin:nsubsmax), omega(nsubsmin:nsubsmax, 1:nconfigmax)
        real(dp), intent(in) :: ene(nsubsmin:nsubsmax, 1:nconfigmax)
        type(ingc_data_t), intent(in) :: ingc
        real(dp), intent(in) :: x_target, mu_input
        real(dp), intent(out) :: enemun(nsubsmin:nsubsmax, 1:nconfigmax)
        real(dp), intent(out) :: enemunrel(nsubsmin:nsubsmax, 1:nconfigmax)
        real(dp), intent(out) :: p(nsubsmin:nsubsmax, 1:nconfigmax)
        real(dp), intent(out) :: pn(nsubsmin:nsubsmax)
        real(dp), intent(out) :: z, e, f, s, mu_out, xeq_out
        real(dp) :: x, mu, nx, emin
        real(dp) :: qn(nsubsmin:nsubsmax), c(nsubsmin:nsubsmax), evsc(nsubsmin:nsubsmax)
        real(dp) :: naver, eneaver, epsilona, epsilonb, epsilon, e0
        integer :: nsubs, m, memin, nsubsemin

        if (trim(ingc%xormu) == 'x') then
            x = x_target
            call compute_weighted_line_fit(nsubsmin, nsubsmax, mm, omega, ene, naver, eneaver, epsilona, epsilonb, epsilon, e0)
            call compute_evsc(nsubsmin, nsubsmax, npos, x, ingc, evsc)
            call compute_qn(nsubsmin, nsubsmax, mm, omega, ene, evsc, epsilon, e0, temp, qn)
            do nsubs = nsubsmin, nsubsmax
                c(nsubs) = (real(nsubs, dp) / real(npos, dp) - x) * qn(nsubs)
            end do
            mu = solve_mu_from_x(nsubsmin, nsubsmax, temp, x, c, epsilon)
            do nsubs = nsubsmin, nsubsmax
                do m = 1, mm(nsubs)
                    enemun(nsubs, m) = ene(nsubs, m) + evsc(nsubs) - real(nsubs, dp) * mu
                end do
            end do
        else
            mu = mu_input
            do nsubs = nsubsmin, nsubsmax
                evsc(nsubs) = 0.0_dp
                do m = 1, mm(nsubs)
                    enemun(nsubs, m) = ene(nsubs, m) - real(nsubs, dp) * mu
                end do
            end do
        end if

        emin = enemun(nsubsmin, 1)
        memin = 1
        nsubsemin = nsubsmin
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                if (emin > enemun(nsubs, m)) then
                    emin = enemun(nsubs, m)
                    memin = m
                    nsubsemin = nsubs
                end if
            end do
        end do

        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                enemunrel(nsubs, m) = enemun(nsubs, m) - emin
            end do
        end do

        z = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                z = z + real(omega(nsubs, m), dp) * exp(-enemunrel(nsubs, m) / (kB_eVk * temp))
            end do
        end do

        nx = 0.0_dp
        pn = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                p(nsubs, m) = real(omega(nsubs, m), dp) * exp(-enemunrel(nsubs, m) / (kB_eVk * temp)) / z
                pn(nsubs) = pn(nsubs) + p(nsubs, m)
                nx = nx + real(nsubs, dp) * p(nsubs, m)
            end do
        end do
        xeq_out = nx / real(npos, dp)
        mu_out = mu

        e = 0.0_dp
        s = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                e = e + ene(nsubs, m) * p(nsubs, m)
                if (p(nsubs, m) > 0.0_dp) then
                    s = s - kB_eVk * p(nsubs, m) * log(p(nsubs, m) / real(omega(nsubs, m), dp))
                end if
            end do
        end do
        f = e - temp * s
    end subroutine compute_gc_for_temperature

    ! Computes the approximate linear coefficients used by the original gcstatsod solver.
    ! Calcula los coeficientes lineales aproximados usados por el solver original de gcstatsod.
    subroutine compute_weighted_line_fit(nsubsmin, nsubsmax, mm, omega, ene, naver, eneaver, epsilona, epsilonb, epsilon, e0)
        integer, intent(in) :: nsubsmin, nsubsmax
        integer, intent(in) :: mm(nsubsmin:), omega(nsubsmin:,:)
        real(dp), intent(in) :: ene(nsubsmin:,:)
        real(dp), intent(out) :: naver, eneaver, epsilona, epsilonb, epsilon, e0
        integer :: nsubs, m
        real(dp) :: omegasum

        omegasum = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                omegasum = omegasum + real(omega(nsubs, m), dp)
            end do
        end do

        naver = 0.0_dp
        eneaver = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                naver = naver + real(omega(nsubs, m), dp) * real(nsubs, dp)
                eneaver = eneaver + real(omega(nsubs, m), dp) * ene(nsubs, m)
            end do
        end do
        naver = naver / omegasum
        eneaver = eneaver / omegasum

        epsilona = 0.0_dp
        epsilonb = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                epsilona = epsilona + real(omega(nsubs, m), dp) * (real(nsubs, dp) - naver) * (ene(nsubs, m) - eneaver)
                epsilonb = epsilonb + real(omega(nsubs, m), dp) * (real(nsubs, dp) - naver) * (real(nsubs, dp) - naver)
            end do
        end do
        epsilon = epsilona / epsilonb
        e0 = eneaver - epsilon * naver
    end subroutine compute_weighted_line_fit

    ! Computes the legacy volume-stress correction EVSC(n).
    ! Calcula la corrección histórica de volumen-esfuerzo EVSC(n).
    subroutine compute_evsc(nsubsmin, nsubsmax, npos, x, ingc, evsc)
        integer, intent(in) :: nsubsmin, nsubsmax, npos
        real(dp), intent(in) :: x
        type(ingc_data_t), intent(in) :: ingc
        real(dp), intent(out) :: evsc(nsubsmin:nsubsmax)
        integer :: nsubs
        real(dp) :: vol, bm, xn, voln, eta

        if (ingc%lambda == 0.0_dp) then
            evsc = 0.0_dp
            return
        end if

        vol = ingc%v0 * (1.0_dp - x) + ingc%v1 * x + ingc%bv * x * (1.0_dp - x)
        bm = (ingc%bm0 * (1.0_dp - x) + ingc%bm1 * x + ingc%bb * x * (1.0_dp - x)) / eva3togpa
        do nsubs = nsubsmin, nsubsmax
            xn = real(nsubs, dp) / real(npos, dp)
            voln = ingc%v0 * (1.0_dp - xn) + ingc%v1 * xn + ingc%bv * xn * (1.0_dp - xn)
            eta = (vol / voln) ** (2.0_dp / 3.0_dp)
            evsc(nsubs) = (9.0_dp / 8.0_dp) * vol * bm * (eta - 1.0_dp) ** 2
        end do
    end subroutine compute_evsc

    ! Computes Q_n coefficients used in the original x-to-mu conversion.
    ! Calcula los coeficientes Q_n usados en la conversión original de x a mu.
    subroutine compute_qn(nsubsmin, nsubsmax, mm, omega, ene, evsc, epsilon, e0, temp, qn)
        integer, intent(in) :: nsubsmin, nsubsmax
        integer, intent(in) :: mm(nsubsmin:), omega(nsubsmin:,:)
        real(dp), intent(in) :: ene(nsubsmin:,:), evsc(nsubsmin:)
        real(dp), intent(in) :: epsilon, e0, temp
        real(dp), intent(out) :: qn(nsubsmin:nsubsmax)
        integer :: nsubs, m

        qn = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                qn(nsubs) = qn(nsubs) + real(omega(nsubs, m), dp) * &
                    exp(-(ene(nsubs, m) + evsc(nsubs) - e0 - real(nsubs, dp) * epsilon) / (kB_eVk * temp))
            end do
        end do
    end subroutine compute_qn

    ! Solves for mu at fixed x using a bracketing+bisection scheme over q = exp((mu-epsilon)/kT).
    ! Resuelve mu a x fijo usando un esquema de acotado+bisección sobre q = exp((mu-epsilon)/kT).
    real(dp) function solve_mu_from_x(nsubsmin, nsubsmax, temp, x, c, epsilon)
        integer, intent(in) :: nsubsmin, nsubsmax
        real(dp), intent(in) :: temp, x, c(nsubsmin:nsubsmax), epsilon
        real(dp) :: qa, qb, q, va, vb, vq, oldmu
        integer :: iter

        q = max(tiny_q, x / max(1.0_dp - x, tiny_q))
        qa = max(tiny_q, q * 0.5_dp)
        qb = max(qa * 2.0_dp, q * 2.0_dp)
        va = evaluate_polynomial(nsubsmin, nsubsmax, c, qa)
        vb = evaluate_polynomial(nsubsmin, nsubsmax, c, qb)

        iter = 0
        do while (va * vb > 0.0_dp .and. iter < 200)
            qa = max(tiny_q, qa * 0.5_dp)
            qb = qb * 2.0_dp
            va = evaluate_polynomial(nsubsmin, nsubsmax, c, qa)
            vb = evaluate_polynomial(nsubsmin, nsubsmax, c, qb)
            iter = iter + 1
        end do
        if (va * vb > 0.0_dp) then
            write(error_unit,'(A)') 'Error: failed to bracket the q root needed to convert x into mu.'
            stop 1
        end if

        oldmu = epsilon + kB_eVk * temp * log(q)
        do iter = 1, 1000
            q = 0.5_dp * (qa + qb)
            vq = evaluate_polynomial(nsubsmin, nsubsmax, c, q)
            solve_mu_from_x = epsilon + kB_eVk * temp * log(q)
            if (vq * va > 0.0_dp) then
                qa = q
                va = vq
            else
                qb = q
                vb = vq
            end if
            if (abs(qa - qb) < tolq .and. abs(solve_mu_from_x - oldmu) < tolmu) exit
            oldmu = solve_mu_from_x
        end do
    end function solve_mu_from_x

    ! Evaluates the polynomial sum_n c_n q^n.
    ! Evalúa el polinomio suma_n c_n q^n.
    real(dp) function evaluate_polynomial(nsubsmin, nsubsmax, c, q)
        integer, intent(in) :: nsubsmin, nsubsmax
        real(dp), intent(in) :: c(nsubsmin:nsubsmax), q
        integer :: nsubs

        evaluate_polynomial = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            evaluate_polynomial = evaluate_polynomial + c(nsubs) * q ** nsubs
        end do
    end function evaluate_polynomial

    ! Computes the ideal-disorder-limit probabilities using the original residual-momenta correction.
    ! Calcula las probabilidades del límite de desorden ideal usando la corrección original de momentos residuales.
    subroutine compute_infinite_limit(nsubsmin, nsubsmax, npos, nconfigmax, mm, omega, ene, x, pinf, tau, xeq, einf, sinf)
        integer, intent(in) :: nsubsmin, nsubsmax, npos, nconfigmax
        integer, intent(in) :: mm(nsubsmin:nsubsmax), omega(nsubsmin:nsubsmax, 1:nconfigmax)
        real(dp), intent(in) :: ene(nsubsmin:nsubsmax, 1:nconfigmax), x
        real(dp), intent(out) :: pinf(nsubsmin:nsubsmax, 1:nconfigmax), tau(nsubsmin:nsubsmax)
        real(dp), intent(out) :: xeq, einf, sinf
        real(dp) :: a0, a1, a2, b0, b1, b2
        real(dp) :: a11, a12, a21, a22, bb1, bb2, alpha, beta
        real(dp) :: nx
        integer :: nsubs, m

        a0 = momentaa0(nsubsmax, npos, x)
        a1 = momentaa1(nsubsmax, npos, x)
        a2 = momentaa2(nsubsmax, npos, x)
        b0 = momentab0(nsubsmin, npos, x)
        b1 = momentab1(nsubsmin, npos, x)
        b2 = momentab2(nsubsmin, npos, x)

        a11 = 1.0_dp - b0 - a0
        a12 = x * real(npos, dp) - b1 - a1
        a21 = a12
        a22 = (x * real(npos, dp) * (1.0_dp + x * real(npos - 1, dp))) - b2 - a2
        bb1 = 1.0_dp
        bb2 = x * real(npos, dp)

        alpha = ((bb1 * a22) - (a12 * bb2)) / ((a11 * a22) - (a12 * a21))
        beta = ((a11 * bb2) - (a21 * bb1)) / ((a11 * a22) - (a12 * a21))

        do nsubs = nsubsmin, nsubsmax
            tau(nsubs) = alpha + beta * real(nsubs, dp)
        end do

        pinf = 0.0_dp
        nx = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                pinf(nsubs, m) = tau(nsubs) * real(omega(nsubs, m), dp) * (x ** nsubs) * ((1.0_dp - x) ** (npos - nsubs))
                nx = nx + real(nsubs, dp) * pinf(nsubs, m)
            end do
        end do
        xeq = nx / real(npos, dp)

        einf = 0.0_dp
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                einf = einf + ene(nsubs, m) * pinf(nsubs, m)
            end do
        end do

        if (x == 0.0_dp .or. x == 1.0_dp) then
            sinf = 0.0_dp
        else
            sinf = -kB_eVk * (x * log(x) + (1.0_dp - x) * log(1.0_dp - x)) * real(npos, dp)
        end if
    end subroutine compute_infinite_limit

    ! Opens the legacy gc output files.
    ! Abre los archivos de salida heredados del modo gc.
    subroutine open_gc_outputs(unit_prob, unit_thermo)
        integer, intent(out) :: unit_prob, unit_thermo
        integer :: ios

        open(newunit=unit_prob, file='probabilities.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to create probabilities.dat.'
            stop 1
        end if

        open(newunit=unit_thermo, file='thermodynamics.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: failed to create thermodynamics.dat.'
            stop 1
        end if
        write(unit_thermo,'(A)') '        T             E               F          S             '
    end subroutine open_gc_outputs

    ! Writes one finite-temperature block to probabilities.dat.
    ! Escribe un bloque a temperatura finita en probabilities.dat.
    subroutine write_probability_block(unit_prob, temp, nsubsmin, nsubsmax, mm, omega, ene, enemun, p, pn, mu, xeq)
        integer, intent(in) :: unit_prob, nsubsmin, nsubsmax
        integer, intent(in) :: mm(nsubsmin:), omega(nsubsmin:,:)
        real(dp), intent(in) :: temp, ene(nsubsmin:,:), enemun(nsubsmin:,:), p(nsubsmin:,:)
        real(dp), intent(in) :: pn(nsubsmin:nsubsmax), mu, xeq
        integer :: nsubs, m

        write(unit_prob,'(A)') ''
        write(unit_prob,'(A)') ' ____________________________________________________________________________________________'
        write(unit_prob,'(A,F24.10,A)') ' Temperature: T = ', temp, ' K'
        write(unit_prob,'(A,F24.10,A)') ' Chemical potential: mu = ', mu, ' eV'
        write(unit_prob,'(A,F24.10)') ' Composition: x = ', xeq
        write(unit_prob,'(A)') ''
        write(unit_prob,'(A)') '  n      m     omega(n,m)       E(n,m)          E-n*mu         prob(n,m,T)   prob(n,m,T)/omega    '
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                write(unit_prob,'(I3,1X,I6,1X,I8,5X,2(4X,F14.5),2(4X,E12.6))') &
                    nsubs, m, omega(nsubs, m), ene(nsubs, m), enemun(nsubs, m), p(nsubs, m), p(nsubs, m)/real(omega(nsubs, m), dp)
            end do
        end do
        write(unit_prob,'(A)') ' ------------------------------'
        write(unit_prob,'(A)') '  n      Cumulative prob(n,T) '
        do nsubs = nsubsmin, nsubsmax
            write(unit_prob,'(I6,5X,F14.5)') nsubs, pn(nsubs)
        end do
        write(unit_prob,'(A)') ' ------------------------------'
    end subroutine write_probability_block

    ! Writes the infinite-temperature limit block to probabilities.dat.
    ! Escribe en probabilities.dat el bloque del límite de temperatura infinita.
    subroutine write_infinite_probability_block(unit_prob, nsubsmin, nsubsmax, mm, omega, ene, pinf, xeq)
        integer, intent(in) :: unit_prob, nsubsmin, nsubsmax
        integer, intent(in) :: mm(nsubsmin:), omega(nsubsmin:,:)
        real(dp), intent(in) :: ene(nsubsmin:,:), pinf(nsubsmin:,:), xeq
        integer :: nsubs, m

        write(unit_prob,'(A)') ''
        write(unit_prob,'(A)') ' ____________________________________________________________________________________________'
        write(unit_prob,'(A)') ' Ideal disorder limit'
        write(unit_prob,'(A,F24.10)') ' Composition: x = ', xeq
        write(unit_prob,'(A)') '  n      m     omega(n,m)     E(n,m)                           prob(n,m,T)   prob(n,m,T)/omega    '
        do nsubs = nsubsmin, nsubsmax
            do m = 1, mm(nsubs)
                write(unit_prob,'(I3,1X,I6,1X,I8,5X,4X,F14.5,16X,2(4X,E12.6))') &
                    nsubs, m, omega(nsubs, m), ene(nsubs, m), pinf(nsubs, m), pinf(nsubs, m)/real(omega(nsubs, m), dp)
            end do
        end do
    end subroutine write_infinite_probability_block

    ! Reads the next non-empty, non-comment line from a text file.
    ! Lee la siguiente línea no vacía y no comentada de un archivo de texto.
    subroutine read_next_data_line(unit_in, line, ios)
        integer, intent(in) :: unit_in
        character(len=*), intent(out) :: line
        integer, intent(out) :: ios

        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) return
            line = strip_comment(line)
            if (len_trim(line) > 0) return
        end do
    end subroutine read_next_data_line

    ! Strips inline comments introduced by # and trims the remaining text.
    ! Elimina comentarios en línea introducidos por # y recorta el texto restante.
    pure function strip_comment(line) result(clean)
        character(len=*), intent(in) :: line
        character(len=len(line)) :: clean
        integer :: hash_pos

        clean = line
        hash_pos = index(clean, '#')
        if (hash_pos > 0) clean(hash_pos:) = ' '
        clean = adjustl(clean)
    end function strip_comment

    ! Extracts all numeric tokens from a text line.
    ! Extrae todos los tokens numéricos de una línea de texto.
    subroutine extract_numeric_tokens(line, values, count)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: count
        integer :: start_pos, end_pos, ios, line_len
        character(len=128) :: token
        character(len=1024) :: text

        values = 0.0_dp
        count = 0
        text = trim(adjustl(line))
        line_len = len_trim(text)
        start_pos = 1
        do while (start_pos <= line_len)
            do while (start_pos <= line_len .and. (text(start_pos:start_pos) == ' ' .or. text(start_pos:start_pos) == achar(9) .or. text(start_pos:start_pos) == ','))
                start_pos = start_pos + 1
            end do
            if (start_pos > line_len) exit
            end_pos = start_pos
            do while (end_pos <= line_len .and. text(end_pos:end_pos) /= ' ' .and. text(end_pos:end_pos) /= achar(9) .and. text(end_pos:end_pos) /= ',')
                end_pos = end_pos + 1
            end do
            token = ' '
            token = text(start_pos:end_pos-1)
            read(token, *, iostat=ios) values(min(size(values), count + 1))
            if (ios == 0) then
                count = min(size(values), count + 1)
            end if
            start_pos = end_pos + 1
        end do
    end subroutine extract_numeric_tokens

    ! Original binomial probability helper from legacy gcstatsod dependencies.
    ! Función auxiliar original de probabilidad binomial de las dependencias históricas de gcstatsod.
    real(dp) function pbinomial(nsubs, npos, x)
        integer, intent(in) :: nsubs, npos
        real(dp), intent(in) :: x
        pbinomial = real(combinations(nsubs, npos), dp) * (x ** nsubs) * ((1.0_dp - x) ** (npos - nsubs))
    end function pbinomial

    ! Original combinations helper from the legacy factorials.f90 code.
    ! Función auxiliar original de combinaciones del código histórico factorials.f90.
    integer(ip) function combinations(nsubs, npos)
        integer, intent(in) :: nsubs, npos
        integer :: k, p
        integer(ip) :: ratio

        if (nsubs <= npos / 2) then
            k = nsubs
        else
            k = npos - nsubs
        end if
        ratio = 1_ip
        do p = npos - k + 1, npos
            ratio = ratio * int(p, ip)
        end do
        combinations = ratio / factorial_int64(k)
    end function combinations

    ! Original recursive factorial helper from legacy code, adapted to int64.
    ! Función auxiliar recursiva original del factorial, adaptada a int64.
    recursive integer(ip) function factorial_int64(n) result(aux)
        integer, intent(in) :: n
        if (n == 0) then
            aux = 1_ip
        else
            aux = int(n, ip) * factorial_int64(n - 1)
        end if
    end function factorial_int64

    ! Original residual momenta helper a0.
    ! Función original del momento residual a0.
    real(dp) function momentaa0(nsubsmax, npos, x)
        integer, intent(in) :: nsubsmax, npos
        real(dp), intent(in) :: x
        integer :: nsubs

        momentaa0 = 0.0_dp
        do nsubs = nsubsmax + 1, npos
            momentaa0 = momentaa0 + pbinomial(nsubs, npos, x)
        end do
    end function momentaa0

    ! Original residual momenta helper a1.
    ! Función original del momento residual a1.
    real(dp) function momentaa1(nsubsmax, npos, x)
        integer, intent(in) :: nsubsmax, npos
        real(dp), intent(in) :: x
        integer :: nsubs

        momentaa1 = 0.0_dp
        do nsubs = nsubsmax + 1, npos
            momentaa1 = momentaa1 + pbinomial(nsubs, npos, x) * real(nsubs, dp)
        end do
    end function momentaa1

    ! Original residual momenta helper a2.
    ! Función original del momento residual a2.
    real(dp) function momentaa2(nsubsmax, npos, x)
        integer, intent(in) :: nsubsmax, npos
        real(dp), intent(in) :: x
        integer :: nsubs

        momentaa2 = 0.0_dp
        do nsubs = nsubsmax + 1, npos
            momentaa2 = momentaa2 + pbinomial(nsubs, npos, x) * (real(nsubs, dp) ** 2)
        end do
    end function momentaa2

    ! Original residual momenta helper b0.
    ! Función original del momento residual b0.
    real(dp) function momentab0(nsubsmin, npos, x)
        integer, intent(in) :: nsubsmin, npos
        real(dp), intent(in) :: x
        integer :: nsubs

        momentab0 = 0.0_dp
        do nsubs = 0, nsubsmin - 1
            momentab0 = momentab0 + pbinomial(nsubs, npos, x)
        end do
    end function momentab0

    ! Original residual momenta helper b1.
    ! Función original del momento residual b1.
    real(dp) function momentab1(nsubsmin, npos, x)
        integer, intent(in) :: nsubsmin, npos
        real(dp), intent(in) :: x
        integer :: nsubs

        momentab1 = 0.0_dp
        do nsubs = 0, nsubsmin - 1
            momentab1 = momentab1 + pbinomial(nsubs, npos, x) * real(nsubs, dp)
        end do
    end function momentab1

    ! Original residual momenta helper b2.
    ! Función original del momento residual b2.
    real(dp) function momentab2(nsubsmin, npos, x)
        integer, intent(in) :: nsubsmin, npos
        real(dp), intent(in) :: x
        integer :: nsubs

        momentab2 = 0.0_dp
        do nsubs = 0, nsubsmin - 1
            momentab2 = momentab2 + pbinomial(nsubs, npos, x) * (real(nsubs, dp) ** 2)
        end do
    end function momentab2

end module gc
