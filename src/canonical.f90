!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!*******************************************************************************

!*******************************************************************************
! Canonical thermodynamic post-processing for one substitution level.
! This module is a modernized port of the original statsod workflow and keeps
! the same outputs: probabilities.dat, thermodynamics.dat, ave_data.dat, and
! ave_spectra.dat when the corresponding inputs are present.
!
! Postproceso termodinámico canónico para un único nivel de sustitución.
! Este módulo es un port modernizado del flujo original de statsod y conserva
! las mismas salidas: probabilities.dat, thermodynamics.dat, ave_data.dat y
! ave_spectra.dat cuando existen las entradas correspondientes.
!*******************************************************************************
module canonical
    use consts, only: dp, error_unit
    use cli, only: compose_mode_command, is_help_token, to_lower_inplace
    use utils, only: format_level_directory, join_path, print_sod_banner
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    private

    real(dp), parameter :: kb = 8.61734d-5
    real(dp), parameter :: tolprob = 1.0d-12
    real(dp), parameter :: tolminspec = 1.0d-6

    type :: canonical_options_t
        integer :: level = -1
        character(len=16) :: source_mode = 'auto'
    end type canonical_options_t

    public :: run_sod_canonical

contains

    ! Runs the canonical thermodynamics workflow for one level directory or for the current directory.
    ! Ejecuta el flujo de termodinámica canónica para un directorio de nivel o para el directorio actual.
    subroutine run_sod_canonical(arg_offset)
        integer, intent(in), optional :: arg_offset
        type(canonical_options_t) :: options
        character(len=512) :: data_dir
        character(len=16) :: source_label
        integer :: nsubs, mm, ncol, npoints, ntt
        integer :: m, col, point, tt
        integer :: memin, memax
        integer, allocatable :: omega(:)
        real(dp) :: emin, emax, maxspec
        real(dp) :: zinf, einf, sinf
        real(dp), allocatable :: z(:), e(:), f(:), s(:), t(:)
        real(dp), allocatable :: ene(:), erel(:)
        real(dp), allocatable :: p(:,:), pinf(:)
        real(dp), allocatable :: data(:,:), avedata(:,:)
        real(dp), allocatable :: spec(:,:), avespec(:,:)
        real(dp), allocatable :: avedatainf(:), xspec(:), avespecinf(:)
        logical :: temperatures_exists, data_exists, spectra_exists
        character(len=64) :: fmtemplist
        character(len=512) :: outsod_path, energies_path, temperatures_path
        character(len=512) :: data_path, spectra_path, xspec_path
        character(len=512) :: probabilities_path, thermodynamics_path
        character(len=512) :: ave_data_path, ave_spectra_path

        call parse_arguments_canonical(options, arg_offset)
        call locate_canonical_dataset(options, data_dir, source_label)

        outsod_path = join_path(data_dir, 'OUTSOD')
        energies_path = join_path(data_dir, 'ENERGIES')
        temperatures_path = join_path(data_dir, 'TEMPERATURES')
        data_path = join_path(data_dir, 'DATA')
        spectra_path = join_path(data_dir, 'SPECTRA')
        xspec_path = join_path(data_dir, 'XSPEC')
        probabilities_path = join_path(data_dir, 'probabilities.dat')
        thermodynamics_path = join_path(data_dir, 'thermodynamics.dat')
        ave_data_path = join_path(data_dir, 'ave_data.dat')
        ave_spectra_path = join_path(data_dir, 'ave_spectra.dat')

        call print_sod_banner()
        write(*,'(A)') '--- Canonical thermodynamic analysis ---'
        write(*,'(A,A)') 'Dataset: ', trim(data_dir)
        write(*,'(A,A)') 'Source: ', trim(source_label)
        flush(output_unit)

        inquire(file=trim(temperatures_path), exist=temperatures_exists)
        inquire(file=trim(data_path), exist=data_exists)
        inquire(file=trim(spectra_path), exist=spectra_exists)

        call read_temperatures_file(trim(temperatures_path), temperatures_exists, t, ntt)
        call read_outsod_file(trim(outsod_path), nsubs, mm, omega)
        call read_energies_file(trim(energies_path), mm, ene)

        if (data_exists) then
            write(*,'(A)') 'DATA file found'
            write(*,*)
            flush(output_unit)
            call read_data_file(trim(data_path), mm, ncol, data)
            allocate(avedata(ncol, ntt))
            allocate(avedatainf(ncol))
        else
            write(*,'(A)') 'DATA file not found. No averaging of DATA will be performed.'
            flush(output_unit)
            ncol = 0
        end if

        if (spectra_exists) then
            write(*,'(A)') 'SPECTRA file found'
            write(*,*)
            write(*,'(A)') 'Reading SPECTRA file...'
            write(*,*)
            flush(output_unit)
            call read_spectra_files(trim(spectra_path), trim(xspec_path), mm, npoints, spec, xspec)
            maxspec = maxval(spec)
            allocate(avespec(npoints, ntt))
            allocate(avespecinf(npoints))
        else
            write(*,'(A)') 'SPECTRA file not found. No averaging of SPECTRA will be performed.'
            flush(output_unit)
            npoints = 0
            maxspec = 0.0_dp
        end if

        allocate(z(ntt), e(ntt), f(ntt), s(ntt))
        allocate(erel(mm))
        allocate(p(mm, ntt))
        allocate(pinf(mm))

        emin = ene(1)
        emax = ene(1)
        memin = 1
        memax = 1
        do m = 2, mm
            if (emin > ene(m)) then
                emin = ene(m)
                memin = m
            end if
            if (emax < ene(m)) then
                emax = ene(m)
                memax = m
            end if
        end do

        do m = 1, mm
            erel(m) = ene(m) - emin
        end do

        call open_probability_output(trim(probabilities_path), memin, memax)
        call open_thermodynamics_output(trim(thermodynamics_path))
        if (data_exists) call open_ave_data_output(trim(ave_data_path))
        if (spectra_exists) call open_ave_spectra_output(trim(ave_spectra_path), t, ntt)

        do tt = 1, ntt
            z(tt) = 0.0_dp
            do m = 1, mm
                z(tt) = z(tt) + real(omega(m), dp) * exp(-erel(m) / (kb * t(tt)))
            end do

            do m = 1, mm
                p(m, tt) = real(omega(m), dp) * exp(-erel(m) / (kb * t(tt))) / z(tt)
            end do

            e(tt) = 0.0_dp
            do m = 1, mm
                e(tt) = e(tt) + ene(m) * p(m, tt)
            end do
            f(tt) = emin - kb * t(tt) * log(z(tt))
            s(tt) = (e(tt) - f(tt)) / t(tt)

            if (data_exists) then
                avedata(:, tt) = 0.0_dp
                do col = 1, ncol
                    do m = 1, mm
                        avedata(col, tt) = avedata(col, tt) + data(col, m) * p(m, tt)
                    end do
                end do
            end if

            if (spectra_exists) then
                avespec(:, tt) = 0.0_dp
                do point = 1, npoints
                    do m = 1, mm
                        avespec(point, tt) = avespec(point, tt) + spec(point, m) * p(m, tt)
                    end do
                    if (maxspec > 0.0_dp) then
                        if (avespec(point, tt) / maxspec < tolminspec) avespec(point, tt) = 0.0_dp
                    end if
                end do
            end if

            call write_probability_block(trim(probabilities_path), t(tt), mm, omega, erel, p(:, tt))
            call append_thermodynamics(trim(thermodynamics_path), t(tt), e(tt), f(tt), s(tt))
            if (data_exists) call append_ave_data(trim(ave_data_path), t(tt), avedata(:, tt))
        end do

        zinf = 0.0_dp
        do m = 1, mm
            zinf = zinf + real(omega(m), dp)
        end do
        do m = 1, mm
            pinf(m) = real(omega(m), dp) / zinf
        end do
        einf = 0.0_dp
        do m = 1, mm
            einf = einf + ene(m) * pinf(m)
        end do
        sinf = kb * log(zinf)

        call write_infinite_probability_block(trim(probabilities_path), mm, omega, erel, pinf)
        call append_infinite_thermodynamics(trim(thermodynamics_path), einf, sinf)

        if (data_exists) then
            avedatainf = 0.0_dp
            do col = 1, ncol
                do m = 1, mm
                    avedatainf(col) = avedatainf(col) + data(col, m) * pinf(m)
                end do
            end do
            call append_infinite_ave_data(trim(ave_data_path), avedatainf)
        end if

        if (spectra_exists) then
            avespecinf = 0.0_dp
            do point = 1, npoints
                do m = 1, mm
                    avespecinf(point) = avespecinf(point) + spec(point, m) * pinf(m)
                end do
                if (maxspec > 0.0_dp) then
                    if (avespecinf(point) / maxspec < tolminspec) avespecinf(point) = 0.0_dp
                end if
            end do
            write(fmtemplist, '(A,I0,A)') '(F10.3,2X,', ntt + 1, '(E12.6,2X))'
            call append_ave_spectra(trim(ave_spectra_path), trim(fmtemplist), xspec, avespec, avespecinf)
        end if

        write(*,'(A)') 'Canonical analysis completed.'
        write(*,'(A,A)') 'probabilities.dat: ', trim(probabilities_path)
        write(*,'(A,A)') 'thermodynamics.dat: ', trim(thermodynamics_path)
        if (data_exists) write(*,'(A,A)') 'ave_data.dat: ', trim(ave_data_path)
        if (spectra_exists) write(*,'(A,A)') 'ave_spectra.dat: ', trim(ave_spectra_path)
        flush(output_unit)
    end subroutine run_sod_canonical

    ! Parses the canonical-mode arguments.
    ! Analiza los argumentos del modo canónico.
    subroutine parse_arguments_canonical(options, arg_offset)
        type(canonical_options_t), intent(out) :: options
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset, eq_pos
        character(len=256) :: arg
        character(len=256) :: lowered

        options%level = -1
        options%source_mode = 'auto'
        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            if (is_help_token(arg)) then
                call print_usage_canonical()
                stop
            end if

            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--source=') == 1) then
                eq_pos = index(arg, '=')
                options%source_mode = adjustl(arg(eq_pos+1:))
                call normalize_source_mode(options%source_mode)
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-n','-N')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing level after -N.'
                    call print_usage_canonical()
                    stop 1
                end if
                call get_command_argument(iarg + 1, arg)
                read(arg, *, err=900) options%level
                iarg = iarg + 2
            case ('--source')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --source.'
                    call print_usage_canonical()
                    stop 1
                end if
                call get_command_argument(iarg + 1, options%source_mode)
                call normalize_source_mode(options%source_mode)
                iarg = iarg + 2
            case default
                write(error_unit,'(A)') 'Error: unrecognized argument in canonical mode.'
                call print_usage_canonical()
                stop 1
            end select
        end do
        return
900     continue
        write(error_unit,'(A)') 'Error: invalid integer value after -N.'
        call print_usage_canonical()
        stop 1
    end subroutine parse_arguments_canonical

    ! Prints the command-line help for canonical thermodynamics.
    ! Imprime la ayuda de línea de comandos para la termodinámica canónica.
    subroutine print_usage_canonical()
        character(len=256) :: command_name

        command_name = compose_mode_command('canonical')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [-N <level>] [--source auto|legacy|n|x]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: canonical, thermo, stats, c'
        write(*,'(A)') ''
        write(*,'(A)') 'The mode is a modernized port of statsod and processes one level at a time.'
        write(*,'(A)') 'It reads OUTSOD and ENERGIES, plus optional DATA, SPECTRA, XSPEC, and TEMPERATURES.'
        write(*,'(A)') ''
        write(*,'(A)') 'If -N is omitted and OUTSOD/ENERGIES exist in the current directory, that dataset is used.'
        write(*,'(A)') 'If -N is provided, the mode reads from nNN or xNN according to --source.'
        write(*,'(A)') ''
        write(*,'(A)') 'Outputs: probabilities.dat, thermodynamics.dat, ave_data.dat, ave_spectra.dat'
    end subroutine print_usage_canonical

    ! Normalizes the canonical data source selector.
    ! Normaliza el selector de fuente de datos canónica.
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

    ! Locates the input directory that contains the canonical dataset files.
    ! Localiza el directorio de entrada que contiene los ficheros del conjunto canónico.
    subroutine locate_canonical_dataset(options, data_dir, source_label)
        type(canonical_options_t), intent(in) :: options
        character(len=*), intent(out) :: data_dir
        character(len=*), intent(out) :: source_label
        logical :: has_current
        character(len=32) :: level_dir

        has_current = dataset_exists('.')
        data_dir = ''
        source_label = ''

        if (options%level < 0) then
            if (.not. has_current) then
                write(error_unit,'(A)') 'Error: no OUTSOD/ENERGIES dataset found in the current directory. Use -N to select a level.'
                stop 1
            end if
            data_dir = '.'
            source_label = 'legacy'
            return
        end if

        select case (trim(options%source_mode))
        case ('legacy')
            if (.not. has_current) then
                write(error_unit,'(A)') 'Error: --source legacy requires OUTSOD and ENERGIES in the current directory.'
                stop 1
            end if
            data_dir = '.'
            source_label = 'legacy'
        case ('n')
            level_dir = trim(format_level_directory('n', options%level))
            if (.not. dataset_exists(trim(level_dir))) then
                write(error_unit,'(A,A)') 'Error: canonical dataset not found in ', trim(level_dir)
                stop 1
            end if
            data_dir = trim(level_dir)
            source_label = 'n'
        case ('x')
            level_dir = trim(format_level_directory('x', options%level))
            if (.not. dataset_exists(trim(level_dir))) then
                write(error_unit,'(A,A)') 'Error: canonical dataset not found in ', trim(level_dir)
                stop 1
            end if
            data_dir = trim(level_dir)
            source_label = 'x'
        case default
            level_dir = trim(format_level_directory('n', options%level))
            if (dataset_exists(trim(level_dir))) then
                data_dir = trim(level_dir)
                source_label = 'n'
                return
            end if
            level_dir = trim(format_level_directory('x', options%level))
            if (dataset_exists(trim(level_dir))) then
                data_dir = trim(level_dir)
                source_label = 'x'
                return
            end if
            if (has_current) then
                data_dir = '.'
                source_label = 'legacy'
                return
            end if
            write(error_unit,'(A,I0,A)') 'Error: failed to locate data for level ', options%level, '.'
            stop 1
        end select
    end subroutine locate_canonical_dataset

    ! Returns whether the requested directory contains both OUTSOD and ENERGIES.
    ! Devuelve si el directorio solicitado contiene tanto OUTSOD como ENERGIES.
    logical function dataset_exists(dir)
        character(len=*), intent(in) :: dir
        logical :: has_outsod, has_energies

        inquire(file=trim(join_path(dir, 'OUTSOD')), exist=has_outsod)
        inquire(file=trim(join_path(dir, 'ENERGIES')), exist=has_energies)
        dataset_exists = has_outsod .and. has_energies
    end function dataset_exists

    ! Reads the TEMPERATURES file or falls back to the original default values.
    ! Lee el fichero TEMPERATURES o usa los valores por defecto originales.
    subroutine read_temperatures_file(path, exists_flag, t, ntt)
        character(len=*), intent(in) :: path
        logical, intent(in) :: exists_flag
        real(dp), allocatable, intent(out) :: t(:)
        integer, intent(out) :: ntt
        integer :: unit_in, ios, count
        real(dp) :: value

        if (.not. exists_flag) then
            allocate(t(3))
            t = [1.0_dp, 300.0_dp, 1000.0_dp]
            ntt = 3
            return
        end if

        open(newunit=unit_in, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to open ', trim(path)
            stop 1
        end if

        count = 0
        do
            read(unit_in, *, iostat=ios) value
            if (ios /= 0) exit
            count = count + 1
        end do
        rewind(unit_in)

        if (count <= 0) then
            close(unit_in)
            allocate(t(3))
            t = [1.0_dp, 300.0_dp, 1000.0_dp]
            ntt = 3
            return
        end if

        allocate(t(count))
        do ntt = 1, count
            read(unit_in, *, iostat=ios) t(ntt)
            if (ios /= 0) then
                close(unit_in)
                write(error_unit,'(A,A)') 'Error: invalid entry in ', trim(path)
                stop 1
            end if
        end do
        close(unit_in)
        ntt = count
    end subroutine read_temperatures_file

    ! Reads OUTSOD and extracts the substitution count, the number of configurations, and the degeneracies.
    ! Lee OUTSOD y extrae el número de sustituciones, el número de configuraciones y las degeneraciones.
    subroutine read_outsod_file(path, nsubs, mm, omega)
        character(len=*), intent(in) :: path
        integer, intent(out) :: nsubs, mm
        integer, allocatable, intent(out) :: omega(:)
        integer :: unit_in, ios, auxm, m
        character(len=512) :: line
        integer :: total_sites
        character(len=32) :: word1, word2

        open(newunit=unit_in, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to open ', trim(path)
            stop 1
        end if

        read(unit_in, '(A)', iostat=ios) line
        if (ios /= 0) then
            close(unit_in)
            write(error_unit,'(A,A)') 'Error: failed to read the header of ', trim(path)
            stop 1
        end if
        read(line, *, iostat=ios) nsubs, word1, word2, total_sites
        if (ios /= 0) then
            close(unit_in)
            write(error_unit,'(A,A)') 'Error: invalid OUTSOD header in ', trim(path)
            stop 1
        end if

        read(unit_in, *, iostat=ios) mm
        if (ios /= 0 .or. mm <= 0) then
            close(unit_in)
            write(error_unit,'(A,A)') 'Error: invalid number of configurations in ', trim(path)
            stop 1
        end if

        allocate(omega(mm))
        do m = 1, mm
            read(unit_in, *, iostat=ios) auxm, omega(m)
            if (ios /= 0) then
                close(unit_in)
                write(error_unit,'(A,A)') 'Error: invalid degeneracy entry in ', trim(path)
                stop 1
            end if
        end do
        close(unit_in)
    end subroutine read_outsod_file

    ! Reads ENERGIES, accepting the historical one-column format and newer commented variants.
    ! Lee ENERGIES, aceptando el formato histórico de una columna y variantes modernas con comentarios.
    subroutine read_energies_file(path, mm, ene)
        character(len=*), intent(in) :: path
        integer, intent(in) :: mm
        real(dp), allocatable, intent(out) :: ene(:)
        integer :: unit_in, ios, count, numeric_count
        character(len=1024) :: line
        real(dp) :: numeric_tokens(32)

        allocate(ene(mm))
        open(newunit=unit_in, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to open ', trim(path)
            stop 1
        end if

        count = 0
        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            call extract_numeric_tokens(line, numeric_tokens, numeric_count)
            if (numeric_count <= 0) cycle
            count = count + 1
            if (count > mm) exit
            ene(count) = numeric_tokens(1)
        end do
        close(unit_in)

        if (count < mm) then
            write(error_unit,'(A,A)') 'Error: not enough energy entries in ', trim(path)
            stop 1
        end if
    end subroutine read_energies_file

    ! Reads DATA and returns the ncol x mm table used for thermal averages.
    ! Lee DATA y devuelve la tabla ncol x mm usada para promedios térmicos.
    subroutine read_data_file(path, mm, ncol, data)
        character(len=*), intent(in) :: path
        integer, intent(in) :: mm
        integer, intent(out) :: ncol
        real(dp), allocatable, intent(out) :: data(:,:)
        integer :: unit_in, ios, m

        open(newunit=unit_in, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to open ', trim(path)
            stop 1
        end if
        read(unit_in, *, iostat=ios) ncol
        if (ios /= 0 .or. ncol <= 0) then
            close(unit_in)
            write(error_unit,'(A,A)') 'Error: invalid DATA header in ', trim(path)
            stop 1
        end if
        allocate(data(ncol, mm))
        do m = 1, mm
            read(unit_in, *, iostat=ios) data(:, m)
            if (ios /= 0) then
                close(unit_in)
                write(error_unit,'(A,A)') 'Error: invalid DATA row in ', trim(path)
                stop 1
            end if
        end do
        close(unit_in)
    end subroutine read_data_file

    ! Reads SPECTRA and XSPEC using the original statsod layout.
    ! Lee SPECTRA y XSPEC usando el formato original de statsod.
    subroutine read_spectra_files(spectra_path, xspec_path, mm, npoints, spec, xspec)
        character(len=*), intent(in) :: spectra_path, xspec_path
        integer, intent(in) :: mm
        integer, intent(out) :: npoints
        real(dp), allocatable, intent(out) :: spec(:,:), xspec(:)
        integer :: unit_spec, unit_xspec, ios, m, point

        open(newunit=unit_spec, file=trim(spectra_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to open ', trim(spectra_path)
            stop 1
        end if
        read(unit_spec, *, iostat=ios) npoints
        if (ios /= 0 .or. npoints <= 0) then
            close(unit_spec)
            write(error_unit,'(A,A)') 'Error: invalid SPECTRA header in ', trim(spectra_path)
            stop 1
        end if

        allocate(spec(npoints, mm))
        do m = 1, mm
            read(unit_spec, *, iostat=ios) spec(:, m)
            if (ios /= 0) then
                close(unit_spec)
                write(error_unit,'(A,A)') 'Error: invalid SPECTRA row in ', trim(spectra_path)
                stop 1
            end if
        end do
        close(unit_spec)

        open(newunit=unit_xspec, file=trim(xspec_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to open ', trim(xspec_path)
            stop 1
        end if
        allocate(xspec(npoints))
        do point = 1, npoints
            read(unit_xspec, *, iostat=ios) xspec(point)
            if (ios /= 0) then
                close(unit_xspec)
                write(error_unit,'(A,A)') 'Error: invalid XSPEC entry in ', trim(xspec_path)
                stop 1
            end if
        end do
        close(unit_xspec)
    end subroutine read_spectra_files

    ! Creates probabilities.dat and writes the classical minimum/maximum summary.
    ! Crea probabilities.dat y escribe el resumen clásico de mínimo/máximo.
    subroutine open_probability_output(path, memin, memax)
        character(len=*), intent(in) :: path
        integer, intent(in) :: memin, memax
        integer :: unit_out, ios

        open(newunit=unit_out, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to create ', trim(path)
            stop 1
        end if
        write(unit_out, *) 'Configuration with minimum energy: ', memin
        write(unit_out, *) 'Configuration with maximum energy: ', memax
        write(unit_out, *)
        close(unit_out)
    end subroutine open_probability_output

    ! Creates thermodynamics.dat with the classical statsod header.
    ! Crea thermodynamics.dat con la cabecera clásica de statsod.
    subroutine open_thermodynamics_output(path)
        character(len=*), intent(in) :: path
        integer :: unit_out, ios

        open(newunit=unit_out, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to create ', trim(path)
            stop 1
        end if
        write(unit_out, *) '       T             E               F          S             '
        close(unit_out)
    end subroutine open_thermodynamics_output

    ! Creates ave_data.dat with the original header.
    ! Crea ave_data.dat con la cabecera original.
    subroutine open_ave_data_output(path)
        character(len=*), intent(in) :: path
        integer :: unit_out, ios

        open(newunit=unit_out, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to create ', trim(path)
            stop 1
        end if
        write(unit_out, *) '       T    Average data'
        close(unit_out)
    end subroutine open_ave_data_output

    ! Creates ave_spectra.dat with the same temperature header used by statsod.
    ! Crea ave_spectra.dat con la misma cabecera de temperaturas usada por statsod.
    subroutine open_ave_spectra_output(path, t, ntt)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: t(:)
        integer, intent(in) :: ntt
        integer :: unit_out, ios

        open(newunit=unit_out, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to create ', trim(path)
            stop 1
        end if
        write(unit_out, *) 'x   ', t(1:ntt)
        close(unit_out)
    end subroutine open_ave_spectra_output

    ! Appends the probability table for one temperature to probabilities.dat.
    ! Añade la tabla de probabilidades para una temperatura a probabilities.dat.
    subroutine write_probability_block(path, temp, mm, omega, erel, pcol)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: temp
        integer, intent(in) :: mm, omega(:)
        real(dp), intent(in) :: erel(:), pcol(:)
        integer :: unit_out, ios, m

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        write(unit_out, 100) 'Temperature=', temp
100     format(a12, 2x, f10.4)
        write(unit_out, *) '        m  omega(m)   Erel(m)/eV       p(m)      p(m)/omega(m)'
        do m = 1, mm
            write(unit_out, 101) m, omega(m), erel(m), pcol(m), pcol(m) / real(omega(m), dp)
101         format(i10, 2x, i8, 2x, e12.6, 2(2x, f10.4))
        end do
        write(unit_out, *)
        write(unit_out, *)
        close(unit_out)
    end subroutine write_probability_block

    ! Appends one finite-temperature line to thermodynamics.dat.
    ! Añade una línea de temperatura finita a thermodynamics.dat.
    subroutine append_thermodynamics(path, temp, e, f, s)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: temp, e, f, s
        integer :: unit_out, ios

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        write(unit_out, 200) temp, e, f, s
200     format(f10.1, 2x, 2(f14.4, 2x), e12.6)
        close(unit_out)
    end subroutine append_thermodynamics

    ! Appends one finite-temperature line to ave_data.dat.
    ! Añade una línea de temperatura finita a ave_data.dat.
    subroutine append_ave_data(path, temp, avedata_col)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: temp, avedata_col(:)
        integer :: unit_out, ios, col
        character(len=64) :: fmt

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        write(fmt, '(A,I0,A)') '(F10.1,2X,', size(avedata_col), '(F10.4,2X))'
        write(unit_out, fmt) temp, avedata_col
        close(unit_out)
    end subroutine append_ave_data

    ! Appends the infinite-temperature probability block.
    ! Añade el bloque de probabilidades a temperatura infinita.
    subroutine write_infinite_probability_block(path, mm, omega, erel, pinf)
        character(len=*), intent(in) :: path
        integer, intent(in) :: mm, omega(:)
        real(dp), intent(in) :: erel(:), pinf(:)
        integer :: unit_out, ios, m

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        write(unit_out, *) 'Infinite Temperature Limit'
        write(unit_out, *) '        m  omega(m)   Erel(m)/eV       p(m)      p(m)/omega(m)'
        do m = 1, mm
            write(unit_out, 101) m, omega(m), erel(m), pinf(m), pinf(m) / real(omega(m), dp)
101         format(i10, 2x, i8, 2x, e12.6, 2(2x, f10.4))
        end do
        close(unit_out)
    end subroutine write_infinite_probability_block

    ! Appends the infinite-temperature line to thermodynamics.dat.
    ! Añade la línea de temperatura infinita a thermodynamics.dat.
    subroutine append_infinite_thermodynamics(path, einf, sinf)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: einf, sinf
        integer :: unit_out, ios

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        write(unit_out, 300) 'Infinite', einf, ' - ', sinf
300     format(a10, 2x, f14.4, 8x, a3, 7x, e12.6)
        close(unit_out)
    end subroutine append_infinite_thermodynamics

    ! Appends the infinite-temperature averages to ave_data.dat.
    ! Añade los promedios de temperatura infinita a ave_data.dat.
    subroutine append_infinite_ave_data(path, avedatainf)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: avedatainf(:)
        integer :: unit_out, ios
        character(len=64) :: fmt

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        write(fmt, '(A,I0,A)') '(A10,2X,', size(avedatainf), '(F10.4,2X))'
        write(unit_out, fmt) adjustr('Infinite'), avedatainf
        close(unit_out)
    end subroutine append_infinite_ave_data

    ! Appends the averaged spectra for all temperatures plus the infinite-temperature limit.
    ! Añade los espectros promediados para todas las temperaturas más el límite a temperatura infinita.
    subroutine append_ave_spectra(path, fmtemplist, xspec, avespec, avespecinf)
        character(len=*), intent(in) :: path, fmtemplist
        real(dp), intent(in) :: xspec(:), avespec(:,:), avespecinf(:)
        integer :: unit_out, ios, point

        open(newunit=unit_out, file=trim(path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,A)') 'Error: failed to append to ', trim(path)
            stop 1
        end if
        do point = 1, size(xspec)
            write(unit_out, fmtemplist) xspec(point), avespec(point, :), avespecinf(point)
        end do
        close(unit_out)
    end subroutine append_ave_spectra

    ! Extracts all numeric tokens from a free-format text line.
    ! Extrae todos los tokens numéricos de una línea de texto en formato libre.
    subroutine extract_numeric_tokens(line, values, count)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: count
        integer :: idx, len_line, start_pos, ios
        character(len=64) :: token

        count = 0
        values = 0.0_dp
        len_line = len_trim(line)
        idx = 1

        do while (idx <= len_line)
            do while (idx <= len_line .and. is_separator(line(idx:idx)))
                idx = idx + 1
            end do
            if (idx > len_line) exit
            if (line(idx:idx) == '#') exit

            start_pos = idx
            do while (idx <= len_line .and. .not. is_separator(line(idx:idx)) .and. line(idx:idx) /= '#')
                idx = idx + 1
            end do

            token = ''
            token(1:min(len(token), idx - start_pos)) = line(start_pos:idx-1)
            read(token, *, iostat=ios) values(min(count + 1, size(values)))
            if (ios == 0) then
                count = count + 1
                if (count >= size(values)) exit
            end if
            if (idx <= len_line .and. line(idx:idx) == '#') exit
        end do
    end subroutine extract_numeric_tokens

    ! Returns whether a character acts as a separator in a free-format numeric line.
    ! Devuelve si un carácter actúa como separador en una línea numérica de formato libre.
    logical function is_separator(ch)
        character(len=1), intent(in) :: ch
        is_separator = (ch == ' ' .or. ch == char(9) .or. ch == ',' .or. ch == ';')
    end function is_separator

end module canonical
