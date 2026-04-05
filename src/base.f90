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
! Shared ensemble Monte Carlo utilities: common constants and helper routines.
! Utilidades compartidas del ensemble Monte Carlo: constantes comunes y rutinas auxiliares.
!*******************************************************************************

! Module `consts` defines shared numerical kinds and physical constants.
! El módulo `consts` define los tipos numéricos y las constantes físicas compartidas.
module consts
    use iso_fortran_env, only: real64, int64, error_unit
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: ip = int64
    real(dp), parameter :: kB_eVk = 8.617333262145d-5
end module consts

! Module `utils` provides general helper routines reused across the code base.
! El módulo `utils` proporciona rutinas auxiliares generales reutilizadas en todo el código.
module utils
    use consts
    implicit none
contains

    ! Prints the standard SOD banner shared by the active executable modes.
    ! Imprime la cabecera estándar de SOD compartida por los modos activos del ejecutable.
    subroutine print_sod_banner()
        write(*,*) '**************************************************************************** '
        write(*,*) '         SOD (Site Occupancy Disorder) version X.YY'
        write(*,*) ' '
        write(*,*) '         Authors: R. Grau-Crespo and S. Hamad'
        write(*,*) '       (modified version by: S. R. G. Balestra)'
        write(*,*) '        Contact:  <r.grau-crespo@reading.ac.uk>'
        write(*,*) ' '
        write(*,*) '**************************************************************************** '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
    end subroutine print_sod_banner

    pure function format_level_directory(prefix, level) result(dir)
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: level
        character(len=32) :: dir

        write(dir,'(A,"0",I0)') trim(prefix), level
    end function format_level_directory

    pure function join_path(dir, filename) result(path)
        character(len=*), intent(in) :: dir, filename
        character(len=512) :: path

        if (len_trim(dir) == 0) then
            path = trim(filename)
        else
            path = trim(dir)//'/'//trim(filename)
        end if
    end function join_path

    subroutine ensure_directory_exists(dir)
        character(len=*), intent(in) :: dir
        integer :: exit_code

        if (len_trim(dir) == 0) return
        call execute_command_line('mkdir -p ' // trim(dir), exitstat=exit_code)
        if (exit_code /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: failed to create directory', trim(dir)
            stop 1
        end if
    end subroutine ensure_directory_exists

    ! Returns "n choose k" using a cumulative product in 64-bit integers.
! Devuelve "n sobre k" usando un producto acumulado con enteros de 64 bits.
    function binomial_int64(n, k) result(val)
        integer, intent(in) :: n, k
        integer(ip) :: val
        integer :: i

        if (k < 0 .or. k > n) then
            val = 0_ip
            return
        end if
        if (k == 0 .or. k == n) then
            val = 1_ip
            return
        end if
        val = 1_ip
        do i = 1, k
            val = val * int(n - k + i, ip) / int(i, ip)
        end do
    end function binomial_int64

    ! Advances an ordered combination to the next lexicographic tuple.
! Avanza una combinación ordenada a la siguiente tupla lexicográfica.
    logical function next_combination(comb, n)
        integer, intent(inout) :: comb(:)
        integer, intent(in) :: n
        integer :: k, i, j

        k = size(comb)
        if (k == 0) then
            next_combination = .false.
            return
        end if

        i = k
        do while (i >= 1 .and. comb(i) == n - k + i)
            i = i - 1
        end do
        if (i == 0) then
            next_combination = .false.
            return
        end if

        comb(i) = comb(i) + 1
        do j = i + 1, k
            comb(j) = comb(j-1) + 1
        end do
        next_combination = .true.
    end function next_combination

    ! Draws an ordered random subset of size k from 1..n without replacement.
! Extrae un subconjunto aleatorio ordenado de tamaño k de 1..n sin reemplazo.
    subroutine random_subset(n, k, subset)
        integer, intent(in) :: n, k
        integer, intent(out) :: subset(:)
        integer :: pool(n)
        integer :: i, j, tmp
        real(dp) :: r

        if (k == 0) then
            return
        end if

        if (size(subset) < k) then
            write(error_unit,'(A)') 'Internal error: insufficient buffer size in random_subset.'
            stop 1
        end if

        do i = 1, n
            pool(i) = i
        end do

        do i = 1, k
            call random_number(r)
            j = i + int(r * real(n - i + 1, dp))
            if (j < i) j = i
            if (j > n) j = n
            tmp = pool(i)
            pool(i) = pool(j)
            pool(j) = tmp
        end do

        subset(1:k) = pool(1:k)

        call sort_int_ascending(subset, k)
    end subroutine random_subset

    ! Sorts the first "length" elements of an integer array using insertion sort.
! Ordena los primeros "length" elementos de un vector entero usando ordenación por inserción.
    subroutine sort_int_ascending(arr, length)
        integer, intent(inout) :: arr(:)
        integer, intent(in) :: length
        integer :: i, j, key

        do i = 2, length
            key = arr(i)
            j = i - 1
            do while (j >= 1 .and. arr(j) > key)
                arr(j+1) = arr(j)
                j = j - 1
            end do
            arr(j+1) = key
        end do
    end subroutine sort_int_ascending

end module utils

! Module `cli` provides command-line parsing and normalization helpers.
! El módulo `cli` proporciona utilidades para analizar y normalizar la línea de comandos.
module cli
    use consts
    implicit none
    private

    public :: append_int_value, basename_from_path, build_level_sequence, compose_mode_command
    public :: classify_template_file, engine_name_to_filer
    public :: is_help_token, parse_level_spec, sanitize_level_list, to_lower_inplace

contains

    subroutine to_lower_inplace(text)
        character(len=*), intent(inout) :: text
        integer :: idx, code

        do idx = 1, len(text)
            code = iachar(text(idx:idx))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                text(idx:idx) = achar(code + 32)
            end if
        end do
    end subroutine to_lower_inplace

    logical function is_help_token(raw)
        character(len=*), intent(in) :: raw
        character(len=len(raw)) :: token

        token = adjustl(raw)
        call to_lower_inplace(token)

        select case (trim(token))
        case ('--help','-h','help','--ayuda','-ayuda','ayuda','/h','/?')
            is_help_token = .true.
        case default
            is_help_token = .false.
        end select
    end function is_help_token

    subroutine append_int_value(vec, count, value)
        integer, allocatable, intent(inout) :: vec(:)
        integer, intent(inout) :: count
        integer, intent(in) :: value
        integer, allocatable :: tmp(:)

        count = count + 1
        allocate(tmp(count))
        if (count > 1 .and. allocated(vec)) tmp(1:count-1) = vec(1:count-1)
        tmp(count) = value
        if (allocated(vec)) deallocate(vec)
        call move_alloc(tmp, vec)
    end subroutine append_int_value

    subroutine parse_level_spec(spec, level_min, level_max, level_list, has_list)
        character(len=*), intent(in) :: spec
        integer, intent(inout) :: level_min, level_max
        integer, allocatable, intent(inout) :: level_list(:)
        logical, intent(inout) :: has_list
        integer :: colon_pos, comma_pos
        integer :: value, ios
        character(len=256) :: item
        integer, allocatable :: temp(:)
        integer :: count, start_idx, delim_pos

        colon_pos = index(spec, ':')
        comma_pos = index(spec, ',')

        if (colon_pos > 0 .and. comma_pos > 0) then
            write(error_unit,'(A)') 'Error: do not mix ranges and lists in the same -N specification.'
            stop 1
        end if

        if (comma_pos > 0) then
            has_list = .true.
            count = 0
            allocate(temp(0))
            start_idx = 1
            do
                delim_pos = index(spec(start_idx:), ',')
                if (delim_pos == 0) then
                    item = adjustl(spec(start_idx:))
                else
                    item = adjustl(spec(start_idx:start_idx + delim_pos - 2))
                end if
                read(item, *, iostat=ios) value
                if (ios /= 0) then
                    write(error_unit,'(A)') 'Error: invalid value in the -N list.'
                    stop 1
                end if
                call append_int_value(temp, count, value)
                if (delim_pos == 0) exit
                start_idx = start_idx + delim_pos
            end do
            if (allocated(level_list)) deallocate(level_list)
            call move_alloc(temp, level_list)
        else if (colon_pos > 0) then
            has_list = .false.
            read(spec(1:colon_pos-1), *, iostat=ios) level_min
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid lower bound in -N.'
                stop 1
            end if
            read(spec(colon_pos+1:), *, iostat=ios) level_max
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid upper bound in -N.'
                stop 1
            end if
        else
            read(spec, *, iostat=ios) value
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid -N specification.'
                stop 1
            end if
            if (value < 0) then
                level_min = 0
                level_max = -1
                has_list = .false.
            else
                level_min = value
                level_max = value
                has_list = .false.
            end if
        end if
    end subroutine parse_level_spec

    subroutine sanitize_level_list(level_list, total_sites)
        integer, allocatable, intent(inout) :: level_list(:)
        integer, intent(in) :: total_sites
        integer, allocatable :: filtered(:)
        integer :: count, idx

        count = 0
        allocate(filtered(0))

        if (allocated(level_list)) then
            do idx = 1, size(level_list)
                if (level_list(idx) < 0 .or. level_list(idx) > total_sites) cycle
                if (count > 0) then
                    if (any(filtered(1:count) == level_list(idx))) cycle
                end if
                call append_int_value(filtered, count, level_list(idx))
            end do
            deallocate(level_list)
        end if

        call move_alloc(filtered, level_list)
    end subroutine sanitize_level_list

    subroutine build_level_sequence(level_min, level_max, level_list, has_list, total_sites, levels)
        integer, intent(in) :: level_min, level_max, total_sites
        integer, allocatable, intent(inout) :: level_list(:)
        logical, intent(in) :: has_list
        integer, allocatable, intent(out) :: levels(:)
        integer :: min_level, max_level, idx

        if (has_list) then
            call sanitize_level_list(level_list, total_sites)
            if (.not. allocated(level_list)) then
                allocate(levels(0))
            else
                allocate(levels(size(level_list)))
                levels = level_list
            end if
            return
        end if

        min_level = level_min
        max_level = level_max
        if (max_level < 0) max_level = total_sites
        if (min_level < 0) min_level = 0
        min_level = min(min_level, total_sites)
        max_level = min(max_level, total_sites)

        if (min_level > max_level) then
            allocate(levels(0))
            return
        end if

        allocate(levels(max_level - min_level + 1))
        do idx = 1, size(levels)
            levels(idx) = min_level + idx - 1
        end do
    end subroutine build_level_sequence

    function basename_from_path(path) result(name)
        character(len=*), intent(in) :: path
        character(len=256) :: name
        integer :: idx

        name = adjustl(path)
        do idx = len_trim(name), 1, -1
            if (name(idx:idx) == '/') then
                name = name(idx+1:)
                exit
            end if
        end do
        name = trim(adjustl(name))
    end function basename_from_path

    function compose_mode_command(mode_name) result(command)
        character(len=*), intent(in) :: mode_name
        character(len=256) :: command
        character(len=256) :: exe_name
        character(len=512) :: raw_exe

        raw_exe = ''
        call get_command_argument(0, raw_exe)
        exe_name = basename_from_path(trim(raw_exe))
        if (len_trim(exe_name) == 0) exe_name = 'sod'

        if (trim(exe_name) == 'sod' .or. trim(exe_name) == 'sod_ensemble') then
            command = trim(exe_name) // ' ' // trim(mode_name)
        else
            command = trim(exe_name)
        end if
    end function compose_mode_command

    ! Maps a human-readable engine name to the INSOD FILER integer.
    ! Devuelve el código FILER a partir del nombre legible del motor.
    integer function engine_name_to_filer(name) result(filer)
        character(len=*), intent(in) :: name
        character(len=32) :: low
        low = adjustl(name)
        call to_lower_inplace(low)
        select case (trim(low))
        case ('gulp')
            filer = 1
        case ('lammps')
            filer = 2
        case ('vasp')
            filer = 11
        case ('castep')
            filer = 12
        case ('qe', 'quantum-espresso')
            filer = 13
        case ('ase', 'mace', 'chgnet')
            filer = 14
        case default
            filer = -1
        end select
    end function engine_name_to_filer

    ! Classifies a --template file by extension into a category string.
    ! Clasifica un archivo --template por extensión en una categoría.
    ! Returns: 'lib', 'gin', 'lammps', 'py', or 'unknown'.
    function classify_template_file(path) result(category)
        character(len=*), intent(in) :: path
        character(len=8) :: category
        character(len=512) :: low
        integer :: dot_pos

        low = adjustl(path)
        call to_lower_inplace(low)
        dot_pos = index(trim(low), '.', back=.true.)

        if (dot_pos <= 0) then
            category = 'unknown'
            return
        end if

        select case (trim(low(dot_pos:)))
        case ('.lib')
            category = 'lib'
        case ('.gin', '.include')
            category = 'gin'
        case ('.lammps')
            category = 'lammps'
        case ('.py')
            category = 'py'
        case default
            category = 'unknown'
        end select
    end function classify_template_file

end module cli
