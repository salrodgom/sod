program energy_stats
    use iso_fortran_env, only: dp => real64, error_unit
    implicit none

    integer :: argc, iarg
    character(len=512) :: arg
    logical :: any_success

    argc = command_argument_count()
    if (argc == 0) then
        write(error_unit,'(A)') 'Uso: energy_stats ruta1 [ruta2 ...]'
        write(error_unit,'(A)') '  Cada ruta puede ser un fichero ENERGIES o un directorio que lo contenga.'
        stop 1
    end if

    any_success = .false.
    do iarg = 1, argc
        call get_command_argument(iarg, arg)
        call process_path(trim(arg), any_success)
    end do

    if (.not. any_success) then
        write(error_unit,'(A)') 'No se pudieron procesar los archivos solicitados.'
        stop 2
    end if

contains

    subroutine process_path(path, any_success)
        character(len=*), intent(in) :: path
        logical, intent(inout) :: any_success
        character(len=512) :: candidate
        logical :: exists_file

        call resolve_target(path, candidate, exists_file)
        if (.not. exists_file) then
            write(error_unit,'(A)') 'Aviso: no se encontró ENERGIES en '//trim(path)
            return
        end if
        call analyze_file(candidate, any_success)
    end subroutine process_path

    subroutine resolve_target(path, out_path, exists_file)
        character(len=*), intent(in) :: path
        character(len=*), intent(out) :: out_path
        logical, intent(out) :: exists_file
        character(len=512) :: candidate
        character(len=512) :: trimmed
        logical :: is_dir, exists_raw
        integer :: last

        out_path = ''
        exists_file = .false.
        is_dir = .false.
        exists_raw = .false.
        trimmed = trim(path)

        if (len_trim(trimmed) == 0) return

        inquire(file=trimmed, exist=exists_raw)
        if (exists_raw) then
            last = len_trim(trimmed)
            if (trimmed(last:last) == '/') then
                is_dir = .true.
            else
                inquire(file=trimmed//'/.', exist=is_dir)
            end if
        else
            inquire(file=trimmed//'/.', exist=is_dir)
        end if

        if (is_dir) then
            candidate = trimmed
            last = len_trim(candidate)
            if (candidate(last:last) /= '/') candidate = candidate(1:last)//'/'
            candidate = trim(candidate)//'ENERGIES'
        else
            candidate = trimmed
        end if

        inquire(file=candidate, exist=exists_file)
        if (exists_file) out_path = candidate
    end subroutine resolve_target

    subroutine analyze_file(filename, any_success)
        character(len=*), intent(in) :: filename
        logical, intent(inout) :: any_success
        integer :: unit, ios, ios_read
        character(len=512) :: line, buffer
        real(dp) :: value
        real(dp) :: sum_e, sum_sq, mean, variance, stddev
        integer :: count

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Aviso: no se pudo abrir '//trim(filename)
            return
        end if

        sum_e = 0.0_dp
        sum_sq = 0.0_dp
        count = 0

        do
            read(unit,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            buffer = line
            call sanitize_line(buffer)
            read(buffer,*, iostat=ios_read) value
            if (ios_read /= 0) cycle
            sum_e = sum_e + value
            sum_sq = sum_sq + value * value
            count = count + 1
        end do
        close(unit)

        if (count == 0) then
            write(error_unit,'(A)') 'Aviso: archivo sin valores -> '//trim(filename)
            return
        end if

        mean = sum_e / real(count, dp)
        variance = max(0.0_dp, (sum_sq / real(count, dp)) - mean * mean)
        stddev = sqrt(variance)

        write(*,'(A)') 'Archivo: '//trim(filename)
        write(*,'(A,I0)') '  Entradas: ', count
        write(*,'(A,F18.8)') '  <E> (media): ', mean
    write(*,'(A,F18.8)') '  Desviacion estandar: ', stddev
        any_success = .true.
    end subroutine analyze_file

    subroutine sanitize_line(text)
        character(len=*), intent(inout) :: text
        integer :: pos
        pos = index(text, '#')
        if (pos > 0) text(pos:) = ' '
    end subroutine sanitize_line

end program energy_stats
