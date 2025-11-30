program energy_stats
    use iso_fortran_env, only: dp => real64, error_unit
    implicit none

    integer :: argc, arg_idx, n_paths
    character(len=512) :: arg, value
    character(len=512), allocatable :: paths(:)
    logical :: any_success, use_weighting, expecting_temp
    real(dp) :: temperature

    argc = command_argument_count()
    if (argc == 0) then
        call print_usage()
        stop 1
    end if

    allocate(character(len=512) :: paths(argc))
    n_paths = 0
    temperature = 0.0_dp
    use_weighting = .false.
    expecting_temp = .false.

    do arg_idx = 1, argc
        call get_command_argument(arg_idx, arg)
        arg = adjustl(arg)

        if (expecting_temp) then
            call parse_temperature(trim(arg), temperature)
            use_weighting = .true.
            expecting_temp = .false.
            cycle
        end if

        select case (trim(arg))
        case ('--help', '-h')
            call print_usage()
            stop 0
        case ('--temperature', '-T')
            expecting_temp = .true.
        case default
            if (index(trim(arg), '--temperature=') == 1) then
                value = trim(arg(len('--temperature=')+1:))
                call parse_temperature(value, temperature)
                use_weighting = .true.
            else if (len_trim(arg) > 0 .and. arg(1:1) == '-') then
                write(error_unit,'("Opción no reconocida: ",A)') trim(arg)
                stop 1
            else
                n_paths = n_paths + 1
                paths(n_paths) = trim(arg)
            end if
        end select
    end do

    if (expecting_temp) then
        write(error_unit,'(A)') 'Error: falta el valor de temperatura tras --temperature.'
        stop 1
    end if

    if (use_weighting .and. temperature <= 0.0_dp) then
        write(error_unit,'(A)') 'Error: la temperatura debe ser positiva.'
        stop 1
    end if

    if (n_paths == 0) then
        call print_usage()
        stop 1
    end if

    any_success = .false.
    do arg_idx = 1, n_paths
        call process_path(paths(arg_idx), temperature, use_weighting, any_success)
    end do

    if (allocated(paths)) deallocate(paths)

    if (.not. any_success) then
        write(error_unit,'(A)') 'No se pudieron procesar los archivos solicitados.'
        stop 2
    end if

contains

    subroutine process_path(path, temperature, use_weighting, any_success)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: temperature
        logical, intent(in) :: use_weighting
        logical, intent(inout) :: any_success
        character(len=512) :: candidate
        logical :: exists_file

        call resolve_target(path, candidate, exists_file)
        if (.not. exists_file) then
            write(error_unit,'(A)') 'Aviso: no se encontró ENERGIES en '//trim(path)
            return
        end if
        call analyze_file(candidate, temperature, use_weighting, any_success)
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

    subroutine analyze_file(filename, temperature, use_weighting, any_success)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: temperature
        logical, intent(in) :: use_weighting
        logical, intent(inout) :: any_success
        integer :: unit, ios, ios_read, count, idx
        character(len=512) :: line, buffer
        real(dp) :: value, sum_e, sum_sq, mean, variance, stddev
        real(dp) :: sum_cub, sum_quart, skew, kurt
        real(dp) :: emin, emax
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: kB, ref_energy, w_sum, w_mean, w_variance, w_std, lnz
        real(dp) :: w_m3, w_m4, w_skew, w_kurt

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Aviso: no se pudo abrir '//trim(filename)
            return
        end if

        count = 0
        do
            read(unit,'(A)', iostat=ios) line
            if (ios /= 0) exit
            buffer = line
            call sanitize_line(buffer)
            if (len_trim(buffer) == 0) cycle
            read(buffer,*, iostat=ios_read) value
            if (ios_read /= 0) cycle
            count = count + 1
        end do

        if (count == 0) then
            close(unit)
            write(error_unit,'(A)') 'Aviso: archivo sin valores -> '//trim(filename)
            return
        end if

        allocate(values(count))
        if (use_weighting) allocate(weights(count))

        rewind(unit)
        sum_e = 0.0_dp
        sum_sq = 0.0_dp
        idx = 0
        do
            read(unit,'(A)', iostat=ios) line
            if (ios /= 0) exit
            buffer = line
            call sanitize_line(buffer)
            if (len_trim(buffer) == 0) cycle
            read(buffer,*, iostat=ios_read) value
            if (ios_read /= 0) cycle
            idx = idx + 1
            values(idx) = value
            sum_e = sum_e + value
            sum_sq = sum_sq + value * value
        end do
        close(unit)

        mean = sum_e / real(count, dp)
        variance = max(0.0_dp, (sum_sq / real(count, dp)) - mean * mean)
        stddev = sqrt(variance)
        ! Momentos centrados para skewness y curtosis (equiprobable)
        sum_cub = 0.0_dp
        sum_quart = 0.0_dp
        do idx = 1, count
            sum_cub = sum_cub + (values(idx) - mean)**3
            sum_quart = sum_quart + (values(idx) - mean)**4
        end do
        if (variance > 0.0_dp) then
            skew = sum_cub / real(count, dp) / (stddev**3)
            kurt = sum_quart / real(count, dp) / (stddev**4)
        else
            skew = 0.0_dp
            kurt = 0.0_dp
        end if
        emin = minval(values)
        emax = maxval(values)

        write(*,'(A)') 'Archivo: '//trim(filename)
        write(*,'(A,I0)') '  Entradas: ', count
        write(*,'(A,F18.8)') '  E_min: ', emin
        write(*,'(A,F18.8)') '  E_max: ', emax
        write(*,'(A,F18.8)') '  <E> (media): ', mean
        write(*,'(A,F18.8)') '  Desviacion estandar: ', stddev
        write(*,'(A,F18.8)') '  Asimetria (skewness): ', skew
        write(*,'(A,F18.8)') '  Curtosis: ', kurt

        if (use_weighting) then
            kB = 8.617333262145e-5_dp
            ref_energy = emin
            w_sum = 0.0_dp
            ! Shift energies by the minimum to avoid overflow when computing weights.
            do idx = 1, count
                weights(idx) = exp(-(values(idx) - ref_energy) / (kB * temperature))
                w_sum = w_sum + weights(idx)
            end do

            if (w_sum <= 0.0_dp) then
                write(error_unit,'(A)') 'Aviso: los pesos fueron nulos en '//trim(filename)
            else
                w_mean = sum(weights * values) / w_sum
                w_variance = sum(weights * (values - w_mean)**2) / w_sum
                w_variance = max(0.0_dp, w_variance)
                w_std = sqrt(w_variance)
                ! Momentos ponderados para skewness y curtosis Boltzmann
                w_m3 = sum(weights * (values - w_mean)**3) / w_sum
                w_m4 = sum(weights * (values - w_mean)**4) / w_sum
                if (w_variance > 0.0_dp) then
                    w_skew = w_m3 / (w_std**3)
                    w_kurt = w_m4 / (w_std**4)
                else
                    w_skew = 0.0_dp
                    w_kurt = 0.0_dp
                end if
                lnz = log(w_sum) - ref_energy / (kB * temperature)
                write(*,'(A,F12.2)') '  Temperatura (K): ', temperature
                write(*,'(A,F18.8)') '  <E>_Boltzmann: ', w_mean
                write(*,'(A,F18.8)') '  Sigma_Boltzmann: ', w_std
                write(*,'(A,F18.8)') '  Asimetria_Boltzmann: ', w_skew
                write(*,'(A,F18.8)') '  Curtosis_Boltzmann: ', w_kurt
                write(*,'(A,F18.8)') '  ln(Z): ', lnz
            end if
        end if

        any_success = .true.

        if (allocated(weights)) deallocate(weights)
        if (allocated(values)) deallocate(values)
    end subroutine analyze_file

    subroutine sanitize_line(text)
        character(len=*), intent(inout) :: text
        integer :: pos
        pos = index(text, '#')
        if (pos > 0) text(pos:) = ' '
    end subroutine sanitize_line

    subroutine parse_temperature(text, temperature)
        character(len=*), intent(in) :: text
        real(dp), intent(out) :: temperature
        integer :: ios

        read(text, *, iostat=ios) temperature
        if (ios /= 0 .or. temperature <= 0.0_dp) then
            write(error_unit,'(A)') 'Error: la temperatura debe ser un número positivo.'
            stop 1
        end if
    end subroutine parse_temperature

    subroutine print_usage()
        write(error_unit,'(A)') 'Uso: energy_stats [--temperature <K>] ruta1 [ruta2 ...]'
        write(error_unit,'(A)') '  Cada ruta puede ser un fichero ENERGIES o un directorio que lo contenga.'
        write(error_unit,'(A)') '  Con --temperature se aplican pesos de Boltzmann a las energías.'
    end subroutine print_usage

end program energy_stats
