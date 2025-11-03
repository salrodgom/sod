!*******************************************************************************
!    Copyright (c) 2025, Salvador R.G. Balestra
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!******************************************************************************

module sod_boltzmann_consts
    use iso_fortran_env, only: real64, int64, error_unit
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: ip = int64
    real(dp), parameter :: kB_eVk = 8.617333262145d-5
end module sod_boltzmann_consts

module sod_boltzmann_utils
    use sod_boltzmann_consts
    implicit none
contains

    ! Returns n choose k using a running product in 64-bit integers.
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

    ! Advances a sorted combination array to the next lexicographic tuple.
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

    ! Draws a random sorted subset of size k from 1..n without replacement.
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
            write(error_unit,'(A)') 'ERROR interno: buffer insuficiente en random_subset'
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

    ! Sorts the first length entries of an integer array via insertion sort.
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

end module sod_boltzmann_utils

program sod_boltzmann_mc
    use sod_boltzmann_consts
    use sod_boltzmann_utils
    use energy_calc
    implicit none

    integer, parameter :: max_exact_combos = 200000
    integer, parameter :: mix_n = 6
    integer, parameter :: mix_m = 12
    real(dp), parameter :: mix_x0 = 0.5_dp
    real(dp), parameter :: mix_d0 = 0.01_dp
    real(dp), parameter :: conv_ev_to_kjmol = 96.48533212331_dp
    real(dp), parameter :: quartz_ge_energy = -121.812200998519_dp
    real(dp), parameter :: quartz_si_energy = -128.799143746482_dp
    character(len=*), parameter :: summary_filename = 'sod_boltzmann_summary.csv'
    character(len=*), parameter :: summary_txt_filename = 'sod_boltzmann_summary.txt'

    real(dp) :: temperature
    integer :: max_substitutions, samples_per_level
    integer :: seed_value
    character(len=16) :: sampling_mode
    integer :: summary_unit, summary_txt_unit
    logical :: use_parallel
    logical :: omp_available
    logical :: force_restart_accept

    integer, allocatable :: eqmatrix(:,:), config(:), local_config(:)
    integer :: nop, total_sites
    integer :: level, effective_max

    use_parallel = .false.
    omp_available = .false.
!$  use_parallel = .true.
!$  omp_available = .true.

    call parse_arguments(temperature, max_substitutions, samples_per_level, seed_value, sampling_mode, use_parallel, omp_available)
    call configure_random_seed(seed_value)
    call configure_restart_mode(force_restart_accept)

    call init_energy_calc()
    call get_eqmatrix(eqmatrix, nop, total_sites)
    if (.not. allocated(eqmatrix) .or. nop <= 0 .or. total_sites <= 0) then
        write(*,'(A)') 'Error: no se pudo obtener EQMATRIX o el número de posiciones/operadores es inválido.'
        stop 1
    end if
    if (allocated(eqmatrix)) deallocate(eqmatrix)

    call init_summary_files(summary_unit, summary_txt_unit)

    if (max_substitutions < 0) then
        effective_max = total_sites
    else
        effective_max = min(max_substitutions, total_sites)
    end if

    write(*,'(A)') '--- Parámetros del cálculo ---'
    write(*,'(A,F10.2)') 'Temperatura (K): ', temperature
    write(*,'(A,I6)') 'Sitios sustituibles (npos): ', total_sites
    write(*,'(A,I6)') 'Max sustituciones evaluadas: ', effective_max
    write(*,'(A,I8)') 'Umbral enumeración exacta: ', max_exact_combos
    write(*,'(A,I8)') 'Muestras aleatorias (si se supera el umbral): ', samples_per_level
    write(*,'(A)') 'Resultados por nivel guardados en: '//trim(summary_filename)
    write(*,'(A)') '                               y: '//trim(summary_txt_filename)
    write(*,'(A)') 'Método de muestreo MC (si aplica): '//trim(sampling_mode)
    write(*,'(A)') 'OpenMP paralelo: '//merge('Si','No',use_parallel)
    if (use_parallel) then
        write(*,'(A)') 'Nota: las salidas por nivel pueden imprimirse en orden no secuencial durante el cálculo paralelo.'
        write(*,'(A)') '      Los archivos de resumen se reordenan al finalizar para dejar los niveles crecientes.'
    end if
    write(*,*)

    if (use_parallel) then
!$omp parallel default(shared) private(local_config)
        allocate(local_config(total_sites))
!$omp do schedule(dynamic)
        do level = 0, effective_max
            call process_level(level, total_sites, local_config, temperature, samples_per_level, &
                max_exact_combos, sampling_mode, use_parallel, summary_unit, summary_txt_unit)
        end do
!$omp end do
        deallocate(local_config)
!$omp end parallel
    else
        allocate(config(total_sites))
        do level = 0, effective_max
            call process_level(level, total_sites, config, temperature, samples_per_level, &
                max_exact_combos, sampling_mode, use_parallel, summary_unit, summary_txt_unit)
        end do
        deallocate(config)
    end if

    call close_summary_files(summary_unit, summary_txt_unit)
    call reorder_summary_outputs()
    call cleanup_energy_calc()

contains

    ! Parses optional command-line arguments and populates runtime parameters.
    subroutine parse_arguments(temp, max_subs, samples_level, seed, sampler, use_parallel, omp_available)
        real(dp), intent(out) :: temp
        integer, intent(out) :: max_subs, samples_level, seed
        character(len=*), intent(out) :: sampler
        logical, intent(inout) :: use_parallel
        logical, intent(in) :: omp_available
        integer :: argc, ios
        character(len=256) :: carg

        temp = 1000.0_dp
        max_subs = -1
        samples_level = 5000
        seed = -1
    sampler = 'uniform'

        argc = command_argument_count()
        if (argc >= 1) then
            call get_command_argument(1, carg)
            if (is_help_argument(carg)) then
                call print_usage(omp_available)
                stop 0
            end if
            read(carg,*,iostat=ios) temp
            if (ios /= 0) then
                write(error_unit,'(A)') 'Advertencia: argumento 1 inválido, se usa T=1000K'
                temp = 1000.0_dp
            end if
        end if
        if (argc >= 2) then
            call get_command_argument(2, carg)
            if (is_help_argument(carg)) then
                call print_usage(omp_available)
                stop 0
            end if
            read(carg,*,iostat=ios) max_subs
            if (ios /= 0) then
                write(error_unit,'(A)') 'Advertencia: argumento 2 inválido, se evalúan todos los valores'
                max_subs = -1
            end if
        end if
        if (argc >= 3) then
            call get_command_argument(3, carg)
            if (is_help_argument(carg)) then
                call print_usage(omp_available)
                stop 0
            end if
            read(carg,*,iostat=ios) samples_level
            if (ios /= 0 .or. samples_level <= 0) then
                write(error_unit,'(A)') 'Advertencia: argumento 3 inválido, se usan 5000 muestras'
                samples_level = 5000
            end if
        end if
        if (argc >= 4) then
            call get_command_argument(4, carg)
            if (is_help_argument(carg)) then
                call print_usage(omp_available)
                stop 0
            end if
            read(carg,*,iostat=ios) seed
            if (ios /= 0) seed = -1
        end if
        if (argc >= 5) then
            call get_command_argument(5, sampler)
            if (is_help_argument(sampler)) then
                call print_usage(omp_available)
                stop 0
            end if
            sampler = adjustl(sampler)
            if (sampler /= 'uniform' .and. sampler /= 'metropolis') then
                write(error_unit,'(A)') 'Advertencia: método MC desconocido, se usa "uniform"'
                sampler = 'uniform'
            end if
        end if
        if (argc >= 6) then
            call get_command_argument(6, carg)
            if (is_help_argument(carg)) then
                call print_usage(omp_available)
                stop 0
            end if
            carg = adjustl(carg)
            select case (trim(carg))
            case ('omp', 'OMP', 'O', '1', 'true', 'TRUE')
                if (omp_available) then
                    use_parallel = .true.
                else
                    write(error_unit,'(A)') 'Advertencia: OpenMP no disponible en esta compilacion, se ignora argumento 6.'
                    use_parallel = .false.
                end if
            case ('noomp', 'NOOMP', 'no', 'NO', '0', 'false', 'FALSE')
                use_parallel = .false.
            case default
                write(error_unit,'(A)') 'Advertencia: argumento 6 desconocido; opciones validas: omp/noomp.'
            end select
        end if
    end subroutine parse_arguments

    ! Emits program usage information including optional OpenMP note.
    subroutine print_usage(omp_available)
        logical, intent(in) :: omp_available

        write(*,'(A)') 'Uso: sod_boltzmann_mc [T_K] [Nmax] [Nsamples] [seed] [sampler] [omp|noomp]'
        write(*,'(A)') '       sod_boltzmann_mc --help'
        write(*,'(A)') ''
    write(*,'(A)') 'Argumentos opcionales (por defecto entre corchetes):'
    write(*,'(A)') '  T_K        Temperatura en Kelvin para los pesos de Boltzmann [1000].'
    write(*,'(A)') '  Nmax       Número máximo de sustituciones evaluadas por nivel [-1 -> todos].'
    write(*,'(A)') '  Nsamples   Muestras MC por nivel cuando el número de combinaciones'// &
             ' supera el umbral [5000].'
    write(*,'(A)') '  seed       Semilla del generador aleatorio [-1 -> semilla desde system_clock].'
    write(*,'(A)') '  sampler    "uniform" para muestreo sin sesgo, "metropolis" para Metropolis-Hastings ['// &
             'uniform].'
        if (omp_available) then
            write(*,'(A)') '  omp|noomp  Control del uso de OpenMP si el binario lo soporta [noomp].'
        else
            write(*,'(A)') '  omp|noomp  Control del uso de OpenMP (esta compilación no lo habilita).'
        end if
        write(*,'(A)') ''
        write(*,'(A)') 'Otros detalles:'
        write(*,'(A)') '  - El programa evalúa todos los niveles de sustitución desde N=0 hasta Nmax.'
        write(*,'(A)') '  - Si C(N,npos) <= 200000 se enumeran todas las configuraciones; en caso contrario'// &
                     ' se toma un muestreo Monte Carlo.'
        write(*,'(A)') '  - Los resultados agregados de cada nivel se guardan en '//trim(summary_filename)//'.'
        write(*,'(A)') '  - El resumen en texto plano se escribe en '//trim(summary_txt_filename)//'.'
        write(*,'(A)') ''
        write(*,'(A)') 'Ejemplos:'
        write(*,'(A)') '  sod_boltzmann_mc                     # Ejecuta con valores por defecto.'
        write(*,'(A)') '  sod_boltzmann_mc 800 6 2000 1234 metropolis omp'
        write(*,'(A)') '                                        # 800 K, hasta 6 sustituciones,'// &
                     ' Metropolis y OpenMP.'
    end subroutine print_usage

    ! Checks whether a token corresponds to any supported help flag.
    logical function is_help_argument(arg)
        character(len=*), intent(in) :: arg
        character(len=len(arg)) :: token
        integer :: lt, i

        token = adjustl(arg)
        lt = len_trim(token)
        if (lt == 0) then
            is_help_argument = .false.
            return
        end if
        do i = 1, lt
            if (token(i:i) >= 'A' .and. token(i:i) <= 'Z') then
                token(i:i) = achar(iachar(token(i:i)) + 32)
            end if
        end do
        token = token(1:lt)
        if (token == '--help' .or. token == '-h' .or. token == 'help' .or. token == '--ayuda' .or. &
            token == '-ayuda' .or. token == 'ayuda' .or. token == '/h' .or. token == '/?') then
            is_help_argument = .true.
        else
            is_help_argument = .false.
        end if
    end function is_help_argument

    ! Initializes the intrinsic random generator with a deterministic or clock seed.
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
    subroutine configure_restart_mode(force_restart)
        logical, intent(out) :: force_restart
        character(len=32) :: env_value
        integer :: status, len_env, i
        character(len=32) :: token

        force_restart = .true.
        call get_environment_variable('SOD_FORCE_RESTART_ACCEPT', env_value, length=len_env, status=status)
        if (status /= 0 .or. len_env <= 0) then
            write(*,'(A)') 'Modo reinicio forzado: activo (por defecto)'
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
            write(*,'(A)') 'Modo reinicio forzado: desactivado (SOD_FORCE_RESTART_ACCEPT='//trim(token)//')'
        case ('1', 'yes', 'true', 'on', 'enable', 'enabled')
            force_restart = .true.
            write(*,'(A)') 'Modo reinicio forzado: activo (SOD_FORCE_RESTART_ACCEPT='//trim(token)//')'
        case default
            force_restart = .true.
            write(*,'(A)') 'Modo reinicio forzado: activo (valor desconocido, se usa por defecto)'
        end select
    end subroutine configure_restart_mode

    ! Opens the summary CSV and text files and writes their headers.
    subroutine init_summary_files(unit_csv, unit_txt)
        integer, intent(out) :: unit_csv, unit_txt
        integer :: ios

        open(newunit=unit_csv, file=summary_filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: no se pudo crear el archivo de resumen '//trim(summary_filename)
            stop 1
        end if
        open(newunit=unit_txt, file=summary_txt_filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: no se pudo crear el archivo de resumen '//trim(summary_txt_filename)
            stop 1
        end if
        write(unit_csv,'(A)') '#N;FracGe;E_exp_total;E_min_total;Var_total;E_exp_ladoSi;E_min_ladoSi;'// &
            'E_exp_ladoGe;E_min_ladoGe;E_exp_combinada;Delta_exp_total;Delta_min_total;'// &
            'Delta_exp_ladoSi;Delta_min_ladoSi;Delta_exp_ladoGe;Delta_min_ladoGe;Delta_exp_combinada;'// &
            'Ratio_aceptacion'
        write(unit_txt,'(A)') '#N FracGe E_exp_total E_min_total Var_total E_exp_ladoSi E_min_ladoSi '// &
            'E_exp_ladoGe E_min_ladoGe E_exp_combinada Delta_exp_total Delta_min_total '// &
            'Delta_exp_ladoSi Delta_min_ladoSi Delta_exp_ladoGe Delta_min_ladoGe Delta_exp_combinada '// &
            'Ratio_aceptacion'
    end subroutine init_summary_files

    ! Closes the summary file units if they are valid.
    subroutine close_summary_files(unit_csv, unit_txt)
        integer, intent(in) :: unit_csv, unit_txt

        if (unit_csv /= 0) close(unit_csv)
        if (unit_txt /= 0) close(unit_txt)
    end subroutine close_summary_files

    ! Chooses exhaustive or stochastic evaluation for a substitution level and dispatches it.
    subroutine process_level(level, total_sites, config, temperature, samples_level, max_exact, sampler, &
                             use_parallel, summary_unit, summary_txt_unit)
        integer, intent(in) :: level, total_sites, samples_level, max_exact
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: temperature
        character(len=*), intent(in) :: sampler
        logical, intent(in) :: use_parallel
        integer, intent(in) :: summary_unit, summary_txt_unit

        integer(ip) :: total_comb
        integer :: ncomb_int
        logical :: use_sampling

        if (level > total_sites) return

        total_comb = binomial_int64(total_sites, level)
        if (total_comb == 0_ip) then
            write(*,'(A,I3)') 'Nivel ', level, ': sin combinaciones válidas.'
            return
        end if

        use_sampling = (total_comb > max_exact)

        if (use_sampling) then
            if (sampler == 'metropolis') then
                write(*,'(A,I0,A)') 'Nivel ', level, ': se usa muestreo Monte Carlo tipo Metropolis-Hastings.'
                call metropolis_level(level, total_sites, config, temperature, samples_level, total_comb, &
                                      use_parallel, summary_unit, summary_txt_unit)
            else
                write(*,'(A,I0,A)') 'Nivel ', level, ': se usa muestreo Monte Carlo uniforme.'
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
    subroutine exhaustive_level(level, total_sites, config, temperature, total_comb, ncomb_int, use_parallel, &
                                summary_unit, summary_txt_unit)
        integer, intent(in) :: level, total_sites, ncomb_int
        integer(ip), intent(in) :: total_comb
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: temperature
        logical, intent(in) :: use_parallel
        integer, intent(in) :: summary_unit, summary_txt_unit

        real(dp), allocatable :: energies(:), energies_low(:), energies_high(:)
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
            write(*,'(A,I3)') 'Nivel ', level, ': no se pudieron generar configuraciones válidas.'
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
    subroutine monte_carlo_level(level, total_sites, config, temperature, samples_level, total_comb, use_parallel, &
                                 summary_unit, summary_txt_unit)
        integer, intent(in) :: level, total_sites, samples_level
        integer(ip), intent(in) :: total_comb
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: temperature
        logical, intent(in) :: use_parallel
        integer, intent(in) :: summary_unit, summary_txt_unit

        integer :: subset(max(1, level))
        integer :: valid_samples
        real(dp), allocatable :: energies(:), energies_low(:), energies_high(:)
        real(dp) :: best_energy
        integer :: best_subset(max(1, level))
        integer :: best_count
        integer :: trace_unit

        if (level == 0) then
            call exhaustive_level(level, total_sites, config, temperature, total_comb, 1, use_parallel, &
                                  summary_unit, summary_txt_unit)
            return
        end if

        allocate(energies(samples_level))
        allocate(energies_low(samples_level))
        allocate(energies_high(samples_level))
        valid_samples = 0
        best_energy = huge(1.0_dp)
        best_subset = 0
        best_count = 0
        call open_mc_trace_file(level, trace_unit)

        do while (valid_samples < samples_level)
            call random_subset(total_sites, level, subset)
            config = 1
            config(subset(1:level)) = 2
            valid_samples = valid_samples + 1
            call calculate_structure_energy(config, total_sites, energies(valid_samples), &
                                            energies_low(valid_samples), energies_high(valid_samples))
            if (energies(valid_samples) < best_energy) then
                best_energy = energies(valid_samples)
                best_subset(1:level) = subset(1:level)
                best_count = level
            end if
            if (trace_unit /= 0) then
                call write_mc_trace_step(trace_unit, valid_samples, valid_samples, energies(valid_samples), &
                                         energies_low(valid_samples), energies_high(valid_samples))
                call flush(trace_unit)
            end if
        end do

    call summarize_level(level, total_sites, temperature, total_comb, real(valid_samples, dp), energies(1:valid_samples), &
                 best_energy, best_subset, best_count, 0, &
                 energies_low(1:valid_samples), energies_high(1:valid_samples), &
                 use_parallel, summary_unit, summary_txt_unit)

        if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
        deallocate(energies, energies_low, energies_high)
    end subroutine monte_carlo_level

    ! Runs a Metropolis-Hastings walk with optional restart moves for a substitution level.
    subroutine metropolis_level(level, total_sites, config, temperature, samples_level, total_comb, use_parallel, &
                                 summary_unit, summary_txt_unit)
        integer, intent(in) :: level, total_sites, samples_level
        integer(ip), intent(in) :: total_comb
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: temperature
        logical, intent(in) :: use_parallel
        integer, intent(in) :: summary_unit, summary_txt_unit

        integer :: current_subset(max(1, level))
        integer :: trial_subset(max(1, level))
        integer :: best_subset(max(1, level))
        integer :: best_count
        integer :: accept_count, max_trials, remove_idx, add_site
        integer :: i_attempt
        integer :: skipped
        integer :: tries
        real(dp) :: beta, delta_e, rand_num
        real(dp) :: current_energy, current_low, current_high
        real(dp) :: trial_energy, trial_low, trial_high
        real(dp) :: best_energy
        real(dp), allocatable :: energies(:), energies_low(:), energies_high(:)
        logical :: accepted, valid_trial
        integer :: trace_unit
        integer :: best_step
        integer :: burn_keep, burn_start
        logical :: best_in_subset
        real(dp), parameter :: restart_move_prob = 0.01_dp
        logical :: use_restart_move
        integer :: restart_attempts, restart_accepts
        integer :: flip_attempts, flip_accepts

        if (level == 0) then
            call exhaustive_level(level, total_sites, config, temperature, total_comb, 1, use_parallel, &
                                  summary_unit, summary_txt_unit)
            return
        end if

        allocate(energies(samples_level))
        allocate(energies_low(samples_level))
        allocate(energies_high(samples_level))

        beta = 1.0_dp / (kB_eVk * temperature)
        skipped = 0

        call random_subset(total_sites, level, current_subset)
        config = 1
        config(current_subset(1:level)) = 2
        call calculate_structure_energy(config, total_sites, current_energy, current_low, current_high)
        best_energy = current_energy
        best_subset = current_subset
        best_count = level
        accept_count = 1
        energies(accept_count) = current_energy
        energies_low(accept_count) = current_low
        energies_high(accept_count) = current_high
        best_step = accept_count
        call open_mc_trace_file(level, trace_unit)
        if (trace_unit /= 0) then
            call write_mc_trace_step(trace_unit, accept_count, 0, current_energy, current_low, current_high)
            call flush(trace_unit)
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
            call calculate_structure_energy(config, total_sites, trial_energy, trial_low, trial_high)

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
                if (current_energy < best_energy) then
                    best_energy = current_energy
                    best_subset = current_subset
                    best_count = level
                    best_step = accept_count
                end if
                if (trace_unit /= 0) then
                    call write_mc_trace_step(trace_unit, accept_count, i_attempt, current_energy, current_low, current_high)
                    call flush(trace_unit)
                end if
            else
                skipped = skipped + 1
            end if
        end do

        burn_keep = max(1, ceiling(0.8_dp * real(accept_count, dp)))
        burn_start = max(1, accept_count - burn_keep + 1)
        best_in_subset = (best_step >= burn_start)
        if (burn_start > 1) then
            write(*,'(A,I0)') 'Metropolis: configuraciones descartadas en equilibrado = ', burn_start - 1
        end if

        call summarize_level(level, total_sites, temperature, total_comb, real(accept_count, dp), &
                 energies(burn_start:accept_count), best_energy, best_subset, best_count, skipped, &
                 energies_low(burn_start:accept_count), energies_high(burn_start:accept_count), &
                 use_parallel, summary_unit, summary_txt_unit, best_in_subset, &
                 restart_attempts, restart_accepts, flip_attempts, flip_accepts)

        if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
        deallocate(energies, energies_low, energies_high)
    end subroutine metropolis_level

    ! Computes Boltzmann statistics for a level and records them to screen and summaries.
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
        real(dp) :: low_ref, high_ref, expected_mix, weight_low, weight_high
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
    character(len=256) :: csv_line, txt_line
    character(len=32) :: accept_ratio_str
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
!$omp parallel do default(shared) private(i) reduction(+:wsum, expected_sum) if(use_parallel)
        do i = 1, ncomb
            weights(i) = exp(-beta * (energies(i) - emin_ref))
            wsum = wsum + weights(i)
            expected_sum = expected_sum + weights(i) * energies(i)
        end do
!$omp end parallel do
        if (wsum > 0.0_dp) then
            expected = expected_sum / wsum
        else
            expected = emin
        end if
        variance_sum = 0.0_dp
        if (wsum > 0.0_dp) then
!$omp parallel do default(shared) private(i) reduction(+:variance_sum) if(use_parallel)
            do i = 1, ncomb
                variance_sum = variance_sum + weights(i) * (energies(i) - expected)**2
            end do
!$omp end parallel do
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
            write(label_low,'(F12.6)') low_min
            valid_low = .true.
            low_ref = low_min
            allocate(weights_low(ncomb))
            wsum_low = 0.0_dp
            expected_low_sum = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:wsum_low, expected_low_sum) if(use_parallel)
            do i = 1, ncomb
                weights_low(i) = exp(-beta * (energies_low(i) - low_ref))
                wsum_low = wsum_low + weights_low(i)
                expected_low_sum = expected_low_sum + weights_low(i) * energies_low(i)
            end do
!$omp end parallel do
            if (wsum_low > 0.0_dp) then
                expected_low = expected_low_sum / wsum_low
            else
                expected_low = low_ref
            end if
            write(exp_low_str,'(F12.6)') expected_low
            deallocate(weights_low)
        end if

        if (have_high) then
            if (.not. all(energies_high > 0.5_dp * huge_marker)) then
                high_min = minval(energies_high)
                write(label_high,'(F12.6)') high_min
                valid_high = .true.
                high_ref = high_min
                allocate(weights_high(ncomb))
                wsum_high = 0.0_dp
                expected_high_sum = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:wsum_high, expected_high_sum) if(use_parallel)
                do i = 1, ncomb
                    weights_high(i) = exp(-beta * (energies_high(i) - high_ref))
                    wsum_high = wsum_high + weights_high(i)
                    expected_high_sum = expected_high_sum + weights_high(i) * energies_high(i)
                end do
!$omp end parallel do
                if (wsum_high > 0.0_dp) then
                    expected_high = expected_high_sum / wsum_high
                else
                    expected_high = high_ref
                end if
                write(exp_high_str,'(F12.6)') expected_high
                deallocate(weights_high)
            end if
        end if

        ge_fraction = 0.0_dp
        if (total_sites > 0) ge_fraction = real(level, dp) / real(total_sites, dp)
        write(frac_str,'(F7.4)') ge_fraction

        if (valid_low .and. valid_high) then
            weight_low = mixing_weight(ge_fraction)
            weight_low = max(0.0_dp, min(1.0_dp, weight_low))
            weight_high = 1.0_dp - weight_low
            expected_mix = weight_low * expected_low + weight_high * expected_high
        else if (valid_low) then
            expected_mix = expected_low
        else if (valid_high) then
            expected_mix = expected_high
        else
            expected_mix = expected
        end if
        if (valid_low .or. valid_high) then
            write(mix_exp_str,'(F12.6)') expected_mix
        end if

        write(exp_total_str,'(F12.6)') expected
        write(min_total_str,'(F12.6)') best_energy
        write(variance_str,'(F12.6)') variance
        rel_exp_total = quartz_relative(level, total_sites, expected)
        rel_min_total = quartz_relative(level, total_sites, best_energy)
        write(delta_exp_total_str,'(F12.6)') rel_exp_total
        write(delta_min_total_str,'(F12.6)') rel_min_total
        if (valid_low) then
            rel_exp_low = quartz_relative(level, total_sites, expected_low)
            rel_min_low = quartz_relative(level, total_sites, low_min)
            write(delta_exp_low_str,'(F12.6)') rel_exp_low
            write(delta_min_low_str,'(F12.6)') rel_min_low
        end if
        if (valid_high) then
            rel_exp_high = quartz_relative(level, total_sites, expected_high)
            rel_min_high = quartz_relative(level, total_sites, high_min)
            write(delta_exp_high_str,'(F12.6)') rel_exp_high
            write(delta_min_high_str,'(F12.6)') rel_min_high
        end if
        if (valid_low .or. valid_high) then
            rel_exp_mix = quartz_relative(level, total_sites, expected_mix)
            write(delta_mix_str,'(F12.6)') rel_exp_mix
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
!$omp critical(summary_io)
        write(*,'(A)') separator
        write(*,'(A,I3)') 'Sustituciones (N): ', level
        write(*,'(A)') 'Combinaciones totales: '//trim(total_str)
        write(*,'(A,F10.0)') 'Configuraciones procesadas: ', processed
        write(*,'(A,I8)') 'Intentos descartados: ', max(skipped, 0)
        write(*,'(A,F12.6)') 'Energía mínima (eV): ', best_energy
        write(*,'(A,F12.6)') 'Energía esperada <E> (eV): ', expected
        write(*,'(A,F12.6)') 'Varianza respecto a <E> (eV^2): ', variance
        write(*,'(A,F12.6)') 'Desviación estándar (eV): ', sqrt(variance)
        write(*,'(A,A)') 'Energía mínima lado Si (eV): ', trim(label_low)
        write(*,'(A,A)') 'Energía esperada lado Si (eV): ', trim(exp_low_str)
        write(*,'(A,A)') 'Energía mínima lado Ge (eV): ', trim(label_high)
        write(*,'(A,A)') 'Energía esperada lado Ge (eV): ', trim(exp_high_str)
        write(*,'(A,A)') 'Energía esperada combinada (eV): ', trim(mix_exp_str)
        write(*,'(A,F10.6)') 'Probabilidad Boltzmann del mínimo: ', prob_best
        write(*,'(A,F8.4)') 'Ratio de aceptación: ', accept_ratio
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
            write(*,'(A,I8,A,I8,A,F8.4)') 'Restart aceptados: ', restart_accepts_val, ' / ', restart_attempts_val, &
                '  ratio: ', restart_ratio
            write(*,'(A,I8,A,I8,A,F8.4)') 'Flip aceptados:    ', flip_accepts_val, ' / ', flip_attempts_val, &
                '  ratio: ', flip_ratio
        end if
        call print_best_positions(best_positions, best_count)
        call save_best_structure_poscar(level, total_sites, best_positions, best_count)
        if (processed < real(total_comb, dp)) then
            write(*,'(A)') 'Aviso: no se cubrieron todas las combinaciones (se empleó muestreo).'
        end if
        if (level == 0) then
            write(*,'(A)') 'FracGe; E_esperada(ladoSi); E_min(ladoSi); E_esperada(ladoGe); E_min(ladoGe)'
        end if
        write(*,'(A)') trim(frac_str)//'; '//trim(exp_low_str)//'; '//trim(label_low)//'; ' &
            //trim(exp_high_str)//'; '//trim(label_high)
        if (summary_unit /= 0) then
            csv_line = trim(adjustl(level_str))//';'//trim(adjustl(frac_str))//';'//trim(adjustl(exp_total_str))//';'// &
                trim(adjustl(min_total_str))//';'//trim(adjustl(variance_str))//';'//trim(adjustl(exp_low_str))//';'//trim(adjustl(label_low))//';'// &
                trim(adjustl(exp_high_str))//';'//trim(adjustl(label_high))//';'//trim(adjustl(mix_exp_str))//';'// &
                trim(adjustl(delta_exp_total_str))//';'//trim(adjustl(delta_min_total_str))//';'// &
                trim(adjustl(delta_exp_low_str))//';'//trim(adjustl(delta_min_low_str))//';'// &
                trim(adjustl(delta_exp_high_str))//';'//trim(adjustl(delta_min_high_str))//';'// &
                trim(adjustl(delta_mix_str))//';'//trim(adjustl(accept_ratio_str))
            write(summary_unit,'(A)') trim(csv_line)
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
        end if
!$omp end critical(summary_io)
        deallocate(weights)
    end subroutine summarize_level

    ! Sorts the summary files by substitution level to keep outputs monotonic.
    subroutine reorder_summary_outputs()
        call reorder_single_summary(summary_filename, ';')
        call reorder_single_summary(summary_txt_filename, ' ')
    end subroutine reorder_summary_outputs

    ! Loads a summary file, sorts its lines by level, and writes it back.
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

    ! Converts an energy to a relative enthalpy per atom with quartz references.
    real(dp) function quartz_relative(level, total_sites, energy) result(rel_val)
        integer, intent(in) :: level, total_sites
        real(dp), intent(in) :: energy
        real(dp) :: ge_atoms, si_atoms, denom

        if (total_sites <= 0) then
            rel_val = 0.0_dp
            return
        end if
        ge_atoms = real(level, dp)
        si_atoms = real(total_sites - level, dp)
        denom = real(total_sites, dp)
        rel_val = conv_ev_to_kjmol * (energy - ge_atoms * quartz_ge_energy - si_atoms * quartz_si_energy) / denom
    end function quartz_relative

    ! Evaluates the empirical mixing weight based on the Ge fraction.
    real(dp) function mixing_weight(ge_fraction) result(weight)
        real(dp), intent(in) :: ge_fraction
        real(dp) :: arg, numerator, denom_base, denominator
        real(dp), parameter :: eps = 1.0e-8_dp

        arg = (ge_fraction - mix_d0) / mix_x0
        denom_base = 1.0_dp - arg
        if (abs(denom_base) < eps) then
            denom_base = merge(eps, -eps, denom_base >= 0.0_dp)
        end if
        denominator = denom_base**mix_m
        numerator = 1.0_dp - arg**mix_n
        if (abs(denominator) < eps) then
            weight = merge(1.0_dp, 0.0_dp, numerator >= 0.0_dp)
        else
            weight = numerator / denominator
        end if
        weight = max(0.0_dp, min(1.0_dp, weight))
    end function mixing_weight

    ! Writes the minimum-energy configuration for a level as a VASP POSCAR file.
    subroutine save_best_structure_poscar(level, total_sites, best_positions, best_count)
        integer, intent(in) :: level, total_sites, best_count
        integer, intent(in) :: best_positions(:)

        integer, allocatable :: best_config(:)
        character(len=64) :: filename
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
            write(filename,'("POSCAR_N",I4.4,".vasp")') level
            call write_vasp_file(best_config, total_sites, trim(filename))
            write(*,'(A)') 'POSCAR minimo guardado en '//trim(filename)
        else
            write(*,'(A)') 'Aviso: no se pudo generar el POSCAR minimo por posiciones invalidas.'
        end if

        deallocate(best_config)
    end subroutine save_best_structure_poscar

    ! Creates a trace CSV for Monte Carlo steps at the specified level.
    subroutine open_mc_trace_file(level, unit_id)
        integer, intent(in) :: level
        integer, intent(out) :: unit_id
        character(len=64) :: filename
        integer :: ios

        write(filename,'("MC_TRACE_N",I4.4,".csv")') level
        open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,I0)') 'Aviso: no se pudo abrir archivo de traza MC para N=', level
            unit_id = 0
        else
            write(unit_id,'(A)') '#step;attempt;energy_eV;energy_low_eV;energy_high_eV'
        end if
    end subroutine open_mc_trace_file

    ! Appends one Monte Carlo sample entry to the trace CSV if the unit is valid.
    subroutine write_mc_trace_step(unit_id, step_idx, attempt_idx, energy, energy_low, energy_high)
        integer, intent(in) :: unit_id
        integer, intent(in) :: step_idx, attempt_idx
        real(dp), intent(in) :: energy, energy_low, energy_high

        if (unit_id == 0) return
        write(unit_id,'(I0,";",I0,";",F18.10,";",F18.10,";",F18.10)') step_idx, attempt_idx, energy, energy_low, energy_high
    end subroutine write_mc_trace_step

    ! Closes the trace CSV unit when tracing has finished.
    subroutine close_mc_trace_file(unit_id)
        integer, intent(in) :: unit_id

        if (unit_id /= 0) close(unit_id)
    end subroutine close_mc_trace_file

    ! Prints the indices of the lowest-energy configuration sampled for a level.
    subroutine print_best_positions(best_positions, best_count)
        integer, intent(in) :: best_positions(:)
        integer, intent(in) :: best_count

        if (best_count <= 0) then
            write(*,'(A)') 'Configuración de mínima energía: configuración base'
        else
            write(*,'(A)', advance='no') 'Configuración de mínima energía (posiciones 1..n): '
            write(*,'(*(1X,I3))') best_positions(1:best_count)
        end if
    end subroutine print_best_positions
end program sod_boltzmann_mc
