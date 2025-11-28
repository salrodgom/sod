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

! Monte Carlo driver that samples substitution levels, calls the SOD energy
! calculator, and aggregates Boltzmann-weighted statistics per level.
program sod_boltzmann_mc
    use sod_boltzmann_consts
    use sod_boltzmann_utils
    use sod_calibration
    use energy_calc
    use omp_lib, only: omp_in_parallel
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    integer, parameter :: max_exact_combos = 200000
    integer, parameter :: mix_n = 6
    integer, parameter :: mix_m = 12
    integer, parameter :: uniform_unique_cap = 100
    integer, parameter :: uniform_unique_min_cap = 25
    real(dp), parameter :: uniform_cap_shrink = 0.75_dp
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
    logical :: force_mc_sampling
    integer :: level_min, level_max
    integer :: level_start, level_end
    
    integer, allocatable :: eqmatrix(:,:), config(:), local_config(:)
    integer, allocatable :: level_overrides(:)
    integer, allocatable :: level_targets(:)
    integer :: nop, total_sites
    integer :: level, effective_max
    integer :: level_idx
    logical :: has_level_overrides
    logical :: allow_parallel_levels
    logical :: effective_use_parallel
    logical :: force_parallel_lists
    character(len=512) :: osda_gin_option
    
    use_parallel = .false.
    omp_available = .false.
    !$  use_parallel = .true.
    !$  omp_available = .true.
    
    osda_gin_option = 'default'
    call parse_arguments(temperature, level_min, level_max, max_substitutions, samples_per_level, seed_value, sampling_mode, use_parallel, omp_available, force_mc_sampling, level_overrides, has_level_overrides, force_parallel_lists, osda_gin_option)
    call set_calibration_osda_gin(trim(osda_gin_option))
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
    allow_parallel_levels = use_parallel
    effective_use_parallel = use_parallel
    
    if (max_substitutions < 0) then
        effective_max = total_sites
    else
        effective_max = min(max_substitutions, total_sites)
    end if
    
    if (has_level_overrides) then
        call finalize_level_overrides(level_overrides, level_targets, total_sites)
        if (.not. allocated(level_targets)) then
            write(*,'(A)') 'Error: la lista de niveles especificada en -N no contiene valores válidos.'
            stop 1
        end if
        if (allow_parallel_levels .and. .not. force_parallel_lists) then
            write(*,'(A)') 'Aviso: se desactiva el paralelismo externo para respetar el orden de -N.'
            allow_parallel_levels = .false.
            if (effective_use_parallel) then
                write(*,'(A)') 'Aviso: el muestreo por nivel se ejecutará en modo secuencial para mayor estabilidad.'
                effective_use_parallel = .false.
            end if
        else if (force_parallel_lists .and. allow_parallel_levels) then
            write(*,'(A)') 'Aviso: se mantiene el paralelismo con lista explícita (--parallel-lists); los resultados pueden llegar fuera de orden.'
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
    
    write(*,'(A)') '--- Parámetros del cálculo ---'
    write(*,'(A,F10.2)') 'Temperatura (K): ', temperature
    write(*,'(A,I6)') 'Sitios sustituibles (npos): ', total_sites
    write(*,'(A,I6)') 'Max sustituciones evaluadas: ', effective_max
    if (has_level_overrides) then
        call print_level_overrides(level_targets)
    else
        write(*,'(A,I0,A,I0)') 'Niveles evaluados: ', level_start, ' .. ', level_end
    end if
    write(*,'(A,I8)') 'Umbral enumeración exacta: ', max_exact_combos
    write(*,'(A,I8)') 'Muestras aleatorias (si se supera el umbral): ', samples_per_level
    write(*,'(A)') 'Resultados por nivel guardados en: '//trim(summary_filename)
    write(*,'(A)') '                               y: '//trim(summary_txt_filename)
    write(*,'(A)') 'Método de muestreo MC (si aplica): '//trim(sampling_mode)
    write(*,'(A)') 'OpenMP paralelo: '//merge('Si','No',effective_use_parallel)
    write(*,'(A)') 'Forzar muestreo MC: '//merge('Si','No',force_mc_sampling)
    if (use_parallel) then
        write(*,'(A)') 'Nota: las salidas por nivel pueden imprimirse en orden no secuencial durante el cálculo paralelo.'
        write(*,'(A)') '      Los archivos de resumen se reordenan al finalizar para dejar los niveles crecientes.'
    end if
    write(*,*)
    
    if (allow_parallel_levels) then
        !$omp parallel default(shared) private(local_config, level, level_idx)
        allocate(local_config(total_sites))
        if (has_level_overrides) then
            !$omp do schedule(dynamic)
            do level_idx = 1, size(level_targets)
                level = level_targets(level_idx)
                call process_level(level, total_sites, local_config, temperature, samples_per_level, &
                max_exact_combos, sampling_mode, force_mc_sampling, effective_use_parallel, summary_unit, summary_txt_unit)
            end do
            !$omp end do
        else
            !$omp do schedule(dynamic)
            do level = level_start, level_end
                call process_level(level, total_sites, local_config, temperature, samples_per_level, &
                max_exact_combos, sampling_mode, force_mc_sampling, effective_use_parallel, summary_unit, summary_txt_unit)
            end do
            !$omp end do
        end if
        deallocate(local_config)
        !$omp end parallel
    else
        allocate(config(total_sites))
        if (has_level_overrides) then
            do level_idx = 1, size(level_targets)
                level = level_targets(level_idx)
                call process_level(level, total_sites, config, temperature, samples_per_level, &
                max_exact_combos, sampling_mode, force_mc_sampling, effective_use_parallel, summary_unit, summary_txt_unit)
            end do
        else
            do level = level_start, level_end
                call process_level(level, total_sites, config, temperature, samples_per_level, &
                max_exact_combos, sampling_mode, force_mc_sampling, effective_use_parallel, summary_unit, summary_txt_unit)
            end do
        end if
        deallocate(config)
    end if
    
    call close_summary_files(summary_unit, summary_txt_unit)
    call reorder_summary_outputs()
    call cleanup_energy_calc()
    
contains
    
    ! Parses optional command-line arguments and populates runtime parameters.
    subroutine parse_arguments(temp, level_min, level_max, max_subs, samples_level, seed, sampler, use_parallel, omp_available, force_mc, level_list, has_level_list, force_parallel_lists, osda_gin_option)
        real(dp), intent(out) :: temp
        integer, intent(out) :: level_min, level_max, max_subs, samples_level, seed
        character(len=*), intent(out) :: sampler
        logical, intent(inout) :: use_parallel
        logical, intent(in) :: omp_available
        logical, intent(out) :: force_mc
        integer, allocatable, intent(out) :: level_list(:)
        logical, intent(out) :: has_level_list
        logical, intent(out) :: force_parallel_lists
        character(len=*), intent(out) :: osda_gin_option
        integer :: argc, ios, i, j, colon_pos
        character(len=256) :: carg, lowered, spec
        character(len=256), allocatable :: args(:)
        logical, allocatable :: skip(:)
        logical :: found_range, handled
        logical :: temp_set, max_subs_set, samples_set, seed_set, sampler_set, omp_set, osda_set
        character(len=512) :: spec_trim
        integer :: list_status
        integer :: eq_pos
        
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
        osda_gin_option = 'default'
        osda_set = .false.
        
        argc = command_argument_count()
        if (argc <= 0) return
        
        allocate(args(argc))
        allocate(skip(argc))
        skip = .false.
        
        do i = 1, argc
            call get_command_argument(i, args(i))
            if (is_help_argument(args(i))) then
                call print_usage(omp_available)
                stop 0
            end if
        end do
        
        do i = 1, argc
            lowered = adjustl(args(i))
            do j = 1, len_trim(lowered)
                if (lowered(j:j) >= 'A' .and. lowered(j:j) <= 'Z') lowered(j:j) = achar(iachar(lowered(j:j)) + 32)
            end do
            if (trim(lowered) == '--force-mc' .or. trim(lowered) == '--force_mc' .or. &
            trim(lowered) == 'force-mc' .or. trim(lowered) == 'forcemc') then
                force_mc = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--parallel-lists' .or. trim(lowered) == '--parallel_lists' .or. &
            trim(lowered) == 'parallel-lists' .or. trim(lowered) == 'parallellists') then
                force_parallel_lists = .true.
                skip(i) = .true.
            else if (index(trim(lowered), '--osda-gin=') == 1) then
                eq_pos = index(args(i), '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(args(i))) then
                    write(error_unit,'(A)') 'Error: argumento inválido para --osda-gin.'
                    call print_usage(omp_available)
                    stop 1
                end if
                osda_gin_option = adjustl(args(i)(eq_pos+1:))
                osda_set = .true.
                skip(i) = .true.
            else if (trim(lowered) == '--osda-gin' .or. trim(lowered) == '--osda_gin') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: falta ruta después de --osda-gin.'
                    call print_usage(omp_available)
                    stop 1
                end if
                osda_gin_option = adjustl(args(i+1))
                osda_set = .true.
                skip(i) = .true.
                skip(i+1) = .true.
            else if (trim(lowered) == '--no-osda-gin' .or. trim(lowered) == '--skip-osda' .or. &
                 trim(lowered) == '--skip_osda') then
                osda_gin_option = 'none'
                osda_set = .true.
                skip(i) = .true.
            end if
        end do
        
        found_range = .false.
        do i = 1, argc
            if (skip(i)) cycle
            lowered = adjustl(args(i))
            do j = 1, len_trim(lowered)
                if (lowered(j:j) >= 'A' .and. lowered(j:j) <= 'Z') lowered(j:j) = achar(iachar(lowered(j:j)) + 32)
            end do
            if (trim(lowered) == '-n') then
                if (i == argc) then
                    write(error_unit,'(A)') 'Error: falta especificación después de -N.'
                    call print_usage(omp_available)
                    stop 1
                end if
                spec = adjustl(args(i+1))
                spec_trim = adjustl(spec)
                colon_pos = index(spec_trim, ':')
                if (index(spec_trim, ',') > 0) then
                    call parse_level_list(spec_trim, level_list, list_status)
                    if (list_status /= 0) then
                        write(error_unit,'(A)') 'Error: especificación de lista inválida en -N.'
                        call print_usage(omp_available)
                        stop 1
                    end if
                    has_level_list = .true.
                    level_min = 0
                    level_max = -1
                else if (colon_pos > 0) then
                    read(spec_trim(1:colon_pos-1),*,iostat=ios) level_min
                    if (ios /= 0) then
                        write(error_unit,'(A)') 'Error: límite inferior inválido en -N.'
                        call print_usage(omp_available)
                        stop 1
                    end if
                    read(spec_trim(colon_pos+1:),*,iostat=ios) level_max
                    if (ios /= 0) then
                        write(error_unit,'(A)') 'Error: límite superior inválido en -N.'
                        call print_usage(omp_available)
                        stop 1
                    end if
                else
                    read(spec_trim,*,iostat=ios) level_min
                    if (ios /= 0) then
                        write(error_unit,'(A)') 'Error: especificación de nivel inválida en -N.'
                        call print_usage(omp_available)
                        stop 1
                    end if
                    if (level_min < 0) then
                        level_min = 0
                        level_max = -1
                    else
                        level_max = level_min
                    end if
                end if
                skip(i) = .true.
                if (i < argc) skip(i+1) = .true.
                found_range = .true.
            end if
        end do
        temp_set = .false.
        max_subs_set = .false.
        samples_set = .false.
        seed_set = .false.
        sampler_set = .false.
        omp_set = .false.
        if (.not. osda_set) osda_gin_option = 'default'
        do i = 1, argc
            if (skip(i)) cycle
            carg = args(i)
            handled = .false.
            
            lowered = adjustl(carg)
            do j = 1, len_trim(lowered)
                if (lowered(j:j) >= 'A' .and. lowered(j:j) <= 'Z') lowered(j:j) = achar(iachar(lowered(j:j)) + 32)
            end do
            
            if (.not. sampler_set) then
                if (trim(lowered) == 'uniform' .or. trim(lowered) == 'metropolis') then
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
                        write(error_unit,'(A)') 'Advertencia: Nsamples debe ser positivo, se usan 5000 muestras'
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
                write(error_unit,'(A)') 'Advertencia: argumento ignorado -> '//trim(carg)
            end if
        end do
        
        if (.not. found_range .and. level_max >= 0) then
            if (level_min < 0) level_min = 0
        end if
        
        if (allocated(args)) deallocate(args)
        if (allocated(skip)) deallocate(skip)
        return
    end subroutine parse_arguments

    subroutine parse_level_list(spec_in, out_levels, status)
        character(len=*), intent(in) :: spec_in
        integer, allocatable, intent(out) :: out_levels(:)
        integer, intent(out) :: status
        character(len=512) :: buffer
        character(len=256) :: token
        integer :: length, count, pos, start_idx, idx, val

        status = 0
        buffer = trim(spec_in)
        length = len_trim(buffer)
        if (length <= 0) then
            status = 1
            return
        end if
        if (buffer(1:1) == ',' .or. buffer(length:length) == ',') then
            status = 1
            return
        end if

        count = 1
        do pos = 1, length
            if (buffer(pos:pos) == ',') then
                if (pos < length .and. buffer(pos+1:pos+1) == ',') then
                    status = 1
                    return
                end if
                count = count + 1
            end if
        end do

        allocate(out_levels(count))
        idx = 0
        start_idx = 1
        do pos = 1, length + 1
            if (pos > length .or. buffer(pos:pos) == ',') then
                idx = idx + 1
                token = buffer(start_idx:pos-1)
                token = adjustl(token)
                read(token,*,iostat=status) val
                if (status /= 0) then
                    if (allocated(out_levels)) deallocate(out_levels)
                    return
                end if
                out_levels(idx) = val
                start_idx = pos + 1
            end if
        end do
    end subroutine parse_level_list

! Interprets an OpenMP flag token and updates the parallel execution mode.
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
            write(error_unit,'(A)') 'Advertencia: OpenMP no disponible en esta compilación, se ignora argumento.'
            use_parallel = .false.
        end if
    case ('noomp', 'no', '0', 'false')
        use_parallel = .false.
    case default
        write(error_unit,'(A)') 'Advertencia: argumento desconocido para modo paralelo; se mantiene configuración previa.'
    end select
end subroutine parse_omp_flag

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
                write(*,'(A,I0,A)') 'Aviso: nivel ', val, ' fuera de rango; se omite.'
                call flush(output_unit)
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

    subroutine print_level_overrides(levels)
        integer, intent(in) :: levels(:)
        integer :: idx
        write(*,'(A)') 'Niveles evaluados: lista específica'
        write(*,'(A)', advance='no') 'Valores: '
        do idx = 1, size(levels)
            if (idx > 1) write(*,'(A)', advance='no') ', '
            write(*,'(I0)', advance='no') levels(idx)
        end do
        write(*,*)
    end subroutine print_level_overrides

! Emits program usage information including optional OpenMP note.
subroutine print_usage(omp_available)
    logical, intent(in) :: omp_available
    
    character(len=32) :: cap_str
    character(len=32) :: shrink_str
    real(dp) :: shrink_percent
    
    write(cap_str,'(I0)') uniform_unique_min_cap
    shrink_percent = uniform_cap_shrink * 100.0_dp
    write(shrink_str,'(F5.1)') shrink_percent
    
    write(*,'(A)') 'Uso: sod_boltzmann_mc [T_K] [Nmax] [Nsamples] [seed] [sampler] [omp|noomp] [--force-mc] [-N rango]'
    write(*,'(A)') '       sod_boltzmann_mc --help'
    write(*,'(A)') ''
    write(*,'(A)') 'Argumentos opcionales (por defecto entre corchetes):'
    write(*,'(A)') '  -N espec   Rango o lista: -N 5 (solo 5), -N 3:8 (3 a 8), -N 12,30,45 (lista puntual).'
    write(*,'(A)') '  --parallel-lists   Mantiene OpenMP incluso con listas de -N (puede alterar el orden).' 
    write(*,'(A)') '  --osda-gin <fichero>  Usa el archivo OSDA indicado al generar .gin (default: OSDA_ITW.gin).'
    write(*,'(A)') '  --no-osda-gin        Evita añadir fragmentos OSDA a los .gin creados.'
    write(*,'(A)') '  T_K        Temperatura en Kelvin para los pesos de Boltzmann [1000].'
    write(*,'(A)') '  Nmax       Número máximo de sustituciones evaluadas cuando no se usa -N [-1 -> todos].'
    write(*,'(A)') '  Nsamples   Muestras MC por nivel cuando C(N,npos) supera el umbral [5000].'
    write(*,'(A)') '  seed       Semilla del generador aleatorio [-1 -> semilla desde system_clock].'
    write(*,'(A)') '  sampler    "uniform" (muestreo sin sesgo) o "metropolis" (Metropolis-Hastings) ['// &
    'uniform].'
    write(*,'(A)') '  --force-mc Fuerza muestreo Monte Carlo incluso si C(N,npos) <= umbral exacto.'
    if (omp_available) then
        write(*,'(A)') '  omp|noomp  Activa/desactiva OpenMP si el binario lo soporta [noomp].'
    else
        write(*,'(A)') '  omp|noomp  Control del uso de OpenMP (esta compilación no lo habilita).'
    end if
    write(*,'(A)') ''
    write(*,'(A)') 'Otros detalles:'
    write(*,'(A)') '  - Por defecto se evalúan todos los niveles desde N=0 hasta Nmax salvo que -N limite el conjunto.'
    write(*,'(A)') '  - Si C(N,npos) <= 200000 se enumeran todas las configuraciones; en caso contrario se muestrea MC.'
    write(*,'(A)') '  - El muestreo uniforme aplica control adaptativo del cupo de configuraciones únicas ' // &
    '(mínimo '//trim(cap_str)//', factor '//trim(adjustl(shrink_str))//'%).'
    write(*,'(A)') '  - Los resultados agregados de cada nivel se guardan en '//trim(summary_filename)//'.'
    write(*,'(A)') '  - El resumen en texto plano se escribe en '//trim(summary_txt_filename)//'.'
    write(*,'(A)') ''
    write(*,'(A)') 'Ejemplos:'
    write(*,'(A)') '  sod_boltzmann_mc                     # Ejecuta con valores por defecto.'
    write(*,'(A)') '  sod_boltzmann_mc 800 6 2000 1234 metropolis omp'
    write(*,'(A)') '                                        # 800 K, hasta 6 sustituciones,'// &
    ' Metropolis y OpenMP.'
    write(*,'(A)') '  sod_boltzmann_mc -N 12,30,45         # Evalúa solo los niveles 12, 30 y 45.'
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
    force_sampling, use_parallel, summary_unit, summary_txt_unit)
    integer, intent(in) :: level, total_sites, samples_level, max_exact
    integer, intent(inout) :: config(:)
    real(dp), intent(in) :: temperature
    character(len=*), intent(in) :: sampler
    logical, intent(in) :: force_sampling
    logical, intent(in) :: use_parallel
    integer, intent(in) :: summary_unit, summary_txt_unit
    
    integer(ip) :: total_comb
    integer :: ncomb_int
    logical :: use_sampling
    logical :: force_this_level
    
    if (level > total_sites) return
    
    total_comb = binomial_int64(total_sites, level)
    if (total_comb == 0_ip) then
        write(*,'(A,I3)') 'Nivel ', level, ': sin combinaciones válidas.'
        return
    end if
    
    force_this_level = force_sampling .and. (level > 0 .and. level < total_sites)
    
    if (level == 0 .or. level == total_sites) then
        use_sampling = .false.
    else
        use_sampling = force_this_level .or. (total_comb > max_exact)
    end if
    
    if (use_sampling) then
        if (force_this_level .and. total_comb <= max_exact) then
            write(*,'(A,I0,A,I0,A,I0,A)') 'Nivel ', level, ': muestreo MC forzado (combinaciones=', int(total_comb, kind=kind(max_exact)), ', umbral exacto=', max_exact, ').'
        else if (sampler == 'metropolis') then
            write(*,'(A,I0,A)') 'Nivel ', level, ': se usa muestreo Monte Carlo tipo Metropolis-Hastings.'
        else
            write(*,'(A,I0,A)') 'Nivel ', level, ': se usa muestreo Monte Carlo uniforme.'
        end if
        if (sampler == 'metropolis') then
            call metropolis_level(level, total_sites, config, temperature, samples_level, total_comb, &
            use_parallel, summary_unit, summary_txt_unit)
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
    integer :: unique_count
    integer :: unique_cap, initial_unique_cap
    integer :: samples_target
    integer :: attempt_count, max_attempts
    integer :: adaptive_window, adaptive_threshold, stall_counter, prev_unique_snapshot
    integer :: new_cap
    real(dp), allocatable :: energies(:), energies_low(:), energies_high(:)
    real(dp) :: best_energy
    integer :: best_subset(max(1, level))
    integer :: best_count
    integer :: best_idx
    integer :: trace_unit
    integer :: calib_best_idx
    integer, allocatable :: unique_subsets(:,:)
    integer, allocatable :: accept_attempt(:)
    integer, allocatable :: config_local(:)
    real(dp), allocatable :: low_contribs(:,:), high_contribs(:,:)
    real(dp) :: low_contrib_tmp(4), high_contrib_tmp(4)
    integer, allocatable :: eqmatrix(:,:)
    integer, allocatable :: canonical_subset(:)
    integer :: existing
    integer :: nop, npos
    integer(ip) :: cap_ip
    real(dp), allocatable :: gulp_energies(:)
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
    allocate(energies(unique_cap))
    allocate(energies_low(unique_cap))
    allocate(energies_high(unique_cap))
    allocate(unique_subsets(level, unique_cap))
    allocate(accept_attempt(unique_cap))
    allocate(low_contribs(4, unique_cap))
    allocate(high_contribs(4, unique_cap))
    allocate(gulp_energies(unique_cap))
    allocate(canonical_subset(level))
    call get_eqmatrix(eqmatrix, nop, npos)
    if (.not. allocated(eqmatrix)) then
        write(*,'(A)') 'Aviso: no se pudo obtener EQMATRIX para la deduplicacion de muestras.'
        call flush(output_unit)
        deallocate(energies, energies_low, energies_high, unique_subsets, low_contribs, high_contribs, gulp_energies, canonical_subset, accept_attempt)
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
    
    ! Deduplica muestras uniformes hasta reunir configuraciones simetricamente unicas.
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
                            write(*,'(A,I0,A,I0,A)') 'Nivel ', level, ': control adaptativo reduce objetivo a ', unique_cap, ' configuraciones unicas.'
                            call flush(output_unit)
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
        write(*,'(A,I0,A,I0,A)') 'Nivel ', level, ': muestreo aleatorio alcanzó ', attempt_count, ' intentos sin cubrir el cupo unico.'
        call flush(output_unit)
    end if
    
    if (unique_count == 0) then
        write(*,'(A,I0)') 'Nivel ', level, ': no se pudo obtener ninguna configuracion unica.'
        call flush(output_unit)
        if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
        if (allocated(eqmatrix)) deallocate(eqmatrix)
        deallocate(energies, energies_low, energies_high, unique_subsets, low_contribs, high_contribs, gulp_energies, canonical_subset, accept_attempt)
        return
    else
        write(*,'(A,I0,A,I0)') 'Nivel ', level, ': configuraciones unicas acumuladas: ', unique_count
        call flush(output_unit)
        if (cap_reduced) then
            write(*,'(A,I0,A,I0)') 'Nivel ', level, ': objetivo adaptativo final: ', unique_cap
            call flush(output_unit)
        else if (unique_count < initial_unique_cap) then
            write(*,'(A,I0,A,I0,A)') 'Nivel ', level, ': objetivo de ', initial_unique_cap, ' muestras unicas no alcanzado.'
            call flush(output_unit)
        end if
    end if
    
    if (unique_count > 0) then
        allow_inner_parallel = use_parallel
        if (allow_inner_parallel) then
            if (omp_in_parallel()) allow_inner_parallel = .false.
        end if
        
        if (allow_inner_parallel) then
            !$omp parallel default(shared) private(idx, low_contrib_tmp, high_contrib_tmp, config_local)
            allocate(config_local(total_sites))
            config_local = 1
            !$omp do schedule(dynamic)
            do idx = 1, unique_count
                config_local = 1
                config_local(unique_subsets(1:level, idx)) = 2
                call calculate_structure_energy(config_local, total_sites, energies(idx), &
                energy_low_side=energies_low(idx), energy_high_side=energies_high(idx), &
                low_contrib=low_contrib_tmp, high_contrib=high_contrib_tmp)
                low_contribs(:, idx) = low_contrib_tmp
                high_contribs(:, idx) = high_contrib_tmp
            end do
            !$omp end do
            deallocate(config_local)
            !$omp end parallel
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
            call flush(trace_unit)
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
    
    call attempt_calibration_from_samples(level, total_sites, 1, unique_count, unique_subsets, &
    energies, energies_low, energies_high, low_contribs, high_contribs, calib_best_idx)
    if (calib_best_idx > 0) then
        best_energy = energies(calib_best_idx)
        best_subset(1:level) = unique_subsets(:, calib_best_idx)
        best_count = level
    end if
    
    gulp_success = .false.
    gulp_energies = 0.0_dp
    ! Evalua energias precisas con GULP para las configuraciones unicas guardadas.
    call evaluate_subsets_with_gulp(level, total_sites, unique_count, unique_subsets(:,1:unique_count), config, gulp_energies, gulp_success)
    if (gulp_success) then
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
        write(*,'(A,I0,A,F16.6,A,F16.6)') 'Nivel ', level, ': media GULP = ', mean_energy, ' eV, desviacion = ', std_energy
        call flush(output_unit)
    else
        write(*,'(A,I0,A)') 'Nivel ', level, ': no se pudieron evaluar todas las muestras con GULP; se mantienen energias internas.'
        call flush(output_unit)
    end if
    config = 1
    
    call summarize_level(level, total_sites, temperature, total_comb, real(unique_count, dp), energies(1:unique_count), &
    best_energy, best_subset, best_count, 0, &
    energies_low(1:unique_count), energies_high(1:unique_count), &
    use_parallel, summary_unit, summary_txt_unit)
    
    if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
    if (allocated(eqmatrix)) deallocate(eqmatrix)
    deallocate(energies, energies_low, energies_high, unique_subsets, low_contribs, high_contribs, gulp_energies, canonical_subset, accept_attempt)
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
    integer :: calib_best_idx
    integer, allocatable :: sampled_subsets(:,:)
    real(dp), allocatable :: low_contribs(:,:), high_contribs(:,:)
    real(dp) :: low_contrib_tmp(4), high_contrib_tmp(4)
    
    if (level == 0) then
        call exhaustive_level(level, total_sites, config, temperature, total_comb, 1, use_parallel, &
        summary_unit, summary_txt_unit)
        return
    end if
    
    allocate(energies(samples_level))
    allocate(energies_low(samples_level))
    allocate(energies_high(samples_level))
    calib_best_idx = 0
    if (level > 0) then
        allocate(sampled_subsets(level, samples_level))
        allocate(low_contribs(4, samples_level))
        allocate(high_contribs(4, samples_level))
    end if
    
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
    
    if (level > 0 .and. accept_count >= burn_start) then
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
    
    call summarize_level(level, total_sites, temperature, total_comb, real(accept_count, dp), &
    energies(burn_start:accept_count), best_energy, best_subset, best_count, skipped, &
    energies_low(burn_start:accept_count), energies_high(burn_start:accept_count), &
    use_parallel, summary_unit, summary_txt_unit, best_in_subset, &
    restart_attempts, restart_accepts, flip_attempts, flip_accepts)
    
    if (trace_unit /= 0) call close_mc_trace_file(trace_unit)
    if (level > 0) then
        if (allocated(sampled_subsets)) deallocate(sampled_subsets)
        if (allocated(low_contribs)) deallocate(low_contribs)
        if (allocated(high_contribs)) deallocate(high_contribs)
    end if
    deallocate(energies, energies_low, energies_high)
end subroutine metropolis_level

subroutine attempt_calibration_from_samples(level, total_sites, sample_start, sample_end, subsets, energies, energies_low, energies_high, low_contribs, high_contribs, best_index)
    integer, intent(in) :: level, total_sites, sample_start, sample_end
    integer, intent(in) :: subsets(:,:)
    real(dp), intent(inout) :: energies(:), energies_low(:), energies_high(:)
    real(dp), intent(in) :: low_contribs(:,:), high_contribs(:,:)
    integer, intent(out), optional :: best_index
    
    integer :: sample_count, idx, actual_idx, best_global_idx
    integer :: hole_count
    integer :: nop, npos
    integer, allocatable :: eqmatrix(:,:), unique_subsets(:,:), unique_deg(:)
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
    real(dp) :: w_low, w_high, mix_beta, x_ge
    
    if (present(best_index)) best_index = 0
    
    if (level <= 0) return
    if (sample_end < sample_start) return
    
    sample_count = sample_end - sample_start + 1
    if (sample_count <= 0) return
    
    hole_count = total_sites - level
    need_low_calib = (level > get_max_low_order())
    need_high_calib = (hole_count > get_max_high_order())
    if (.not. (need_low_calib .or. need_high_calib)) return
    
    call get_eqmatrix(eqmatrix, nop, npos)
    if (.not. allocated(eqmatrix)) return
    if (level > size(subsets,1)) then
        if (allocated(eqmatrix)) deallocate(eqmatrix)
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
        write(*,'(A,I0,A,F16.6,A,F16.6)') 'Nivel ', level, ': energia lado Si min=', min_low, ', max=', max_low
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
        write(*,'(A,I0,A,F16.6,A,F16.6)') 'Nivel ', level, ': energia lado Ge min=', min_high, ', max=', max_high
    end if
    
    if (do_calibrate .and. unique_count >= 5) then
        allocate(config(total_sites))
        config = 1
        call calibrate_level_with_gulp(level, total_sites, unique_count, unique_subsets(:,1:unique_count), &
        unique_low_contrib(:,1:unique_count), unique_high_contrib(:,1:unique_count), unique_deg(1:unique_count), &
        unique_low(1:unique_count), unique_high(1:unique_count), config, need_low_calib, need_high_calib)
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
        
        ! Recompose total energies using calibrated low/high estimates.
        mix_beta = 8.0_dp
        x_ge = real(level, dp) / real(total_sites, dp)
        if (have_high) then
            w_high = 0.5_dp * (1.0_dp + tanh(mix_beta * (x_ge - 0.5_dp)))
        else
            w_high = 0.0_dp
        end if
        w_low = 1.0_dp - w_high
        
        do idx = 1, sample_count
            actual_idx = sample_start + idx - 1
            if (have_high .and. energies_high(actual_idx) < huge_marker) then
                energies(actual_idx) = w_low * energies_low(actual_idx) + w_high * energies_high(actual_idx)
            else
                energies(actual_idx) = energies_low(actual_idx)
            end if
        end do
        
        ! Propagate the updated best-sample index to the caller if needed.
        if (present(best_index)) then
            best_loc = minloc(energies(sample_start:sample_end))
            best_global_idx = sample_start + best_loc(1) - 1
            best_index = best_global_idx
        end if
        
        write(*,'(A,I0,A,I0)') 'Calibracion MC completada en nivel ', level, ' con ', unique_count
        call flush(output_unit)
    else
        if (do_calibrate) then
            write(*,'(A,I0)') 'Calibracion MC omitida por configuraciones unicas insuficientes en nivel ', level
        else if (need_low_calib .neqv. need_high_calib) then
            write(*,'(A,I0)') 'Calibracion MC omitida porque solo un lado requiere ajuste en nivel ', level
        end if
    end if
    
    if (allocated(eqmatrix)) deallocate(eqmatrix)
    deallocate(unique_subsets, unique_low, unique_high, unique_low_contrib, unique_high_contrib)
    deallocate(unique_low_max, unique_high_max)
    deallocate(unique_deg, sample_map, subset_buf, canonical)
    
end subroutine attempt_calibration_from_samples

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
