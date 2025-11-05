!*******************************************************************************
!  Exhaustive enumeration driver that prints all minimum-energy configurations
!  for each substitution level found via the SOD energy calculator.
!*******************************************************************************

program sod_boltzmann_exact
    use sod_boltzmann_consts
    use sod_boltzmann_utils
    use energy_calc
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none

    integer :: level_min, level_max
    real(dp) :: energy_tolerance
    integer :: nop, total_sites
    integer, allocatable :: eqmatrix(:,:)
    integer, allocatable :: scratch_config(:)
    integer :: level
    real(dp), parameter :: conv_ev_to_kjmol = 96.48533212331_dp
    real(dp), parameter :: quartz_ge_energy = -121.812200998519_dp
    real(dp), parameter :: quartz_si_energy = -128.799143746482_dp

    call parse_arguments_exact(level_min, level_max, energy_tolerance)

    call init_energy_calc()
    call get_eqmatrix(eqmatrix, nop, total_sites)
    if (.not. allocated(eqmatrix) .or. total_sites <= 0) then
        write(*,'(A)') 'Error: no se pudo obtener EQMATRIX o no hay sitios sustituibles.'
        stop 1
    end if
    if (allocated(eqmatrix)) deallocate(eqmatrix)

    ! Adjust range based on actual number of sites
    if (level_min < 0) level_min = 0
    if (level_max < 0) level_max = total_sites
    level_min = min(level_min, total_sites)
    level_max = min(level_max, total_sites)
    if (level_min > level_max) then
        write(*,'(A)') 'Error: rango de niveles inválido (min > max).'
        stop 1
    end if

    allocate(scratch_config(total_sites))

    write(*,'(A)') '--- Enumeración exhaustiva de configuraciones ---'
    write(*,'(A,I0)') 'Sitios sustituibles (npos): ', total_sites
    write(*,'(A,I0,A,I0)') 'Niveles evaluados: ', level_min, ' .. ', level_max
    write(*,'(A,ES12.5)') 'Tolerancia para agrupar mínimos (eV): ', energy_tolerance
    write(*,*)
    call flush(output_unit)

    do level = level_min, level_max
        call process_level_exact(level, total_sites, scratch_config, energy_tolerance)
    end do

    deallocate(scratch_config)
    call cleanup_energy_calc()
contains

    subroutine parse_arguments_exact(level_min, level_max, tol)
        integer, intent(out) :: level_min, level_max
        real(dp), intent(out) :: tol
        integer :: argc, ios, colon_pos, iarg
        character(len=256) :: carg, range_spec, tol_arg
        logical :: found_n_flag

        ! Defaults: all levels (0 to max), tolerance 1e-6
        level_min = 0
        level_max = -1
        tol = 1.0e-6_dp
        found_n_flag = .false.

        argc = command_argument_count()
        
        ! Check for help
        if (argc >= 1) then
            call get_command_argument(1, carg)
            if (is_help_token(carg)) then
                call print_usage_exact()
                stop 0
            end if
        end if

        ! Parse arguments looking for -N flag
        do iarg = 1, argc
            call get_command_argument(iarg, carg)
            if (trim(carg) == '-N' .or. trim(carg) == '-n') then
                if (iarg + 1 > argc) then
                    write(*,'(A)') 'Error: falta especificación después de -N.'
                    call print_usage_exact()
                    stop 1
                end if
                call get_command_argument(iarg + 1, range_spec)
                found_n_flag = .true.
                
                ! Parse range_spec: could be "N", "N:M", or "-1"
                colon_pos = index(range_spec, ':')
                if (colon_pos > 0) then
                    ! Range format "N:M"
                    read(range_spec(1:colon_pos-1), *, iostat=ios) level_min
                    read(range_spec(colon_pos+1:), *, iostat=ios) level_max
                    if (ios /= 0) then
                        write(*,'(A)') 'Error: formato de rango inválido. Use N:M (ej: 1:12).'
                        call print_usage_exact()
                        stop 1
                    end if
                else
                    ! Single value or -1
                    read(range_spec, *, iostat=ios) level_min
                    if (ios /= 0) then
                        write(*,'(A)') 'Error: especificación de nivel inválida.'
                        call print_usage_exact()
                        stop 1
                    end if
                    if (level_min < 0) then
                        ! -1 means all levels
                        level_min = 0
                        level_max = -1
                    else
                        ! Single level
                        level_max = level_min
                    end if
                end if
                exit
            end if
        end do

        ! Look for tolerance argument (last numeric argument not part of -N)
        if (argc >= 1) then
            if (found_n_flag) then
                ! If -N was used, tolerance would be after the range spec
                if (argc >= 3) then
                    call get_command_argument(3, tol_arg)
                    if (trim(tol_arg) /= '-N' .and. trim(tol_arg) /= '-n') then
                        read(tol_arg, *, iostat=ios) tol
                        if (ios /= 0 .or. tol <= 0.0_dp) then
                            write(*,'(A)') 'Aviso: tolerancia inválida, se usa 1e-6.'
                            tol = 1.0e-6_dp
                        end if
                    end if
                end if
            else
                ! Old-style positional arguments: Nmax [tol]
                call get_command_argument(1, carg)
                read(carg, *, iostat=ios) level_max
                if (ios == 0) then
                    if (level_max < 0) then
                        level_min = 0
                        level_max = -1
                    else
                        level_min = 0
                    end if
                else
                    ! Could not parse, default to all
                    level_min = 0
                    level_max = -1
                end if
                
                if (argc >= 2) then
                    call get_command_argument(2, carg)
                    read(carg, *, iostat=ios) tol
                    if (ios /= 0 .or. tol <= 0.0_dp) then
                        write(*,'(A)') 'Aviso: tolerancia inválida, se usa 1e-6.'
                        tol = 1.0e-6_dp
                    end if
                end if
            end if
        end if
    end subroutine parse_arguments_exact

    subroutine print_usage_exact()
        write(*,'(A)') 'Uso: sod_boltzmann_exact [-N <especificación>] [tol_eV]'
        write(*,'(A)') '     sod_boltzmann_exact [Nmax] [tol_eV]  (modo compatibilidad)'
        write(*,'(A)') ''
        write(*,'(A)') 'Argumentos opcionales:'
        write(*,'(A)') '  -N <spec>  Especifica niveles de sustitución a evaluar:'
        write(*,'(A)') '             -N -1      : Todos los niveles (0..npos)'
        write(*,'(A)') '             -N 12      : Sólo el nivel 12'
        write(*,'(A)') '             -N 1:12    : Rango del nivel 1 al 12 (inclusive)'
        write(*,'(A)') ''
        write(*,'(A)') '  tol_eV     Tolerancia en energía (eV) para considerar configuraciones'
        write(*,'(A)') '             equivalentes [por defecto: 1e-6].'
        write(*,'(A)') ''
        write(*,'(A)') 'Ejemplos:'
        write(*,'(A)') '  sod_boltzmann_exact -N 5:10 1e-5    # Niveles 5-10, tolerancia 1e-5'
        write(*,'(A)') '  sod_boltzmann_exact -N 8             # Sólo nivel 8'
        write(*,'(A)') '  sod_boltzmann_exact -N -1            # Todos los niveles'
        write(*,'(A)') '  sod_boltzmann_exact 12 1e-6          # Niveles 0-12 (modo antiguo)'
        write(*,'(A)') ''
        write(*,'(A)') 'El programa enumera exhaustivamente las combinaciones en cada nivel,'
        write(*,'(A)') 'guardando las configuraciones con energía mínima como archivos POSCAR_EXACT.'
    end subroutine print_usage_exact

    logical function is_help_token(raw)
        character(len=*), intent(in) :: raw
        character(len=len(raw)) :: token
        integer :: i, lt, offset

        token = adjustl(raw)
        lt = len_trim(token)
        if (lt == 0) then
            is_help_token = .false.
            return
        end if
        do i = 1, lt
            offset = iachar(token(i:i))
            if (offset >= iachar('A') .and. offset <= iachar('Z')) then
                token(i:i) = achar(offset + 32)
            end if
        end do
        token = token(1:lt)
        select case (token)
        case ('--help', '-h', 'help', '--ayuda', '-ayuda', 'ayuda', '/h', '/?')
            is_help_token = .true.
        case default
            is_help_token = .false.
        end select
    end function is_help_token

    subroutine process_level_exact(level, total_sites, config, tol)
        integer, intent(in) :: level, total_sites
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: tol

        integer(ip) :: total_comb
        real(dp) :: energy, low_estimate, high_estimate
        real(dp) :: best_energy
    integer :: best_count
    integer :: capacity, idx
    integer :: i, nop_local
    integer :: unique_capacity, unique_count, subset_idx
    integer :: total_min_degeneracy
    integer :: total_degeneracy_weight
    integer, allocatable :: subset(:)
    integer, allocatable :: best_subsets(:,:)
    integer, allocatable :: canonical_subset(:)
    integer, allocatable :: unique_subsets(:,:)
    integer, allocatable :: unique_deg(:)
    integer, allocatable :: best_deg(:)
    real(dp), allocatable :: best_low(:), best_high(:)
    real(dp), allocatable :: unique_energy(:), unique_low(:), unique_high(:)
    real(dp) :: degeneracy_energy_sum, degeneracy_mean_energy
    real(dp) :: boltzmann_weight_sum, boltzmann_energy_sum, boltzmann_energy_sq_sum
    real(dp) :: boltzmann_mean_energy, boltzmann_variance_energy
    real(dp) :: inverse_boltzmann_temperature, boltzmann_reference_energy
    integer :: unit_summary, ios_summary
    logical :: summary_exists
    character(len=64) :: summary_filename
    integer, allocatable :: eqmatrix_local(:,:)
    character(len=80), parameter :: separator = '------------------------------------------------------------------------'
    real(dp), parameter :: sentinel_energy = huge(1.0_dp)
    real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
    real(dp), parameter :: temp_targets(3) = (/300.0_dp, 800.0_dp, 1200.0_dp/)
    real(dp), parameter :: boltzmann_temperature = 300.0_dp
    real(dp) :: entropy_total
    real(dp) :: ts_terms(3)
    real(dp) :: free_energy_est(3)
    real(dp) :: deltaF_quartz_300
    real(dp) :: max_entropy, ideal_entropy, x, min_entropy
    real(dp) :: mean_low_all, mean_high_all
    real(dp) :: weighted_low_sum, weighted_high_sum
    integer :: total_low_weight, total_high_weight
    real(dp) :: min_low_energy, min_high_energy
    real(dp) :: lowest_energy_in_level
    real(dp) :: degeneracy_weight, boltzmann_factor
    real(dp) :: boltzmann_low_weight_sum, boltzmann_high_weight_sum
    real(dp) :: boltzmann_low_energy_sum, boltzmann_high_energy_sum

        ! Get EQMATRIX for symmetry checking
        call get_eqmatrix(eqmatrix_local, nop_local, i)

        total_comb = binomial_int64(total_sites, level)

        write(*,'(A)') separator
        write(*,'(A,I0)') 'Nivel: ', level
        write(*,'(A,I0)') 'Combinaciones totales: ', total_comb

        if (total_comb > 0_int64) then
            ! kB_eVk: Boltzmann constant in eV/K; dp: double precision kind
            max_entropy = kB_eVk * log(real(total_comb, dp))
        else
            max_entropy = 0.0_dp
        end if
        write(*,'(A,F12.6)') 'Entropia máxima (k_B ln W) [eV/K]: ', max_entropy

        if (total_sites > 0) then
            x = real(level, dp) / real(total_sites, dp)
            if (x <= 0.0_dp .or. x >= 1.0_dp) then
                ideal_entropy = 0.0_dp
            else
                ! Binary mixing entropy: S = -N k_B [x ln x + (1-x) ln(1-x)]
                ! Reference: Kittel, "Introduction to Solid State Physics", mixing entropy equation
                ideal_entropy = -real(total_sites, dp) * kB_eVk * (x * log(x) + (1.0_dp - x) * log(1.0_dp - x))
            end if
        else
            ideal_entropy = 0.0_dp
        end if
        write(*,'(A,F12.6)') 'Entropia ideal (mezcla binaria) [eV/K]: ', ideal_entropy
    call flush(output_unit)

        if (level == 0) then
            config = 1
            call calculate_structure_energy(config, total_sites, energy, low_estimate, high_estimate)
            total_degeneracy_weight = 1
            unique_count = 1
            degeneracy_mean_energy = energy
            boltzmann_mean_energy = energy
            boltzmann_variance_energy = 0.0_dp
            entropy_total = 0.0_dp
            ts_terms = 0.0_dp
            mean_low_all = low_estimate
            mean_high_all = high_estimate
            min_low_energy = low_estimate
            min_high_energy = high_estimate
            free_energy_est = boltzmann_mean_energy
            deltaF_quartz_300 = quartz_relative(level, total_sites, free_energy_est(1))

            write(*,'(A,F12.6)') 'Energía mínima (eV): ', energy
            write(*,'(A,F12.6)') 'Energía media (todos los inequivalentes) [eV]: ', boltzmann_mean_energy
            write(*,'(A,F12.6)') 'Varianza de la energía [eV^2]: ', boltzmann_variance_energy
            call flush(output_unit)
            call emit_level_header(1)
            min_entropy = 0.0_dp
            write(*,'(A,F12.6)') 'Entropia de los mínimos (g = 1) [eV/K]: ', min_entropy
            call emit_configuration_info(level, 1, energy, low_estimate, high_estimate)
            call write_exact_poscar(level, 1, config, total_sites)
            call flush(output_unit)

            summary_filename = 'sod_boltzmann_exact.txt'
            inquire(file=summary_filename, exist=summary_exists)
            open(newunit=unit_summary, file=summary_filename, status='unknown', position='append', action='write', iostat=ios_summary)
            if (ios_summary == 0) then
                if (.not. summary_exists) then
                    write(unit_summary,'(A)') '#Nivel DeltaF_300_kJmol Media_eV Media_ladoSi_eV Media_ladoGe_eV Min_ladoSi_eV Min_ladoGe_eV Varianza_eV2 Degeneracion_total Configs_inequivalentes Entropia_eV_perK TS_300_eV F_300_eV TS_800_eV F_800_eV TS_1200_eV F_1200_eV'
                end if
                write(unit_summary,'(I6,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,I0,1X,I0,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8)') level, deltaF_quartz_300, boltzmann_mean_energy, mean_low_all, mean_high_all, min_low_energy, min_high_energy, boltzmann_variance_energy, total_degeneracy_weight, unique_count, entropy_total, ts_terms(1), free_energy_est(1), ts_terms(2), free_energy_est(2), ts_terms(3), free_energy_est(3)
                close(unit_summary)
            else
                write(*,'(A,I0,A)') 'Aviso: no se pudo escribir resumen en txt (iostat=', ios_summary, ')'
                call flush(output_unit)
            end if
            return
        end if

        allocate(subset(level))
        allocate(canonical_subset(level))
        do i = 1, level
            subset(i) = i
        end do

        capacity = max(8, level)
        allocate(best_subsets(level, capacity))
        allocate(best_low(capacity))
        allocate(best_high(capacity))
        allocate(best_deg(capacity))
    best_low = 0.0_dp
    best_high = 0.0_dp
    best_deg = 0
        best_count = 0
        best_energy = huge(1.0_dp)

        unique_capacity = max(8, level)
        allocate(unique_subsets(level, unique_capacity))
        allocate(unique_energy(unique_capacity))
        allocate(unique_low(unique_capacity))
        allocate(unique_high(unique_capacity))
        allocate(unique_deg(unique_capacity))
    unique_energy = sentinel_energy
    unique_low = sentinel_energy
    unique_high = sentinel_energy
        unique_deg = 0
        unique_count = 0
        config = 1

        do
            call canonicalize_subset(subset, level, eqmatrix_local, nop_local, canonical_subset)
            subset_idx = find_subset_index(canonical_subset, level, unique_subsets, unique_count)
            if (subset_idx > 0) then
                unique_deg(subset_idx) = unique_deg(subset_idx) + 1
            else
                unique_count = unique_count + 1
                call ensure_unique_capacity(level, unique_subsets, unique_energy, unique_low, unique_high, unique_deg, unique_capacity, unique_count)
                unique_subsets(1:level, unique_count) = canonical_subset(1:level)
                unique_energy(unique_count) = sentinel_energy
                unique_low(unique_count) = sentinel_energy
                unique_high(unique_count) = sentinel_energy
                unique_deg(unique_count) = 1
            end if

            if (.not. next_combination(subset, total_sites)) exit
        end do

        write(*,'(A,I0)') 'Configuraciones unicas evaluadas: ', unique_count
        call flush(output_unit)

        degeneracy_energy_sum = 0.0_dp
        weighted_low_sum = 0.0_dp
        weighted_high_sum = 0.0_dp
        boltzmann_weight_sum = 0.0_dp
    boltzmann_energy_sum = 0.0_dp
    boltzmann_energy_sq_sum = 0.0_dp
    boltzmann_low_weight_sum = 0.0_dp
    boltzmann_high_weight_sum = 0.0_dp
    boltzmann_low_energy_sum = 0.0_dp
    boltzmann_high_energy_sum = 0.0_dp
        total_degeneracy_weight = 0
        total_low_weight = 0
        total_high_weight = 0
        min_low_energy = huge_marker
        min_high_energy = huge_marker
        lowest_energy_in_level = huge_marker
        do idx = 1, unique_count
            if (unique_energy(idx) == sentinel_energy) then
                config = 1
                config(unique_subsets(1:level, idx)) = 2
                call calculate_structure_energy(config, total_sites, energy, low_estimate, high_estimate)
                unique_energy(idx) = energy
                unique_low(idx) = low_estimate
                unique_high(idx) = high_estimate
            end if
            total_degeneracy_weight = total_degeneracy_weight + unique_deg(idx)
            degeneracy_weight = real(unique_deg(idx), dp)
            degeneracy_energy_sum = degeneracy_energy_sum + degeneracy_weight * unique_energy(idx)
            if (unique_energy(idx) < lowest_energy_in_level) lowest_energy_in_level = unique_energy(idx)
            if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                total_low_weight = total_low_weight + unique_deg(idx)
                weighted_low_sum = weighted_low_sum + degeneracy_weight * unique_low(idx)
                if (unique_low(idx) < min_low_energy) min_low_energy = unique_low(idx)
            end if
            if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                total_high_weight = total_high_weight + unique_deg(idx)
                weighted_high_sum = weighted_high_sum + degeneracy_weight * unique_high(idx)
                if (unique_high(idx) < min_high_energy) min_high_energy = unique_high(idx)
            end if
        end do

        if (total_degeneracy_weight > 0) then
            degeneracy_mean_energy = degeneracy_energy_sum / real(total_degeneracy_weight, dp)
        else
            degeneracy_mean_energy = 0.0_dp
        end if

        if (total_low_weight > 0) then
            mean_low_all = weighted_low_sum / real(total_low_weight, dp)
        else
            mean_low_all = degeneracy_mean_energy
        end if
        if (total_high_weight > 0) then
            mean_high_all = weighted_high_sum / real(total_high_weight, dp)
        else
            mean_high_all = degeneracy_mean_energy
        end if

        if (total_degeneracy_weight > 0) then
            inverse_boltzmann_temperature = 1.0_dp / (kB_eVk * boltzmann_temperature)
            boltzmann_reference_energy = min(lowest_energy_in_level, degeneracy_mean_energy)
            do idx = 1, unique_count
                energy = unique_energy(idx)
                degeneracy_weight = real(unique_deg(idx), dp)
                boltzmann_factor = exp(-inverse_boltzmann_temperature * (energy - boltzmann_reference_energy))
                boltzmann_weight_sum = boltzmann_weight_sum + degeneracy_weight * boltzmann_factor
                boltzmann_energy_sum = boltzmann_energy_sum + degeneracy_weight * boltzmann_factor * energy
                boltzmann_energy_sq_sum = boltzmann_energy_sq_sum + degeneracy_weight * boltzmann_factor * energy * energy
                if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                    boltzmann_low_weight_sum = boltzmann_low_weight_sum + degeneracy_weight * boltzmann_factor
                    boltzmann_low_energy_sum = boltzmann_low_energy_sum + degeneracy_weight * boltzmann_factor * unique_low(idx)
                end if
                if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                    boltzmann_high_weight_sum = boltzmann_high_weight_sum + degeneracy_weight * boltzmann_factor
                    boltzmann_high_energy_sum = boltzmann_high_energy_sum + degeneracy_weight * boltzmann_factor * unique_high(idx)
                end if
            end do
            if (boltzmann_weight_sum > 0.0_dp) then
                boltzmann_mean_energy = boltzmann_energy_sum / boltzmann_weight_sum
                boltzmann_variance_energy = max(0.0_dp, boltzmann_energy_sq_sum / boltzmann_weight_sum - boltzmann_mean_energy * boltzmann_mean_energy)
            else
                boltzmann_mean_energy = degeneracy_mean_energy
                boltzmann_variance_energy = 0.0_dp
            end if
            if (boltzmann_low_weight_sum > 0.0_dp) mean_low_all = boltzmann_low_energy_sum / boltzmann_low_weight_sum
            if (boltzmann_high_weight_sum > 0.0_dp) mean_high_all = boltzmann_high_energy_sum / boltzmann_high_weight_sum
        else
            boltzmann_mean_energy = 0.0_dp
            boltzmann_variance_energy = 0.0_dp
        end if

        if (min_low_energy >= huge_marker) min_low_energy = mean_low_all
        if (min_high_energy >= huge_marker) min_high_energy = mean_high_all

        do idx = 1, unique_count
            call update_best_structures(level, unique_subsets(1:level, idx), unique_energy(idx), unique_low(idx), unique_high(idx), unique_deg(idx), tol, best_energy, best_subsets, best_low, best_high, best_deg, best_count, capacity, eqmatrix_local, nop_local)
        end do

    write(*,'(A,F12.6)') 'Energía mínima (eV): ', best_energy
    write(*,'(A,F12.6)') 'Energía media (todos los inequivalentes) [eV]: ', boltzmann_mean_energy
    write(*,'(A,F12.6)') 'Varianza de la energía [eV^2]: ', boltzmann_variance_energy
        call flush(output_unit)
        call emit_level_header(best_count)
        total_min_degeneracy = 0
        do idx = 1, best_count
            total_min_degeneracy = total_min_degeneracy + best_deg(idx)
        end do
        write(*,'(A,I0)') 'Degeneracion total de los minimos: ', total_min_degeneracy
        if (total_min_degeneracy > 0) then
            min_entropy = kB_eVk * log(real(total_min_degeneracy, dp))
        else
            min_entropy = 0.0_dp
        end if
        if (total_min_degeneracy <= 1) then
            write(*,'(A,F12.6)') 'Entropia de los mínimos (g = 1) [eV/K]: ', 0.0_dp
        else
            write(*,'(A,F12.6)') 'Entropia de los mínimos (g = g_min) [eV/K]: ', min_entropy
        end if
        call flush(output_unit)

        do idx = 1, best_count
            config = 1
            config(best_subsets(1:level, idx)) = 2
            call emit_configuration_info(level, idx, best_energy, best_low(idx), best_high(idx), best_subsets(1:level, idx), best_deg(idx))
            call write_exact_poscar(level, idx, config, total_sites)
        end do

        if (total_degeneracy_weight > 0) then
            entropy_total = kB_eVk * log(real(total_degeneracy_weight, dp))
        else
            entropy_total = 0.0_dp
        end if
        ts_terms = entropy_total * temp_targets
        free_energy_est = boltzmann_mean_energy - ts_terms
        deltaF_quartz_300 = quartz_relative(level, total_sites, free_energy_est(1))

        summary_filename = 'sod_boltzmann_exact.txt'
        inquire(file=summary_filename, exist=summary_exists)
        open(newunit=unit_summary, file=summary_filename, status='unknown', position='append', action='write', iostat=ios_summary)
        if (ios_summary == 0) then
            if (.not. summary_exists) then
                write(unit_summary,'(A)') '#Nivel DeltaF_300_kJmol Media_eV Media_ladoSi_eV Media_ladoGe_eV Min_ladoSi_eV Min_ladoGe_eV Varianza_eV2 Degeneracion_total Configs_inequivalentes Entropia_eV_perK TS_300_eV F_300_eV TS_800_eV F_800_eV TS_1200_eV F_1200_eV'
            end if
            write(unit_summary,'(I6,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,I0,1X,I0,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8,1X,F18.8)') level, deltaF_quartz_300, boltzmann_mean_energy, mean_low_all, mean_high_all, min_low_energy, min_high_energy, boltzmann_variance_energy, total_degeneracy_weight, unique_count, entropy_total, ts_terms(1), free_energy_est(1), ts_terms(2), free_energy_est(2), ts_terms(3), free_energy_est(3)
            close(unit_summary)
        else
            write(*,'(A,I0,A)') 'Aviso: no se pudo escribir resumen en txt (iostat=', ios_summary, ')'
            call flush(output_unit)
        end if

    deallocate(unique_subsets)
    deallocate(unique_energy)
    deallocate(unique_low)
    deallocate(unique_high)
    deallocate(unique_deg)
    deallocate(best_subsets)
    deallocate(best_low)
    deallocate(best_high)
    deallocate(best_deg)
    deallocate(subset)
    deallocate(canonical_subset)
    if (allocated(eqmatrix_local)) deallocate(eqmatrix_local)
    end subroutine process_level_exact

    subroutine emit_level_header(count)
        integer, intent(in) :: count
        write(*,'(A,I0)') 'Configuraciones con energía mínima: ', count
        call flush(output_unit)
    end subroutine emit_level_header

    subroutine emit_configuration_info(level, index, energy, low_val, high_val, subset, degeneracy)
        integer, intent(in) :: level, index
        real(dp), intent(in) :: energy, low_val, high_val
        integer, intent(in), optional :: subset(:)
        integer, intent(in), optional :: degeneracy
        real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
        write(*,'(A,I0)') '  Configuración ', index

        if (level > 0 .and. present(subset)) then
            write(*,'(A)', advance='no') '    Sitios Ge: '
            write(*,'( *(1X,I0) )') subset
        else
            write(*,'(A)') '    Configuración base sin sustituciones.'
        end if
        if (present(degeneracy)) then
            write(*,'(A,I0)') '    Degeneración (g): ', degeneracy
        end if
        write(*,'(A,F12.6)') '    Energía estimada (eV): ', energy
        if (abs(low_val) < huge_marker) then
            write(*,'(A,F12.6)') '    Estimación lado Si (eV): ', low_val
        end if
        if (abs(high_val) < huge_marker) then
            write(*,'(A,F12.6)') '    Estimación lado Ge (eV): ', high_val
        end if
        call flush(output_unit)
    end subroutine emit_configuration_info

    subroutine write_exact_poscar(level, index, config, total_sites)
        integer, intent(in) :: level, index, total_sites
        integer, intent(in) :: config(:)
        character(len=64) :: filename

        write(filename,'("POSCAR_EXACT_N",I4.4,"_",I4.4,".vasp")') level, index
        call write_vasp_file(config, total_sites, trim(filename))
        write(*,'(A)') '    POSCAR guardado en '//trim(filename)
        call flush(output_unit)
    end subroutine write_exact_poscar

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

    subroutine update_best_structures(level, subset, energy, low_estimate, high_estimate, degeneracy, tol, best_energy, best_subsets, best_low, best_high, best_deg, best_count, capacity, eqmatrix, nop)
        integer, intent(in) :: level, degeneracy, nop
        integer, intent(in) :: subset(:)
        integer, intent(in) :: eqmatrix(:,:)
        real(dp), intent(in) :: energy, low_estimate, high_estimate, tol
        real(dp), intent(inout) :: best_energy
        integer, intent(inout) :: best_count, capacity
        integer, allocatable, intent(inout) :: best_subsets(:,:)
        real(dp), allocatable, intent(inout) :: best_low(:), best_high(:)
        integer, allocatable, intent(inout) :: best_deg(:)
        real(dp) :: diff
        logical :: is_equivalent

        diff = energy - best_energy
        if (best_count == 0 .or. diff < -tol) then
            best_energy = energy
            best_count = 1
            call ensure_capacity(level, best_subsets, best_low, best_high, best_deg, capacity, best_count)
            best_subsets(1:level,1) = subset(1:level)
            best_low(1) = low_estimate
            best_high(1) = high_estimate
            best_deg(1) = degeneracy
        else if (abs(diff) <= tol) then
            ! Check if this configuration is symmetry-equivalent to any already stored
            is_equivalent = is_symmetry_equivalent(subset, level, best_subsets, best_count, eqmatrix, nop)
            if (.not. is_equivalent) then
                best_count = best_count + 1
                call ensure_capacity(level, best_subsets, best_low, best_high, best_deg, capacity, best_count)
                best_subsets(1:level,best_count) = subset(1:level)
                best_low(best_count) = low_estimate
                best_high(best_count) = high_estimate
                best_deg(best_count) = degeneracy
            end if
        end if
    end subroutine update_best_structures

    logical function is_symmetry_equivalent(subset, level, best_subsets, best_count, eqmatrix, nop)
        ! Check if subset is equivalent to any of the stored best_subsets under symmetry operations
        integer, intent(in) :: level, best_count, nop
        integer, intent(in) :: subset(:)
        integer, intent(in) :: best_subsets(:,:)
        integer, intent(in) :: eqmatrix(:,:)
        integer :: idx, op, i, j
        integer :: mapped_subset(level), sorted_mapped(level), sorted_stored(level)
        logical :: match

        is_symmetry_equivalent = .false.

        ! Check against each stored configuration
        do idx = 1, best_count
            ! Try all symmetry operations
            do op = 1, nop
                ! Apply symmetry operation to current subset
                do i = 1, level
                    mapped_subset(i) = eqmatrix(op, subset(i))
                end do

                ! Sort both arrays for comparison
                sorted_mapped = mapped_subset
                call sort_int_ascending(sorted_mapped, level)
                sorted_stored(1:level) = best_subsets(1:level, idx)
                call sort_int_ascending(sorted_stored, level)

                ! Compare
                match = .true.
                do j = 1, level
                    if (sorted_mapped(j) /= sorted_stored(j)) then
                        match = .false.
                        exit
                    end if
                end do

                if (match) then
                    is_symmetry_equivalent = .true.
                    return
                end if
            end do
        end do
    end function is_symmetry_equivalent

    subroutine ensure_capacity(level, best_subsets, best_low, best_high, best_deg, capacity, required)
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: best_subsets(:,:)
        real(dp), allocatable, intent(inout) :: best_low(:), best_high(:)
        integer, allocatable, intent(inout) :: best_deg(:)
        integer :: new_capacity, old_capacity
        integer, allocatable :: tmp_subsets(:,:)
        real(dp), allocatable :: tmp(:)
        integer, allocatable :: tmp_deg(:)

        if (required <= capacity) return
        old_capacity = capacity
        new_capacity = max(required, max(2*max(1,old_capacity), 8))
        if (level > 0) then
            allocate(tmp_subsets(level, new_capacity))
            if (old_capacity > 0) tmp_subsets(:,1:old_capacity) = best_subsets(:,1:old_capacity)
            if (allocated(best_subsets)) deallocate(best_subsets)
            call move_alloc(tmp_subsets, best_subsets)
        end if
        allocate(tmp(new_capacity))
        if (old_capacity > 0) tmp(1:old_capacity) = best_low(1:old_capacity)
        if (allocated(best_low)) deallocate(best_low)
        call move_alloc(tmp, best_low)
        allocate(tmp(new_capacity))
        if (old_capacity > 0) tmp(1:old_capacity) = best_high(1:old_capacity)
        if (allocated(best_high)) deallocate(best_high)
        call move_alloc(tmp, best_high)
        allocate(tmp_deg(new_capacity))
        if (old_capacity > 0) tmp_deg(1:old_capacity) = best_deg(1:old_capacity)
        if (allocated(best_deg)) deallocate(best_deg)
        call move_alloc(tmp_deg, best_deg)
        capacity = new_capacity
    end subroutine ensure_capacity

    subroutine ensure_unique_capacity(level, unique_subsets, unique_energy, unique_low, unique_high, unique_deg, capacity, required)
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: unique_subsets(:,:)
        real(dp), allocatable, intent(inout) :: unique_energy(:), unique_low(:), unique_high(:)
        integer, allocatable, intent(inout) :: unique_deg(:)
        integer :: old_capacity, new_capacity
        integer, allocatable :: tmp_subsets(:,:)
        real(dp), allocatable :: tmp_real(:)
        integer, allocatable :: tmp_int(:)

        if (required <= capacity) return
        old_capacity = capacity
        new_capacity = max(required, max(2*max(1,old_capacity), 8))

        allocate(tmp_subsets(level, new_capacity))
        if (old_capacity > 0) tmp_subsets(:,1:old_capacity) = unique_subsets(:,1:old_capacity)
        if (allocated(unique_subsets)) deallocate(unique_subsets)
        call move_alloc(tmp_subsets, unique_subsets)

        allocate(tmp_real(new_capacity))
        if (old_capacity > 0) tmp_real(1:old_capacity) = unique_energy(1:old_capacity)
        if (allocated(unique_energy)) deallocate(unique_energy)
        call move_alloc(tmp_real, unique_energy)

        allocate(tmp_real(new_capacity))
        if (old_capacity > 0) tmp_real(1:old_capacity) = unique_low(1:old_capacity)
        if (allocated(unique_low)) deallocate(unique_low)
        call move_alloc(tmp_real, unique_low)

        allocate(tmp_real(new_capacity))
        if (old_capacity > 0) tmp_real(1:old_capacity) = unique_high(1:old_capacity)
        if (allocated(unique_high)) deallocate(unique_high)
        call move_alloc(tmp_real, unique_high)

        allocate(tmp_int(new_capacity))
        if (old_capacity > 0) tmp_int(1:old_capacity) = unique_deg(1:old_capacity)
        if (allocated(unique_deg)) deallocate(unique_deg)
        call move_alloc(tmp_int, unique_deg)

        capacity = new_capacity
    end subroutine ensure_unique_capacity

    subroutine canonicalize_subset(subset, level, eqmatrix, nop, canonical)
        integer, intent(in) :: subset(:)
        integer, intent(in) :: level, nop
        integer, intent(in) :: eqmatrix(:,:)
        integer, intent(out) :: canonical(:)
        integer :: temp(level)
        integer :: op

        canonical(1:level) = subset(1:level)
        call sort_int_ascending(canonical, level)

        do op = 1, nop
            temp(1:level) = eqmatrix(op, subset(1:level))
            call sort_int_ascending(temp, level)
            if (lexicographically_less(temp, canonical, level)) then
                canonical(1:level) = temp(1:level)
            end if
        end do
    end subroutine canonicalize_subset

    logical function lexicographically_less(a, b, level)
        integer, intent(in) :: a(:), b(:)
        integer, intent(in) :: level
        integer :: i

        do i = 1, level
            if (a(i) < b(i)) then
                lexicographically_less = .true.
                return
            else if (a(i) > b(i)) then
                lexicographically_less = .false.
                return
            end if
        end do
        lexicographically_less = .false.
    end function lexicographically_less

    integer function find_subset_index(subset, level, stored_subsets, stored_count)
        integer, intent(in) :: subset(:)
        integer, intent(in) :: level, stored_count
        integer, intent(in) :: stored_subsets(:,:)
        integer :: idx, j
        logical :: match

        find_subset_index = 0
        do idx = 1, stored_count
            match = .true.
            do j = 1, level
                if (stored_subsets(j, idx) /= subset(j)) then
                    match = .false.
                    exit
                end if
            end do
            if (match) then
                find_subset_index = idx
                return
            end if
        end do
    end function find_subset_index

end program sod_boltzmann_exact
