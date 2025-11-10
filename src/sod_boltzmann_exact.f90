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
        implicit none
        integer, intent(out) :: level_min, level_max
        real(dp), intent(out) :: tol
        integer :: argc, ios, colon_pos, iarg
        character(len=256) :: carg, range_spec, tol_arg
        logical :: found_n_flag

        level_min = 0
        level_max = -1
        tol = 1.0e-6_dp
        found_n_flag = .false.

        argc = command_argument_count()

        if (argc >= 1) then
            call get_command_argument(1, carg)
            if (is_help_token(carg)) then
                call print_usage_exact()
                stop 0
            end if
        end if

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

                colon_pos = index(range_spec, ':')
                if (colon_pos > 0) then
                    read(range_spec(1:colon_pos-1), *, iostat=ios) level_min
                    read(range_spec(colon_pos+1:), *, iostat=ios) level_max
                    if (ios /= 0) then
                        write(*,'(A)') 'Error: formato de rango inválido. Use N:M (ej: 1:12).'
                        call print_usage_exact()
                        stop 1
                    end if
                else
                    read(range_spec, *, iostat=ios) level_min
                    if (ios /= 0) then
                        write(*,'(A)') 'Error: especificación de nivel inválida.'
                        call print_usage_exact()
                        stop 1
                    end if
                    if (level_min < 0) then
                        level_min = 0
                        level_max = -1
                    else
                        level_max = level_min
                    end if
                end if
                exit
            end if
        end do

        if (argc >= 1) then
            if (found_n_flag) then
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
                if (argc >= 1) then
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
        end if
    end subroutine parse_arguments_exact

    subroutine print_usage_exact()
        implicit none
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

    pure logical function is_help_token(raw)
        implicit none
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
        implicit none
        integer, intent(in) :: level, total_sites
        integer, intent(inout) :: config(:)
        real(dp), intent(in) :: tol

        integer(ip) :: total_comb
        real(dp) :: energy, low_estimate, high_estimate
        integer :: idx
        integer :: i, nop_local
        integer :: unique_capacity, unique_count, subset_idx
        integer :: total_degeneracy_weight
    integer, allocatable :: subset(:)
    integer, allocatable :: best_subsets_si(:,:), best_subsets_ge(:,:)
    integer, allocatable :: canonical_subset(:)
    integer, allocatable :: unique_subsets(:,:)
    integer, allocatable :: unique_deg(:)
    integer, allocatable :: best_deg_si(:), best_deg_ge(:)
    real(dp), allocatable :: best_values_si(:), best_values_ge(:)
    real(dp), allocatable :: unique_low(:), unique_high(:)
    real(dp), allocatable :: unique_low_contrib(:,:), unique_high_contrib(:,:)
    integer, allocatable :: tmp_subsets(:,:)
    integer, allocatable :: tmp_deg(:)
    real(dp), allocatable :: tmp_low(:), tmp_high(:)
    integer :: unit_summary, ios_summary
    logical :: summary_exists
    logical :: allow_low_estimate, allow_high_estimate
    logical :: has_low_data, has_high_data
    character(len=64) :: summary_filename
    integer, allocatable :: eqmatrix_local(:,:)
    character(len=80), parameter :: separator = '------------------------------------------------------------------------'
    real(dp), parameter :: sentinel_energy = huge(1.0_dp)
    real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
    real(dp), parameter :: temp_targets(3) = (/300.0_dp, 800.0_dp, 1200.0_dp/)
    real(dp), parameter :: boltzmann_temperature = 300.0_dp
        real(dp) :: entropy_total
        real(dp) :: max_entropy, ideal_entropy, x, min_entropy
    real(dp) :: mean_low_all, mean_high_all
    real(dp) :: variance_low_all, variance_high_all
    real(dp) :: weighted_low_sum, weighted_high_sum
    real(dp) :: weighted_low_sq_sum, weighted_high_sq_sum
        integer :: total_low_weight, total_high_weight
        real(dp) :: min_low_energy, min_high_energy
        real(dp) :: degeneracy_weight
        real(dp) :: boltzmann_low_weight_sum, boltzmann_high_weight_sum
        real(dp) :: boltzmann_low_energy_sum, boltzmann_high_energy_sum
        real(dp) :: boltzmann_low_energy_sq_sum, boltzmann_high_energy_sq_sum
        real(dp) :: boltzmann_reference_low, boltzmann_reference_high
        real(dp) :: boltzmann_factor
        real(dp) :: boltzmann_mean_low, boltzmann_variance_low
        real(dp) :: boltzmann_mean_high, boltzmann_variance_high
        real(dp) :: entropy_low, entropy_high
        real(dp) :: ts_low(3), ts_high(3)
        real(dp) :: free_energy_low(3), free_energy_high(3)
        real(dp) :: deltaF_quartz_si, deltaF_quartz_ge
        integer :: capacity_si, capacity_ge
    integer :: best_count_si, best_count_ge
        real(dp) :: best_energy_si, best_energy_ge
        integer :: hole_count
        logical :: need_low_calibration, need_high_calibration

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
        write(*,'(A,F16.6)') 'Entropia máxima (k_B ln W) [eV/K]: ', max_entropy

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
        write(*,'(A,F16.6)') 'Entropia ideal (mezcla binaria) [eV/K]: ', ideal_entropy
        call flush(output_unit)

        hole_count = total_sites - level
        allow_low_estimate = .true.
        allow_high_estimate = .true.
        if (hole_count <= get_max_high_order()) allow_low_estimate = .false.
        if (level <= get_max_low_order()) allow_high_estimate = .false.

        if (level == 0) then
            config = 1
            call calculate_structure_energy(config, total_sites, energy, low_estimate, high_estimate)
            if (.not. allow_low_estimate) low_estimate = sentinel_energy
            if (.not. allow_high_estimate) high_estimate = sentinel_energy

            total_degeneracy_weight = 1

            total_low_weight = 0
            total_high_weight = 0
            mean_low_all = 0.0_dp
            mean_high_all = 0.0_dp
            variance_low_all = 0.0_dp
            variance_high_all = 0.0_dp
            min_low_energy = sentinel_energy
            min_high_energy = sentinel_energy
            boltzmann_mean_low = 0.0_dp
            boltzmann_mean_high = 0.0_dp
            boltzmann_variance_low = 0.0_dp
            boltzmann_variance_high = 0.0_dp
            entropy_total = 0.0_dp
            entropy_low = 0.0_dp
            entropy_high = 0.0_dp
            ts_low = 0.0_dp
            ts_high = 0.0_dp
            free_energy_low = 0.0_dp
            free_energy_high = 0.0_dp
            deltaF_quartz_si = 0.0_dp
            deltaF_quartz_ge = 0.0_dp

            if (allow_low_estimate) then
                total_low_weight = 1
                mean_low_all = low_estimate
                min_low_energy = low_estimate
                boltzmann_mean_low = low_estimate
                free_energy_low = (/low_estimate, low_estimate, low_estimate/)
                deltaF_quartz_si = quartz_relative(level, total_sites, free_energy_low(1))
            end if

            if (allow_high_estimate) then
                total_high_weight = 1
                mean_high_all = high_estimate
                min_high_energy = high_estimate
                boltzmann_mean_high = high_estimate
                free_energy_high = (/high_estimate, high_estimate, high_estimate/)
                deltaF_quartz_ge = quartz_relative(level, total_sites, free_energy_high(1))
            end if

            has_low_data = allow_low_estimate
            has_high_data = allow_high_estimate

            allocate(tmp_subsets(0,1))
            allocate(tmp_deg(1))
            allocate(tmp_low(1))
            allocate(tmp_high(1))
            tmp_deg(1) = 1
            if (allow_low_estimate) then
                tmp_low(1) = low_estimate
            else
                tmp_low(1) = sentinel_energy
            end if
            if (allow_high_estimate) then
                tmp_high(1) = high_estimate
            else
                tmp_high(1) = sentinel_energy
            end if
            call write_level_outputs(level, total_sites, 1, tmp_subsets, tmp_deg, tmp_low, tmp_high)
            deallocate(tmp_subsets)
            deallocate(tmp_deg)
            deallocate(tmp_low)
            deallocate(tmp_high)

            allocate(best_subsets_si(0,1))
            allocate(best_subsets_ge(0,1))
            allocate(best_values_si(1))
            allocate(best_values_ge(1))
            allocate(best_deg_si(1))
            allocate(best_deg_ge(1))

            if (allow_low_estimate) then
                best_values_si(1) = low_estimate
                best_deg_si(1) = 1
                best_energy_si = low_estimate
                best_count_si = 1
            else
                best_values_si(1) = 0.0_dp
                best_deg_si(1) = 0
                best_energy_si = sentinel_energy
                best_count_si = 0
            end if

            if (allow_high_estimate) then
                best_values_ge(1) = high_estimate
                best_deg_ge(1) = 1
                best_energy_ge = high_estimate
                best_count_ge = 1
            else
                best_values_ge(1) = 0.0_dp
                best_deg_ge(1) = 0
                best_energy_ge = sentinel_energy
                best_count_ge = 0
            end if

            call emit_side_statistics('lado Si', 'SI', level, best_count_si, best_energy_si, best_subsets_si, best_values_si, best_deg_si, total_sites, config, min_low_energy, mean_low_all, variance_low_all, entropy_low, has_low_data)
            call emit_side_statistics('lado Ge', 'GE', level, best_count_ge, best_energy_ge, best_subsets_ge, best_values_ge, best_deg_ge, total_sites, config, min_high_energy, mean_high_all, variance_high_all, entropy_high, has_high_data)

            summary_filename = 'sod_boltzmann_exact.txt'
            inquire(file=summary_filename, exist=summary_exists)
            open(newunit=unit_summary, file=summary_filename, status='unknown', position='append', action='write', iostat=ios_summary)
            if (ios_summary == 0) then
                if (.not. summary_exists) then
                    write(unit_summary,'(A)') '#Nivel DeltaF_Si_300_kJmol DeltaF_Ge_300_kJmol Media_Si_eV Media_Ge_eV Min_Si_eV Min_Ge_eV Varianza_Si_eV2 Varianza_Ge_eV2 Degeneracion_total Configs_inequivalentes Entropia_total_eV_perK Entropia_Si_eV_perK Entropia_Ge_eV_perK TS_Si_300_eV F_Si_300_eV TS_Ge_300_eV F_Ge_300_eV TS_Si_800_eV F_Si_800_eV TS_Ge_800_eV F_Ge_800_eV TS_Si_1200_eV F_Si_1200_eV TS_Ge_1200_eV F_Ge_1200_eV'
                end if
                write(unit_summary,'(I6,1X,8F18.8,1X,I0,1X,I0,1X,3F18.8,1X,12F18.8)') &
                    level, deltaF_quartz_si, deltaF_quartz_ge, mean_low_all, mean_high_all, min_low_energy, min_high_energy, variance_low_all, variance_high_all, total_degeneracy_weight, 1, entropy_total, entropy_low, entropy_high, &
                    ts_low(1), free_energy_low(1), ts_high(1), free_energy_high(1), ts_low(2), free_energy_low(2), ts_high(2), free_energy_high(2), ts_low(3), free_energy_low(3), ts_high(3), free_energy_high(3)
                close(unit_summary)
            else
                write(*,'(A,I0,A)') 'Aviso: no se pudo escribir resumen en txt (iostat=', ios_summary, ')'
                call flush(output_unit)
            end if

            deallocate(best_subsets_si)
            deallocate(best_subsets_ge)
            deallocate(best_values_si)
            deallocate(best_values_ge)
            deallocate(best_deg_si)
            deallocate(best_deg_ge)
            return
        end if

        allocate(subset(level))
        allocate(canonical_subset(level))
        do i = 1, level
            subset(i) = i
        end do

        capacity_si = max(8, level)
        capacity_ge = max(8, level)
        allocate(best_subsets_si(level, capacity_si))
        allocate(best_subsets_ge(level, capacity_ge))
        allocate(best_values_si(capacity_si))
        allocate(best_values_ge(capacity_ge))
        allocate(best_deg_si(capacity_si))
        allocate(best_deg_ge(capacity_ge))
        best_values_si = 0.0_dp
        best_values_ge = 0.0_dp
        best_deg_si = 0
        best_deg_ge = 0
        best_count_si = 0
        best_count_ge = 0
        best_energy_si = huge(1.0_dp)
        best_energy_ge = huge(1.0_dp)

        unique_capacity = max(8, level)
        allocate(unique_subsets(level, unique_capacity))
        allocate(unique_low(unique_capacity))
        allocate(unique_high(unique_capacity))
        allocate(unique_deg(unique_capacity))
        allocate(unique_low_contrib(4, unique_capacity))
        allocate(unique_high_contrib(4, unique_capacity))
        unique_low = sentinel_energy
        unique_high = sentinel_energy
        unique_low_contrib = 0.0_dp
        unique_high_contrib = 0.0_dp
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
                call ensure_unique_capacity(level, unique_subsets, unique_low, unique_high, unique_deg, unique_low_contrib, unique_high_contrib, unique_capacity, unique_count)
                unique_subsets(1:level, unique_count) = canonical_subset(1:level)
                unique_low(unique_count) = sentinel_energy
                unique_high(unique_count) = sentinel_energy
                unique_deg(unique_count) = 1
                unique_low_contrib(:, unique_count) = 0.0_dp
                unique_high_contrib(:, unique_count) = 0.0_dp
            end if

            if (.not. next_combination(subset, total_sites)) exit
        end do

        write(*,'(A,I0)') 'Configuraciones unicas evaluadas: ', unique_count
        call flush(output_unit)

        do idx = 1, unique_count
            if (unique_low(idx) == sentinel_energy) then
                config = 1
                config(unique_subsets(1:level, idx)) = 2
                call calculate_structure_energy(config, total_sites, energy, energy_low_side=low_estimate, &
                 energy_high_side=high_estimate, low_contrib=unique_low_contrib(:, idx), &
                 high_contrib=unique_high_contrib(:, idx))
                if (.not. allow_low_estimate) then
                    low_estimate = sentinel_energy
                    unique_low_contrib(:, idx) = 0.0_dp
                end if
                if (.not. allow_high_estimate) then
                    high_estimate = sentinel_energy
                    unique_high_contrib(:, idx) = 0.0_dp
                end if
                unique_low(idx) = low_estimate
                unique_high(idx) = high_estimate
            end if
        end do

        hole_count = total_sites - level
        need_low_calibration = allow_low_estimate .and. (level > get_max_low_order())
        need_high_calibration = allow_high_estimate .and. (hole_count > get_max_high_order())

        if (need_low_calibration .and. need_high_calibration) then
            call calibrate_level_with_gulp(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, unique_deg, unique_low, unique_high, config, need_low_calibration, need_high_calibration)
        end if

        weighted_low_sum = 0.0_dp
        weighted_high_sum = 0.0_dp
        weighted_low_sq_sum = 0.0_dp
        weighted_high_sq_sum = 0.0_dp
        boltzmann_low_weight_sum = 0.0_dp
        boltzmann_high_weight_sum = 0.0_dp
        boltzmann_low_energy_sum = 0.0_dp
        boltzmann_high_energy_sum = 0.0_dp
        boltzmann_low_energy_sq_sum = 0.0_dp
        boltzmann_high_energy_sq_sum = 0.0_dp
        total_degeneracy_weight = 0
        total_low_weight = 0
        total_high_weight = 0
        min_low_energy = huge_marker
        min_high_energy = huge_marker

        do idx = 1, unique_count
            degeneracy_weight = real(unique_deg(idx), dp)
            total_degeneracy_weight = total_degeneracy_weight + unique_deg(idx)

            if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                total_low_weight = total_low_weight + unique_deg(idx)
                weighted_low_sum = weighted_low_sum + degeneracy_weight * unique_low(idx)
                weighted_low_sq_sum = weighted_low_sq_sum + degeneracy_weight * unique_low(idx) * unique_low(idx)
                if (unique_low(idx) < min_low_energy) min_low_energy = unique_low(idx)
            end if

            if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                total_high_weight = total_high_weight + unique_deg(idx)
                weighted_high_sum = weighted_high_sum + degeneracy_weight * unique_high(idx)
                weighted_high_sq_sum = weighted_high_sq_sum + degeneracy_weight * unique_high(idx) * unique_high(idx)
                if (unique_high(idx) < min_high_energy) min_high_energy = unique_high(idx)
            end if
        end do

        if (total_low_weight > 0) then
            mean_low_all = weighted_low_sum / real(total_low_weight, dp)
            variance_low_all = max(0.0_dp, weighted_low_sq_sum / real(total_low_weight, dp) - mean_low_all * mean_low_all)
        else
            mean_low_all = 0.0_dp
            variance_low_all = 0.0_dp
            min_low_energy = huge_marker
        end if

        if (total_high_weight > 0) then
            mean_high_all = weighted_high_sum / real(total_high_weight, dp)
            variance_high_all = max(0.0_dp, weighted_high_sq_sum / real(total_high_weight, dp) - mean_high_all * mean_high_all)
        else
            mean_high_all = 0.0_dp
            variance_high_all = 0.0_dp
            min_high_energy = huge_marker
        end if

        call write_level_outputs(level, total_sites, unique_count, unique_subsets, unique_deg, unique_low, unique_high)

        boltzmann_mean_low = mean_low_all
        boltzmann_variance_low = variance_low_all
        boltzmann_mean_high = mean_high_all
        boltzmann_variance_high = variance_high_all

        if (total_low_weight > 0) then
            boltzmann_reference_low = min_low_energy
            if (boltzmann_reference_low >= huge_marker) boltzmann_reference_low = mean_low_all
            boltzmann_low_weight_sum = 0.0_dp
            boltzmann_low_energy_sum = 0.0_dp
            boltzmann_low_energy_sq_sum = 0.0_dp
            do idx = 1, unique_count
                if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                    degeneracy_weight = real(unique_deg(idx), dp)
                    boltzmann_factor = exp(-(unique_low(idx) - boltzmann_reference_low) / (kB_eVk * boltzmann_temperature))
                    boltzmann_low_weight_sum = boltzmann_low_weight_sum + degeneracy_weight * boltzmann_factor
                    boltzmann_low_energy_sum = boltzmann_low_energy_sum + degeneracy_weight * boltzmann_factor * unique_low(idx)
                    boltzmann_low_energy_sq_sum = boltzmann_low_energy_sq_sum + degeneracy_weight * boltzmann_factor * unique_low(idx) * unique_low(idx)
                end if
            end do
            if (boltzmann_low_weight_sum > 0.0_dp) then
                boltzmann_mean_low = boltzmann_low_energy_sum / boltzmann_low_weight_sum
                boltzmann_variance_low = max(0.0_dp, boltzmann_low_energy_sq_sum / boltzmann_low_weight_sum - boltzmann_mean_low * boltzmann_mean_low)
            end if
        end if

        if (total_high_weight > 0) then
            boltzmann_reference_high = min_high_energy
            if (boltzmann_reference_high >= huge_marker) boltzmann_reference_high = mean_high_all
            boltzmann_high_weight_sum = 0.0_dp
            boltzmann_high_energy_sum = 0.0_dp
            boltzmann_high_energy_sq_sum = 0.0_dp
            do idx = 1, unique_count
                if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                    degeneracy_weight = real(unique_deg(idx), dp)
                    boltzmann_factor = exp(-(unique_high(idx) - boltzmann_reference_high) / (kB_eVk * boltzmann_temperature))
                    boltzmann_high_weight_sum = boltzmann_high_weight_sum + degeneracy_weight * boltzmann_factor
                    boltzmann_high_energy_sum = boltzmann_high_energy_sum + degeneracy_weight * boltzmann_factor * unique_high(idx)
                    boltzmann_high_energy_sq_sum = boltzmann_high_energy_sq_sum + degeneracy_weight * boltzmann_factor * unique_high(idx) * unique_high(idx)
                end if
            end do
            if (boltzmann_high_weight_sum > 0.0_dp) then
                boltzmann_mean_high = boltzmann_high_energy_sum / boltzmann_high_weight_sum
                boltzmann_variance_high = max(0.0_dp, boltzmann_high_energy_sq_sum / boltzmann_high_weight_sum - boltzmann_mean_high * boltzmann_mean_high)
            end if
        end if

        if (total_low_weight > 0) then
            entropy_low = kB_eVk * log(real(total_low_weight, dp))
            ts_low = entropy_low * temp_targets
            free_energy_low = boltzmann_mean_low - ts_low
            deltaF_quartz_si = quartz_relative(level, total_sites, free_energy_low(1))
        else
            entropy_low = 0.0_dp
            ts_low = 0.0_dp
            free_energy_low = 0.0_dp
            deltaF_quartz_si = 0.0_dp
        end if

        if (total_high_weight > 0) then
            entropy_high = kB_eVk * log(real(total_high_weight, dp))
            ts_high = entropy_high * temp_targets
            free_energy_high = boltzmann_mean_high - ts_high
            deltaF_quartz_ge = quartz_relative(level, total_sites, free_energy_high(1))
        else
            entropy_high = 0.0_dp
            ts_high = 0.0_dp
            free_energy_high = 0.0_dp
            deltaF_quartz_ge = 0.0_dp
        end if

        if (total_degeneracy_weight > 0) then
            entropy_total = kB_eVk * log(real(total_degeneracy_weight, dp))
        else
            entropy_total = 0.0_dp
        end if

        if (min_low_energy >= huge_marker) min_low_energy = mean_low_all
        if (min_high_energy >= huge_marker) min_high_energy = mean_high_all

        do idx = 1, unique_count
            if (unique_low(idx) /= sentinel_energy .and. abs(unique_low(idx)) < huge_marker) then
                call update_best_structures_side(level, unique_subsets(1:level, idx), unique_low(idx), unique_deg(idx), tol, best_energy_si, best_subsets_si, best_values_si, best_deg_si, best_count_si, capacity_si, eqmatrix_local, nop_local)
            end if
            if (unique_high(idx) /= sentinel_energy .and. abs(unique_high(idx)) < huge_marker) then
                call update_best_structures_side(level, unique_subsets(1:level, idx), unique_high(idx), unique_deg(idx), tol, best_energy_ge, best_subsets_ge, best_values_ge, best_deg_ge, best_count_ge, capacity_ge, eqmatrix_local, nop_local)
            end if
        end do

        has_low_data = (total_low_weight > 0)
        has_high_data = (total_high_weight > 0)

        call emit_side_statistics('lado Si', 'SI', level, best_count_si, best_energy_si, best_subsets_si, best_values_si, best_deg_si, total_sites, config, min_low_energy, mean_low_all, variance_low_all, entropy_low, has_low_data)
        call emit_side_statistics('lado Ge', 'GE', level, best_count_ge, best_energy_ge, best_subsets_ge, best_values_ge, best_deg_ge, total_sites, config, min_high_energy, mean_high_all, variance_high_all, entropy_high, has_high_data)

        summary_filename = 'sod_boltzmann_exact.txt'
        inquire(file=summary_filename, exist=summary_exists)
        open(newunit=unit_summary, file=summary_filename, status='unknown', position='append', action='write', iostat=ios_summary)
        if (ios_summary == 0) then
            if (.not. summary_exists) then
                write(unit_summary,'(A)') '#Nivel DeltaF_Si_300_kJmol DeltaF_Ge_300_kJmol Media_Si_eV Media_Ge_eV Min_Si_eV Min_Ge_eV Varianza_Si_eV2 Varianza_Ge_eV2 Degeneracion_total Configs_inequivalentes Entropia_total_eV_perK Entropia_Si_eV_perK Entropia_Ge_eV_perK TS_Si_300_eV F_Si_300_eV TS_Ge_300_eV F_Ge_300_eV TS_Si_800_eV F_Si_800_eV TS_Ge_800_eV F_Ge_800_eV TS_Si_1200_eV F_Si_1200_eV TS_Ge_1200_eV F_Ge_1200_eV'
            end if
            write(unit_summary,'(I6,1X,8F18.8,1X,I0,1X,I0,1X,3F18.8,1X,12F18.8)') &
                level, deltaF_quartz_si, deltaF_quartz_ge, mean_low_all, mean_high_all, min_low_energy, min_high_energy, variance_low_all, variance_high_all, total_degeneracy_weight, unique_count, entropy_total, entropy_low, entropy_high, &
                ts_low(1), free_energy_low(1), ts_high(1), free_energy_high(1), ts_low(2), free_energy_low(2), ts_high(2), free_energy_high(2), ts_low(3), free_energy_low(3), ts_high(3), free_energy_high(3)
            close(unit_summary)
        else
            write(*,'(A,I0,A)') 'Aviso: no se pudo escribir resumen en txt (iostat=', ios_summary, ')'
            call flush(output_unit)
        end if

    deallocate(unique_subsets)
    deallocate(unique_low)
    deallocate(unique_high)
    deallocate(unique_deg)
    deallocate(unique_low_contrib)
    deallocate(unique_high_contrib)
    deallocate(best_subsets_si)
    deallocate(best_subsets_ge)
    deallocate(best_values_si)
    deallocate(best_values_ge)
    deallocate(best_deg_si)
    deallocate(best_deg_ge)
    deallocate(subset)
    deallocate(canonical_subset)
    if (allocated(eqmatrix_local)) deallocate(eqmatrix_local)
    end subroutine process_level_exact

    subroutine emit_side_statistics(label, tag, level, count, best_energy, best_subsets, best_values, best_deg, total_sites, config, min_energy, mean_energy, variance_energy, entropy_side, has_data)
        implicit none
        character(len=*), intent(in) :: label
        character(len=*), intent(in) :: tag
        integer, intent(in) :: level, count, total_sites
        real(dp), intent(in) :: best_energy, min_energy, mean_energy, variance_energy, entropy_side
        integer, intent(inout) :: config(:)
        integer, intent(in) :: best_subsets(:,:)
        real(dp), intent(in) :: best_values(:)
        integer, intent(in) :: best_deg(:)
        logical, intent(in) :: has_data
        integer :: idx
        integer :: total_deg
        real(dp) :: entropy_min

        write(*,'(A)') trim(label)//' — Estadísticas'
        if (.not. has_data) then
            write(*,'(A)') '  Energía mínima (eV): -'
            write(*,'(A)') '  Energía media (eV): -'
            write(*,'(A)') '  Varianza de la energía [eV^2]: -'
            write(*,'(A)') '  Energía mínima usada para referencia (eV): -'
            write(*,'(A)') '  Configuraciones con energía mínima: 0'
            write(*,'(A)') '  Degeneración total de los mínimos: 0'
            write(*,'(A)') '  Entropia de los mínimos (eV/K): -'
            write(*,'(A)') '  Entropia total lado ('//trim(tag)//') (eV/K): -'
            call flush(output_unit)
            return
        end if

        if (count <= 0) then
            write(*,'(A)') '  No hay configuraciones medibles para este lado.'
            call flush(output_unit)
            return
        end if

        write(*,'(A,F16.6)') '  Energía mínima (eV): ', best_energy
        write(*,'(A,F16.6)') '  Energía media (eV): ', mean_energy
        write(*,'(A,F16.6)') '  Varianza de la energía [eV^2]: ', variance_energy
        write(*,'(A,F16.6)') '  Energía mínima usada para referencia (eV): ', min_energy

        total_deg = 0
        do idx = 1, count
            total_deg = total_deg + best_deg(idx)
        end do
        write(*,'(A,I0)') '  Configuraciones con energía mínima: ', count
        write(*,'(A,I0)') '  Degeneración total de los mínimos: ', total_deg
        if (total_deg > 0) then
            entropy_min = kB_eVk * log(real(total_deg, dp))
        else
            entropy_min = 0.0_dp
        end if
        write(*,'(A,F16.6)') '  Entropia de los mínimos (eV/K): ', entropy_min
        write(*,'(A,F16.6)') '  Entropia total lado ('//trim(tag)//') (eV/K): ', entropy_side
        call flush(output_unit)

        do idx = 1, count
            config = 1
            if (level > 0) config(best_subsets(1:level, idx)) = 2
            call emit_configuration_info_side(trim(label), idx, best_values(idx), best_subsets(1:level, idx), best_deg(idx))
            call write_exact_poscar_with_tag(level, idx, config, total_sites, trim(tag))
        end do
    end subroutine emit_side_statistics

    subroutine emit_configuration_info_side(label, index, energy, subset, degeneracy)
        implicit none
        character(len=*), intent(in) :: label
        integer, intent(in) :: index
        real(dp), intent(in) :: energy
        integer, intent(in) :: subset(:)
        integer, intent(in) :: degeneracy

        write(*,'(A,I0)') '  '//trim(label)//' — Configuración ', index
        if (size(subset) > 0) then
            write(*,'(A)', advance='no') '    Sitios Ge: '
            if (size(subset) > 0) write(*,'( *(1X,I0) )') subset
        else
            write(*,'(A)') '    Configuración base sin sustituciones.'
        end if
        write(*,'(A,I0)') '    Degeneración (g): ', degeneracy
        write(*,'(A,F16.6)') '    Energía corregida (eV): ', energy
        call flush(output_unit)
    end subroutine emit_configuration_info_side

    subroutine write_exact_poscar_with_tag(level, index, config, total_sites, tag)
        implicit none
        integer, intent(in) :: level, index, total_sites
        integer, intent(in) :: config(:)
        character(len=*), intent(in) :: tag
        character(len=64) :: filename

        write(filename,'("POSCAR_",A,"_N",I4.4,"_",I4.4,".vasp")') trim(tag), level, index
        call write_vasp_file(config, total_sites, trim(filename))
        write(*,'(A)') '    POSCAR guardado en '//trim(filename)
        call flush(output_unit)
    end subroutine write_exact_poscar_with_tag

    subroutine write_level_outputs(level, total_sites, unique_count, unique_subsets, unique_deg, unique_low, unique_high)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_subsets(:,:), unique_deg(:)
        real(dp), intent(in) :: unique_low(:), unique_high(:)
        character(len=64) :: outsod_name, energies_name
        character(len=32) :: low_field, high_field
        integer :: unit_outsod, unit_energy
        integer :: idx, site
        integer, allocatable :: subset(:)
        real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp

        write(outsod_name,'("OUTSOD_N",I4.4)') level
        write(energies_name,'("ENERGIES_N",I4.4)') level

        open(newunit=unit_outsod, file=trim(outsod_name), status='replace', action='write')
        open(newunit=unit_energy, file=trim(energies_name), status='replace', action='write')

        write(unit_outsod,'(I12,"  substitutions in",I12," sites")') level, total_sites
        write(unit_outsod,'(I12,"  configurations")') unique_count
        write(unit_energy,'(A)') '# idx degeneracy energy_Si_eV energy_Ge_eV'

        if (unique_count <= 0) then
            close(unit_outsod)
            close(unit_energy)
            return
        end if

        if (level > 0) allocate(subset(level))

        do idx = 1, unique_count
            if (level > 0) subset(1:level) = unique_subsets(1:level, idx)
            write(unit_outsod,'(I6,1X,I12)', advance='no') idx, unique_deg(idx)
            if (level > 0) then
                do site = 1, level
                    write(unit_outsod,'(1X,I6)', advance='no') subset(site)
                end do
            end if
            write(unit_outsod,*)

            if (abs(unique_low(idx)) < huge_marker) then
                write(low_field,'(F18.8)') unique_low(idx)
            else
                low_field = 'NaN'
            end if
            if (abs(unique_high(idx)) < huge_marker) then
                write(high_field,'(F18.8)') unique_high(idx)
            else
                high_field = 'NaN'
            end if

            write(unit_energy,'(I6,1X,I12,2(1X,A))') idx, unique_deg(idx), adjustl(low_field), adjustl(high_field)
        end do

        if (allocated(subset)) deallocate(subset)
        close(unit_outsod)
        close(unit_energy)
    end subroutine write_level_outputs

    pure real(dp) function quartz_relative(level, total_sites, energy) result(rel_val)
        implicit none
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
        implicit none
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

    subroutine update_best_structures_side(level, subset, energy, degeneracy, tol, best_energy, best_subsets, best_values, best_deg, best_count, capacity, eqmatrix, nop)
        implicit none
        integer, intent(in) :: level, degeneracy, nop
        integer, intent(in) :: subset(:)
        real(dp), intent(in) :: energy, tol
        real(dp), intent(inout) :: best_energy
        integer, allocatable, intent(inout) :: best_subsets(:,:)
        real(dp), allocatable, intent(inout) :: best_values(:)
        integer, allocatable, intent(inout) :: best_deg(:)
        integer, intent(inout) :: best_count, capacity
        integer, intent(in) :: eqmatrix(:,:)
        real(dp) :: diff
        logical :: equivalent

        diff = energy - best_energy

        if (best_count == 0 .or. diff < -tol) then
            best_energy = energy
            best_count = 1
            call ensure_side_capacity(level, best_subsets, best_values, best_deg, capacity, best_count)
            if (level > 0) best_subsets(1:level,1) = subset(1:level)
            best_values(1) = energy
            best_deg(1) = degeneracy
            return
        end if

        if (abs(diff) <= tol) then
            if (best_count > 0) then
                equivalent = is_symmetry_equivalent(subset, level, best_subsets, best_count, eqmatrix, nop)
                if (equivalent) return
            end if

            best_count = best_count + 1
            call ensure_side_capacity(level, best_subsets, best_values, best_deg, capacity, best_count)
            if (level > 0) best_subsets(1:level, best_count) = subset(1:level)
            best_values(best_count) = energy
            best_deg(best_count) = degeneracy
        end if
    end subroutine update_best_structures_side

    subroutine ensure_side_capacity(level, best_subsets, best_values, best_deg, capacity, required)
        implicit none
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: best_subsets(:,:)
        real(dp), allocatable, intent(inout) :: best_values(:)
        integer, allocatable, intent(inout) :: best_deg(:)
        integer :: old_capacity, new_capacity
        integer, allocatable :: tmp_subsets(:,:)
        real(dp), allocatable :: tmp_values(:)
        integer, allocatable :: tmp_int(:)

        if (required <= capacity) return

        old_capacity = capacity
        new_capacity = max(required, max(2*max(1, old_capacity), 8))

        if (level > 0) then
            allocate(tmp_subsets(level, new_capacity))
            if (old_capacity > 0) tmp_subsets(:,1:old_capacity) = best_subsets(:,1:old_capacity)
            if (allocated(best_subsets)) deallocate(best_subsets)
            call move_alloc(tmp_subsets, best_subsets)
        end if

        allocate(tmp_values(new_capacity))
        if (old_capacity > 0) tmp_values(1:old_capacity) = best_values(1:old_capacity)
        if (allocated(best_values)) deallocate(best_values)
        call move_alloc(tmp_values, best_values)

        allocate(tmp_int(new_capacity))
        if (old_capacity > 0) tmp_int(1:old_capacity) = best_deg(1:old_capacity)
        if (allocated(best_deg)) deallocate(best_deg)
        call move_alloc(tmp_int, best_deg)

        capacity = new_capacity
    end subroutine ensure_side_capacity

    logical function is_symmetry_equivalent(subset, level, best_subsets, best_count, eqmatrix, nop)
        implicit none
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
        implicit none
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

    subroutine ensure_unique_capacity(level, unique_subsets, unique_low, unique_high, unique_deg, unique_low_contrib, unique_high_contrib, capacity, required)
        implicit none
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: unique_subsets(:,:)
        real(dp), allocatable, intent(inout) :: unique_low(:), unique_high(:)
        integer, allocatable, intent(inout) :: unique_deg(:)
        real(dp), allocatable, intent(inout) :: unique_low_contrib(:,:), unique_high_contrib(:,:)
        integer :: old_capacity, new_capacity
        integer, allocatable :: tmp_subsets(:,:)
        real(dp), allocatable :: tmp_real(:)
        integer, allocatable :: tmp_int(:)
        real(dp), allocatable :: tmp_contrib(:,:)
        integer :: order_dim

        if (required <= capacity) return
        old_capacity = capacity
        new_capacity = max(required, max(2*max(1,old_capacity), 8))

        allocate(tmp_subsets(level, new_capacity))
        if (old_capacity > 0) tmp_subsets(:,1:old_capacity) = unique_subsets(:,1:old_capacity)
        if (allocated(unique_subsets)) deallocate(unique_subsets)
        call move_alloc(tmp_subsets, unique_subsets)

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

        order_dim = 4
        if (allocated(unique_low_contrib)) order_dim = ubound(unique_low_contrib,1)
        allocate(tmp_contrib(order_dim, new_capacity))
        if (allocated(unique_low_contrib) .and. old_capacity > 0) tmp_contrib(:,1:old_capacity) = unique_low_contrib(:,1:old_capacity)
        if (allocated(unique_low_contrib)) deallocate(unique_low_contrib)
        call move_alloc(tmp_contrib, unique_low_contrib)

        order_dim = 4
        if (allocated(unique_high_contrib)) order_dim = ubound(unique_high_contrib,1)
        allocate(tmp_contrib(order_dim, new_capacity))
        if (allocated(unique_high_contrib) .and. old_capacity > 0) tmp_contrib(:,1:old_capacity) = unique_high_contrib(:,1:old_capacity)
        if (allocated(unique_high_contrib)) deallocate(unique_high_contrib)
        call move_alloc(tmp_contrib, unique_high_contrib)

        capacity = new_capacity
    end subroutine ensure_unique_capacity

    logical function find_scripts_directory(out_dir)
        implicit none
        character(len=*), intent(out) :: out_dir
        character(len=512) :: env_dir
        integer :: lenv, i
        character(len=512), dimension(6) :: candidates

        out_dir = ''
        env_dir = ''
        lenv = 0
        call get_environment_variable('SOD_SCRIPTS', env_dir, length=lenv)
        if (lenv > 0) then
            env_dir = env_dir(1:lenv)
            if (directory_has_scripts(trim(env_dir))) then
                out_dir = trim(env_dir)
                find_scripts_directory = .true.
                return
            end if
        end if

    candidates = ''
    candidates(1) = 'scripts'
    candidates(2) = '../scripts'
    candidates(3) = '../../scripts'
    candidates(4) = '../../../scripts'
    candidates(5) = '../../../../scripts'
    candidates(6) = '/home/salvador/coding/sod/sod/scripts'
        do i = 1, size(candidates)
            if (directory_has_scripts(trim(candidates(i)))) then
                out_dir = trim(candidates(i))
                find_scripts_directory = .true.
                return
            end if
        end do
        find_scripts_directory = .false.
    end function find_scripts_directory

    logical function directory_has_scripts(dir)
        implicit none
        character(len=*), intent(in) :: dir
        logical :: exist_vasp, exist_run, exist_extract

        inquire(file=trim(dir)//'/vasp2gin.sh', exist=exist_vasp)
        inquire(file=trim(dir)//'/run_jobs.sh', exist=exist_run)
        inquire(file=trim(dir)//'/extract.sh', exist=exist_extract)
        directory_has_scripts = exist_vasp .AND. exist_run .AND. exist_extract
    end function directory_has_scripts

    subroutine select_calibration_indices(values, n, indices, ok)
        implicit none
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: n
        integer, intent(out) :: indices(:)
        logical, intent(out) :: ok
        integer, allocatable :: sorted_idx(:)
        integer :: s, m, pick
        real(dp) :: position

        if (n < size(indices)) then
            ok = .false.
            return
        end if

        allocate(sorted_idx(n))
        call sort_indices_by_values(values, sorted_idx, n)

        m = size(indices)
        do s = 1, m
            position = 1.0_dp + floor(real(n - 1, dp) * real(s - 1, dp) / real(m - 1, dp))
            pick = int(position)
            pick = max(1, min(n, pick))
            indices(s) = sorted_idx(pick)
        end do

        deallocate(sorted_idx)
        ok = .true.
    end subroutine select_calibration_indices

    subroutine sort_indices_by_values(values, indices, n)
        implicit none
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: indices(:)
        integer, intent(in) :: n
        integer :: i, j, tmp

        do i = 1, n
            indices(i) = i
        end do

        do i = 1, n - 1
            do j = i + 1, n
                if (values(indices(j)) < values(indices(i))) then
                    tmp = indices(i)
                    indices(i) = indices(j)
                    indices(j) = tmp
                end if
            end do
        end do
    end subroutine sort_indices_by_values

    subroutine copy_calibration_scripts(script_dir, calib_dir)
        implicit none
        character(len=*), intent(in) :: script_dir
        character(len=*), intent(in) :: calib_dir
        integer :: exit_code

        call execute_command_line('cp ' // trim(script_dir)//'/vasp2gin.sh ' // trim(calib_dir), exitstat=exit_code)
        call execute_command_line('cp ' // trim(script_dir)//'/run_jobs.sh ' // trim(calib_dir), exitstat=exit_code)
        call execute_command_line('cp ' // trim(script_dir)//'/extract.sh ' // trim(calib_dir), exitstat=exit_code)
        call execute_command_line('chmod +x ' // trim(calib_dir)//'/vasp2gin.sh ' // trim(calib_dir)//'/run_jobs.sh ' // trim(calib_dir)//'/extract.sh', exitstat=exit_code)
    end subroutine copy_calibration_scripts

    subroutine solve_normal_equations(mat, rhs, coeff, ok)
        implicit none
        real(dp), intent(inout) :: mat(:,:)
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: coeff(:)
        logical, intent(out) :: ok
        real(dp) :: aug(4,5)
        real(dp) :: pivot, factor
        integer :: i, j, k, pivot_row
        integer, parameter :: n = 4
        real(dp) :: temp_row(5)

        aug = 0.0_dp
        do i = 1, n
            aug(i,1:n) = mat(i,1:n)
            aug(i,n+1) = rhs(i)
        end do

        do k = 1, n
            pivot_row = k
            pivot = abs(aug(k,k))
            do i = k+1, n
                if (abs(aug(i,k)) > pivot) then
                    pivot = abs(aug(i,k))
                    pivot_row = i
                end if
            end do

            if (pivot < 1.0e-10_dp) then
                ok = .false.
                return
            end if

            if (pivot_row /= k) then
                temp_row = aug(k,:)
                aug(k,:) = aug(pivot_row,:)
                aug(pivot_row,:) = temp_row
            end if

            pivot = aug(k,k)
            aug(k,:) = aug(k,:) / pivot

            do i = 1, n
                if (i == k) cycle
                factor = aug(i,k)
                aug(i,:) = aug(i,:) - factor * aug(k,:)
            end do
        end do

        do i = 1, n
            coeff(i) = aug(i,n+1)
        end do
        ok = .true.
    end subroutine solve_normal_equations

    subroutine calibrate_level_with_gulp(level, total_sites, unique_count, unique_subsets, unique_low_contrib, unique_high_contrib, unique_deg, unique_low, unique_high, config, do_low_calibration, do_high_calibration)
        implicit none
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_subsets(:,:)
        real(dp), intent(in) :: unique_low_contrib(:,:), unique_high_contrib(:,:)
        integer, intent(in) :: unique_deg(:)
        real(dp), intent(inout) :: unique_low(:)
        real(dp), intent(inout) :: unique_high(:)
        integer, intent(inout) :: config(:)
        logical, intent(in) :: do_low_calibration, do_high_calibration

        integer, parameter :: n_orders = 4
        integer, parameter :: n_targets = n_orders + 1
        integer, parameter :: max_union = 2 * n_targets
        real(dp), parameter :: huge_marker = huge(1.0_dp) * 0.5_dp
        character(len=256) :: calib_dir, script_dir
        character(len=64) :: base_name
        character(len=16) :: level_tag
        character(len=128) :: gout_names(max_union)
        integer :: selected_low_idx(n_targets)
        integer :: selected_high_idx(n_targets)
        integer :: union_idx(max_union)
        real(dp) :: measured_energy(max_union)
        logical :: filled_union(max_union)
        integer :: i, j, k, ios, unit_file, exit_code
        integer :: union_count, pos_union
        integer :: n_valid_high
        integer :: idx_global
        integer, allocatable :: valid_high_idx(:), tmp_selected(:)
        real(dp), allocatable :: tmp_values(:)
        logical :: ok_scripts, ok_low_selection, ok_high_selection, calibrate_high, ok_system, already
        real(dp) :: base_energy_low, base_energy_high
        real(dp) :: mat_low(n_orders, n_orders), rhs_low(n_orders), coeff_low(n_orders)
        real(dp) :: mat_high(n_orders, n_orders), rhs_high(n_orders), coeff_high(n_orders)
        real(dp) :: energy_val, corrected, contrib_j
        character(len=512) :: command, line
        character(len=256) :: gout_label
        integer :: pos_hash, match_idx

    if (unique_count < n_targets) return

    if (.not. (do_low_calibration .and. do_high_calibration)) return

        selected_low_idx = 0
        selected_high_idx = 0

        ok_scripts = find_scripts_directory(script_dir)
        if (.not. ok_scripts) then
            write(*,'(A,I0)') 'Aviso: no se encontro directorio de scripts GULP para el nivel ', level
            return
        end if

        ok_low_selection = .false.
        if (do_low_calibration) then
            call select_calibration_indices(unique_low, unique_count, selected_low_idx, ok_low_selection)
            if (.not. ok_low_selection) return
        end if

        allocate(valid_high_idx(unique_count))
        n_valid_high = 0
        do i = 1, unique_count
            if (abs(unique_high(i)) < huge_marker) then
                n_valid_high = n_valid_high + 1
                valid_high_idx(n_valid_high) = i
            end if
        end do

        calibrate_high = .false.
        ok_high_selection = .false.
        if (do_high_calibration .and. n_valid_high >= n_targets) then
            allocate(tmp_values(n_valid_high))
            allocate(tmp_selected(n_targets))
            do i = 1, n_valid_high
                tmp_values(i) = unique_high(valid_high_idx(i))
            end do
            call select_calibration_indices(tmp_values, n_valid_high, tmp_selected, ok_high_selection)
            if (ok_high_selection) then
                calibrate_high = .true.
                do i = 1, n_targets
                    selected_high_idx(i) = valid_high_idx(tmp_selected(i))
                end do
            end if
            deallocate(tmp_values)
            deallocate(tmp_selected)
        end if
        if (allocated(valid_high_idx)) deallocate(valid_high_idx)

        write(level_tag,'(I4.4)') level
        write(calib_dir,'("gulp_calib_N",A)') adjustl(level_tag)
        calib_dir = trim(calib_dir)

        union_idx = 0
        union_count = 0
        if (do_low_calibration) then
            do i = 1, n_targets
                idx_global = selected_low_idx(i)
                if (idx_global <= 0) cycle
                already = .false.
                do j = 1, union_count
                    if (union_idx(j) == idx_global) then
                        already = .true.
                        exit
                    end if
                end do
                if (.not. already) then
                    if (union_count >= max_union) then
                        write(*,'(A,I0)') 'Aviso: capacidad insuficiente para union de calibracion en nivel ', level
                        return
                    end if
                    union_count = union_count + 1
                    union_idx(union_count) = idx_global
                end if
            end do
        end if
        if (calibrate_high) then
            do i = 1, n_targets
                idx_global = selected_high_idx(i)
                if (idx_global <= 0) cycle
                already = .false.
                do j = 1, union_count
                    if (union_idx(j) == idx_global) then
                        already = .true.
                        exit
                    end if
                end do
                if (.not. already) then
                    if (union_count >= max_union) then
                        write(*,'(A,I0)') 'Aviso: capacidad insuficiente para union de calibracion en nivel ', level
                        return
                    end if
                    union_count = union_count + 1
                    union_idx(union_count) = idx_global
                end if
            end do
        end if

        if (union_count == 0) return

        call execute_command_line('rm -rf ' // trim(calib_dir), exitstat=exit_code)
        call execute_command_line('mkdir -p ' // trim(calib_dir), exitstat=exit_code)

        gout_names = ''
        measured_energy = 0.0_dp
        filled_union = .false.

        do i = 1, union_count
            config = 1
            if (level > 0) config(unique_subsets(1:level, union_idx(i))) = 2
            write(base_name,'("calib_N",I4.4,"_c",I2.2)') level, i
            call write_vasp_file(config, total_sites, trim(calib_dir)//'/'//trim(base_name)//'.vasp')
            gout_names(i) = trim(base_name)//'.vasp.gout'
        end do

        call copy_calibration_scripts(script_dir, calib_dir)

        command = 'cd ' // trim(calib_dir) // ' && bash run_jobs.sh'
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A,I0)') 'Aviso: fallo run_jobs.sh durante calibracion en nivel ', level
            return
        end if

        command = 'cd ' // trim(calib_dir) // ' && bash extract.sh'
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A,I0)') 'Aviso: fallo extract.sh durante calibracion en nivel ', level
            return
        end if

        open(newunit=unit_file, file=trim(calib_dir)//'/ENERGIES', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,I0)') 'Aviso: no se pudo leer ENERGIES para el nivel ', level
            return
        end if

        do
            read(unit_file,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            pos_hash = index(line, '#')
            if (pos_hash <= 0) cycle
            read(line(1:pos_hash-1),*,iostat=ios) energy_val
            if (ios /= 0) cycle
            gout_label = adjustl(trim(line(pos_hash+1:)))
            gout_label = trim(gout_label)
            match_idx = 0
            do i = 1, union_count
                if (trim(gout_label) == trim(gout_names(i))) then
                    match_idx = i
                    exit
                end if
            end do
            if (match_idx > 0) then
                measured_energy(match_idx) = energy_val
                filled_union(match_idx) = .true.
            end if
        end do
        close(unit_file)

        do i = 1, union_count
            if (.not. filled_union(i)) then
                write(*,'(A,I0)') 'Aviso: energias incompletas en calibracion del nivel ', level
                return
            end if
        end do

        if (do_low_calibration) then
            base_energy_low = get_base_energy()
            mat_low = 0.0_dp
            rhs_low = 0.0_dp
            do i = 1, n_targets
                idx_global = selected_low_idx(i)
                if (idx_global <= 0) cycle
                pos_union = 0
                do j = 1, union_count
                    if (union_idx(j) == idx_global) then
                        pos_union = j
                        exit
                    end if
                end do
                if (pos_union == 0) cycle
                energy_val = measured_energy(pos_union)
                do j = 1, n_orders
                    contrib_j = unique_low_contrib(j, idx_global)
                    rhs_low(j) = rhs_low(j) + contrib_j * (energy_val - base_energy_low)
                    do k = 1, n_orders
                        mat_low(j,k) = mat_low(j,k) + contrib_j * unique_low_contrib(k, idx_global)
                    end do
                end do
            end do

            call solve_normal_equations(mat_low, rhs_low, coeff_low, ok_system)
            if (.not. ok_system) then
                write(*,'(A,I0)') 'Aviso: sistema singular en calibracion (lado Si) del nivel ', level
                return
            end if
        end if

        coeff_high = 0.0_dp
        if (calibrate_high) then
            base_energy_high = get_high_base_energy()
            mat_high = 0.0_dp
            rhs_high = 0.0_dp
            do i = 1, n_targets
                idx_global = selected_high_idx(i)
                if (idx_global <= 0) cycle
                pos_union = 0
                do j = 1, union_count
                    if (union_idx(j) == idx_global) then
                        pos_union = j
                        exit
                    end if
                end do
                if (pos_union == 0) cycle
                energy_val = measured_energy(pos_union)
                do j = 1, n_orders
                    contrib_j = unique_high_contrib(j, idx_global)
                    rhs_high(j) = rhs_high(j) + contrib_j * (energy_val - base_energy_high)
                    do k = 1, n_orders
                        mat_high(j,k) = mat_high(j,k) + contrib_j * unique_high_contrib(k, idx_global)
                    end do
                end do
            end do

            call solve_normal_equations(mat_high, rhs_high, coeff_high, ok_system)
            if (.not. ok_system) then
                write(*,'(A,I0)') 'Aviso: sistema singular en calibracion (lado Ge) del nivel ', level
                calibrate_high = .false.
            end if
        end if

        if (do_low_calibration) then
            do i = 1, unique_count
                corrected = base_energy_low + sum(coeff_low(:) * unique_low_contrib(:, i))
                unique_low(i) = corrected
            end do
        end if

        if (calibrate_high) then
            do i = 1, unique_count
                corrected = base_energy_high + sum(coeff_high(:) * unique_high_contrib(:, i))
                unique_high(i) = corrected
            end do
        end if

        do i = 1, union_count
            idx_global = union_idx(i)
            energy_val = measured_energy(i)
            unique_low(idx_global) = energy_val
            unique_high(idx_global) = energy_val
        end do
    end subroutine calibrate_level_with_gulp

    subroutine canonicalize_subset(subset, level, eqmatrix, nop, canonical)
        implicit none
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

    ! Compares two integer arrays of length 'level' and returns .true. if array 'a' is lexicographically less than array 'b'.
    ! Inputs: a(:), b(:) - integer arrays; level - number of elements to compare.
    ! Output: logical - .true. if a < b lexicographically, .false. otherwise.
    pure logical function lexicographically_less(a, b, level)
        implicit none
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

    pure integer function find_subset_index(subset, level, stored_subsets, stored_count)
        implicit none
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
