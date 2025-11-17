!*******************************************************************************
!  Common calibration routines shared between SOD Boltzmann drivers.
!*******************************************************************************
module sod_calibration
    use sod_boltzmann_consts
    use sod_boltzmann_utils
    use energy_calc
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    private
    public :: calibrate_level_with_gulp, canonicalize_subset, find_subset_index, evaluate_subsets_with_gulp

contains

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

    subroutine evaluate_subsets_with_gulp(level, total_sites, subset_count, subsets, config, energies, success)
        implicit none
        integer, intent(in) :: level, total_sites, subset_count
        integer, intent(in) :: subsets(:,:)
        integer, intent(inout) :: config(:)
        real(dp), intent(out) :: energies(:)
        logical, intent(out) :: success

        character(len=256) :: script_dir, sample_dir
        character(len=64) :: base_name
        character(len=256) :: gout_label
        character(len=512) :: command, line
        character(len=128), allocatable :: gout_files(:)
        logical, allocatable :: filled(:)
        integer :: exit_code, ios, unit_energy
        integer :: i, pos_hash, match_idx
        real(dp) :: energy_val

        success = .false.

        if (subset_count <= 0) return
        if (size(energies) < subset_count) then
            write(*,'(A)') 'Aviso: buffer de energias insuficiente en evaluate_subsets_with_gulp.'
            call flush(output_unit)
            return
        end if

        if (.not. find_scripts_directory(script_dir)) then
            write(*,'(A)') 'Aviso: no se encontro directorio de scripts para GULP.'
            call flush(output_unit)
            return
        end if

        write(sample_dir,'("mc_samples_N",I4.4)') level
        call execute_command_line('rm -rf ' // trim(sample_dir), exitstat=exit_code)
        call execute_command_line('mkdir -p ' // trim(sample_dir), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A)') 'Aviso: no se pudo crear el directorio de evaluacion GULP.'
            call flush(output_unit)
            return
        end if

        allocate(gout_files(subset_count))
        allocate(filled(subset_count))
        gout_files = ''
        filled = .false.

        do i = 1, subset_count
            config = 1
            if (level > 0 .and. size(subsets,1) >= level) then
                config(subsets(1:level, i)) = 2
            end if
            write(base_name,'("sample_N",I4.4,"_c",I5.5)') level, i
            call write_vasp_file(config, total_sites, trim(sample_dir)//'/'//trim(base_name)//'.vasp')
            gout_files(i) = trim(base_name)//'.vasp.gout'
        end do
        config = 1

        call copy_calibration_scripts(script_dir, sample_dir)
    write(*,'(A,I0,A,A)') 'Nivel ', level, ': configuraciones GULP en ', trim(sample_dir)
    call flush(output_unit)

        command = 'cd ' // trim(sample_dir) // ' && bash run_jobs.sh'
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A)') 'Aviso: run_jobs.sh fallo al evaluar las muestras con GULP.'
            call flush(output_unit)
            call cleanup()
            return
        end if

        command = 'cd ' // trim(sample_dir) // ' && bash extract.sh'
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(*,'(A)') 'Aviso: extract.sh fallo al evaluar las muestras con GULP.'
            call flush(output_unit)
            call cleanup()
            return
        end if

        open(newunit=unit_energy, file=trim(sample_dir)//'/ENERGIES', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') 'Aviso: no se pudo abrir el archivo ENERGIES generado por GULP.'
            call flush(output_unit)
            call cleanup()
            return
        end if

        do
            read(unit_energy,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            pos_hash = index(line, '#')
            if (pos_hash <= 0) cycle
            read(line(1:pos_hash-1), *, iostat=ios) energy_val
            if (ios /= 0) cycle
            gout_label = adjustl(trim(line(pos_hash+1:)))
            gout_label = trim(gout_label)
            match_idx = 0
            do i = 1, subset_count
                if (gout_label == trim(gout_files(i))) then
                    match_idx = i
                    exit
                end if
            end do
            if (match_idx > 0) then
                energies(match_idx) = energy_val
                filled(match_idx) = .true.
            end if
        end do
        close(unit_energy)

        if (.not. all(filled)) then
            write(*,'(A)') 'Aviso: faltan energias en la salida de GULP.'
            call flush(output_unit)
            call cleanup()
            return
        end if

        success = .true.
        call cleanup()
        return

    contains
        subroutine cleanup()
            if (allocated(gout_files)) deallocate(gout_files)
            if (allocated(filled)) deallocate(filled)
        end subroutine cleanup
    end subroutine evaluate_subsets_with_gulp

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
        directory_has_scripts = exist_vasp .and. exist_run .and. exist_extract
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

end module sod_calibration
