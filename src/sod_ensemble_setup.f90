!*******************************************************************************
!  Setup mode: enumerates symmetry-inequivalent configurations, prepares POSCAR
!  files per level, runs the external job pipeline, and collects energies.
!*******************************************************************************
module sod_ensemble_setup_mod
    use sod_ensemble_consts
    use sod_ensemble_utils, only: binomial_int64, next_combination, sort_int_ascending
    use sod_ensemble_calibration, only: canonicalize_subset, find_subset_index
    use sod_ensemble_energy_calculations, only: init_energy_calc, cleanup_energy_calc, get_eqmatrix, write_vasp_file
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    character(len=512), save :: scripts_directory = ''
    logical, save :: scripts_ready = .false.
    contains

    subroutine run_sod_ensemble_setup(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: level_min, level_max
        logical :: has_list
        integer, allocatable :: level_list(:)
        integer, allocatable :: levels(:)
        integer :: argc, offset
        integer :: nop, total_sites
        integer, allocatable :: eqmatrix(:,:)
        integer :: i
        logical :: success

        level_min = 0
        level_max = -1
        has_list = .false.
        allocate(level_list(0))

        call parse_setup_arguments(level_min, level_max, level_list, has_list, arg_offset)
        call verify_support_scripts()

        call init_energy_calc(skip_energy_files=.true.)
        call get_eqmatrix(eqmatrix, nop, total_sites)
        if (.not. allocated(eqmatrix) .or. total_sites <= 0) then
            write(error_unit,'(A)') 'Error: no se pudo obtener EQMATRIX o no hay sitios sustituibles.'
            call cleanup_energy_calc()
            stop 1
        end if

        call build_level_sequence(level_min, level_max, level_list, has_list, total_sites, levels)
        if (.not. allocated(levels) .or. size(levels) == 0) then
            write(error_unit,'(A)') 'Error: la especificacion -N no produjo niveles validos.'
            deallocate(eqmatrix)
            call cleanup_energy_calc()
            stop 1
        end if

        write(*,'(A)') '--- Setup de configuraciones exactas ---'
        write(*,'(A,I0)') 'Sitios sustituibles (npos): ', total_sites
        write(*,'(A)', advance='no') 'Niveles a procesar: '
        do i = 1, size(levels)
            write(*,'(I0)', advance='no') levels(i)
            if (i < size(levels)) write(*,'(A)', advance='no') ', '
        end do
        write(*,*)
        call flush(output_unit)

        success = .true.
        do i = 1, size(levels)
            if (.not. process_level_setup(levels(i), total_sites, eqmatrix, nop)) then
                success = .false.
                exit
            end if
        end do

        deallocate(levels)
        deallocate(eqmatrix)
        call cleanup_energy_calc()

        if (.not. success) then
            stop 1
        end if
    end subroutine run_sod_ensemble_setup

    subroutine parse_setup_arguments(level_min, level_max, level_list, has_list, arg_offset)
        integer, intent(inout) :: level_min, level_max
        integer, allocatable, intent(inout) :: level_list(:)
        logical, intent(inout) :: has_list
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset
        character(len=256) :: arg, spec
        character(len=256) :: lowered

        if (allocated(level_list)) deallocate(level_list)
        allocate(level_list(0))
        has_list = .false.
        level_min = 0
        level_max = -1

        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)
        if (argc <= offset) return

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            select case (trim(lowered))
            case ('-n')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: falta especificacion despues de -N.'
                    call print_setup_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, spec)
                call parse_level_spec(spec, level_min, level_max, level_list, has_list)
                iarg = iarg + 2
            case ('--help','-h','--ayuda','-ayuda','ayuda')
                call print_setup_usage()
                stop 0
            case default
                write(error_unit,'(A)') 'Error: argumento no reconocido en modo setup. Use --help para mas informacion.'
                call print_setup_usage()
                stop 1
            end select
        end do
    end subroutine parse_setup_arguments

    subroutine print_setup_usage()
        write(*,'(A)') 'Uso: sod_ensemble setup -N <especificacion>'
        write(*,'(A)') ''
        write(*,'(A)') '  -N, -n    Define los niveles (átomos Ge) a preparar. Acepta:'
        write(*,'(A)') '            * un valor único (ej. 4)'
        write(*,'(A)') '            * un rango con dos puntos (ej. 2:6)'
        write(*,'(A)') '            * una lista separada por comas (ej. 0,3,5,7)'
        write(*,'(A)') ''
        write(*,'(A)') 'El modo setup genera carpetas n0X con OUTSOD, POSCARs, ejecuta run_jobs.sh y'
        write(*,'(A)') 'extract.sh para construir el archivo ENERGIES por nivel.'
    end subroutine print_setup_usage

    subroutine parse_level_spec(spec, level_min, level_max, level_list, has_list)
        character(len=*), intent(in) :: spec
        integer, intent(inout) :: level_min, level_max
        integer, allocatable, intent(inout) :: level_list(:)
        logical, intent(inout) :: has_list
        integer :: colon_pos, comma_pos
        character(len=len(spec)) :: token
        integer :: value, ios
        integer, allocatable :: temp(:)
        integer :: count, start, finish

        colon_pos = index(spec, ':')
        comma_pos = index(spec, ',')

        if (colon_pos > 0 .and. comma_pos > 0) then
            write(error_unit,'(A)') 'Error: no mezcle rangos y listas en la misma especificacion -N.'
            stop 1
        end if

        if (comma_pos > 0) then
            has_list = .true.
            count = 0
            allocate(temp(0))
            start = 1
            do
                finish = index(spec(start:), ',')
                if (finish == 0) then
                    token = adjustl(spec(start:))
                else
                    token = adjustl(spec(start:start+finish-2))
                end if
                read(token, *, iostat=ios) value
                if (ios /= 0) then
                    write(error_unit,'(A)') 'Error: valor invalido en la lista -N.'
                    stop 1
                end if
                call append_value(temp, count, value)
                if (finish == 0) exit
                start = start + finish
            end do
            if (allocated(level_list)) deallocate(level_list)
            call move_alloc(temp, level_list)
        else if (colon_pos > 0) then
            has_list = .false.
            read(spec(1:colon_pos-1), *, iostat=ios) level_min
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: limite inferior invalido en -N.'
                stop 1
            end if
            read(spec(colon_pos+1:), *, iostat=ios) level_max
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: limite superior invalido en -N.'
                stop 1
            end if
        else
            has_list = .false.
            read(spec, *, iostat=ios) value
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: especificacion -N invalida.'
                stop 1
            end if
            level_min = value
            level_max = value
        end if
    end subroutine parse_level_spec

    subroutine append_value(vec, count, value)
        integer, allocatable, intent(inout) :: vec(:)
        integer, intent(inout) :: count
        integer, intent(in) :: value
        integer, allocatable :: tmp(:)

        count = count + 1
        if (.not. allocated(vec)) then
            allocate(vec(1))
            vec(1) = value
            return
        end if

        allocate(tmp(count))
        if (count > 1) tmp(1:count-1) = vec(1:count-1)
        tmp(count) = value
        call move_alloc(tmp, vec)
    end subroutine append_value

    subroutine build_level_sequence(level_min, level_max, level_list, has_list, total_sites, levels)
        integer, intent(in) :: level_min, level_max, total_sites
        integer, allocatable, intent(in) :: level_list(:)
        logical, intent(in) :: has_list
        integer, allocatable, intent(out) :: levels(:)
        integer :: min_level, max_level
        integer :: i, n_valid, unique_count
        integer, allocatable :: work(:)

        if (has_list) then
            if (.not. allocated(level_list) .or. size(level_list) == 0) then
                allocate(levels(0))
                return
            end if
            allocate(work(size(level_list)))
            n_valid = 0
            do i = 1, size(level_list)
                if (level_list(i) < 0 .or. level_list(i) > total_sites) cycle
                n_valid = n_valid + 1
                work(n_valid) = level_list(i)
            end do
            if (n_valid == 0) then
                allocate(levels(0))
                deallocate(work)
                return
            end if
            call sort_int_ascending(work, n_valid)
            unique_count = 0
            do i = 1, n_valid
                if (unique_count == 0 .or. work(i) /= work(unique_count)) then
                    unique_count = unique_count + 1
                    work(unique_count) = work(i)
                end if
            end do
            allocate(levels(unique_count))
            levels = work(1:unique_count)
            deallocate(work)
        else
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
            do i = 1, size(levels)
                levels(i) = min_level + i - 1
            end do
        end if
    end subroutine build_level_sequence

    logical function process_level_setup(level, total_sites, eqmatrix, nop) result(success)
        integer, intent(in) :: level, total_sites, nop
        integer, intent(in) :: eqmatrix(:,:)
        integer :: unique_capacity, unique_count
        integer, allocatable :: unique_subsets(:,:)
        integer, allocatable :: unique_deg(:)
        integer, allocatable :: subset(:), canonical(:)
        integer :: idx
        integer(ip) :: total_comb
        character(len=16) :: level_dir
        logical :: ok
        character(len=512) :: command
        integer, allocatable :: config(:)

        success = .false.

        total_comb = binomial_int64(total_sites, level)
        write(*,'(A,I0)') 'Nivel: ', level
        write(*,'(A,I0)') 'Combinaciones totales: ', total_comb
        call flush(output_unit)

        level_dir = format_level_directory(level)
        command = 'mkdir -p ' // trim(level_dir)
        if (.not. run_shell_command(command, 'crear directorio ' // trim(level_dir))) return

        command = 'bash -c "rm -f ' // trim(level_dir)//'/c*.vasp ' // trim(level_dir)//'/*.vasp.gin ' // trim(level_dir)//'/*.vasp.gout ' // trim(level_dir)//'/*.vasp.grs ' // trim(level_dir)//'/OUTSOD ' // trim(level_dir)//'/ENERGIES"'
        if (.not. run_shell_command(command, 'limpiar carpeta ' // trim(level_dir))) return

        if (level == 0) then
            unique_count = 1
            allocate(unique_deg(1))
            unique_deg = 1
        else
            unique_capacity = max(16, level)
            allocate(unique_subsets(level, unique_capacity))
            allocate(unique_deg(unique_capacity))
            unique_deg = 0
            unique_count = 0
            allocate(subset(level))
            allocate(canonical(level))
            do idx = 1, level
                subset(idx) = idx
            end do
            do
                call canonicalize_subset(subset, level, eqmatrix, nop, canonical)
                idx = find_subset_index(canonical, level, unique_subsets, unique_count)
                if (idx > 0) then
                    unique_deg(idx) = unique_deg(idx) + 1
                else
                    unique_count = unique_count + 1
                    call ensure_unique_capacity(level, unique_subsets, unique_deg, unique_capacity, unique_count)
                    unique_subsets(1:level, unique_count) = canonical(1:level)
                    unique_deg(unique_count) = 1
                end if
                if (.not. next_combination(subset, total_sites)) exit
            end do
            deallocate(subset)
            deallocate(canonical)
        end if

        if (level > 0 .and. unique_count <= 0) then
            write(error_unit,'(A)') 'Error: no se encontraron configuraciones unicas.'
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        write(*,'(A,I0)') 'Configuraciones unicas: ', max(1, unique_count)
        call flush(output_unit)

        if (level > 0) then
            ok = write_outsod_file(level_dir, level, total_sites, unique_count, unique_deg, unique_subsets)
        else
            ok = write_outsod_file(level_dir, level, total_sites, 1, unique_deg)
        end if
        if (.not. ok) then
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        allocate(config(total_sites))
        config = 1
        if (level > 0) then
            do idx = 1, unique_count
                config = 1
                config(unique_subsets(1:level, idx)) = 2
                if (.not. write_poscar(level_dir, idx, config, total_sites)) then
                    deallocate(config)
                    if (allocated(unique_subsets)) deallocate(unique_subsets)
                    deallocate(unique_deg)
                    return
                end if
            end do
        else
            if (.not. write_poscar(level_dir, 1, config, total_sites)) then
                deallocate(config)
                deallocate(unique_deg)
                return
            end if
        end if
        deallocate(config)

        if (.not. copy_support_scripts(level_dir)) then
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        command = 'bash -c "cd ' // trim(level_dir) // ' && ./run_jobs.sh"'
        if (.not. run_shell_command(command, 'run_jobs.sh en '//trim(level_dir))) then
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        command = 'bash -c "cd ' // trim(level_dir) // ' && ./extract.sh"'
        if (.not. run_shell_command(command, 'extract.sh en '//trim(level_dir))) then
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        if (.not. energies_generated(level_dir)) then
            write(error_unit,'(A)') 'Error: no se encontro el archivo ENERGIES despues de extract.sh.'
            if (allocated(unique_subsets)) deallocate(unique_subsets)
            deallocate(unique_deg)
            return
        end if

        if (allocated(unique_subsets)) deallocate(unique_subsets)
        deallocate(unique_deg)

        success = .true.
    end function process_level_setup

    logical function write_outsod_file(dir, level, total_sites, unique_count, unique_deg, unique_subsets) result(ok)
        character(len=*), intent(in) :: dir
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_deg(:)
        integer, intent(in), optional :: unique_subsets(:,:)
        integer :: unit_outsod
        integer :: idx, site
        character(len=256) :: filename

        ok = .false.
        filename = trim(dir)//'/OUTSOD'
        open(newunit=unit_outsod, file=trim(filename), status='replace', action='write')
        write(unit_outsod,'(I12,"  substitutions in",I12," sites")') level, total_sites
        write(unit_outsod,'(I12,"  configurations")') unique_count
        if (unique_count > 0) then
            do idx = 1, unique_count
                write(unit_outsod,'(I6,1X,I12)', advance='no') idx, unique_deg(idx)
                if (level > 0 .and. present(unique_subsets)) then
                    do site = 1, level
                        write(unit_outsod,'(1X,I6)', advance='no') unique_subsets(site, idx)
                    end do
                end if
                write(unit_outsod,*)
            end do
        end if
        close(unit_outsod)
        ok = .true.
    end function write_outsod_file

    logical function write_poscar(dir, index, config, total_sites) result(ok)
        character(len=*), intent(in) :: dir
        integer, intent(in) :: index, total_sites
        integer, intent(in) :: config(:)
        character(len=256) :: filename
        character(len=32) :: index_str
        integer :: pad_len

        ok = .false.
        index_str = ''
        write(index_str,'(I0)') index
        pad_len = max(0, 5 - len_trim(index_str))
        filename = trim(dir)//'/c'//repeat('0', pad_len)//trim(index_str)//'.vasp'
        call write_vasp_file(config, total_sites, trim(filename))
        ok = .true.
    end function write_poscar

    logical function copy_support_scripts(dir) result(ok)
        character(len=*), intent(in) :: dir
        character(len=512) :: command

        if (.not. scripts_ready) call verify_support_scripts()

        command = 'cp ' // trim(scripts_directory)//'/run_jobs.sh ' // trim(scripts_directory)//'/extract.sh ' // trim(scripts_directory)//'/vasp2gin.sh ' // trim(dir)
        if (.not. run_shell_command(command, 'copiar scripts a '//trim(dir))) then
            ok = .false.
            return
        end if

        command = 'chmod +x ' // trim(dir)//'/run_jobs.sh ' // trim(dir)//'/extract.sh ' // trim(dir)//'/vasp2gin.sh'
        if (.not. run_shell_command(command, 'otorgar permisos en '//trim(dir))) then
            ok = .false.
            return
        end if

        ok = .true.
    end function copy_support_scripts

    logical function run_shell_command(command, context) result(ok)
        character(len=*), intent(in) :: command
        character(len=*), intent(in) :: context
        integer :: cmdstat, exitstat

        call execute_command_line(trim(command), wait=.true., cmdstat=cmdstat, exitstat=exitstat)
        if (cmdstat /= 0) then
            write(error_unit,'(A,": fallo al lanzar comando (cmdstat=",I0,")")') trim(context), cmdstat
            ok = .false.
            return
        end if
        if (exitstat /= 0) then
            write(error_unit,'(A,": el comando devolvio codigo ",I0)') trim(context), exitstat
            ok = .false.
            return
        end if
        ok = .true.
    end function run_shell_command

    subroutine ensure_unique_capacity(level, subsets, deg, capacity, required)
        integer, intent(in) :: level, required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: subsets(:,:)
        integer, allocatable, intent(inout) :: deg(:)
        integer :: new_capacity
        integer, allocatable :: tmp_subsets(:,:)
        integer, allocatable :: tmp_deg(:)

        if (required <= capacity) return
        new_capacity = max(required, int(capacity * 1.5_dp) + 4)
        if (new_capacity <= capacity) new_capacity = capacity + 4
        allocate(tmp_subsets(level, new_capacity))
        allocate(tmp_deg(new_capacity))
        if (capacity > 0) then
            tmp_subsets(:,1:capacity) = subsets(:,1:capacity)
            tmp_deg(1:capacity) = deg(1:capacity)
        end if
        if (new_capacity > capacity) then
            tmp_subsets(:,capacity+1:new_capacity) = 0
            tmp_deg(capacity+1:new_capacity) = 0
        end if
        call move_alloc(tmp_subsets, subsets)
        call move_alloc(tmp_deg, deg)
        capacity = new_capacity
    end subroutine ensure_unique_capacity

    logical function energies_generated(dir) result(ok)
        character(len=*), intent(in) :: dir
        logical :: exists
        inquire(file=trim(dir)//'/ENERGIES', exist=exists)
        ok = exists
    end function energies_generated

    subroutine verify_support_scripts()
        logical :: ok

        if (scripts_ready) return

        ok = locate_scripts_directory(scripts_directory)
        if (.not. ok) then
            write(error_unit,'(A)') 'Error: no se encontro el directorio scripts con run_jobs.sh, extract.sh y vasp2gin.sh.'
            stop 1
        end if

        scripts_ready = .true.
    end subroutine verify_support_scripts

    logical function locate_scripts_directory(out_dir) result(found)
        character(len=*), intent(out) :: out_dir
        integer :: depth, step
        character(len=512) :: candidate
        logical :: exists_run, exists_extract, exists_vasp2gin

        out_dir = ''
        found = .false.

        do depth = 0, 6
            if (depth == 0) then
                candidate = 'scripts'
            else
                candidate = ''
                do step = 1, depth
                    candidate = trim(candidate)//'../'
                end do
                candidate = trim(candidate)//'scripts'
            end if

            inquire(file=trim(candidate)//'/run_jobs.sh', exist=exists_run)
            if (.not. exists_run) cycle
            inquire(file=trim(candidate)//'/extract.sh', exist=exists_extract)
            if (.not. exists_extract) cycle
            inquire(file=trim(candidate)//'/vasp2gin.sh', exist=exists_vasp2gin)
            if (.not. exists_vasp2gin) cycle

            out_dir = trim(candidate)
            found = .true.
            exit
        end do
    end function locate_scripts_directory

    pure function format_level_directory(level) result(dir)
        integer, intent(in) :: level
        character(len=16) :: dir
        write(dir,'("n0",I0)') level
    end function format_level_directory

    subroutine to_lower_inplace(str)
        character(len=*), intent(inout) :: str
        integer :: i, code
        do i = 1, len(str)
            code = iachar(str(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                str(i:i) = achar(code + 32)
            end if
        end do
    end subroutine to_lower_inplace

end module sod_ensemble_setup_mod
