program sod_ensemble
    use sod_ensemble_mc_mod, only: run_sod_ensemble_mc
    use sod_ensemble_exact_mod, only: run_sod_ensemble_exact
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none

    integer :: argc, offset
    character(len=256) :: first_arg
    character(len=256) :: mode_token
    logical :: mode_selected

    argc = command_argument_count()
    if (argc == 0) then
        call print_combined_usage()
        stop 1
    end if

    call get_command_argument(1, first_arg)
    if (is_help_token(first_arg)) then
        call print_combined_usage()
        stop 0
    end if

    mode_selected = .false.
    mode_token = ''
    offset = 0

    call classify_first_argument(first_arg, argc, mode_selected, mode_token, offset)
    if (.not. mode_selected) then
        write(error_unit,'(A)') 'Error: se debe indicar el modo "mc" o "exact" como primer argumento o mediante --mode.'
        call print_combined_usage()
        stop 1
    end if

    if (.not. normalize_mode(mode_token)) then
        write(error_unit,'(A)') 'Error: modo desconocido "'//trim(mode_token)//'". Use "mc" o "exact".'
        call print_combined_usage()
        stop 1
    end if

    select case (trim(mode_token))
    case ('mc')
        call run_sod_ensemble_mc(arg_offset=offset)
    case ('exact')
        call run_sod_ensemble_exact(arg_offset=offset)
    case default
        write(error_unit,'(A)') 'Error interno: modo no reconocido tras normalización.'
        stop 1
    end select

contains

    subroutine classify_first_argument(raw_arg, argc, mode_selected, mode_token, offset)
        character(len=*), intent(in) :: raw_arg
        integer, intent(in) :: argc
        logical, intent(out) :: mode_selected
        character(len=*), intent(inout) :: mode_token
        integer, intent(out) :: offset
        character(len=256) :: lowered
        character(len=256) :: lowered_trim
        integer :: eq_pos

        lowered = adjustl(raw_arg)
        call to_lower_inplace(lowered)
        lowered_trim = trim(lowered)
        mode_selected = .false.
        mode_token = ''
        offset = 0

        select case (trim(lowered_trim))
        case ('mc', 'montecarlo', 'monte-carlo')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('exact', 'enumeracion', 'enumeration', 'exhaustiva', 'exhaustive')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('--mode')
            if (argc < 2) then
                write(error_unit,'(A)') 'Error: se esperaba un valor después de --mode.'
                call print_combined_usage()
                stop 1
            end if
            call get_command_argument(2, mode_token)
            mode_token = adjustl(mode_token)
            call to_lower_inplace(mode_token)
            if (len_trim(mode_token) == 0) then
                write(error_unit,'(A)') 'Error: se esperaba un valor después de --mode.'
                call print_combined_usage()
                stop 1
            end if
            mode_selected = .true.
            offset = 2
        case default
            if (index(lowered_trim, '--mode=') == 1) then
                eq_pos = index(raw_arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(raw_arg)) then
                    write(error_unit,'(A)') 'Error: se esperaba un valor después de --mode=.'
                    call print_combined_usage()
                    stop 1
                end if
                mode_token = raw_arg(eq_pos+1:)
                mode_token = adjustl(mode_token)
                call to_lower_inplace(mode_token)
                if (len_trim(mode_token) == 0) then
                    write(error_unit,'(A)') 'Error: se esperaba un valor después de --mode=.'
                    call print_combined_usage()
                    stop 1
                end if
                mode_selected = .true.
                offset = 1
            else if (is_help_token(raw_arg)) then
                call print_combined_usage()
                stop 0
            else
                mode_selected = .false.
            end if
        end select
    end subroutine classify_first_argument

    logical function normalize_mode(token)
        character(len=*), intent(inout) :: token
        character(len=256) :: lowered
        character(len=256) :: trimmed_token

        lowered = adjustl(token)
        call to_lower_inplace(lowered)
        trimmed_token = trim(lowered)

        select case (trim(trimmed_token))
        case ('mc', 'montecarlo', 'monte-carlo', 'stochastic', 'mcmc')
            token = 'mc'
            normalize_mode = .true.
        case ('exact', 'enumeracion', 'enumeration', 'enumerate', 'exhaustive', 'exhaustiva')
            token = 'exact'
            normalize_mode = .true.
        case default
            normalize_mode = .false.
        end select
    end function normalize_mode

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

    logical function is_help_token(raw)
        character(len=*), intent(in) :: raw
        character(len=len(raw)) :: token
        integer :: i, lt, code

        token = adjustl(raw)
        lt = len_trim(token)
        if (lt == 0) then
            is_help_token = .false.
            return
        end if

        do i = 1, lt
            code = iachar(token(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) token(i:i) = achar(code + 32)
        end do
        token = token(1:lt)

        select case (trim(token))
        case ('--help', '-h', 'help', '--ayuda', '-ayuda', 'ayuda', '/h', '/?')
            is_help_token = .true.
        case default
            is_help_token = .false.
        end select
    end function is_help_token

    subroutine print_combined_usage()
        write(*,'(A)') 'Uso: sod_ensemble <mc|exact> [argumentos del modo]'
        write(*,'(A)') '     sod_ensemble --mode <mc|exact> [argumentos del modo]'
        write(*,'(A)') '     sod_ensemble --mode=<mc|exact> [argumentos del modo]'
        write(*,'(A)') '     sod_ensemble --help'
        write(*,'(A)') ''
        write(*,'(A)') '  mc    -> ejecuta el muestreo Monte Carlo (equivalente a sod_ensemble_mc).'
        write(*,'(A)') '  exact -> ejecuta la enumeración exhaustiva (equivalente a sod_ensemble_exact).'
        write(*,'(A)') ''
        write(*,'(A)') 'Ejemplos:'
        write(*,'(A)') '  sod_ensemble mc -T 800 -M 6 -C 2000 -s 1234 -a metropolis --omp'
        write(*,'(A)') '  sod_ensemble exact -N 5:10 -t 1e-5'
        write(*,'(A)') ''
        write(*,'(A)') 'Los ejecutables tradicionales sod_ensemble_mc y sod_ensemble_exact siguen disponibles como'
        write(*,'(A)') 'envoltorios que delegan en este binario unificado.'
    end subroutine print_combined_usage

end program sod_ensemble
