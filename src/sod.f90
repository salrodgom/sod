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

! Program `sod` dispatches the unified executable to the selected operating mode.
! Entry point for the unified SOD executable; selects MC, exact, setup, clean, comb, canonical, SPBE, entropy, gc, compare, or eqmatrix modes.
! Punto de entrada del ejecutable unificado de SOD; selecciona los modos MC, exact, setup, clean, comb, canonical, SPBE, entropy, gc, compare o eqmatrix.
! El programa `sod` redirige el ejecutable unificado al modo de operación seleccionado.
program sod
    use canonical, only: run_sod_canonical
    use cli, only: is_help_token, to_lower_inplace
    use clean, only: run_sod_clean
    use comb, only: run_sod_ensemble_comb
    use compare, only: run_sod_ensemble_compare
    use eqmatrix, only: run_sod_ensemble_eqmatrix
    use entropy, only: run_sod_ensemble_entropy
    use gc, only: run_sod_gc
    use mc, only: run_sod_ensemble_mc
    use exact, only: run_sod_ensemble_exact
    use spbe, only: run_sod_ensemble_spbe
    use setup, only: run_sod_ensemble_setup
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
        stop
    end if

    mode_selected = .false.
    mode_token = ''
    offset = 0

    call classify_first_argument(first_arg, argc, mode_selected, mode_token, offset)
    if (.not. mode_selected) then
        write(error_unit,'(A)') 'Error: you must specify "mc", "exact", "setup", "clean", "comb", "canonical", "spbe", "entropy", "gc", "compare", or "eqmatrix" as the mode.'
        call print_combined_usage()
        stop 1
    end if
    
    if (.not. normalize_mode(mode_token)) then
        write(error_unit,'(A)') 'Error: unknown mode "'//trim(mode_token)//'". Use "mc", "exact", "setup", "clean", "comb", "canonical", "spbe", "entropy", "gc", "compare", or "eqmatrix".'
        call print_combined_usage()
        stop 1
    end if

    select case (trim(mode_token))
    case ('mc')
        call run_sod_ensemble_mc(arg_offset=offset)
    case ('exact')
        call run_sod_ensemble_exact(arg_offset=offset)
    case ('setup')
        call run_sod_ensemble_setup(arg_offset=offset)
    case ('clean')
        call run_sod_clean(arg_offset=offset)
    case ('comb')
        call run_sod_ensemble_comb(arg_offset=offset)
    case ('canonical')
        call run_sod_canonical(arg_offset=offset)
    case ('spbe')
        call run_sod_ensemble_spbe(arg_offset=offset)
    case ('entropy')
        call run_sod_ensemble_entropy(arg_offset=offset)
    case ('gc')
        call run_sod_gc(arg_offset=offset)
    case ('compare')
        call run_sod_ensemble_compare(arg_offset=offset)
    case ('eqmatrix')
        call run_sod_ensemble_eqmatrix(arg_offset=offset)
    case default
        write(error_unit,'(A)') 'Internal error: mode not recognized after normalization.'
        stop 1
    end select

contains

    ! Classifies the first CLI argument, extracting the desired execution mode and argument offset.
! Clasifica el primer argumento de la CLI y extrae el modo de ejecución deseado y el desplazamiento de argumentos.
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
        case ('mc', 'montecarlo', 'monte-carlo', 'sampling', 'sample', 'spbemc')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('exact', 'enumeracion', 'enumeration', 'enumerate', 'exhaustiva', 'exhaustive', 'full')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('setup', 'preparacion', 'prepare', 'inputs', 'energy')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('clean', 'cleanup', 'prune', 'tidy')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('comb', 'combsod', 'configs', 'combinatorial', 'list', 'unique')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('canonical', 'thermo', 'thermal', 'stats', 'statsod', 'c')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('spbe', 'pair', 'pairs', 'pair-based', 'pair_based')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('entropy', 'config-entropy', 'config_entropy', 'entropia', 'sconf')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('gc', 'grand', 'grandcanonical', 'grand-canonical', 'gcstat', 'populations')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('compare', 'comparison', 'compare-config', 'compare_config', 'diff', 'reference')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('eqmatrix', 'eq-matrix', 'matrix')
            mode_selected = .true.
            mode_token = trim(lowered_trim)
            offset = 1
        case ('--mode')
            if (argc < 2) then
                write(error_unit,'(A)') 'Error: expected a value after --mode.'
                call print_combined_usage()
                stop 1
            end if
            call get_command_argument(2, mode_token)
            mode_token = adjustl(mode_token)
            call to_lower_inplace(mode_token)
            if (len_trim(mode_token) == 0) then
                write(error_unit,'(A)') 'Error: expected a value after --mode.'
                call print_combined_usage()
                stop 1
            end if
            mode_selected = .true.
            offset = 2
        case default
            if (index(lowered_trim, '--mode=') == 1) then
                eq_pos = index(raw_arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(raw_arg)) then
                    write(error_unit,'(A)') 'Error: expected a value after --mode=.'
                    call print_combined_usage()
                    stop 1
                end if
                mode_token = raw_arg(eq_pos+1:)
                mode_token = adjustl(mode_token)
                call to_lower_inplace(mode_token)
                if (len_trim(mode_token) == 0) then
                    write(error_unit,'(A)') 'Error: expected a value after --mode=.'
                    call print_combined_usage()
                    stop 1
                end if
                mode_selected = .true.
                offset = 1
            else if (is_help_token(raw_arg)) then
                call print_combined_usage()
                stop
            else
                mode_selected = .false.
            end if
        end select
    end subroutine classify_first_argument

    ! Normalizes a mode token to one of the canonical identifiers.
! Normaliza un token de modo a uno de los identificadores canónicos.
    logical function normalize_mode(token)
        character(len=*), intent(inout) :: token
        character(len=256) :: lowered
        character(len=256) :: trimmed_token

        lowered = adjustl(token)
        call to_lower_inplace(lowered)
        trimmed_token = trim(lowered)

        select case (trim(trimmed_token))
        case ('mc', 'montecarlo', 'monte-carlo', 'stochastic', 'mcmc', 'sampling', 'sample', 'spbemc')
            token = 'mc'
            normalize_mode = .true.
        case ('exact', 'enumeracion', 'enumeration', 'enumerate', 'exhaustive', 'exhaustiva', 'full')
            token = 'exact'
            normalize_mode = .true.
        case ('setup', 'preparacion', 'prepare', 'inputs', 'energy')
            token = 'setup'
            normalize_mode = .true.
        case ('clean', 'cleanup', 'prune', 'tidy')
            token = 'clean'
            normalize_mode = .true.
        case ('comb', 'combsod', 'configs', 'combinatorial', 'list', 'unique')
            token = 'comb'
            normalize_mode = .true.
        case ('canonical', 'thermo', 'thermal', 'stats', 'statsod', 'c')
            token = 'canonical'
            normalize_mode = .true.
        case ('spbe', 'pair', 'pairs', 'pair-based', 'pair_based')
            token = 'spbe'
            normalize_mode = .true.
        case ('entropy', 'config-entropy', 'config_entropy', 'entropia', 'sconf')
            token = 'entropy'
            normalize_mode = .true.
        case ('gc', 'grand', 'grandcanonical', 'grand-canonical', 'gcstat', 'populations')
            token = 'gc'
            normalize_mode = .true.
        case ('compare', 'comparison', 'compare-config', 'compare_config', 'diff', 'reference')
            token = 'compare'
            normalize_mode = .true.
        case ('eqmatrix', 'eq-matrix', 'matrix')
            token = 'eqmatrix'
            normalize_mode = .true.
        case default
            normalize_mode = .false.
        end select
    end function normalize_mode

    ! Prints usage instructions combining the MC, exact, setup, clean, comb, canonical, SPBE, entropy, gc, compare, and eqmatrix modes.
! Imprime las instrucciones de uso combinando los modos MC, exact, setup, clean, comb, canonical, SPBE, entropy, gc, compare y eqmatrix.
    subroutine print_combined_usage()
        write(*,'(A)') 'Usage: sod <mc|exact|setup|clean|comb|canonical|thermo|spbe|entropy|gc|compare|eqmatrix> [mode arguments]'
        write(*,'(A)') '       sod --mode <mc|exact|setup|clean|comb|canonical|thermo|spbe|entropy|gc|compare|eqmatrix> [mode arguments]'
        write(*,'(A)') '       sod --mode=<mc|exact|setup|clean|comb|canonical|thermo|spbe|entropy|gc|compare|eqmatrix> [mode arguments]'
        write(*,'(A)') '       sod --help'
        write(*,'(A)') ''
        write(*,'(A)') '  mc|sampling|sample|spbemc -> runs the Monte Carlo sampling workflow.'
        write(*,'(A)') '  exact|exhaustive|enumerate|full -> performs exhaustive enumeration.'
        write(*,'(A)') '  setup|prepare|inputs|energy -> prepares n0X folders and triggers the external energy pipeline.'
        write(*,'(A)') '  clean|prune|tidy -> prunes transient helper files from finished n0X setup directories.'
        write(*,'(A)') '  comb|list|unique|configs -> reproduces the classic combinatorial workflow without external calculations.'
        write(*,'(A)') '  canonical|thermo|stats|c -> runs the classic canonical thermodynamic post-processing for one level.'
        write(*,'(A)') '  spbe|pair|pairs -> runs the classic pair-based extrapolation on top of the modern Hamiltonian core.'
        write(*,'(A)') '  entropy|sconf -> computes configurational entropy from OUTSOD files.'
        write(*,'(A)') '  gc|grand|populations -> performs grand-canonical post-processing over several substitution levels.'
        write(*,'(A)') '  compare|diff|reference -> builds DeltaG_config comparison tables and a gnuplot fitting script.'
        write(*,'(A)') '  eqmatrix|matrix -> writes the symmetry mapping EQMATRIX from INSOD and SGO.'
        write(*,'(A)') ''
        write(*,'(A)') 'Examples:'
        write(*,'(A)') '  sod mc -T 800 -M 6 -C 2000 -s 1234 -a metropolis --omp'
        write(*,'(A)') '  sod exact -N 5:10 -t 1e-5'
        write(*,'(A)') '  sod setup -N 3,6,9 --template-gin default'
        write(*,'(A)') '  sod clean -N 3:9'
        write(*,'(A)') '  sod comb'
        write(*,'(A)') '  sod canonical -N 4 --source n'
        write(*,'(A)') '  sod spbe -N 4 --side X'
        write(*,'(A)') '  sod entropy -N -1'
        write(*,'(A)') '  sod gc'
        write(*,'(A)') '  sod compare --system phase_A --reference phase_B -T 1500'
        write(*,'(A)') '  sod eqmatrix -o EQMATRIX'
    end subroutine print_combined_usage

end program sod
