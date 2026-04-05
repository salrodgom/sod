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

!*******************************************************************************
! Utility mode that initializes the symmetry mapping from INSOD/SGO and writes
! the resulting EQMATRIX file to disk for inspection or reuse.
! Modo utilitario que inicializa el mapeo de simetría desde INSOD/SGO y escribe
! el archivo EQMATRIX resultante a disco para inspección o reutilización.
!*******************************************************************************
module eqmatrix
    use cli, only: compose_mode_command, is_help_token, to_lower_inplace
    use symmetry, only: symmetry_finalize, symmetry_initialize, symmetry_write_matrix
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    implicit none
    private

    public :: run_sod_ensemble_eqmatrix

contains

    subroutine run_sod_ensemble_eqmatrix(arg_offset)
        integer, intent(in), optional :: arg_offset
        character(len=512) :: output_name

        call parse_arguments_eqmatrix(output_name, arg_offset)

        write(*,'(A)') '--- EQMATRIX generation ---'
        write(*,'(A)') 'Initializing symmetry data from INSOD and SGO.'
        flush(output_unit)

        call symmetry_initialize(force_generate=.true.)
        call symmetry_write_matrix(trim(output_name))
        call symmetry_finalize()

        write(*,'(A)') 'EQMATRIX file ready: '//trim(output_name)
        flush(output_unit)
    end subroutine run_sod_ensemble_eqmatrix

    subroutine parse_arguments_eqmatrix(output_name, arg_offset)
        character(len=*), intent(out) :: output_name
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset, eq_pos
        character(len=512) :: arg, lowered

        output_name = 'EQMATRIX'

        argc = command_argument_count()
        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)
        if (argc <= offset) return

        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            if (is_help_token(arg)) then
                call print_usage_eqmatrix()
                stop
            end if

            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--output=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: missing file name after --output=.'
                    call print_usage_eqmatrix()
                    stop 1
                end if
                output_name = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-o', '--output')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing file name after --output.'
                    call print_usage_eqmatrix()
                    stop 1
                end if
                call get_command_argument(iarg + 1, output_name)
                output_name = adjustl(output_name)
                iarg = iarg + 2
            case default
                write(error_unit,'(A)') 'Error: unrecognized argument in eqmatrix mode. Use --help for more information.'
                call print_usage_eqmatrix()
                stop 1
            end select
        end do

        if (len_trim(output_name) == 0) then
            write(error_unit,'(A)') 'Error: the EQMATRIX output file name cannot be empty.'
            call print_usage_eqmatrix()
            stop 1
        end if
    end subroutine parse_arguments_eqmatrix

    subroutine print_usage_eqmatrix()
        character(len=256) :: command_name

        command_name = compose_mode_command('eqmatrix')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [-o <file>]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: eqmatrix, matrix'
        write(*,'(A)') ''
        write(*,'(A)') '  -o, --output <file>   Output file name. Default: EQMATRIX'
        write(*,'(A)') ''
        write(*,'(A)') 'The current directory must contain INSOD and SGO. The mode initializes'
        write(*,'(A)') 'the symmetry mapping directly from those files and writes the resulting'
        write(*,'(A)') 'EQMATRIX without reading the level ENERGIES files.'
    end subroutine print_usage_eqmatrix

end module eqmatrix
