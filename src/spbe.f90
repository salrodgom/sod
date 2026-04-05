!*******************************************************************************
! Copyright (c) 2026, Salvador R.G. Balestra
!
! This file is part of the SOD package.
!
! SOD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!*******************************************************************************

!*******************************************************************************
! Unified SPBE frontend built on top of the modern effective-Hamiltonian core.
!*******************************************************************************
! Module `spbe` implements the pair-based extrapolation workflow on top of the modern core.
! El módulo `spbe` implementa el flujo de extrapolación por pares sobre el núcleo moderno.
module spbe
    use consts, only: dp, error_unit
    use cli, only: compose_mode_command, is_help_token, to_lower_inplace
    use utils, only: ensure_directory_exists, format_level_directory, join_path
    use inputs, only: insod_file_data, read_insod_file
    use configurations, only: read_outsod_file
    use energy_calculations, only: calculate_structure_energy, cleanup_energy_calc, get_base_energy, &
        get_high_base_energy, get_max_high_order, get_max_low_order, init_energy_calc
    implicit none
    private

    integer, parameter :: spbe_case_none = -1
    integer, parameter :: spbe_case_default = 0
    integer, parameter :: spbe_case_manual = 1
    integer, parameter :: spbe_case_reference = 2

    public :: run_sod_ensemble_spbe

contains

    subroutine run_sod_ensemble_spbe(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: argc, iarg, offset, eq_pos
        integer :: level, total_sites, unique_count
        integer :: selected_side
        integer, allocatable :: unique_subsets(:,:), unique_deg(:)
        integer, allocatable :: config(:)
        real(dp), allocatable :: energies(:), a1(:), a2(:)
        real(dp) :: mu1, mu2, base_energy
        real(dp) :: low_terms(4), high_terms(4), low_energy, high_energy
        real(dp) :: emin, emax
        integer :: m, mmin, mmax
        logical :: success, level_override, mu1_set, mu2_set, use_inspbe
        logical :: used_manual_mu, wrote_template
        integer :: rescale_case
        character(len=256) :: arg, lowered
        character(len=512) :: level_dir, output_dir, outsod_path, report_path, energies_path
        character(len=512) :: inspbe_path, inspbe_tmp_path
        character(len=8) :: side_name, side_folder
        character(len=128) :: rescale_message
        type(insod_file_data) :: insod

        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)

        level_override = .false.
        level = -1
        selected_side = 0
        mu1 = 1.0_dp
        mu2 = 1.0_dp
        mu1_set = .false.
        mu2_set = .false.
        use_inspbe = .true.

        argc = command_argument_count()
        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--level=') == 1 .or. index(trim(lowered), '-n=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid level specification for spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                read(arg(eq_pos+1:),*,err=901) level
                level_override = .true.
                iarg = iarg + 1
                cycle
            end if
            if (index(trim(lowered), '--side=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid side specification for spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                call parse_spbe_side(arg(eq_pos+1:), selected_side)
                iarg = iarg + 1
                cycle
            end if
            if (index(trim(lowered), '--mu1=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid --mu1 value in spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                read(arg(eq_pos+1:),*,err=902) mu1
                mu1_set = .true.
                iarg = iarg + 1
                cycle
            end if
            if (index(trim(lowered), '--mu2=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid --mu2 value in spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                read(arg(eq_pos+1:),*,err=903) mu2
                mu2_set = .true.
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-n','--level')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after -N/--level in spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, arg)
                read(arg,*,err=901) level
                level_override = .true.
                iarg = iarg + 2
            case ('--side')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --side in spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, arg)
                call parse_spbe_side(arg, selected_side)
                iarg = iarg + 2
            case ('--mu1')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --mu1 in spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, arg)
                read(arg,*,err=902) mu1
                mu1_set = .true.
                iarg = iarg + 2
            case ('--mu2')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after --mu2 in spbe mode.'
                    call print_spbe_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, arg)
                read(arg,*,err=903) mu2
                mu2_set = .true.
                iarg = iarg + 2
            case ('--no-inspbe')
                use_inspbe = .false.
                iarg = iarg + 1
            case ('--help','-h','help','--ayuda','-ayuda','ayuda')
                call print_spbe_usage()
                stop
            case default
                write(error_unit,'(A,1X,A)') 'Error: unrecognized argument in spbe mode:', trim(arg)
                call print_spbe_usage()
                stop 1
            end select
        end do

        if (mu1_set .neqv. mu2_set) then
            write(error_unit,'(A)') 'Error: you must provide both --mu1 and --mu2, or neither.'
            call print_spbe_usage()
            stop 1
        end if
        used_manual_mu = mu1_set .and. mu2_set

        call read_insod_file(insod)
        if (.not. level_override) level = insod%nsubs

        if (selected_side /= 0 .and. selected_side /= 1) then
            write(error_unit,'(A)') 'Error: spbe side must be X/low (0) or Y/high (1).'
            stop 1
        end if

        level_dir = trim(format_level_directory('n', level))
        outsod_path = join_path(trim(level_dir), 'OUTSOD')
        call read_outsod_file(trim(outsod_path), level, total_sites, unique_count, unique_subsets, unique_deg, success)
        if (.not. success) then
            write(error_unit,'(A,1X,A)') 'Error: could not read OUTSOD for spbe mode from', trim(outsod_path)
            stop 1
        end if
        if (level <= 0 .or. level >= total_sites) then
            write(error_unit,'(A)') 'Error: spbe mode requires 1 <= N < total number of sites.'
            stop 1
        end if

        if (selected_side == 0) then
            side_name = 'X'
            side_folder = 'spbe0'
        else
            side_name = 'Y'
            side_folder = 'spbe1'
        end if

        output_dir = join_path(trim(level_dir), trim(side_folder))
        call ensure_directory_exists(trim(output_dir))
        report_path = join_path(trim(output_dir), 'OUTSPBE')
        energies_path = join_path(trim(output_dir), 'ENERGIES')
        inspbe_path = join_path(trim(output_dir), 'INSPBE')
        inspbe_tmp_path = join_path(trim(output_dir), 'INSPBE.tmp')

        call init_energy_calc()
        if (selected_side == 0) then
            if (get_max_low_order() < 2) then
                write(error_unit,'(A)') 'Error: low-side order-2 data are not available for SPBE.'
                call cleanup_energy_calc()
                stop 1
            end if
            base_energy = get_base_energy()
        else
            if (get_max_high_order() < 2) then
                write(error_unit,'(A)') 'Error: high-side order-2 data are not available for SPBE.'
                call cleanup_energy_calc()
                stop 1
            end if
            base_energy = get_high_base_energy()
        end if

        allocate(config(total_sites))
        allocate(energies(unique_count))
        allocate(a1(unique_count))
        allocate(a2(unique_count))

        do m = 1, unique_count
            config = 1
            if (level > 0) config(unique_subsets(1:level, m)) = 2
            call calculate_structure_energy(config, total_sites, energies(m), low_energy, high_energy, low_terms, high_terms)
            if (selected_side == 0) then
                a1(m) = low_terms(1)
                a2(m) = low_terms(2)
            else
                a1(m) = high_terms(1)
                a2(m) = high_terms(2)
            end if
        end do

        rescale_case = spbe_case_none
        rescale_message = 'No rescaling applied (mu1=mu2=1.0).'
        if (used_manual_mu) then
            rescale_case = spbe_case_manual
            write(rescale_message,'(A,F12.6,A,F12.6)') 'Rescaling parameters provided on the CLI: mu1=', mu1, ', mu2=', mu2
        else if (use_inspbe) then
            call resolve_inspbe_parameters(trim(inspbe_path), base_energy, a1, a2, unique_count, mu1, mu2, rescale_case, rescale_message)
        else
            mu1 = 1.0_dp
            mu2 = 1.0_dp
            rescale_case = spbe_case_default
            rescale_message = 'INSPBE ignored by request; using mu1=mu2=1.0.'
        end if

        do m = 1, unique_count
            energies(m) = base_energy + mu1 * a1(m) + mu2 * a2(m)
        end do

        mmin = minloc(energies, dim=1)
        emin = minval(energies)
        mmax = maxloc(energies, dim=1)
        emax = maxval(energies)

        call write_spbe_outputs(trim(report_path), trim(energies_path), side_name, level, total_sites, unique_count, base_energy, &
            mu1, mu2, a1, a2, energies, mmin, emin, mmax, emax, trim(rescale_message))

        wrote_template = .false.
        if (.not. used_manual_mu) then
            if (.not. file_exists(trim(inspbe_path))) then
                call write_inspbe_template(trim(inspbe_tmp_path), mmin, emin, mmax, emax)
                wrote_template = .true.
            end if
        end if

        write(*,'(A)') '--- SPBE summary ---'
        write(*,'(A,I0)') 'Level                : ', level
        write(*,'(A,A)')  'Side                 : ', trim(side_name)
        write(*,'(A,F16.6)') 'Base energy (eV)     : ', base_energy
        write(*,'(A,F12.6)') 'mu1                  : ', mu1
        write(*,'(A,F12.6)') 'mu2                  : ', mu2
        write(*,'(A,F16.6)') 'Minimum energy (eV)  : ', emin
        write(*,'(A,I0)') 'Minimum config index : ', mmin
        write(*,'(A,F16.6)') 'Maximum energy (eV)  : ', emax
        write(*,'(A,I0)') 'Maximum config index : ', mmax
        write(*,'(A,A)') 'SPBE report          : ', trim(report_path)
        write(*,'(A,A)') 'SPBE energies        : ', trim(energies_path)
        if (wrote_template) write(*,'(A,A)') 'INSPBE template      : ', trim(inspbe_tmp_path)

        deallocate(config, energies, a1, a2, unique_subsets, unique_deg)
        call cleanup_energy_calc()
        return

901     continue
        write(error_unit,'(A)') 'Error: invalid integer level specification in spbe mode.'
        call print_spbe_usage()
        stop 1
902     continue
        write(error_unit,'(A)') 'Error: invalid numeric value for --mu1 in spbe mode.'
        call print_spbe_usage()
        stop 1
903     continue
        write(error_unit,'(A)') 'Error: invalid numeric value for --mu2 in spbe mode.'
        call print_spbe_usage()
        stop 1
    end subroutine run_sod_ensemble_spbe

    subroutine parse_spbe_side(raw_side, selected_side)
        character(len=*), intent(in) :: raw_side
        integer, intent(out) :: selected_side
        character(len=256) :: lowered

        lowered = adjustl(raw_side)
        call to_lower_inplace(lowered)

        select case (trim(lowered))
        case ('x','0','low','x-side','low-side')
            selected_side = 0
        case ('y','1','high','y-side','high-side')
            selected_side = 1
        case default
            write(error_unit,'(A,1X,A)') 'Error: invalid SPBE side:', trim(raw_side)
            stop 1
        end select
    end subroutine parse_spbe_side

    subroutine print_spbe_usage()
        character(len=256) :: command_name

        command_name = compose_mode_command('spbe')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [-N <level>] [--side X|Y] [--mu1 <value> --mu2 <value>] [--no-inspbe]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: spbe, pair, pairs'
        write(*,'(A)') ''
        write(*,'(A)') 'SPBE mode reproduces the classic pair-based extrapolation workflow using'
        write(*,'(A)') 'the modern effective-Hamiltonian module only. Outputs are written to'
        write(*,'(A)') 'nNN/spbe0 (X-side) or nNN/spbe1 (Y-side).'
        write(*,'(A)') ''
        write(*,'(A)') '  -N, --level <level>  Target combinatorial level. If omitted, INSOD nsubs is used.'
        write(*,'(A)') '  --side X|Y           Choose the expansion side. Default: X.'
        write(*,'(A)') '  --mu1, --mu2         Override the pair-extrapolation scaling parameters.'
        write(*,'(A)') '  --no-inspbe          Ignore nNN/spbe0|spbe1/INSPBE and use the defaults or CLI mu values.'
    end subroutine print_spbe_usage

    subroutine resolve_inspbe_parameters(inspbe_path, base_energy, a1, a2, unique_count, mu1, mu2, rescale_case, rescale_message)
        character(len=*), intent(in) :: inspbe_path
        real(dp), intent(in) :: base_energy
        real(dp), intent(in) :: a1(:), a2(:)
        integer, intent(in) :: unique_count
        real(dp), intent(out) :: mu1, mu2
        integer, intent(out) :: rescale_case
        character(len=*), intent(out) :: rescale_message
        integer :: unit_in, ios, irescale
        integer :: m1ref, m2ref
        real(dp) :: e1ref, e2ref, denom
        character(len=512) :: line

        mu1 = 1.0_dp
        mu2 = 1.0_dp
        rescale_case = spbe_case_default
        rescale_message = 'No INSPBE file; using mu1=mu2=1.0.'

        if (.not. file_exists(trim(inspbe_path))) return

        open(newunit=unit_in, file=trim(inspbe_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            rescale_message = 'Failed to read INSPBE; using mu1=mu2=1.0.'
            return
        end if

        read(unit_in,'(A)',iostat=ios) line
        if (ios /= 0) then
            close(unit_in)
            rescale_message = 'Malformed INSPBE header; using mu1=mu2=1.0.'
            return
        end if
        read(unit_in,*,iostat=ios) irescale
        if (ios /= 0) then
            close(unit_in)
            rescale_message = 'Malformed INSPBE case selector; using mu1=mu2=1.0.'
            return
        end if
        read(unit_in,'(A)',iostat=ios) line
        if (ios /= 0) then
            close(unit_in)
            rescale_message = 'Malformed INSPBE body; using mu1=mu2=1.0.'
            return
        end if

        select case (irescale)
        case (0)
            read(unit_in,'(A)',iostat=ios) line
            mu1 = 1.0_dp
            mu2 = 1.0_dp
            rescale_case = spbe_case_default
            rescale_message = 'INSPBE requests no rescaling (mu1=mu2=1.0).'
        case (1)
            read(unit_in,*,iostat=ios) mu1, mu2
            if (ios /= 0) then
                mu1 = 1.0_dp
                mu2 = 1.0_dp
                rescale_case = spbe_case_default
                rescale_message = 'Invalid INSPBE mu1/mu2 entry; using defaults.'
            else
                rescale_case = spbe_case_manual
                write(rescale_message,'(A,F12.6,A,F12.6)') 'Rescaling parameters read from INSPBE: mu1=', mu1, ', mu2=', mu2
            end if
        case (2)
            read(unit_in,*,iostat=ios) m1ref, e1ref
            if (ios == 0) read(unit_in,*,iostat=ios) m2ref, e2ref
            if (ios /= 0 .or. m1ref < 1 .or. m1ref > unique_count .or. m2ref < 1 .or. m2ref > unique_count) then
                mu1 = 1.0_dp
                mu2 = 1.0_dp
                rescale_case = spbe_case_default
                rescale_message = 'Invalid INSPBE reference indices; using defaults.'
            else
                denom = a1(m1ref) * a2(m2ref) - a2(m1ref) * a1(m2ref)
                if (abs(denom) <= 1.0e-12_dp) then
                    mu1 = 1.0_dp
                    mu2 = 1.0_dp
                    rescale_case = spbe_case_default
                    rescale_message = 'Degenerate INSPBE reference system; using defaults.'
                else
                    mu1 = (a2(m2ref) * (e1ref - base_energy) - a2(m1ref) * (e2ref - base_energy)) / denom
                    mu2 = (a1(m1ref) * (e2ref - base_energy) - a1(m2ref) * (e1ref - base_energy)) / denom
                    rescale_case = spbe_case_reference
                    write(rescale_message,'(A,I0,A,I0,A,F12.6,A,F12.6)') 'Rescaling from INSPBE reference configurations ', &
                        m1ref, ' and ', m2ref, ': mu1=', mu1, ', mu2=', mu2
                end if
            end if
        case default
            mu1 = 1.0_dp
            mu2 = 1.0_dp
            rescale_case = spbe_case_default
            rescale_message = 'Unknown INSPBE mode; using mu1=mu2=1.0.'
        end select

        close(unit_in)
    end subroutine resolve_inspbe_parameters

    subroutine write_spbe_outputs(report_path, energies_path, side_name, level, total_sites, unique_count, base_energy, mu1, mu2, &
        a1, a2, energies, mmin, emin, mmax, emax, rescale_message)
        character(len=*), intent(in) :: report_path, energies_path, side_name, rescale_message
        integer, intent(in) :: level, total_sites, unique_count, mmin, mmax
        real(dp), intent(in) :: base_energy, mu1, mu2, emin, emax
        real(dp), intent(in) :: a1(:), a2(:), energies(:)
        integer :: unit_report, unit_energy, m

        open(newunit=unit_report, file=trim(report_path), status='replace', action='write')
        open(newunit=unit_energy, file=trim(energies_path), status='replace', action='write')

        write(unit_report,'(A)') '-----------------------------------------------------------------'
        write(unit_report,'(A)') 'Simple pair-based extrapolation from the modern SOD core'
        write(unit_report,'(A)') '-----------------------------------------------------------------'
        write(unit_report,*)
        write(unit_report,'(A,I0)') 'Level: ', level
        write(unit_report,'(A,I0)') 'Total sites: ', total_sites
        write(unit_report,'(A,A)') 'Side: ', trim(side_name)
        write(unit_report,'(A,F20.10)') 'Base energy (eV): ', base_energy
        write(unit_report,'(A,F20.10)') 'mu1: ', mu1
        write(unit_report,'(A,F20.10)') 'mu2: ', mu2
        write(unit_report,'(A,A)') 'Rescaling: ', trim(rescale_message)
        write(unit_report,*)
        write(unit_report,'(A,I0,A)') 'Zeroth-, first-, and second-order contribution by configuration for ', level, ' substitutions:'
        write(unit_report,*)
        write(unit_report,'(A)') '   m     E0(eV)             dE1(eV)         dE2(eV)      E[total] (eV)'
        write(unit_report,'(A)') ' --------------------------------------------------------------------------'

        do m = 1, unique_count
            write(unit_report,'(I5,4F16.6)') m, base_energy, mu1 * a1(m), mu2 * a2(m), energies(m)
            write(unit_energy,'(F20.10)') energies(m)
        end do

        write(unit_report,*)
        write(unit_report,'(A,I0,A,F20.10,A)') 'Minimum-energy configuration: ', mmin, ' with energy: ', emin, ' eV.'
        write(unit_report,'(A,I0,A,F20.10,A)') 'Maximum-energy configuration: ', mmax, ' with energy: ', emax, ' eV.'

        close(unit_report)
        close(unit_energy)
    end subroutine write_spbe_outputs

    subroutine write_inspbe_template(filename, mmin, emin, mmax, emax)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: mmin, mmax
        real(dp), intent(in) :: emin, emax
        integer :: unit_tmp

        open(newunit=unit_tmp, file=trim(filename), status='replace', action='write')
        write(unit_tmp,'(A)') '# irescale case: 0 = no rescaling; 1 = enter mu1 and mu2 manually; 2 = enter two reference energies'
        write(unit_tmp,'(A)') '2'
        write(unit_tmp,'(A)') '# If irescale=1, enter one line with mu1, mu2; if irescale=2, enter two lines (m1, E1), (m2, E2)'
        write(unit_tmp,'(I0,1X,F20.10)') mmin, emin
        write(unit_tmp,'(I0,1X,F20.10)') mmax, emax
        close(unit_tmp)
    end subroutine write_inspbe_template

    logical function file_exists(path)
        character(len=*), intent(in) :: path
        inquire(file=trim(path), exist=file_exists)
    end function file_exists

end module spbe
