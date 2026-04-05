!*******************************************************************************
! Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
!
! This file is part of the SOD package.
!
! SOD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SOD is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!*******************************************************************************
!
! Modern shared-module update of the classic combinatorial SOD workflow.
!
!*******************************************************************************
! Module `comb` implements the pure combinatorial workflow for inequivalent configurations.
! El módulo `comb` implementa el flujo puramente combinatorio para configuraciones inequivalentes.
module comb
    use consts, only: dp, ip, kB_eVk, error_unit
    use cli, only: compose_mode_command, is_help_token, parse_level_spec, build_level_sequence, to_lower_inplace
    use utils, only: format_level_directory, join_path, ensure_directory_exists, print_sod_banner
    use inputs, only: sgo_file_data, insod_file_data, read_sgo_file, read_insod_file
    use configurations, only: collect_subset_stabilizer_operators, enumerate_unique_subsets, write_outsod_unit
    use structure_io, only: write_castep_configuration, write_cif_configuration, write_gulp_configuration, &
        write_lammps_configuration, write_qe_configuration, write_vasp_configuration, &
        motif_atom_type, motif_data_type, read_motif_file
    use symmetry, only: symmetry_restrict_supercell_operators, symmetry_wrap_fractional_vector
    implicit none
    private

    real(dp), parameter :: coord_tol = 1.0e-3_dp

    public :: run_sod_ensemble_comb
    public :: run_combination_workflow

contains

    subroutine run_sod_ensemble_comb(arg_offset)
        integer, intent(in), optional :: arg_offset
        integer :: offset, argc, iarg, eq_pos
        integer :: level_min, level_max
        integer, allocatable :: level_list(:), levels(:)
        logical :: has_list, level_specified
        integer :: level_index
        character(len=256) :: arg, spec, lowered
        character(len=512) :: motif_option
        type(insod_file_data) :: insod

        offset = 0
        if (present(arg_offset)) offset = max(0, arg_offset)

        level_min = 0
        level_max = -1
        has_list = .false.
        level_specified = .false.
        motif_option = ''
        allocate(level_list(0))
        argc = command_argument_count()
        iarg = 1 + offset
        do while (iarg <= argc)
            call get_command_argument(iarg, arg)
            lowered = adjustl(arg)
            call to_lower_inplace(lowered)

            if (index(trim(lowered), '--level=') == 1 .or. index(trim(lowered), '-n=') == 1) then
                eq_pos = index(arg, '=')
                if (eq_pos <= 0 .or. eq_pos == len_trim(arg)) then
                    write(error_unit,'(A)') 'Error: invalid level specification for comb mode.'
                    call print_comb_usage()
                    stop 1
                end if
                spec = adjustl(arg(eq_pos+1:))
                call parse_level_spec(spec, level_min, level_max, level_list, has_list)
                level_specified = .true.
                iarg = iarg + 1
                cycle
            end if
            if (index(trim(lowered), '--motif=') == 1) then
                eq_pos = index(arg, '=')
                motif_option = adjustl(arg(eq_pos+1:))
                iarg = iarg + 1
                cycle
            end if

            select case (trim(lowered))
            case ('-n','--level')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing value after -N/--level in comb mode.'
                    call print_comb_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, spec)
                call parse_level_spec(spec, level_min, level_max, level_list, has_list)
                level_specified = .true.
                iarg = iarg + 2
            case ('--motif')
                if (iarg + 1 > argc) then
                    write(error_unit,'(A)') 'Error: missing path after --motif.'
                    call print_comb_usage()
                    stop 1
                end if
                call get_command_argument(iarg + 1, motif_option)
                motif_option = adjustl(motif_option)
                iarg = iarg + 2
            case ('--help','-h','--ayuda','-ayuda','ayuda')
                call print_comb_usage()
                stop
            case default
                write(error_unit,'(A,1X,A)') 'Error: unrecognized argument in comb mode:', trim(arg)
                call print_comb_usage()
                stop 1
            end select
        end do

        call read_insod_file(insod)
        if (.not. level_specified) then
            level_min = insod%nsubs_min
            level_max = insod%nsubs_max
        end if

        ! Use a generous upper bound; run_combination_workflow validates each level
        call build_level_sequence(level_min, level_max, level_list, has_list, &
            insod%nat0 * insod%na * insod%nb * insod%nc, levels)
        if (allocated(level_list)) deallocate(level_list)
        if (.not. allocated(levels) .or. size(levels) == 0) then
            write(error_unit,'(A)') 'Error: the level selection did not produce any valid levels.'
            stop 1
        end if

        do level_index = 1, size(levels)
            if (len_trim(motif_option) > 0) then
                call run_combination_workflow(levels(level_index), motif_file=trim(motif_option))
            else
                call run_combination_workflow(levels(level_index))
            end if
        end do
        if (allocated(levels)) deallocate(levels)
    end subroutine run_sod_ensemble_comb

    subroutine print_comb_usage()
        character(len=256) :: command_name

        command_name = compose_mode_command('comb')
        write(*,'(A)') 'Usage: '//trim(command_name)//' [-N <spec>] [--motif <file>]'
        write(*,'(A)') ''
        write(*,'(A)') 'Aliases: comb, list, unique, configs'
        write(*,'(A)') ''
        write(*,'(A)') 'Comb mode reproduces the classic combinatorial SOD workflow in a modern'
        write(*,'(A)') 'per-level directory nNN. It writes OUTSOD, SUPERCELL, EQMATRIX, OPERATORS,'
        write(*,'(A)') 'cSGO, filer, and, when FILER requests it, the corresponding calculation'
        write(*,'(A)') 'input files such as cXXXXX.vasp, cXXXXX.cif, cXXXXX.gin, cXXXXX.cell,'
        write(*,'(A)') 'cXXXXX.pw.in, or cXXXXX.out.ase/job_sender stubs for ASE-based workflows.'
        write(*,'(A)') ''
        write(*,'(A)') '  -N, --level <spec>   Override INSOD nsubs for the combinatorial run.'
        write(*,'(A)') '                       Supports a single value (e.g. 3), a range (e.g. 0:5),'
        write(*,'(A)') '                       or a comma-separated list (e.g. 0,1,2,3,4,5).'
        write(*,'(A)') '                       If omitted, the full nsubs_min:nsubs_max range from INSOD is used.'
        write(*,'(A)') '  --motif <file>       Include molecular motif atoms from the specified .include'
        write(*,'(A)') '                       file into each generated structure. The motif file contains'
        write(*,'(A)') '                       atom lines in GULP format (Symbol role x y z) or plain'
        write(*,'(A)') '                       format (Symbol x y z).'
    end subroutine print_comb_usage

    subroutine run_combination_workflow(level, motif_file)
        integer, intent(in) :: level
        character(len=*), intent(in), optional :: motif_file
        type(sgo_file_data) :: sgo
        type(insod_file_data) :: insod
        integer :: nop1, nop, npos, nsp, sptarget, nat0, nat1r, nat1_unit, nat
        integer :: na, nb, nc, op, at, at0, at1, at1i, op1, t, pos, attmp, i, j, k
        integer :: unit_eq, unit_super, unit_ops, unit_csg, unit_filer
        integer :: unique_count, stabilizer_count, sp
        integer(ip) :: total_comb, candidate_count
        logical :: found, used_recursive
        real(dp) :: prod, x, maxentropy, ientropy
        real(dp) :: a, b, c, alpha, beta, gamma
        real(dp) :: coordstemp(3)
        character(len=32) :: level_dir
        character(len=512) :: neighbor_paths(2), chosen_neighbor
        integer, allocatable :: natsp1(:), natsp(:), spat1r(:), spat_unit(:), spat(:)
        integer, allocatable :: pos2coord(:), coord2pos(:), fulleqmatrix(:,:), eqmatrixtarget(:,:)
        integer, allocatable :: unique_subsets(:,:), unique_deg(:), stabilizer_ops(:)
        real(dp), allocatable :: coords1r(:,:), coords1(:,:), coords(:,:)
        real(dp), allocatable :: mgroup(:,:,:), vgroup(:,:), vt(:,:)

        call print_sod_banner()
        write(*,*) 'Reading input files...'
        write(*,*) ' '

        call read_sgo_file(sgo)
        call read_insod_file(insod)

        nsp = insod%nsp
        sptarget = insod%sptarget
        nat0 = insod%nat0
        na = insod%na
        nb = insod%nb
        nc = insod%nc
        nop1 = sgo%nop1

        if (level < 0) then
            write(error_unit,'(A)') 'Error: comb mode received a negative substitution level.'
            stop 1
        end if
        if (sptarget < 1 .or. sptarget > nsp) then
            write(error_unit,'(A)') 'Error: invalid target species index in INSOD.'
            stop 1
        end if

        write(*,*) 'Generating the supercell...'
        write(*,*) ' '

        nat1r = nat0 * nop1
        allocate(coords1r(nat1r,3))
        allocate(spat1r(nat1r))
        at1 = 0
        do at0 = 1, nat0
            do op1 = 1, nop1
                at1 = at1 + 1
                coords1r(at1,1:3) = matmul(sgo%mgroup1(op1,1:3,1:3), insod%coords0(at0,1:3)) + sgo%vgroup1(op1,1:3)
                call symmetry_wrap_fractional_vector(coords1r(at1,1:3))
                spat1r(at1) = insod%spat0(at0)
            end do
        end do

        allocate(coords1(nat1r,3))
        allocate(spat_unit(nat1r))
        coords1(1,1:3) = coords1r(1,1:3)
        spat_unit(1) = spat1r(1)
        at1 = 1
        do at0 = 2, nat1r
            found = .false.
            do at1i = 1, at1
                prod = dot_product(coords1r(at0,1:3) - coords1(at1i,1:3), coords1r(at0,1:3) - coords1(at1i,1:3))
                if (prod <= coord_tol) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                at1 = at1 + 1
                coords1(at1,1:3) = coords1r(at0,1:3)
                spat_unit(at1) = spat1r(at0)
            end if
        end do
        nat1_unit = at1

        allocate(natsp1(nsp))
        natsp1 = 0
        do at1 = 1, nat1_unit
            natsp1(spat_unit(at1)) = natsp1(spat_unit(at1)) + 1
        end do

        allocate(natsp(nsp))
        natsp(1:nsp) = natsp1(1:nsp) * na * nb * nc
        npos = natsp(sptarget)
        if (npos <= 0) then
            write(error_unit,'(A)') 'Error: no substitutable sites found in comb mode.'
            stop 1
        end if
        if (level > npos) then
            write(error_unit,'(A)') 'Error: comb mode level exceeds the number of substitutable sites.'
            stop 1
        end if

        call symmetry_restrict_supercell_operators(sgo%mgroup1, sgo%vgroup1, nop1, na, nb, nc)

        allocate(vt(na*nb*nc,3))
        t = 0
        do i = 0, na-1
            do j = 0, nb-1
                do k = 0, nc-1
                    t = t + 1
                    vt(t,1) = real(i,dp) / real(na,dp)
                    vt(t,2) = real(j,dp) / real(nb,dp)
                    vt(t,3) = real(k,dp) / real(nc,dp)
                end do
            end do
        end do

        nat = nat1_unit * na * nb * nc
        allocate(coords(nat,3))
        allocate(spat(nat))
        at = 0
        do at0 = 1, nat1_unit
            do t = 1, na*nb*nc
                at = at + 1
                coords(at,1) = vt(t,1) + coords1(at0,1) / real(na,dp)
                coords(at,2) = vt(t,2) + coords1(at0,2) / real(nb,dp)
                coords(at,3) = vt(t,3) + coords1(at0,3) / real(nc,dp)
                call symmetry_wrap_fractional_vector(coords(at,1:3))
                spat(at) = spat_unit(at0)
            end do
        end do

        a = real(na,dp) * insod%cell_params(1)
        b = real(nb,dp) * insod%cell_params(2)
        c = real(nc,dp) * insod%cell_params(3)
        alpha = insod%cell_params(4)
        beta = insod%cell_params(5)
        gamma = insod%cell_params(6)

        level_dir = format_level_directory('n', level)
        call ensure_directory_exists(level_dir)
        call clean_comb_level_directory(level_dir)

        open(newunit=unit_super, file=trim(join_path(level_dir,'SUPERCELL')), status='replace', action='write')
        write(unit_super,'(6(f10.4,2x))') a, b, c, alpha, beta, gamma
        write(unit_super,*) natsp(1:nsp)
        do at = 1, nat
            write(unit_super,'(A3,2X,3(f11.7,2X))') insod%symbol(spat(at)), coords(at,1), coords(at,2), coords(at,3)
        end do
        close(unit_super)

        nop = nop1 * na * nb * nc
        allocate(mgroup(nop,3,3))
        allocate(vgroup(nop,3))
        op = 0
        do op1 = 1, nop1
            do t = 1, na*nb*nc
                op = op + 1
                mgroup(op,1:3,1:3) = sgo%mgroup1(op1,1:3,1:3)
                vgroup(op,1) = sgo%vgroup1(op1,1) / real(na,dp) + vt(t,1)
                vgroup(op,2) = sgo%vgroup1(op1,2) / real(nb,dp) + vt(t,2)
                vgroup(op,3) = sgo%vgroup1(op1,3) / real(nc,dp) + vt(t,3)
            end do
        end do

        open(newunit=unit_ops, file=trim(join_path(level_dir,'OPERATORS')), status='replace', action='write')
        write(unit_ops,*) nop
        do op = 1, nop
            write(unit_ops,*) 'Operator number ', op
            do i = 1, 3
                write(unit_ops,*) mgroup(op,i,1:3), vgroup(op,i)
            end do
        end do
        close(unit_ops)

        write(*,*) '       Composition of the original supercell (parent structure):'
        do sp = 1, nsp
            write(*,*) '                                                         ', trim(insod%symbol(sp)), natsp(sp)
        end do
        write(*,*) ' '
        write(*,*) '       Number of symmetry operators in the supercell:       ', nop
        write(*,*) ' '
        write(*,*) '       Composition of the substituted supercell:'
        do sp = 1, nsp
            if (sp == sptarget) then
                write(*,*) '                                                         ', trim(insod%newsymbol(1)), level
                write(*,*) '                                                         ', trim(insod%newsymbol(2)), natsp(sp) - level
            else
                write(*,*) '                                                         ', trim(insod%symbol(sp)), natsp(sp)
            end if
        end do
        write(*,*) ' '
        write(*,*) 'Generating the complete configurational space...'
        write(*,*) ' '

        allocate(pos2coord(npos))
        pos = 0
        do at = 1, nat
            if (spat(at) == sptarget) then
                pos = pos + 1
                pos2coord(pos) = at
            end if
        end do

        allocate(fulleqmatrix(nop,nat))
        allocate(coord2pos(nat))
        coord2pos = 0
        do pos = 1, npos
            coord2pos(pos2coord(pos)) = pos
        end do

        do op = 1, nop
            do at = 1, nat
                coordstemp = matmul(mgroup(op,1:3,1:3), coords(at,1:3)) + vgroup(op,1:3)
                call symmetry_wrap_fractional_vector(coordstemp)
                found = .false.
                do attmp = 1, nat
                    if (abs(coordstemp(1) - coords(attmp,1)) < coord_tol .and. &
                        abs(coordstemp(2) - coords(attmp,2)) < coord_tol .and. &
                        abs(coordstemp(3) - coords(attmp,3)) < coord_tol) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    write(error_unit,'(A,I0,A,I0,A)') 'Error: operator ', op, ' applied to atom ', at, &
                        ' does not map onto a valid supercell atom.'
                    stop 1
                end if
                fulleqmatrix(op,at) = attmp
            end do
        end do

        allocate(eqmatrixtarget(nop,npos))
        do op = 1, nop
            do pos = 1, npos
                eqmatrixtarget(op,pos) = coord2pos(fulleqmatrix(op,pos2coord(pos)))
            end do
        end do

        open(newunit=unit_eq, file=trim(join_path(level_dir,'EQMATRIX')), status='replace', action='write')
        write(unit_eq,*) nop, npos
        do op = 1, nop
            write(unit_eq,*) eqmatrixtarget(op,1:npos)
        end do
        close(unit_eq)

        total_comb = combinations_ip(level, npos)
        maxentropy = kB_eVk * log(real(total_comb,dp))
        write(*,*) ' '
        write(*,*) '       Total number of configurations in the supercell:     ', total_comb
        write(*,*) ' '
        write(*,*) '       Maximum entropy for this composition and supercell:', maxentropy, ' eV/K'
        write(*,*) ' '

        x = real(level,dp) / real(npos,dp)
        write(*,*) '       Fraction of substituted sites:                 x = ', x
        write(*,*) ' '
        if (level == 0 .or. level == npos) then
            ientropy = 0.0_dp
        else
            ientropy = -real(npos,dp) * kB_eVk * (x * log(x) + (1.0_dp - x) * log(1.0_dp - x))
        end if
        write(*,*) '       Ideal entropy (per cell) for this composition:     ', ientropy, ' eV/K'
        write(*,*) ' '

        neighbor_paths = ''
        if (level > 0) neighbor_paths(1) = join_path(format_level_directory('n', level - 1), 'OUTSOD')
        if (level < npos) neighbor_paths(2) = join_path(format_level_directory('n', level + 1), 'OUTSOD')
        call enumerate_unique_subsets(npos, level, eqmatrixtarget, nop, unique_subsets, unique_deg, unique_count, &
            neighbor_paths=neighbor_paths, used_recursive=used_recursive, candidate_count=candidate_count, &
            chosen_neighbor=chosen_neighbor)
        allocate(stabilizer_ops(nop))

        if (used_recursive) then
            write(*,*) 'Using recursive generation from ', trim(chosen_neighbor)
            write(*,*) ' '
            write(*,*) '       Number of candidate configurations:            ', candidate_count
            write(*,*) ' '
            write(*,*) 'Finding the inequivalent configurations (recursive)...'
        else
            write(*,*) 'Generating the complete configurational space...'
            write(*,*) ' '
            write(*,*) 'Finding the inequivalent configurations...'
        end if
        write(*,*) ' '
        write(*,*) '       Found    Completion '
        write(*,*) '       =====    ========== '

        open(newunit=unit_csg, file=trim(join_path(level_dir,'cSGO')), status='replace', action='write')
        do i = 1, unique_count
            call collect_subset_stabilizer_operators(unique_subsets(1:level,i), level, eqmatrixtarget, nop, &
                stabilizer_ops, stabilizer_count)
            write(unit_csg,*) 'List of operators for configuration: ', i
            do j = 1, stabilizer_count
                op = stabilizer_ops(j)
                write(unit_csg,*) j
                do k = 1, 3
                    write(unit_csg,*) mgroup(op,k,1:3), vgroup(op,k)
                end do
            end do
            write(unit_csg,*) 0
        end do
        close(unit_csg)

        write(*,'(4X,I6,7X,F5.1,A2)') unique_count, 100.0_dp, '% '
        write(*,*) ' '
        write(*,*) '       Number of inequivalent configurations:               ', unique_count
        write(*,*) ' '

        call write_legacy_outsod(join_path(level_dir,'OUTSOD'), level, npos, unique_count, unique_deg, unique_subsets)

        open(newunit=unit_filer, file=trim(join_path(level_dir,'filer')), status='replace', action='write')
        write(unit_filer,*) insod%filer
        close(unit_filer)

        if (present(motif_file)) then
            call write_comb_calculation_files(level_dir, insod%filer, (/a,b,c,alpha,beta,gamma/), insod%symbol, spat, coords, &
                sptarget, insod%newsymbol, pos2coord, level, unique_count, unique_subsets, motif_file=motif_file)
        else
            call write_comb_calculation_files(level_dir, insod%filer, (/a,b,c,alpha,beta,gamma/), insod%symbol, spat, coords, &
                sptarget, insod%newsymbol, pos2coord, level, unique_count, unique_subsets)
        end if
        write(*,'(A,1X,A)') 'Level outputs written to:', trim(level_dir)
        write(*,*) 'Done!!!'
        write(*,*) ' '
        write(*,*) ' '
    end subroutine run_combination_workflow

    subroutine write_legacy_outsod(filename, level, total_sites, unique_count, unique_deg, unique_subsets)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: level, total_sites, unique_count
        integer, intent(in) :: unique_deg(:)
        integer, intent(in) :: unique_subsets(:,:)
        integer :: unit_outsod

        open(newunit=unit_outsod, file=trim(filename), status='replace', action='write')
        call write_outsod_unit(unit_outsod, level, total_sites, unique_count, unique_deg, unique_subsets, &
            legacy_format=.true.)
        close(unit_outsod)
    end subroutine write_legacy_outsod

    subroutine clean_comb_level_directory(level_dir)
        character(len=*), intent(in) :: level_dir
        integer :: exit_code
        character(len=2048) :: command

        command = 'rm -f ' // trim(level_dir)//'/OUTSOD ' // trim(level_dir)//'/SUPERCELL ' // trim(level_dir)// &
            '/EQMATRIX ' // trim(level_dir)//'/OPERATORS ' // trim(level_dir)//'/cSGO ' // trim(level_dir)// &
            '/filer ' // trim(level_dir)//'/job_sender ' // trim(level_dir)//'/c*.vasp ' // trim(level_dir)// &
            '/c*.vout ' // trim(level_dir)//'/c*.cif ' // trim(level_dir)//'/c*.cell ' // trim(level_dir)// &
            '/c*.pw.in ' // trim(level_dir)//'/c*.pw.out ' // trim(level_dir)//'/c*.gin ' // trim(level_dir)// &
            '/c*.min ' // trim(level_dir)//'/c*.data ' // trim(level_dir)//'/c*.min.data ' // trim(level_dir)// &
            '/c*.in.lammps ' // trim(level_dir)//'/c*.out.lammps ' // trim(level_dir)//'/c*.log.lammps ' // &
            '/c*.out.ase ' // trim(level_dir)//'/c*.relaxed.cif ' // trim(level_dir)//'/c*.relaxed.vasp ' // &
            '/c*.ase.log ' // trim(level_dir)//'/c*.traj ' // trim(level_dir)//'/RELAXED_STRUCTURES'
        call execute_command_line(trim(command), exitstat=exit_code)
        if (exit_code /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: failed to clean comb output directory', trim(level_dir)
            stop 1
        end if
    end subroutine clean_comb_level_directory

    subroutine write_comb_calculation_files(level_dir, filer, cell_params, symbols, spat, coords, target_species, &
        replacement_symbols, target_positions, level, unique_count, unique_subsets, motif_file)
        character(len=*), intent(in) :: level_dir
        integer, intent(in) :: filer, target_species, level, unique_count
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: unique_subsets(:,:)
        character(len=*), intent(in), optional :: motif_file
        integer :: i, unit_job
        character(len=32) :: basename
        character(len=512) :: path
        logical :: template_exists
        type(motif_data_type) :: motif
        logical :: have_motif

        have_motif = .false.
        if (present(motif_file)) then
            if (len_trim(motif_file) > 0) then
                call read_motif_file(trim(motif_file), motif)
                have_motif = (motif%natoms > 0)
                if (have_motif) then
                    write(*,'(A,I0,A,1X,A)') 'Loaded ', motif%natoms, ' motif atoms from', trim(motif_file)
                    if (motif%ntail > 0) then
                        write(*,'(A,I0,A)') 'Loaded ', motif%ntail, ' motif tail lines (connectivity, etc.)'
                    end if
                end if
            end if
        end if

        select case (filer)
        case (-1)
            write(*,*) 'No calculation files requested (FILER=-1).'
        case (0)
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating CIF files in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                path = join_path(level_dir, trim(basename)//'.cif')
                if (have_motif) then
                    call write_cif_configuration(trim(path), i, cell_params, symbols, spat, coords, target_species, &
                        replacement_symbols, target_positions, unique_subsets(1:level,i), motif%atoms, motif%natoms)
                else
                    call write_cif_configuration(trim(path), i, cell_params, symbols, spat, coords, target_species, &
                        replacement_symbols, target_positions, unique_subsets(1:level,i))
                end if
            end do
        case (11)
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating input files for VASP in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                path = join_path(level_dir, trim(basename)//'.vasp')
                if (have_motif) then
                    call write_vasp_configuration(trim(path), 'vasp', cell_params, symbols, spat, coords, target_species, &
                        replacement_symbols, target_positions, unique_subsets(1:level,i), motif%atoms, motif%natoms)
                else
                    call write_vasp_configuration(trim(path), 'vasp', cell_params, symbols, spat, coords, target_species, &
                        replacement_symbols, target_positions, unique_subsets(1:level,i))
                end if
            end do
            open(newunit=unit_job, file=trim(join_path(level_dir,'job_sender')), status='replace', action='write')
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                write(unit_job,'(A)') 'vasp < ' // trim(basename)//'.vasp > ' // trim(basename)//'.vout'
            end do
            close(unit_job)
        case (1)
            inquire(file='template_input.gin', exist=template_exists)
            if (.not. template_exists) then
                write(error_unit,'(A)') 'Error: template_input.gin not found for FILER=1.'
                stop 1
            end if
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating input files for GULP in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                path = join_path(level_dir, trim(basename)//'.gin')
                if (have_motif) then
                    call write_gulp_configuration(trim(path), 'template_input.gin', trim(basename), cell_params, symbols, &
                        spat, coords, target_species, replacement_symbols, target_positions, unique_subsets(1:level,i), &
                        motif=motif)
                else
                    call write_gulp_configuration(trim(path), 'template_input.gin', trim(basename), cell_params, symbols, &
                        spat, coords, target_species, replacement_symbols, target_positions, unique_subsets(1:level,i))
                end if
            end do
            open(newunit=unit_job, file=trim(join_path(level_dir,'job_sender')), status='replace', action='write')
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                write(unit_job,'(A)') 'gulp < ' // trim(basename)//'.gin > ' // trim(basename)//'.gout'
            end do
            close(unit_job)
        case (14)
            inquire(file='template_ase.py', exist=template_exists)
            if (.not. template_exists) then
                write(error_unit,'(A)') 'Error: template_ase.py not found for FILER=14.'
                stop 1
            end if
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating input files for ASE in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                path = join_path(level_dir, trim(basename)//'.vasp')
                if (have_motif) then
                    call write_vasp_configuration(trim(path), 'vasp', cell_params, symbols, spat, coords, target_species, &
                        replacement_symbols, target_positions, unique_subsets(1:level,i), motif%atoms, motif%natoms)
                else
                    call write_vasp_configuration(trim(path), 'vasp', cell_params, symbols, spat, coords, target_species, &
                        replacement_symbols, target_positions, unique_subsets(1:level,i))
                end if
            end do
            open(newunit=unit_job, file=trim(join_path(level_dir,'job_sender')), status='replace', action='write')
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                write(unit_job,'(A)') 'python3 vasp2ase.py ' // trim(basename)//'.vasp template_ase.py > ' // &
                    trim(basename)//'.out.ase'
            end do
            close(unit_job)
        case (2)
            inquire(file='template_in.lammps', exist=template_exists)
            if (.not. template_exists) then
                write(error_unit,'(A)') 'Error: template_in.lammps not found for FILER=2.'
                stop 1
            end if
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating input files for LAMMPS in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                if (have_motif) then
                    call write_lammps_configuration(join_path(level_dir, trim(basename)//'.data'), &
                        join_path(level_dir, trim(basename)//'.in.lammps'), 'template_in.lammps', trim(basename), &
                        cell_params, symbols, spat, coords, target_species, replacement_symbols, target_positions, &
                        unique_subsets(1:level,i), motif%atoms, motif%natoms)
                else
                    call write_lammps_configuration(join_path(level_dir, trim(basename)//'.data'), &
                        join_path(level_dir, trim(basename)//'.in.lammps'), 'template_in.lammps', trim(basename), &
                        cell_params, symbols, spat, coords, target_species, replacement_symbols, target_positions, &
                        unique_subsets(1:level,i))
                end if
            end do
            open(newunit=unit_job, file=trim(join_path(level_dir,'job_sender')), status='replace', action='write')
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                write(unit_job,'(A)') 'lmp < ' // trim(basename)//'.in.lammps > ' // trim(basename)//'.out.lammps'
            end do
            close(unit_job)
        case (12)
            inquire(file='template_castep.cell', exist=template_exists)
            if (.not. template_exists) then
                write(error_unit,'(A)') 'Error: template_castep.cell not found for FILER=12.'
                stop 1
            end if
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating input files for CASTEP in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                path = join_path(level_dir, trim(basename)//'.cell')
                if (have_motif) then
                    call write_castep_configuration(trim(path), 'template_castep.cell', basename(2:), cell_params, symbols, &
                        spat, coords, target_species, replacement_symbols, target_positions, unique_subsets(1:level,i), &
                        motif%atoms, motif%natoms)
                else
                    call write_castep_configuration(trim(path), 'template_castep.cell', basename(2:), cell_params, symbols, &
                        spat, coords, target_species, replacement_symbols, target_positions, unique_subsets(1:level,i))
                end if
            end do
            open(newunit=unit_job, file=trim(join_path(level_dir,'job_sender')), status='replace', action='write')
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                write(unit_job,'(A)') 'castep ' // trim(basename)
            end do
            close(unit_job)
        case (13)
            inquire(file='template_pw.in', exist=template_exists)
            if (.not. template_exists) then
                write(error_unit,'(A)') 'Error: template_pw.in not found for FILER=13.'
                stop 1
            end if
            write(*,*) ' '
            write(*,'(A,1X,A,A)') 'Creating input files for Quantum ESPRESSO in ./', trim(level_dir), '...'
            write(*,*) ' '
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                path = join_path(level_dir, trim(basename)//'.pw.in')
                if (have_motif) then
                    call write_qe_configuration(trim(path), 'template_pw.in', basename(2:), cell_params, symbols, spat, &
                        coords, target_species, replacement_symbols, target_positions, unique_subsets(1:level,i), &
                        motif%atoms, motif%natoms)
                else
                    call write_qe_configuration(trim(path), 'template_pw.in', basename(2:), cell_params, symbols, spat, &
                        coords, target_species, replacement_symbols, target_positions, unique_subsets(1:level,i))
                end if
            end do
            open(newunit=unit_job, file=trim(join_path(level_dir,'job_sender')), status='replace', action='write')
            do i = 1, unique_count
                write(basename,'("c",I5.5)') i
                write(unit_job,'(A)') 'pw.x < ' // trim(basename)//'.pw.in > ' // trim(basename)//'.pw.out'
            end do
            close(unit_job)
        case default
            write(*,'(A,I0,A)') 'Warning: FILER=', filer, ' is not yet implemented in comb mode.'
            write(*,*) 'Implemented FILER values in comb mode: -1, 0, 1, 2, 11, 12, 13, 14.'
        end select
        if (allocated(motif%atoms)) deallocate(motif%atoms)
        if (allocated(motif%tail_lines)) deallocate(motif%tail_lines)
    end subroutine write_comb_calculation_files

    pure function combinations_ip(k, n) result(val)
        integer, intent(in) :: k, n
        integer(ip) :: val
        integer :: i, kk

        if (k < 0 .or. k > n) then
            val = 0_ip
            return
        end if
        kk = k
        if (kk > n - kk) kk = n - kk
        if (kk == 0) then
            val = 1_ip
            return
        end if

        val = 1_ip
        do i = 1, kk
            val = val * int(n - kk + i, ip) / int(i, ip)
        end do
    end function combinations_ip

end module comb
