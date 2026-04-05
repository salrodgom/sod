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
! Modern shared-module update of the classic VASP/POSCAR writing logic.
!
!*******************************************************************************
! Module `structure_io` writes crystal structures in the supported file formats.
! El módulo `structure_io` escribe estructuras cristalinas en los formatos soportados.
module structure_io
    use consts, only: dp, error_unit
    use cell_mod, only: cell
    implicit none
    private

    type :: motif_atom_type
        character(len=16) :: symbol = ''
        character(len=8)  :: role   = 'core'
        real(dp)          :: x = 0.0_dp, y = 0.0_dp, z = 0.0_dp
    end type motif_atom_type

    type :: motif_data_type
        type(motif_atom_type), allocatable :: atoms(:)
        integer :: natoms = 0
        character(len=1024), allocatable :: tail_lines(:)
        integer :: ntail = 0
        integer :: min_connect_index = 0
    end type motif_data_type

    interface write_vasp_configuration
        module procedure write_vasp_configuration_dp
        module procedure write_vasp_configuration_sp
    end interface

    public :: motif_atom_type
    public :: motif_data_type
    public :: read_motif_file
    public :: write_vasp_configuration
    public :: write_cif_configuration
    public :: write_gulp_configuration
    public :: write_lammps_configuration
    public :: write_castep_configuration
    public :: write_qe_configuration

contains

    subroutine write_vasp_configuration_dp(filename, title, cell_params, symbols, spat, coords, target_species, &
        replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: title
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        integer :: unit_id, ios
        integer :: nsp, nat, npos, sp, at, pos
        integer, allocatable :: counts(:)
        character(len=16), allocatable :: ordered_symbols(:)
        logical, allocatable :: is_selected(:)
        real(dp) :: lattice_vectors(3,3)
        integer :: nmotif, nuniq_motif, m, j
        character(len=16), allocatable :: motif_unique(:)
        integer, allocatable :: motif_ucounts(:)
        logical :: motif_found

        nsp = size(symbols)
        nat = size(spat)
        npos = size(target_positions)

        nmotif = 0
        if (present(motif_count) .and. present(motif_atoms)) nmotif = motif_count

        nuniq_motif = 0
        if (nmotif > 0) then
            allocate(motif_unique(nmotif))
            allocate(motif_ucounts(nmotif))
            motif_ucounts = 0
            do m = 1, nmotif
                motif_found = .false.
                do j = 1, nuniq_motif
                    if (trim(motif_atoms(m)%symbol) == trim(motif_unique(j))) then
                        motif_ucounts(j) = motif_ucounts(j) + 1
                        motif_found = .true.
                        exit
                    end if
                end do
                if (.not. motif_found) then
                    nuniq_motif = nuniq_motif + 1
                    motif_unique(nuniq_motif) = motif_atoms(m)%symbol
                    motif_ucounts(nuniq_motif) = 1
                end if
            end do
        end if

        allocate(counts(nsp + 1 + nuniq_motif))
        allocate(ordered_symbols(nsp + 1 + nuniq_motif))
        allocate(is_selected(npos))
        counts = 0
        ordered_symbols = ''
        is_selected = .false.

        do pos = 1, size(selected_sites)
            if (selected_sites(pos) >= 1 .and. selected_sites(pos) <= npos) then
                is_selected(selected_sites(pos)) = .true.
            end if
        end do

        do sp = 1, nsp
            if (sp < target_species) then
                ordered_symbols(sp) = adjustl(symbols(sp))
                counts(sp) = count(spat == sp)
            else if (sp == target_species) then
                ordered_symbols(sp) = adjustl(replacement_symbols(1))
                ordered_symbols(sp + 1) = adjustl(replacement_symbols(2))
                counts(sp) = count(is_selected)
                counts(sp + 1) = npos - counts(sp)
            else
                ordered_symbols(sp + 1) = adjustl(symbols(sp))
                counts(sp + 1) = count(spat == sp)
            end if
        end do

        do j = 1, nuniq_motif
            ordered_symbols(nsp + 1 + j) = adjustl(motif_unique(j))
            counts(nsp + 1 + j) = motif_ucounts(j)
        end do

        call cell(lattice_vectors, cell_params(1), cell_params(2), cell_params(3), cell_params(4), cell_params(5), &
            cell_params(6))

        open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) stop 1

        write(unit_id,'(A)') trim(title)
        write(unit_id,'(A)') '1.00000000'
        write(unit_id,'(3(f10.6,2x))') lattice_vectors(1,1), lattice_vectors(2,1), lattice_vectors(3,1)
        write(unit_id,'(3(f10.6,2x))') lattice_vectors(1,2), lattice_vectors(2,2), lattice_vectors(3,2)
        write(unit_id,'(3(f10.6,2x))') lattice_vectors(1,3), lattice_vectors(2,3), lattice_vectors(3,3)
        write(unit_id,*) ordered_symbols(1:nsp+1+nuniq_motif)
        write(unit_id,'(20(i4,1x))') counts(1:nsp+1+nuniq_motif)
        write(unit_id,'(A)') 'Direct'

        do at = 1, nat
            if (spat(at) < target_species) then
                write(unit_id,'(3(f10.6,2x))') coords(at,1), coords(at,2), coords(at,3)
            end if
        end do

        do pos = 1, npos
            if (is_selected(pos)) then
                at = target_positions(pos)
                write(unit_id,'(3(f10.6,2x))') coords(at,1), coords(at,2), coords(at,3)
            end if
        end do

        do pos = 1, npos
            if (.not. is_selected(pos)) then
                at = target_positions(pos)
                write(unit_id,'(3(f10.6,2x))') coords(at,1), coords(at,2), coords(at,3)
            end if
        end do

        do at = 1, nat
            if (spat(at) > target_species) then
                write(unit_id,'(3(f10.6,2x))') coords(at,1), coords(at,2), coords(at,3)
            end if
        end do

        do j = 1, nuniq_motif
            do m = 1, nmotif
                if (trim(motif_atoms(m)%symbol) == trim(motif_unique(j))) then
                    write(unit_id,'(3(f10.6,2x))') motif_atoms(m)%x, motif_atoms(m)%y, motif_atoms(m)%z
                end if
            end do
        end do

        close(unit_id)
        if (allocated(motif_ucounts)) deallocate(motif_ucounts)
        if (allocated(motif_unique)) deallocate(motif_unique)
        deallocate(is_selected)
        deallocate(ordered_symbols)
        deallocate(counts)
    end subroutine write_vasp_configuration_dp

    subroutine write_vasp_configuration_sp(filename, title, cell_params, symbols, spat, coords, target_species, &
        replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: title
        real, intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real, intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        real(dp) :: cell_params_dp(6)
        real(dp), allocatable :: coords_dp(:,:)

        cell_params_dp = real(cell_params, dp)
        allocate(coords_dp(size(coords,1), size(coords,2)))
        coords_dp = real(coords, dp)

        call write_vasp_configuration_dp(filename, title, cell_params_dp, symbols, spat, coords_dp, target_species, &
            replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)

        deallocate(coords_dp)
    end subroutine write_vasp_configuration_sp

    subroutine write_cif_configuration(filename, config_index, cell_params, symbols, spat, coords, target_species, &
        replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: config_index
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        integer :: unit_id, ios
        integer :: nat, npos, at, sp
        integer, allocatable :: target_lookup(:)
        logical, allocatable :: is_selected(:)
        character(len=16) :: atom_symbol

        nat = size(spat)
        npos = size(target_positions)
        allocate(target_lookup(nat))
        allocate(is_selected(npos))
        call build_selection_lookup(nat, npos, target_positions, selected_sites, target_lookup, is_selected)

        open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: could not create CIF file', trim(filename)
            stop 1
        end if

        write(unit_id,'("data_c",I0)') config_index
        write(unit_id,'(A)') ' '
        write(unit_id,'(A,F12.6)') '_cell_length_a     ', cell_params(1)
        write(unit_id,'(A,F12.6)') '_cell_length_b     ', cell_params(2)
        write(unit_id,'(A,F12.6)') '_cell_length_c     ', cell_params(3)
        write(unit_id,'(A,F12.6)') '_cell_angle_alpha  ', cell_params(4)
        write(unit_id,'(A,F12.6)') '_cell_angle_beta   ', cell_params(5)
        write(unit_id,'(A,F12.6)') '_cell_angle_gamma  ', cell_params(6)
        write(unit_id,'(A)') ' '
        write(unit_id,'(A)') "_symmetry_space_group_name_H-M  'P 1'"
        write(unit_id,'(A)') '_symmetry_Int_Tables_number      1'
        write(unit_id,'(A)') ' '
        write(unit_id,'(A)') 'loop_'
        write(unit_id,'(A)') '_symmetry_equiv_pos_as_xyz'
        write(unit_id,'(A)') "'x, y, z'"
        write(unit_id,'(A)') ' '
        write(unit_id,'(A)') 'loop_'
        write(unit_id,'(A)') '_atom_site_label'
        write(unit_id,'(A)') '_atom_site_type_symbol'
        write(unit_id,'(A)') '_atom_site_fract_x'
        write(unit_id,'(A)') '_atom_site_fract_y'
        write(unit_id,'(A)') '_atom_site_fract_z'

        do at = 1, nat
            sp = spat(at)
            atom_symbol = configuration_symbol(sp, at, target_species, replacement_symbols, target_lookup, is_selected, symbols)
            write(unit_id,'(A,I0,2X,A,3(2X,F11.7))') trim(atom_symbol), at, trim(atom_symbol), &
                coords(at,1), coords(at,2), coords(at,3)
        end do
        if (present(motif_count) .and. present(motif_atoms)) then
            do at = 1, motif_count
                write(unit_id,'(A,I0,2X,A,3(2X,F11.7))') trim(motif_atoms(at)%symbol), nat + at, &
                    trim(motif_atoms(at)%symbol), motif_atoms(at)%x, motif_atoms(at)%y, motif_atoms(at)%z
            end do
        end if

        close(unit_id)
        deallocate(is_selected)
        deallocate(target_lookup)
    end subroutine write_cif_configuration

    subroutine write_castep_configuration(filename, template_filename, config_tag, cell_params, symbols, spat, coords, &
        target_species, replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: template_filename
        character(len=*), intent(in) :: config_tag
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        integer :: unit_in, unit_out, ios, at, sp
        integer :: nat, npos
        integer, allocatable :: target_lookup(:)
        logical, allocatable :: is_selected(:)
        real(dp) :: lattice_vectors(3,3)
        logical :: inserted
        character(len=1024) :: line
        character(len=16) :: atom_symbol

        nat = size(spat)
        npos = size(target_positions)
        allocate(target_lookup(nat))
        allocate(is_selected(npos))
        call build_selection_lookup(nat, npos, target_positions, selected_sites, target_lookup, is_selected)

        call cell(lattice_vectors, cell_params(1), cell_params(2), cell_params(3), cell_params(4), cell_params(5), &
            cell_params(6))

        open(newunit=unit_in, file=trim(template_filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: template file not found', trim(template_filename)
            stop 1
        end if
        open(newunit=unit_out, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: could not create CASTEP file', trim(filename)
            stop 1
        end if

        inserted = .false.
        do
            read(unit_in,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (trim(adjustl(line)) == '@configuration_structure@') then
                inserted = .true.
                write(unit_out,'(A)') '%BLOCK lattice_cart'
                write(unit_out,'(3(F10.6,2X))') lattice_vectors(1,1), lattice_vectors(2,1), lattice_vectors(3,1)
                write(unit_out,'(3(F10.6,2X))') lattice_vectors(1,2), lattice_vectors(2,2), lattice_vectors(3,2)
                write(unit_out,'(3(F10.6,2X))') lattice_vectors(1,3), lattice_vectors(2,3), lattice_vectors(3,3)
                write(unit_out,'(A)') '%ENDBLOCK lattice_cart'
                write(unit_out,'(A)') '%BLOCK positions_frac'
                do at = 1, nat
                    sp = spat(at)
                    atom_symbol = configuration_symbol(sp, at, target_species, replacement_symbols, target_lookup, is_selected, symbols)
                    write(unit_out,'(A,3(2X,F11.7))') trim(atom_symbol), coords(at,1), coords(at,2), coords(at,3)
                end do
                if (present(motif_count) .and. present(motif_atoms)) then
                    do at = 1, motif_count
                        write(unit_out,'(A,3(2X,F11.7))') trim(motif_atoms(at)%symbol), &
                            motif_atoms(at)%x, motif_atoms(at)%y, motif_atoms(at)%z
                    end do
                end if
                write(unit_out,'(A)') '%ENDBLOCK positions_frac'
            else
                if (index(line, '@configuration_structure@') /= 0) then
                    write(error_unit,'(A,1X,A)') 'Error: @configuration_structure@ must appear alone in template', &
                        trim(template_filename)
                    stop 1
                end if
                call replace_all_tokens(line, '@configuration_number@', trim(config_tag))
                write(unit_out,'(A)') trim(line)
            end if
        end do
        if (.not. inserted) then
            write(error_unit,'(A,1X,A)') 'Error: @configuration_structure@ not found in template', trim(template_filename)
            stop 1
        end if

        close(unit_out)
        close(unit_in)
        deallocate(is_selected)
        deallocate(target_lookup)
    end subroutine write_castep_configuration

    subroutine write_qe_configuration(filename, template_filename, config_tag, cell_params, symbols, spat, coords, &
        target_species, replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: template_filename
        character(len=*), intent(in) :: config_tag
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        integer :: unit_in, unit_out, ios, at, sp
        integer :: nat, npos
        integer, allocatable :: target_lookup(:)
        logical, allocatable :: is_selected(:)
        real(dp) :: lattice_vectors(3,3)
        logical :: inserted
        character(len=1024) :: line
        character(len=16) :: atom_symbol

        nat = size(spat)
        npos = size(target_positions)
        allocate(target_lookup(nat))
        allocate(is_selected(npos))
        call build_selection_lookup(nat, npos, target_positions, selected_sites, target_lookup, is_selected)

        call cell(lattice_vectors, cell_params(1), cell_params(2), cell_params(3), cell_params(4), cell_params(5), &
            cell_params(6))

        open(newunit=unit_in, file=trim(template_filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: template file not found', trim(template_filename)
            stop 1
        end if
        open(newunit=unit_out, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: could not create Quantum ESPRESSO file', trim(filename)
            stop 1
        end if

        inserted = .false.
        do
            read(unit_in,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (trim(adjustl(line)) == '@configuration_structure@') then
                inserted = .true.
                write(unit_out,'(A)') 'CELL_PARAMETERS {angstrom}'
                write(unit_out,'(3(F10.6,2X))') lattice_vectors(1,1), lattice_vectors(2,1), lattice_vectors(3,1)
                write(unit_out,'(3(F10.6,2X))') lattice_vectors(1,2), lattice_vectors(2,2), lattice_vectors(3,2)
                write(unit_out,'(3(F10.6,2X))') lattice_vectors(1,3), lattice_vectors(2,3), lattice_vectors(3,3)
                write(unit_out,'(A)') 'ATOMIC_POSITIONS {crystal}'
                do at = 1, nat
                    sp = spat(at)
                    atom_symbol = configuration_symbol(sp, at, target_species, replacement_symbols, target_lookup, is_selected, symbols)
                    write(unit_out,'(A,3(2X,F11.7))') trim(atom_symbol), coords(at,1), coords(at,2), coords(at,3)
                end do
                if (present(motif_count) .and. present(motif_atoms)) then
                    do at = 1, motif_count
                        write(unit_out,'(A,3(2X,F11.7))') trim(motif_atoms(at)%symbol), &
                            motif_atoms(at)%x, motif_atoms(at)%y, motif_atoms(at)%z
                    end do
                end if
            else
                if (index(line, '@configuration_structure@') /= 0) then
                    write(error_unit,'(A,1X,A)') 'Error: @configuration_structure@ must appear alone in template', &
                        trim(template_filename)
                    stop 1
                end if
                call replace_all_tokens(line, '@configuration_number@', trim(config_tag))
                write(unit_out,'(A)') trim(line)
            end if
        end do
        if (.not. inserted) then
            write(error_unit,'(A,1X,A)') 'Error: @configuration_structure@ not found in template', trim(template_filename)
            stop 1
        end if

        close(unit_out)
        close(unit_in)
        deallocate(is_selected)
        deallocate(target_lookup)
    end subroutine write_qe_configuration

    subroutine write_lammps_configuration(data_filename, input_filename, template_filename, config_tag, cell_params, &
        symbols, spat, coords, target_species, replacement_symbols, target_positions, selected_sites, motif_atoms, motif_count)
        character(len=*), intent(in) :: data_filename
        character(len=*), intent(in) :: input_filename
        character(len=*), intent(in) :: template_filename
        character(len=*), intent(in) :: config_tag
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        integer, parameter :: max_line_len = 1024
        integer :: nat, nsp, npos, nlines, unit_data, unit_input, ios
        integer :: i, at, sp, map_count, angle_map_count, atom_id, bond_id, mol_id, angle_id
        integer :: natoms_lammps, nbonds_lammps, nangles_lammps, natom_types, nbond_types, nangle_types
        integer :: current_species_index
        integer, allocatable :: target_lookup(:)
        logical, allocatable :: is_selected(:)
        character(len=max_line_len), allocatable :: template_lines(:)
        character(len=32) :: atom_style
        character(len=16), allocatable :: all_symbols(:), map_symbols(:), map_roles(:)
        character(len=16), allocatable :: angle_center_symbols(:), angle_outer1_symbols(:), angle_outer2_symbols(:)
        integer, allocatable :: map_types(:), map_bonds(:)
        integer, allocatable :: angle_type_ids(:)
        integer, allocatable :: core_type(:), shell_type(:), shell_bond(:)
        integer, allocatable :: atom_species_index(:), core_id(:), shell_id(:), mol_id_per_atom(:)
        integer, allocatable :: angle_atom1(:), angle_atom2(:), angle_atom3(:), angle_type(:)
        real(dp) :: lattice_vectors(3,3), xy_tilt, xz_tilt, yz_tilt, xc, yc, zc
        real(dp), allocatable :: angle_cutoff1(:), angle_cutoff2(:), angle_cutoff12(:)
        logical :: is_triclinic, has_shells
        character(len=max_line_len) :: line
        integer :: nmotif, m
        integer, allocatable :: motif_core_type_lmp(:), motif_shell_type_lmp(:), motif_shell_bond_lmp(:)
        integer, allocatable :: motif_core_id_lmp(:), motif_shell_id_lmp(:), motif_mol_id_lmp(:)

        nat = size(spat)
        nsp = size(symbols)
        npos = size(target_positions)
        allocate(target_lookup(nat))
        allocate(is_selected(npos))
        call build_selection_lookup(nat, npos, target_positions, selected_sites, target_lookup, is_selected)

        call read_text_template(template_filename, template_lines, nlines)
        call parse_lammps_template(template_lines, nlines, atom_style, map_symbols, map_roles, map_types, map_bonds, map_count)
        call parse_lammps_angle_maps(template_lines, nlines, angle_center_symbols, angle_outer1_symbols, angle_outer2_symbols, &
            angle_type_ids, angle_cutoff1, angle_cutoff2, angle_cutoff12, angle_map_count)

        allocate(all_symbols(nsp + 1))
        allocate(core_type(nsp + 1))
        allocate(shell_type(nsp + 1))
        allocate(shell_bond(nsp + 1))
        call build_output_species(symbols, target_species, replacement_symbols, all_symbols)
        core_type = 0
        shell_type = 0
        shell_bond = 0
        do i = 1, size(all_symbols)
            do sp = 1, map_count
                if (trim(all_symbols(i)) == trim(map_symbols(sp))) then
                    if (trim(map_roles(sp)) == 'core') then
                        core_type(i) = map_types(sp)
                    else if (trim(map_roles(sp)) == 'shell') then
                        shell_type(i) = map_types(sp)
                        shell_bond(i) = map_bonds(sp)
                    end if
                end if
            end do
            if (core_type(i) == 0) then
                write(error_unit,'(A,1X,A,1X,A)') 'Error: missing LAMMPS core type mapping for species', &
                    trim(all_symbols(i)), 'in template_in.lammps.'
                stop 1
            end if
        end do

        has_shells = any(shell_type > 0)
        if (has_shells .and. trim(atom_style) /= 'full') then
            write(error_unit,'(A)') 'Error: shell species in LAMMPS require atom_style full.'
            stop 1
        end if

        natom_types = 0
        nbond_types = 0
        nangle_types = 0
        if (map_count > 0) natom_types = maxval(map_types(1:map_count))
        if (map_count > 0) nbond_types = maxval(map_bonds(1:map_count))
        if (angle_map_count > 0) nangle_types = maxval(angle_type_ids(1:angle_map_count))

        allocate(atom_species_index(nat))
        allocate(core_id(nat))
        allocate(shell_id(nat))
        allocate(mol_id_per_atom(nat))
        atom_id = 0
        mol_id = 0
        do at = 1, nat
            atom_species_index(at) = find_species_index(configuration_symbol(spat(at), at, target_species, replacement_symbols, &
                target_lookup, is_selected, symbols), all_symbols)
            if (atom_species_index(at) <= 0) then
                write(error_unit,'(A,I0)') 'Error: internal failure while mapping LAMMPS species for atom ', at
                stop 1
            end if
            atom_id = atom_id + 1
            core_id(at) = atom_id
            if (shell_type(atom_species_index(at)) > 0) then
                mol_id = mol_id + 1
                mol_id_per_atom(at) = mol_id
                atom_id = atom_id + 1
                shell_id(at) = atom_id
            else
                mol_id_per_atom(at) = 0
                shell_id(at) = 0
            end if
        end do
        natoms_lammps = atom_id
        nbonds_lammps = mol_id

        nmotif = 0
        if (present(motif_count) .and. present(motif_atoms)) nmotif = motif_count
        if (nmotif > 0) then
            allocate(motif_core_type_lmp(nmotif))
            allocate(motif_shell_type_lmp(nmotif))
            allocate(motif_shell_bond_lmp(nmotif))
            allocate(motif_core_id_lmp(nmotif))
            allocate(motif_shell_id_lmp(nmotif))
            allocate(motif_mol_id_lmp(nmotif))
            do m = 1, nmotif
                motif_core_type_lmp(m) = 0
                motif_shell_type_lmp(m) = 0
                motif_shell_bond_lmp(m) = 0
                do sp = 1, map_count
                    if (trim(motif_atoms(m)%symbol) == trim(map_symbols(sp))) then
                        if (trim(map_roles(sp)) == 'core') then
                            motif_core_type_lmp(m) = map_types(sp)
                        else if (trim(map_roles(sp)) == 'shell') then
                            motif_shell_type_lmp(m) = map_types(sp)
                            motif_shell_bond_lmp(m) = map_bonds(sp)
                        end if
                    end if
                end do
                if (motif_core_type_lmp(m) == 0) then
                    write(error_unit,'(A,1X,A,1X,A)') 'Error: missing LAMMPS core type mapping for motif species', &
                        trim(motif_atoms(m)%symbol), 'in template_in.lammps.'
                    stop 1
                end if
                atom_id = atom_id + 1
                motif_core_id_lmp(m) = atom_id
                if (motif_shell_type_lmp(m) > 0) then
                    mol_id = mol_id + 1
                    motif_mol_id_lmp(m) = mol_id
                    atom_id = atom_id + 1
                    motif_shell_id_lmp(m) = atom_id
                else
                    motif_mol_id_lmp(m) = 0
                    motif_shell_id_lmp(m) = 0
                end if
            end do
            natoms_lammps = atom_id
            nbonds_lammps = mol_id
        end if

        call cell(lattice_vectors, cell_params(1), cell_params(2), cell_params(3), cell_params(4), cell_params(5), &
            cell_params(6))
        xy_tilt = lattice_vectors(1,2)
        xz_tilt = lattice_vectors(1,3)
        yz_tilt = lattice_vectors(2,3)
        is_triclinic = (abs(xy_tilt) > 1.0e-6_dp .or. abs(xz_tilt) > 1.0e-6_dp .or. abs(yz_tilt) > 1.0e-6_dp)
        call build_lammps_angles(coords, lattice_vectors, all_symbols, atom_species_index, core_id, shell_id, &
            angle_center_symbols, angle_outer1_symbols, angle_outer2_symbols, angle_type_ids, angle_cutoff1, angle_cutoff2, &
            angle_cutoff12, angle_map_count, angle_atom1, angle_atom2, angle_atom3, angle_type, nangles_lammps)

        open(newunit=unit_data, file=trim(data_filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: could not create LAMMPS data file', trim(data_filename)
            stop 1
        end if
        write(unit_data,'(A,1X,A)') 'LAMMPS data file generated by SOD - configuration', trim(config_tag)
        write(unit_data,'(A)') ' '
        write(unit_data,'(I0,1X,A)') natoms_lammps, 'atoms'
        if (nbonds_lammps > 0) write(unit_data,'(I0,1X,A)') nbonds_lammps, 'bonds'
        if (nangles_lammps > 0) write(unit_data,'(I0,1X,A)') nangles_lammps, 'angles'
        write(unit_data,'(I0,1X,A)') natom_types, 'atom types'
        if (nbond_types > 0) write(unit_data,'(I0,1X,A)') nbond_types, 'bond types'
        if (nangles_lammps > 0) write(unit_data,'(I0,1X,A)') nangle_types, 'angle types'
        write(unit_data,'(A)') ' '
        write(unit_data,'(2(F14.6,2X),A)') 0.0_dp, lattice_vectors(1,1), 'xlo xhi'
        write(unit_data,'(2(F14.6,2X),A)') 0.0_dp, lattice_vectors(2,2), 'ylo yhi'
        write(unit_data,'(2(F14.6,2X),A)') 0.0_dp, lattice_vectors(3,3), 'zlo zhi'
        if (is_triclinic) then
            write(unit_data,'(3(F14.6,2X),A)') xy_tilt, xz_tilt, yz_tilt, 'xy xz yz'
        end if
        write(unit_data,'(A)') ' '

        select case (trim(atom_style))
        case ('atomic')
            write(unit_data,'(A)') 'Atoms  # atomic'
        case ('charge')
            write(unit_data,'(A)') 'Atoms  # charge'
        case default
            write(unit_data,'(A)') 'Atoms  # full'
        end select
        write(unit_data,'(A)') ' '

        do at = 1, nat
            current_species_index = atom_species_index(at)
            xc = lattice_vectors(1,1)*coords(at,1) + lattice_vectors(1,2)*coords(at,2) + lattice_vectors(1,3)*coords(at,3)
            yc = lattice_vectors(2,2)*coords(at,2) + lattice_vectors(2,3)*coords(at,3)
            zc = lattice_vectors(3,3)*coords(at,3)
            select case (trim(atom_style))
            case ('atomic')
                write(unit_data,'(I0,2X,I0,3(2X,F14.6))') core_id(at), core_type(current_species_index), xc, yc, zc
            case ('charge')
                write(unit_data,'(I0,2X,I0,2X,F8.4,3(2X,F14.6))') core_id(at), core_type(current_species_index), 0.0_dp, &
                    xc, yc, zc
            case default
                write(unit_data,'(I0,2X,I0,2X,I0,2X,F8.4,3(2X,F14.6))') core_id(at), mol_id_per_atom(at), &
                    core_type(current_species_index), 0.0_dp, xc, yc, zc
                if (shell_id(at) > 0) then
                    write(unit_data,'(I0,2X,I0,2X,I0,2X,F8.4,3(2X,F14.6))') shell_id(at), mol_id_per_atom(at), &
                        shell_type(current_species_index), 0.0_dp, xc, yc, zc
                end if
            end select
        end do

        do m = 1, nmotif
            xc = lattice_vectors(1,1)*motif_atoms(m)%x + lattice_vectors(1,2)*motif_atoms(m)%y + &
                lattice_vectors(1,3)*motif_atoms(m)%z
            yc = lattice_vectors(2,2)*motif_atoms(m)%y + lattice_vectors(2,3)*motif_atoms(m)%z
            zc = lattice_vectors(3,3)*motif_atoms(m)%z
            select case (trim(atom_style))
            case ('atomic')
                write(unit_data,'(I0,2X,I0,3(2X,F14.6))') motif_core_id_lmp(m), motif_core_type_lmp(m), xc, yc, zc
            case ('charge')
                write(unit_data,'(I0,2X,I0,2X,F8.4,3(2X,F14.6))') motif_core_id_lmp(m), motif_core_type_lmp(m), &
                    0.0_dp, xc, yc, zc
            case default
                write(unit_data,'(I0,2X,I0,2X,I0,2X,F8.4,3(2X,F14.6))') motif_core_id_lmp(m), motif_mol_id_lmp(m), &
                    motif_core_type_lmp(m), 0.0_dp, xc, yc, zc
                if (motif_shell_id_lmp(m) > 0) then
                    write(unit_data,'(I0,2X,I0,2X,I0,2X,F8.4,3(2X,F14.6))') motif_shell_id_lmp(m), motif_mol_id_lmp(m), &
                        motif_shell_type_lmp(m), 0.0_dp, xc, yc, zc
                end if
            end select
        end do

        if (nbonds_lammps > 0) then
            write(unit_data,'(A)') ' '
            write(unit_data,'(A)') 'Bonds'
            write(unit_data,'(A)') ' '
            bond_id = 0
            do at = 1, nat
                current_species_index = atom_species_index(at)
                if (shell_id(at) > 0) then
                    bond_id = bond_id + 1
                    write(unit_data,'(I0,2X,I0,2X,I0,2X,I0)') bond_id, shell_bond(current_species_index), core_id(at), &
                        shell_id(at)
                end if
            end do
            do m = 1, nmotif
                if (motif_shell_id_lmp(m) > 0) then
                    bond_id = bond_id + 1
                    write(unit_data,'(I0,2X,I0,2X,I0,2X,I0)') bond_id, motif_shell_bond_lmp(m), &
                        motif_core_id_lmp(m), motif_shell_id_lmp(m)
                end if
            end do
        end if
        if (nangles_lammps > 0) then
            write(unit_data,'(A)') ' '
            write(unit_data,'(A)') 'Angles'
            write(unit_data,'(A)') ' '
            do angle_id = 1, nangles_lammps
                write(unit_data,'(I0,2X,I0,2X,I0,2X,I0,2X,I0)') angle_id, angle_type(angle_id), angle_atom1(angle_id), &
                    angle_atom2(angle_id), angle_atom3(angle_id)
            end do
        end if
        close(unit_data)

        open(newunit=unit_input, file=trim(input_filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: could not create LAMMPS input file', trim(input_filename)
            stop 1
        end if
        do i = 1, nlines
            line = template_lines(i)
            if (index(adjustl(line), '# sod_type_map ') == 1) cycle
            if (index(adjustl(line), '# sod_angle_map ') == 1) cycle
            call replace_all_tokens(line, '@configuration_number@', trim(config_tag))
            if (index(adjustl(line), 'read_data') == 1) then
                write(unit_input,'(A)') 'read_data       ' // trim(path_basename(data_filename))
            else
                write(unit_input,'(A)') trim(line)
            end if
        end do
        write(unit_input,'(A)') 'variable sod_final_pe equal pe'
        write(unit_input,'(A)') 'print "SOD_FINAL_ENERGY ${sod_final_pe}"'
        write(unit_input,'(A)') 'write_data ' // trim(path_basename(data_filename(1:len_trim(data_filename)-5)//'.min.data'))
        close(unit_input)

        deallocate(mol_id_per_atom)
        deallocate(shell_id)
        deallocate(core_id)
        deallocate(atom_species_index)
        deallocate(shell_bond)
        deallocate(shell_type)
        deallocate(core_type)
        deallocate(all_symbols)
        if (allocated(motif_mol_id_lmp)) deallocate(motif_mol_id_lmp)
        if (allocated(motif_shell_id_lmp)) deallocate(motif_shell_id_lmp)
        if (allocated(motif_core_id_lmp)) deallocate(motif_core_id_lmp)
        if (allocated(motif_shell_bond_lmp)) deallocate(motif_shell_bond_lmp)
        if (allocated(motif_shell_type_lmp)) deallocate(motif_shell_type_lmp)
        if (allocated(motif_core_type_lmp)) deallocate(motif_core_type_lmp)
        deallocate(map_bonds)
        deallocate(map_types)
        deallocate(map_roles)
        deallocate(map_symbols)
        if (allocated(angle_type)) deallocate(angle_type)
        if (allocated(angle_atom3)) deallocate(angle_atom3)
        if (allocated(angle_atom2)) deallocate(angle_atom2)
        if (allocated(angle_atom1)) deallocate(angle_atom1)
        if (allocated(angle_cutoff12)) deallocate(angle_cutoff12)
        if (allocated(angle_cutoff2)) deallocate(angle_cutoff2)
        if (allocated(angle_cutoff1)) deallocate(angle_cutoff1)
        if (allocated(angle_type_ids)) deallocate(angle_type_ids)
        if (allocated(angle_outer2_symbols)) deallocate(angle_outer2_symbols)
        if (allocated(angle_outer1_symbols)) deallocate(angle_outer1_symbols)
        if (allocated(angle_center_symbols)) deallocate(angle_center_symbols)
        deallocate(template_lines)
        deallocate(is_selected)
        deallocate(target_lookup)
    end subroutine write_lammps_configuration

    subroutine read_motif_file(filename, motif)
        character(len=*), intent(in) :: filename
        type(motif_data_type), intent(out) :: motif
        integer :: unit_id, ios, atom_cap, tail_cap, idx1, idx2, min_idx
        character(len=1024) :: line
        character(len=16) :: w1, w2
        real(dp) :: rx, ry, rz
        type(motif_atom_type), allocatable :: atmp(:)
        character(len=1024), allocatable :: ttmp(:)
        logical :: in_tail

        motif%natoms = 0
        motif%ntail = 0
        motif%min_connect_index = 0
        atom_cap = 64
        tail_cap = 64
        allocate(motif%atoms(atom_cap))
        allocate(motif%tail_lines(tail_cap))
        in_tail = .false.
        min_idx = huge(1)

        open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: motif file not found', trim(filename)
            stop 1
        end if

        do
            read(unit_id,'(A)',iostat=ios) line
            if (ios /= 0) exit

            if (.not. in_tail) then
                w1 = ''
                w2 = ''
                ! Check if this is a tail keyword (space, connect, etc.)
                line = adjustl(line)
                read(line,*,iostat=ios) w1
                if (ios == 0) then
                    if (trim(w1) == 'space' .or. trim(w1) == 'connect') then
                        in_tail = .true.
                    end if
                end if
            end if

            if (in_tail) then
                motif%ntail = motif%ntail + 1
                if (motif%ntail > tail_cap) then
                    tail_cap = tail_cap * 2
                    allocate(ttmp(tail_cap))
                    ttmp(1:motif%ntail-1) = motif%tail_lines(1:motif%ntail-1)
                    call move_alloc(ttmp, motif%tail_lines)
                end if
                motif%tail_lines(motif%ntail) = line
                ! Parse connect indices
                if (index(adjustl(line), 'connect') == 1) then
                    read(line,*,iostat=ios) w1, idx1, idx2
                    if (ios == 0) then
                        if (idx1 < min_idx) min_idx = idx1
                        if (idx2 < min_idx) min_idx = idx2
                    end if
                end if
                cycle
            end if

            ! Parse atom line: try GULP format (Symbol role x y z)
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle

            read(line,*,iostat=ios) w1, w2, rx, ry, rz
            if (ios == 0 .and. (trim(w2) == 'core' .or. trim(w2) == 'shel')) then
                motif%natoms = motif%natoms + 1
                if (motif%natoms > atom_cap) then
                    atom_cap = atom_cap * 2
                    allocate(atmp(atom_cap))
                    atmp(1:motif%natoms-1) = motif%atoms(1:motif%natoms-1)
                    call move_alloc(atmp, motif%atoms)
                end if
                motif%atoms(motif%natoms)%symbol = trim(w1)
                motif%atoms(motif%natoms)%role   = trim(w2)
                motif%atoms(motif%natoms)%x = rx
                motif%atoms(motif%natoms)%y = ry
                motif%atoms(motif%natoms)%z = rz
                cycle
            end if

            ! Try plain format (Symbol x y z)
            read(line,*,iostat=ios) w1, rx, ry, rz
            if (ios == 0) then
                motif%natoms = motif%natoms + 1
                if (motif%natoms > atom_cap) then
                    atom_cap = atom_cap * 2
                    allocate(atmp(atom_cap))
                    atmp(1:motif%natoms-1) = motif%atoms(1:motif%natoms-1)
                    call move_alloc(atmp, motif%atoms)
                end if
                motif%atoms(motif%natoms)%symbol = trim(w1)
                motif%atoms(motif%natoms)%role   = 'core'
                motif%atoms(motif%natoms)%x = rx
                motif%atoms(motif%natoms)%y = ry
                motif%atoms(motif%natoms)%z = rz
                cycle
            end if
        end do
        close(unit_id)

        if (min_idx < huge(1)) then
            motif%min_connect_index = min_idx
        end if

        ! Trim arrays to actual size
        if (motif%natoms > 0 .and. motif%natoms < atom_cap) then
            allocate(atmp(motif%natoms))
            atmp = motif%atoms(1:motif%natoms)
            call move_alloc(atmp, motif%atoms)
        end if
        if (motif%ntail > 0 .and. motif%ntail < tail_cap) then
            allocate(ttmp(motif%ntail))
            ttmp = motif%tail_lines(1:motif%ntail)
            call move_alloc(ttmp, motif%tail_lines)
        end if
    end subroutine read_motif_file

    subroutine write_gulp_configuration(filename, template_filename, config_tag, cell_params, symbols, spat, coords, &
        target_species, replacement_symbols, target_positions, selected_sites, motif, motif_atoms, motif_count)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: template_filename
        character(len=*), intent(in) :: config_tag
        real(dp), intent(in) :: cell_params(6)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: spat(:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        type(motif_data_type), intent(in), optional :: motif
        type(motif_atom_type), intent(in), optional :: motif_atoms(:)
        integer, intent(in), optional :: motif_count
        integer :: unit_out, ios, at, sp, i, nlines, nmotif
        integer :: nat, npos
        integer, allocatable :: target_lookup(:)
        logical, allocatable :: is_selected(:)
        character(len=1024), allocatable :: template_lines(:)
        character(len=1024) :: line
        character(len=16) :: atom_symbol
        character(len=16), allocatable :: shell_species(:)
        character(len=16), allocatable :: oxy_retype_symbols(:,:)
        integer :: nshell, noxy_rules
        logical :: inserted
        real(dp) :: lattice_vectors(3,3)

        nat = size(spat)
        npos = size(target_positions)
        nmotif = 0
        if (present(motif)) then
            nmotif = motif%natoms
        else if (present(motif_count)) then
            nmotif = motif_count
        end if
        allocate(target_lookup(nat))
        allocate(is_selected(npos))
        call build_selection_lookup(nat, npos, target_positions, selected_sites, target_lookup, is_selected)

        call read_text_template(template_filename, template_lines, nlines)
        call parse_gulp_shell_directives(template_lines, nlines, shell_species, nshell)
        call parse_gulp_oxygen_retype(template_lines, nlines, oxy_retype_symbols, noxy_rules)

        call cell(lattice_vectors, cell_params(1), cell_params(2), cell_params(3), &
            cell_params(4), cell_params(5), cell_params(6))

        open(newunit=unit_out, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: could not create GULP file', trim(filename)
            stop 1
        end if

        inserted = .false.
        do i = 1, nlines
            line = template_lines(i)
            if (trim(adjustl(line)) == '@configuration_structure@') then
                inserted = .true.
                write(unit_out,'(A)') 'cell'
                write(unit_out,'(6(F10.4,2X))') cell_params(1), cell_params(2), cell_params(3), &
                    cell_params(4), cell_params(5), cell_params(6)
                write(unit_out,'(A)') 'frac'
                do at = 1, nat
                    sp = spat(at)
                    atom_symbol = configuration_symbol(sp, at, target_species, replacement_symbols, &
                        target_lookup, is_selected, symbols)
                    if (noxy_rules > 0) then
                        atom_symbol = retype_oxygen(atom_symbol, at, nat, coords, lattice_vectors, &
                            spat, target_species, target_lookup, is_selected, replacement_symbols, &
                            oxy_retype_symbols, noxy_rules)
                    end if
                    write(unit_out,'(A3,2X,A4,2X,3(F11.7,2X))') trim(atom_symbol), 'core', &
                        coords(at,1), coords(at,2), coords(at,3)
                    if (species_has_shell(atom_symbol, shell_species, nshell) .and. &
                        .not. (present(motif) .and. motif%ntail > 0)) then
                        write(unit_out,'(A3,2X,A4,2X,3(F11.7,2X))') trim(atom_symbol), 'shel', &
                            coords(at,1), coords(at,2), coords(at,3)
                    end if
                end do
                ! Write motif atoms
                if (present(motif) .and. motif%natoms > 0) then
                    do at = 1, motif%natoms
                        write(unit_out,'(A3,2X,A4,2X,3(F11.7,2X))') trim(motif%atoms(at)%symbol), &
                            trim(motif%atoms(at)%role), motif%atoms(at)%x, motif%atoms(at)%y, motif%atoms(at)%z
                    end do
                else if (present(motif_atoms) .and. nmotif > 0) then
                    do at = 1, nmotif
                        write(unit_out,'(A3,2X,A4,2X,3(F11.7,2X))') trim(motif_atoms(at)%symbol), &
                            trim(motif_atoms(at)%role), motif_atoms(at)%x, motif_atoms(at)%y, motif_atoms(at)%z
                    end do
                end if
                ! Write motif tail (space, connect lines) with reindexed atom numbers
                if (present(motif) .and. motif%ntail > 0) then
                    call write_motif_tail(unit_out, motif, nat)
                end if
            else
                if (index(line, '@configuration_structure@') /= 0) then
                    write(error_unit,'(A,1X,A)') &
                        'Error: @configuration_structure@ must appear alone in template', trim(template_filename)
                    stop 1
                end if
                if (index(adjustl(line), '# sod_shell ') == 1) cycle
                if (index(adjustl(line), '# sod_oxygen_retype ') == 1) cycle
                call replace_all_tokens(line, '@configuration_number@', trim(config_tag))
                write(unit_out,'(A)') trim(line)
            end if
        end do
        if (.not. inserted) then
            write(error_unit,'(A,1X,A)') 'Error: @configuration_structure@ not found in template', trim(template_filename)
            stop 1
        end if

        close(unit_out)
        if (allocated(oxy_retype_symbols)) deallocate(oxy_retype_symbols)
        if (allocated(shell_species)) deallocate(shell_species)
        deallocate(template_lines)
        deallocate(is_selected)
        deallocate(target_lookup)
    end subroutine write_gulp_configuration

    subroutine parse_gulp_shell_directives(template_lines, nlines, shell_species, nshell)
        character(len=*), intent(in) :: template_lines(:)
        integer, intent(in) :: nlines
        character(len=16), allocatable, intent(out) :: shell_species(:)
        integer, intent(out) :: nshell
        integer :: i, ios, capacity, pos
        character(len=1024) :: line, rest
        character(len=16) :: word
        character(len=16), allocatable :: tmp(:)

        nshell = 0
        capacity = 8
        allocate(shell_species(capacity))

        do i = 1, nlines
            line = adjustl(template_lines(i))
            if (index(line, '# sod_shell ') /= 1) cycle
            rest = adjustl(line(13:))
            do while (len_trim(rest) > 0)
                read(rest,*,iostat=ios) word
                if (ios /= 0) exit
                nshell = nshell + 1
                if (nshell > capacity) then
                    capacity = capacity * 2
                    allocate(tmp(capacity))
                    tmp(1:nshell-1) = shell_species(1:nshell-1)
                    call move_alloc(tmp, shell_species)
                end if
                shell_species(nshell) = trim(word)
                pos = index(rest, trim(word))
                if (pos > 0) then
                    rest = adjustl(rest(pos + len_trim(word):))
                else
                    exit
                end if
            end do
        end do
    end subroutine parse_gulp_shell_directives

    pure logical function species_has_shell(symbol, shell_species, nshell)
        character(len=*), intent(in) :: symbol
        character(len=16), intent(in) :: shell_species(:)
        integer, intent(in) :: nshell
        integer :: i

        species_has_shell = .false.
        do i = 1, nshell
            if (trim(symbol) == trim(shell_species(i))) then
                species_has_shell = .true.
                return
            end if
        end do
    end function species_has_shell

    subroutine write_motif_tail(unit_out, motif, framework_atom_count)
        integer, intent(in) :: unit_out, framework_atom_count
        type(motif_data_type), intent(in) :: motif
        integer :: i, ios, idx1, idx2, new1, new2, offset
        character(len=1024) :: line
        character(len=16) :: keyword, bond_type
        integer :: ix, iy, iz

        if (motif%min_connect_index > 0) then
            offset = framework_atom_count + 1 - motif%min_connect_index
        else
            offset = 0
        end if

        do i = 1, motif%ntail
            line = adjustl(motif%tail_lines(i))
            if (index(line, 'connect') == 1) then
                read(line,*,iostat=ios) keyword, idx1, idx2, bond_type, ix, iy, iz
                if (ios == 0) then
                    new1 = idx1 + offset
                    new2 = idx2 + offset
                    write(unit_out,'(A,1X,I0,1X,I0,1X,A,2X,I0,1X,I0,1X,I0)') &
                        'connect', new1, new2, trim(bond_type), ix, iy, iz
                else
                    write(unit_out,'(A)') trim(motif%tail_lines(i))
                end if
            else
                write(unit_out,'(A)') trim(motif%tail_lines(i))
            end if
        end do
    end subroutine write_motif_tail

    ! Parse # sod_oxygen_retype directives from GULP template.
    ! Format: # sod_oxygen_retype O cutoff=2.2 O1=Ge,Ge O2=Si,Si O3=Si,Ge
    ! oxy_retype_symbols(rule, 1:3) = [label, neighbor1, neighbor2]
    subroutine parse_gulp_oxygen_retype(template_lines, nlines, oxy_retype_symbols, nrules)
        character(len=*), intent(in) :: template_lines(:)
        integer, intent(in) :: nlines
        character(len=16), allocatable, intent(out) :: oxy_retype_symbols(:,:)
        integer, intent(out) :: nrules
        integer :: i, pos, eq_pos, comma_pos
        character(len=1024) :: line, rest, token
        character(len=16) :: label, neigh1, neigh2
        character(len=16), allocatable :: tmp(:,:)
        integer :: cap

        nrules = 0
        cap = 8
        allocate(oxy_retype_symbols(cap, 3))

        do i = 1, nlines
            line = adjustl(template_lines(i))
            if (index(line, '# sod_oxygen_retype ') /= 1) cycle
            rest = adjustl(line(21:))
            ! Parse tokens like O1=Ge,Ge O2=Si,Si O3=Si,Ge
            ! Note: cannot use list-directed read because Fortran treats ',' as delimiter
            do while (len_trim(rest) > 0)
                pos = index(trim(rest), ' ')
                if (pos > 0) then
                    token = rest(1:pos-1)
                    rest = adjustl(rest(pos+1:))
                else
                    token = trim(rest)
                    rest = ''
                end if
                if (len_trim(token) == 0) exit
                eq_pos = index(token, '=')
                if (eq_pos > 1) then
                    label = token(1:eq_pos-1)
                    comma_pos = index(token(eq_pos+1:), ',')
                    if (comma_pos > 0) then
                        neigh1 = token(eq_pos+1:eq_pos+comma_pos-1)
                        neigh2 = token(eq_pos+comma_pos+1:)
                        nrules = nrules + 1
                        if (nrules > cap) then
                            cap = cap * 2
                            allocate(tmp(cap, 3))
                            tmp(1:nrules-1, :) = oxy_retype_symbols(1:nrules-1, :)
                            call move_alloc(tmp, oxy_retype_symbols)
                        end if
                        oxy_retype_symbols(nrules, 1) = trim(label)
                        oxy_retype_symbols(nrules, 2) = trim(neigh1)
                        oxy_retype_symbols(nrules, 3) = trim(neigh2)
                    end if
                end if
            end do
        end do
    end subroutine parse_gulp_oxygen_retype

    ! Reclassify oxygen label based on its two nearest T-site neighbors.
    function retype_oxygen(base_symbol, oxygen_idx, nat, coords, lattice_vectors, &
        spat, target_species, target_lookup, is_selected, replacement_symbols, &
        oxy_retype_symbols, nrules) result(new_symbol)
        character(len=16), intent(in) :: base_symbol
        integer, intent(in) :: oxygen_idx, nat, target_species, nrules
        real(dp), intent(in) :: coords(:,:), lattice_vectors(3,3)
        integer, intent(in) :: spat(:), target_lookup(:)
        logical, intent(in) :: is_selected(:)
        character(len=*), intent(in) :: replacement_symbols(2)
        character(len=16), intent(in) :: oxy_retype_symbols(:,:)
        character(len=16) :: new_symbol
        real(dp) :: d, d1, d2
        integer :: at, pos, t1_idx, t2_idx, rule
        character(len=16) :: t1_sym, t2_sym, n1, n2

        new_symbol = base_symbol
        ! Only retype atoms of the target species (oxygens in the framework)
        if (spat(oxygen_idx) == target_species) return

        ! Find two nearest T-site atoms (target_species atoms)
        d1 = huge(1.0_dp)
        d2 = huge(1.0_dp)
        t1_idx = 0
        t2_idx = 0
        do at = 1, nat
            if (spat(at) /= target_species) cycle
            d = fractional_distance(coords(oxygen_idx,1:3), coords(at,1:3), lattice_vectors)
            if (d < d1) then
                d2 = d1
                t2_idx = t1_idx
                d1 = d
                t1_idx = at
            else if (d < d2) then
                d2 = d
                t2_idx = at
            end if
        end do

        if (t1_idx == 0 .or. t2_idx == 0) return

        ! Determine T-site species labels
        pos = target_lookup(t1_idx)
        if (pos >= 1 .and. pos <= size(is_selected)) then
            if (is_selected(pos)) then
                t1_sym = trim(replacement_symbols(1))
            else
                t1_sym = trim(replacement_symbols(2))
            end if
        else
            t1_sym = trim(replacement_symbols(2))
        end if

        pos = target_lookup(t2_idx)
        if (pos >= 1 .and. pos <= size(is_selected)) then
            if (is_selected(pos)) then
                t2_sym = trim(replacement_symbols(1))
            else
                t2_sym = trim(replacement_symbols(2))
            end if
        else
            t2_sym = trim(replacement_symbols(2))
        end if

        ! Match against rules (order-independent: Si,Ge matches Ge,Si)
        do rule = 1, nrules
            n1 = trim(oxy_retype_symbols(rule, 2))
            n2 = trim(oxy_retype_symbols(rule, 3))
            if ((trim(t1_sym) == n1 .and. trim(t2_sym) == n2) .or. &
                (trim(t1_sym) == n2 .and. trim(t2_sym) == n1)) then
                new_symbol = trim(oxy_retype_symbols(rule, 1))
                return
            end if
        end do
    end function retype_oxygen

    subroutine build_selection_lookup(nat, npos, target_positions, selected_sites, target_lookup, is_selected)
        integer, intent(in) :: nat, npos
        integer, intent(in) :: target_positions(:)
        integer, intent(in) :: selected_sites(:)
        integer, intent(out) :: target_lookup(nat)
        logical, intent(out) :: is_selected(npos)
        integer :: pos

        target_lookup = 0
        is_selected = .false.
        do pos = 1, npos
            if (target_positions(pos) >= 1 .and. target_positions(pos) <= nat) then
                target_lookup(target_positions(pos)) = pos
            end if
        end do
        do pos = 1, size(selected_sites)
            if (selected_sites(pos) >= 1 .and. selected_sites(pos) <= npos) then
                is_selected(selected_sites(pos)) = .true.
            end if
        end do
    end subroutine build_selection_lookup

    subroutine read_text_template(filename, lines, line_count)
        character(len=*), intent(in) :: filename
        character(len=1024), allocatable, intent(out) :: lines(:)
        integer, intent(out) :: line_count
        integer :: unit_id, ios, i
        character(len=1024) :: line

        line_count = 0
        open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: template file not found', trim(filename)
            stop 1
        end if
        do
            read(unit_id,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line_count = line_count + 1
        end do
        rewind(unit_id)
        allocate(lines(line_count))
        do i = 1, line_count
            read(unit_id,'(A)',iostat=ios) lines(i)
            if (ios /= 0) then
                write(error_unit,'(A,1X,A)') 'Error: failed while reading template file', trim(filename)
                stop 1
            end if
        end do
        close(unit_id)
    end subroutine read_text_template

    subroutine parse_lammps_template(lines, line_count, atom_style, map_symbols, map_roles, map_types, map_bonds, map_count)
        character(len=*), intent(in) :: lines(:)
        integer, intent(in) :: line_count
        character(len=*), intent(out) :: atom_style
        character(len=16), allocatable, intent(out) :: map_symbols(:), map_roles(:)
        integer, allocatable, intent(out) :: map_types(:), map_bonds(:)
        integer, intent(out) :: map_count
        integer :: i, ios, bond_pos, type_value, bond_value
        character(len=1024) :: line, tail, masked
        character(len=32) :: token1, token2
        logical :: has_bond

        atom_style = ''
        map_count = 0
        do i = 1, line_count
            line = adjustl(lines(i))
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') then
                if (index(line, '# sod_type_map ') == 1) map_count = map_count + 1
                cycle
            end if
            read(line,*,iostat=ios) token1, token2
            if (ios == 0) then
                if (trim(token1) == 'atom_style') atom_style = trim(token2)
            end if
        end do
        if (len_trim(atom_style) == 0) then
            write(error_unit,'(A)') 'Error: atom_style not found in template_in.lammps.'
            stop 1
        end if

        allocate(map_symbols(max(1,map_count)))
        allocate(map_roles(max(1,map_count)))
        allocate(map_types(max(1,map_count)))
        allocate(map_bonds(max(1,map_count)))
        map_symbols = ''
        map_roles = ''
        map_types = 0
        map_bonds = 0
        map_count = 0
        do i = 1, line_count
            line = adjustl(lines(i))
            if (index(line, '# sod_type_map ') /= 1) cycle
            map_count = map_count + 1
            tail = adjustl(line(16:))
            read(tail,*,iostat=ios) map_symbols(map_count), map_roles(map_count)
            if (ios /= 0) then
                write(error_unit,'(A,1X,A)') 'Error: invalid sod_type_map line in template_in.lammps:', trim(lines(i))
                stop 1
            end if
            call read_integer_after_key(tail, 'bond_type=', bond_value, found=has_bond)
            masked = tail
            bond_pos = index(masked, 'bond_type=')
            if (has_bond .and. bond_pos > 0) then
                masked(bond_pos:min(len(masked), bond_pos + 9)) = ' '
            end if
            call read_integer_after_key(masked, 'type=', type_value)
            map_types(map_count) = type_value
            map_bonds(map_count) = max(0, bond_value)
        end do
    end subroutine parse_lammps_template

    subroutine parse_lammps_angle_maps(lines, line_count, center_symbols, outer1_symbols, outer2_symbols, type_ids, &
        cutoff1, cutoff2, cutoff12, angle_count)
        character(len=*), intent(in) :: lines(:)
        integer, intent(in) :: line_count
        character(len=16), allocatable, intent(out) :: center_symbols(:), outer1_symbols(:), outer2_symbols(:)
        integer, allocatable, intent(out) :: type_ids(:)
        real(dp), allocatable, intent(out) :: cutoff1(:), cutoff2(:), cutoff12(:)
        integer, intent(out) :: angle_count
        integer :: i, ios
        character(len=1024) :: line, tail
        character(len=16) :: c, o1, o2
        integer :: t_id
        real(dp) :: c1, c2, c12

        angle_count = 0
        do i = 1, line_count
            line = adjustl(lines(i))
            if (index(line, '# sod_angle_map ') == 1) angle_count = angle_count + 1
        end do

        allocate(center_symbols(max(1, angle_count)))
        allocate(outer1_symbols(max(1, angle_count)))
        allocate(outer2_symbols(max(1, angle_count)))
        allocate(type_ids(max(1, angle_count)))
        allocate(cutoff1(max(1, angle_count)))
        allocate(cutoff2(max(1, angle_count)))
        allocate(cutoff12(max(1, angle_count)))
        center_symbols = ''
        outer1_symbols = ''
        outer2_symbols = ''
        type_ids = 0
        cutoff1 = 0.0_dp
        cutoff2 = 0.0_dp
        cutoff12 = 0.0_dp

        angle_count = 0
        do i = 1, line_count
            line = adjustl(lines(i))
            if (index(line, '# sod_angle_map ') /= 1) cycle
            angle_count = angle_count + 1
            tail = adjustl(line(17:))
            read(tail,*,iostat=ios) c, o1, o2
            if (ios /= 0) then
                write(error_unit,'(A,1X,A)') 'Error: invalid sod_angle_map line in template_in.lammps:', trim(lines(i))
                stop 1
            end if
            call read_integer_after_key(tail, 'type=', t_id)
            call read_real_after_key(tail, 'cutoff1=', c1)
            call read_real_after_key(tail, 'cutoff2=', c2)
            call read_real_after_key(tail, 'cutoff12=', c12)
            center_symbols(angle_count) = trim(c)
            outer1_symbols(angle_count) = trim(o1)
            outer2_symbols(angle_count) = trim(o2)
            type_ids(angle_count) = t_id
            cutoff1(angle_count) = c1
            cutoff2(angle_count) = c2
            cutoff12(angle_count) = c12
        end do
    end subroutine parse_lammps_angle_maps

    subroutine build_lammps_angles(coords_frac, lattice_vectors, all_symbols, atom_species_index, core_id, shell_id, &
        center_symbols, outer1_symbols, outer2_symbols, type_ids, cutoff1, cutoff2, cutoff12, angle_map_count, &
        angle_atom1, angle_atom2, angle_atom3, angle_type, angle_count)
        real(dp), intent(in) :: coords_frac(:,:), lattice_vectors(3,3)
        character(len=*), intent(in) :: all_symbols(:)
        integer, intent(in) :: atom_species_index(:), core_id(:), shell_id(:)
        character(len=*), intent(in) :: center_symbols(:), outer1_symbols(:), outer2_symbols(:)
        integer, intent(in) :: type_ids(:), angle_map_count
        real(dp), intent(in) :: cutoff1(:), cutoff2(:), cutoff12(:)
        integer, allocatable, intent(out) :: angle_atom1(:), angle_atom2(:), angle_atom3(:), angle_type(:)
        integer, intent(out) :: angle_count
        integer :: nat, map_idx, center_at, outer_i, outer_j, idx1, idx2, capacity
        real(dp) :: d1, d2, d12

        nat = size(coords_frac, 1)
        angle_count = 0
        capacity = 16
        allocate(angle_atom1(capacity))
        allocate(angle_atom2(capacity))
        allocate(angle_atom3(capacity))
        allocate(angle_type(capacity))

        do map_idx = 1, angle_map_count
            do center_at = 1, nat
                if (trim(all_symbols(atom_species_index(center_at))) /= trim(center_symbols(map_idx))) cycle
                do outer_i = 1, nat
                    if (shell_id(outer_i) <= 0) cycle
                    if (trim(all_symbols(atom_species_index(outer_i))) /= trim(outer1_symbols(map_idx))) cycle
                    d1 = fractional_distance(coords_frac(center_at,1:3), coords_frac(outer_i,1:3), lattice_vectors)
                    if (d1 > cutoff1(map_idx)) cycle
                    do outer_j = 1, nat
                        if (shell_id(outer_j) <= 0) cycle
                        if (outer_j == outer_i) cycle
                        if (trim(all_symbols(atom_species_index(outer_j))) /= trim(outer2_symbols(map_idx))) cycle
                        d2 = fractional_distance(coords_frac(center_at,1:3), coords_frac(outer_j,1:3), lattice_vectors)
                        if (d2 > cutoff2(map_idx)) cycle
                        d12 = fractional_distance(coords_frac(outer_i,1:3), coords_frac(outer_j,1:3), lattice_vectors)
                        if (d12 > cutoff12(map_idx)) cycle
                        if (trim(outer1_symbols(map_idx)) == trim(outer2_symbols(map_idx)) .and. outer_j <= outer_i) cycle
                        if (angle_exists(shell_id(outer_i), core_id(center_at), shell_id(outer_j), type_ids(map_idx), &
                            angle_atom1, angle_atom2, angle_atom3, angle_type, angle_count)) cycle
                        call ensure_angle_capacity(angle_count + 1, capacity, angle_atom1, angle_atom2, angle_atom3, angle_type)
                        angle_count = angle_count + 1
                        angle_atom1(angle_count) = shell_id(outer_i)
                        angle_atom2(angle_count) = core_id(center_at)
                        angle_atom3(angle_count) = shell_id(outer_j)
                        angle_type(angle_count) = type_ids(map_idx)
                    end do
                end do
            end do
        end do

        if (angle_count == 0) then
            deallocate(angle_atom1)
            deallocate(angle_atom2)
            deallocate(angle_atom3)
            deallocate(angle_type)
            allocate(angle_atom1(1))
            allocate(angle_atom2(1))
            allocate(angle_atom3(1))
            allocate(angle_type(1))
            angle_atom1 = 0
            angle_atom2 = 0
            angle_atom3 = 0
            angle_type = 0
        end if
    end subroutine build_lammps_angles

    pure function fractional_distance(frac_a, frac_b, lattice_vectors) result(distance)
        real(dp), intent(in) :: frac_a(3), frac_b(3), lattice_vectors(3,3)
        real(dp) :: distance
        real(dp) :: delta_frac(3), delta_cart(3)

        delta_frac = frac_b - frac_a
        delta_frac = delta_frac - anint(delta_frac)
        delta_cart = matmul(transpose(lattice_vectors), delta_frac)
        distance = sqrt(dot_product(delta_cart, delta_cart))
    end function fractional_distance

    pure logical function angle_exists(atom1, atom2, atom3, atype, angle_atom1, angle_atom2, angle_atom3, angle_type, angle_count)
        integer, intent(in) :: atom1, atom2, atom3, atype, angle_count
        integer, intent(in) :: angle_atom1(:), angle_atom2(:), angle_atom3(:), angle_type(:)
        integer :: i

        angle_exists = .false.
        do i = 1, angle_count
            if (angle_type(i) /= atype) cycle
            if (angle_atom2(i) /= atom2) cycle
            if ((angle_atom1(i) == atom1 .and. angle_atom3(i) == atom3) .or. &
                (angle_atom1(i) == atom3 .and. angle_atom3(i) == atom1)) then
                angle_exists = .true.
                return
            end if
        end do
    end function angle_exists

    subroutine ensure_angle_capacity(required, capacity, angle_atom1, angle_atom2, angle_atom3, angle_type)
        integer, intent(in) :: required
        integer, intent(inout) :: capacity
        integer, allocatable, intent(inout) :: angle_atom1(:), angle_atom2(:), angle_atom3(:), angle_type(:)
        integer, allocatable :: tmp1(:), tmp2(:), tmp3(:), tmpt(:)

        if (required <= capacity) return
        do while (capacity < required)
            capacity = capacity * 2
        end do
        allocate(tmp1(capacity))
        allocate(tmp2(capacity))
        allocate(tmp3(capacity))
        allocate(tmpt(capacity))
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmpt = 0
        tmp1(1:size(angle_atom1)) = angle_atom1
        tmp2(1:size(angle_atom2)) = angle_atom2
        tmp3(1:size(angle_atom3)) = angle_atom3
        tmpt(1:size(angle_type)) = angle_type
        call move_alloc(tmp1, angle_atom1)
        call move_alloc(tmp2, angle_atom2)
        call move_alloc(tmp3, angle_atom3)
        call move_alloc(tmpt, angle_type)
    end subroutine ensure_angle_capacity

    pure function configuration_symbol(sp, atom_index, target_species, replacement_symbols, target_lookup, is_selected, &
        symbols) result(atom_symbol)
        integer, intent(in) :: sp, atom_index, target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        integer, intent(in) :: target_lookup(:)
        logical, intent(in) :: is_selected(:)
        character(len=*), intent(in) :: symbols(:)
        character(len=16) :: atom_symbol
        integer :: pos

        if (sp /= target_species) then
            atom_symbol = adjustl(symbols(sp))
            return
        end if

        pos = target_lookup(atom_index)
        if (pos >= 1 .and. pos <= size(is_selected)) then
            if (is_selected(pos)) then
                atom_symbol = adjustl(replacement_symbols(1))
            else
                atom_symbol = adjustl(replacement_symbols(2))
            end if
        else
            atom_symbol = adjustl(replacement_symbols(2))
        end if
    end function configuration_symbol

    pure function find_species_index(symbol, species_list) result(idx)
        character(len=*), intent(in) :: symbol
        character(len=*), intent(in) :: species_list(:)
        integer :: idx
        integer :: i

        idx = 0
        do i = 1, size(species_list)
            if (trim(symbol) == trim(species_list(i))) then
                idx = i
                return
            end if
        end do
    end function find_species_index

    subroutine build_output_species(symbols, target_species, replacement_symbols, all_symbols)
        character(len=*), intent(in) :: symbols(:)
        integer, intent(in) :: target_species
        character(len=*), intent(in) :: replacement_symbols(2)
        character(len=*), intent(out) :: all_symbols(:)
        integer :: sp

        all_symbols = ''
        do sp = 1, size(symbols)
            if (sp < target_species) then
                all_symbols(sp) = adjustl(symbols(sp))
            else if (sp == target_species) then
                all_symbols(sp) = adjustl(replacement_symbols(1))
                all_symbols(sp + 1) = adjustl(replacement_symbols(2))
            else
                all_symbols(sp + 1) = adjustl(symbols(sp))
            end if
        end do
    end subroutine build_output_species

    subroutine read_integer_after_key(line, key, value, found)
        character(len=*), intent(in) :: line
        character(len=*), intent(in) :: key
        integer, intent(out) :: value
        logical, intent(out), optional :: found
        integer :: pos, ios

        value = 0
        pos = index(line, key)
        if (present(found)) found = (pos > 0)
        if (pos <= 0) return
        read(line(pos + len_trim(key):),*,iostat=ios) value
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: invalid integer keyword in template line:', trim(line)
            stop 1
        end if
    end subroutine read_integer_after_key

    subroutine read_real_after_key(line, key, value, found)
        character(len=*), intent(in) :: line
        character(len=*), intent(in) :: key
        real(dp), intent(out) :: value
        logical, intent(out), optional :: found
        integer :: pos, ios

        value = 0.0_dp
        pos = index(line, key)
        if (present(found)) found = (pos > 0)
        if (pos <= 0) then
            write(error_unit,'(A,1X,A)') 'Error: missing real-valued keyword in template line:', trim(line)
            stop 1
        end if
        read(line(pos + len_trim(key):),*,iostat=ios) value
        if (ios /= 0) then
            write(error_unit,'(A,1X,A)') 'Error: invalid real-valued keyword in template line:', trim(line)
            stop 1
        end if
    end subroutine read_real_after_key

    pure function path_basename(path) result(name)
        character(len=*), intent(in) :: path
        character(len=len(path)) :: name
        integer :: pos

        pos = len_trim(path)
        do while (pos > 0)
            if (path(pos:pos) == '/') exit
            pos = pos - 1
        end do
        if (pos > 0) then
            name = trim(path(pos+1:))
        else
            name = trim(path)
        end if
    end function path_basename

    subroutine replace_all_tokens(buffer, token, replacement)
        character(len=*), intent(inout) :: buffer
        character(len=*), intent(in) :: token
        character(len=*), intent(in) :: replacement
        character(len=len(buffer)) :: updated
        integer :: pos, token_len, replacement_len, current_len

        token_len = len_trim(token)
        replacement_len = len_trim(replacement)
        if (token_len <= 0) return

        do
            pos = index(buffer, token(1:token_len))
            if (pos <= 0) exit
            updated = ''
            if (pos > 1) updated(1:pos-1) = buffer(1:pos-1)
            if (replacement_len > 0) then
                updated(pos:pos+replacement_len-1) = replacement(1:replacement_len)
            end if
            current_len = len_trim(buffer)
            if (pos + token_len <= current_len) then
                updated(pos+replacement_len:) = buffer(pos+token_len:)
            end if
            buffer = updated
        end do
    end subroutine replace_all_tokens

end module structure_io