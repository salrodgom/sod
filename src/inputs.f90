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
! Modern shared-module update of the classic SGO/INSOD parsing logic.
!
!*******************************************************************************
! Module `inputs` reads and validates the SGO and INSOD input files.
! El módulo `inputs` lee y valida los ficheros de entrada SGO e INSOD.
module inputs
    use consts, only: dp, error_unit
    implicit none
    private

    type, public :: sgo_file_data
        integer :: nop1 = 0
        real(dp), allocatable :: mgroup1(:,:,:)
        real(dp), allocatable :: vgroup1(:,:)
    end type sgo_file_data

    type, public :: insod_file_data
        character(len=80) :: runtitle = ''
        real(dp) :: cell_params(6) = 0.0_dp
        integer :: nsp = 0
        integer :: nat0 = 0
        integer :: na = 0
        integer :: nb = 0
        integer :: nc = 0
        integer :: sptarget = 0
        integer :: nsubs = 0
        integer :: nsubs_min = 0
        integer :: nsubs_max = 0
        integer :: filer = 0
        integer :: mapper = 0
        logical :: has_shell_data = .false.
        character(len=3), allocatable :: symbol(:)
        character(len=3) :: newsymbol(2) = ''
        integer, allocatable :: natsp0(:)
        integer, allocatable :: spat0(:)
        integer, allocatable :: ishell(:)
        integer :: newshell(2) = 0
        real(dp), allocatable :: coords0(:,:)
    end type insod_file_data

    public :: read_sgo_file
    public :: read_insod_file

contains

    subroutine read_sgo_file(data)
        type(sgo_file_data), intent(out) :: data
        integer :: ios, i, j, op1, max_op
        character(len=512) :: line

        if (allocated(data%mgroup1)) deallocate(data%mgroup1)
        if (allocated(data%vgroup1)) deallocate(data%vgroup1)
        data%nop1 = 0

        open(unit=21, file='SGO', status='old', iostat=ios)
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: cannot open SGO.'
            stop 1
        end if

        max_op = 0
        read(21,'(A)',iostat=ios) line
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: unexpected end of SGO while reading the header.'
            stop 1
        end if
        read(line,*,iostat=ios) op1
        if (ios /= 0) then
            call next_data_line(21, line, 'SGO operator header')
            read(line,*,iostat=ios) op1
        end if
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: invalid SGO header.'
            stop 1
        end if

        do while (op1 > 0)
            max_op = max(max_op, op1)
            do i = 1, 3
                read(21,*,iostat=ios)
                if (ios /= 0) then
                    write(error_unit,'(A,I0)') 'Error: unexpected end of SGO while counting operator ', op1
                    stop 1
                end if
            end do
            call next_data_line(21, line, 'SGO operator header')
            read(line,*,iostat=ios) op1
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid SGO operator header.'
                stop 1
            end if
        end do

        if (max_op <= 0) then
            write(error_unit,'(A)') 'Error: SGO contains no symmetry operators.'
            stop 1
        end if

        rewind(21)
        read(21,'(A)',iostat=ios) line
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: unexpected end of SGO while restarting the reader.'
            stop 1
        end if
        read(line,*,iostat=ios) op1
        if (ios /= 0) then
            call next_data_line(21, line, 'SGO operator header')
            read(line,*,iostat=ios) op1
        end if
        if (ios /= 0) then
            write(error_unit,'(A)') 'Error: invalid SGO header on second pass.'
            stop 1
        end if

        allocate(data%mgroup1(max_op,3,3))
        allocate(data%vgroup1(max_op,3))
        data%mgroup1 = 0.0_dp
        data%vgroup1 = 0.0_dp

        do while (op1 > 0)
            do i = 1, 3
                read(21,*,iostat=ios) (data%mgroup1(op1,i,j), j=1,3), data%vgroup1(op1,i)
                if (ios /= 0) then
                    write(error_unit,'(A,I0)') 'Error: invalid SGO row for operator ', op1
                    stop 1
                end if
            end do
            data%nop1 = max(data%nop1, op1)
            call next_data_line(21, line, 'SGO operator header')
            read(line,*,iostat=ios) op1
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid SGO operator header.'
                stop 1
            end if
        end do
        close(21)
    end subroutine read_sgo_file

    subroutine read_insod_file(data)
        type(insod_file_data), intent(out) :: data
        integer :: i, ios, pair_values(2), supercell_values(3)
        character(len=512) :: line_buffer

        call reset_insod_data(data)

        open(unit=20, file='INSOD', status='old', action='read')

        call next_data_line(20, data%runtitle, 'INSOD runtitle')
        call next_real_line(20, data%cell_params, 'INSOD cell parameters')
        call next_int_scalar(20, data%nsp, 'INSOD nsp')

        if (data%nsp <= 0) then
            write(error_unit,'(A)') 'Error: INSOD declares an invalid number of species.'
            stop 1
        end if

        allocate(data%symbol(data%nsp))
        allocate(data%natsp0(data%nsp))
        call next_char_line(20, data%symbol, 'INSOD species symbols')
        call next_int_line(20, data%natsp0, 'INSOD species populations')

        data%nat0 = sum(data%natsp0)
        if (data%nat0 <= 0) then
            write(error_unit,'(A)') 'Error: INSOD contains no asymmetric-unit atoms.'
            stop 1
        end if

        allocate(data%coords0(data%nat0,3))
        do i = 1, data%nat0
            call next_real_line(20, data%coords0(i,1:3), 'INSOD coordinates')
        end do

        allocate(data%spat0(data%nat0))
        call build_spat0(data%natsp0, data%spat0)

        call next_int_line(20, supercell_values, 'INSOD supercell')
        data%na = supercell_values(1)
        data%nb = supercell_values(2)
        data%nc = supercell_values(3)

        call next_int_scalar(20, data%sptarget, 'INSOD target species')
        call next_data_line(20, line_buffer, 'INSOD substitution range')
        read(line_buffer,*,iostat=ios) data%nsubs_min, data%nsubs_max
        if (ios /= 0) then
            read(line_buffer,*,iostat=ios) data%nsubs_min
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid substitution count or range in INSOD.'
                stop 1
            end if
            data%nsubs_max = data%nsubs_min
        end if
        if (data%nsubs_min < 0) then
            write(error_unit,'(A)') 'Error: INSOD declares a negative minimum number of substitutions.'
            stop 1
        end if
        if (data%nsubs_max < data%nsubs_min) then
            write(error_unit,'(A)') 'Error: INSOD declares nsubs_max smaller than nsubs_min.'
            stop 1
        end if
        data%nsubs = data%nsubs_min

        call next_char_line(20, data%newsymbol, 'INSOD substitution symbols')
        call next_data_line(20, line_buffer, 'INSOD FILER or FILER/MAPPER')
        pair_values = 0
        read(line_buffer,*,iostat=ios) pair_values(1), pair_values(2)
        if (ios /= 0) then
            read(line_buffer,*,iostat=ios) pair_values(1)
            if (ios /= 0) then
                write(error_unit,'(A)') 'Error: invalid FILER or FILER/MAPPER specification in INSOD.'
                stop 1
            end if
            pair_values(2) = 0
        end if
        data%filer = pair_values(1)
        data%mapper = pair_values(2)

        close(20)
    end subroutine read_insod_file

    subroutine reset_insod_data(data)
        type(insod_file_data), intent(inout) :: data

        if (allocated(data%symbol)) deallocate(data%symbol)
        if (allocated(data%natsp0)) deallocate(data%natsp0)
        if (allocated(data%spat0)) deallocate(data%spat0)
        if (allocated(data%ishell)) deallocate(data%ishell)
        if (allocated(data%coords0)) deallocate(data%coords0)

        data%runtitle = ''
        data%cell_params = 0.0_dp
        data%nsp = 0
        data%nat0 = 0
        data%na = 0
        data%nb = 0
        data%nc = 0
        data%sptarget = 0
        data%nsubs = 0
        data%nsubs_min = 0
        data%nsubs_max = 0
        data%filer = 0
        data%mapper = 0
        data%has_shell_data = .false.
        data%newsymbol = ''
        data%newshell = 0
    end subroutine reset_insod_data

    subroutine build_spat0(natsp0, spat0)
        integer, intent(in) :: natsp0(:)
        integer, intent(out) :: spat0(:)
        integer :: sp, i, cumnatsp

        cumnatsp = 0
        do sp = 1, size(natsp0)
            do i = cumnatsp + 1, cumnatsp + natsp0(sp)
                spat0(i) = sp
            end do
            cumnatsp = cumnatsp + natsp0(sp)
        end do
    end subroutine build_spat0

    subroutine next_data_line(unit_id, line, label)
        integer, intent(in) :: unit_id
        character(len=*), intent(out) :: line
        character(len=*), intent(in) :: label
        integer :: ios

        do
            read(unit_id,'(A)',iostat=ios) line
            if (ios /= 0) then
                write(error_unit,'(A,1X,A)') 'Error: unexpected end of file while reading', trim(label)
                stop 1
            end if
            if (trim(line) == '') cycle
            if (line(1:1) == '#') cycle
            return
        end do
    end subroutine next_data_line

    subroutine next_real_line(unit_id, values, label)
        integer, intent(in) :: unit_id
        real(dp), intent(out) :: values(:)
        character(len=*), intent(in) :: label
        integer :: ios
        character(len=512) :: line

        do
            call next_data_line(unit_id, line, label)
            read(line,*,iostat=ios) values
            if (ios == 0) return
        end do
    end subroutine next_real_line

    subroutine next_int_line(unit_id, values, label)
        integer, intent(in) :: unit_id
        integer, intent(out) :: values(:)
        character(len=*), intent(in) :: label
        integer :: ios
        character(len=512) :: line

        do
            call next_data_line(unit_id, line, label)
            read(line,*,iostat=ios) values
            if (ios == 0) return
        end do
    end subroutine next_int_line

    subroutine next_int_scalar(unit_id, value, label)
        integer, intent(in) :: unit_id
        integer, intent(out) :: value
        character(len=*), intent(in) :: label
        integer :: ios
        character(len=512) :: line

        do
            call next_data_line(unit_id, line, label)
            read(line,*,iostat=ios) value
            if (ios == 0) return
        end do
    end subroutine next_int_scalar

    subroutine next_char_line(unit_id, values, label)
        integer, intent(in) :: unit_id
        character(len=*), intent(out) :: values(:)
        character(len=*), intent(in) :: label
        integer :: ios
        character(len=512) :: line

        do
            call next_data_line(unit_id, line, label)
            read(line,*,iostat=ios) values
            if (ios == 0) return
        end do
    end subroutine next_char_line

end module inputs
