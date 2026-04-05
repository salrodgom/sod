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
! Modern shared-module update of the classic `cell` helper.
! Module `cell_mod` converts lattice parameters into Cartesian cell vectors.
! El módulo `cell_mod` convierte los parámetros de red en vectores cartesianos de celda.
!
!*******************************************************************************
module cell_mod
    use iso_fortran_env, only: real32, real64
    implicit none
    private
    public :: cell

    interface cell
        module procedure cell_sp
        module procedure cell_dp
    end interface cell

contains

    pure subroutine cell_sp(cellvector, a, b, c, alpha, beta, gamma)
        real(real32), intent(out) :: cellvector(3,3)
        real(real32), intent(in) :: a, b, c, alpha, beta, gamma
        real(real64) :: work(3,3)

        call cell_dp(work, real(a, real64), real(b, real64), real(c, real64), &
            real(alpha, real64), real(beta, real64), real(gamma, real64))
        cellvector = real(work, real32)
    end subroutine cell_sp

    pure subroutine cell_dp(cellvector, a, b, c, alpha, beta, gamma)
        real(real64), intent(out) :: cellvector(3,3)
        real(real64), intent(in) :: a, b, c, alpha, beta, gamma
        real(real64), parameter :: degtorad = acos(-1.0_real64) / 180.0_real64
        real(real64) :: alp, bet, gam
        real(real64) :: cosa, cosb, cosg, sing, trm1

        if (alpha == 90.0_real64) then
            cosa = 0.0_real64
        else
            alp = alpha * degtorad
            cosa = cos(alp)
        end if

        if (beta == 90.0_real64) then
            cosb = 0.0_real64
        else
            bet = beta * degtorad
            cosb = cos(bet)
        end if

        if (gamma == 90.0_real64) then
            sing = 1.0_real64
            cosg = 0.0_real64
        else
            gam = gamma * degtorad
            sing = sin(gam)
            cosg = cos(gam)
        end if

        cellvector = 0.0_real64
        cellvector(1,1) = a
        cellvector(1,2) = b * cosg
        cellvector(2,2) = b * sing
        cellvector(1,3) = c * cosb
        cellvector(2,3) = c * (cosa - cosg * cosb) / sing
        trm1 = cellvector(2,3) / c
        cellvector(3,3) = c * sqrt(max(0.0_real64, 1.0_real64 - cosb**2 - trm1**2))
    end subroutine cell_dp

end module cell_mod
