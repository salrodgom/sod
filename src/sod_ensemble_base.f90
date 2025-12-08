!*******************************************************************************
!  Shared ensemble Monte Carlo utilities: common constants and helper routines.
!*******************************************************************************

module sod_ensemble_consts
    use iso_fortran_env, only: real64, int64, error_unit
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: ip = int64
    real(dp), parameter :: kB_eVk = 8.617333262145d-5
end module sod_ensemble_consts

module sod_ensemble_utils
    use sod_ensemble_consts
    implicit none
contains

    ! Devuelve "n sobre k" usando un producto acumulativo en enteros de 64 bits.
    function binomial_int64(n, k) result(val)
        integer, intent(in) :: n, k
        integer(ip) :: val
        integer :: i

        if (k < 0 .or. k > n) then
            val = 0_ip
            return
        end if
        if (k == 0 .or. k == n) then
            val = 1_ip
            return
        end if
        val = 1_ip
        do i = 1, k
            val = val * int(n - k + i, ip) / int(i, ip)
        end do
    end function binomial_int64

    ! Avanza una combinación ordenada a la siguiente tupla lexicográfica.
    logical function next_combination(comb, n)
        integer, intent(inout) :: comb(:)
        integer, intent(in) :: n
        integer :: k, i, j

        k = size(comb)
        if (k == 0) then
            next_combination = .false.
            return
        end if

        i = k
        do while (i >= 1 .and. comb(i) == n - k + i)
            i = i - 1
        end do
        if (i == 0) then
            next_combination = .false.
            return
        end if

        comb(i) = comb(i) + 1
        do j = i + 1, k
            comb(j) = comb(j-1) + 1
        end do
        next_combination = .true.
    end function next_combination

    ! Extrae un subconjunto ordenado aleatorio de tamaño k dentro de 1..n sin reemplazo.
    subroutine random_subset(n, k, subset)
        integer, intent(in) :: n, k
        integer, intent(out) :: subset(:)
        integer :: pool(n)
        integer :: i, j, tmp
        real(dp) :: r

        if (k == 0) then
            return
        end if

        if (size(subset) < k) then
            write(error_unit,'(A)') 'ERROR interno: buffer insuficiente en random_subset'
            stop 1
        end if

        do i = 1, n
            pool(i) = i
        end do

        do i = 1, k
            call random_number(r)
            j = i + int(r * real(n - i + 1, dp))
            if (j < i) j = i
            if (j > n) j = n
            tmp = pool(i)
            pool(i) = pool(j)
            pool(j) = tmp
        end do

        subset(1:k) = pool(1:k)

        call sort_int_ascending(subset, k)
    end subroutine random_subset

    ! Ordena los primeros "length" elementos de un vector entero mediante inserción.
    subroutine sort_int_ascending(arr, length)
        integer, intent(inout) :: arr(:)
        integer, intent(in) :: length
        integer :: i, j, key

        do i = 2, length
            key = arr(i)
            j = i - 1
            do while (j >= 1 .and. arr(j) > key)
                arr(j+1) = arr(j)
                j = j - 1
            end do
            arr(j+1) = key
        end do
    end subroutine sort_int_ascending

end module sod_ensemble_utils
