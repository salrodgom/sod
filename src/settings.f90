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
! Shared configuration constants and helper routines for ensemble drivers.
!*******************************************************************************

! Module `settings` stores shared runtime settings and low/high energy blending rules.
! El módulo `settings` almacena ajustes compartidos de ejecución y reglas de mezcla de energías low/high.
module settings
    use consts
    implicit none
    private

    ! Mixing parameters shared by MC and exact drivers
    integer, parameter, public :: mix_n = 6
    integer, parameter, public :: mix_m = 12
    real(dp), parameter, public :: mix_x0 = 0.5_dp
    real(dp), parameter, public :: mix_d0 = 0.01_dp

    ! Monte Carlo sampling defaults
    integer, parameter, public :: default_max_exact_combos = 200000
    integer, parameter, public :: uniform_unique_cap = 1000
    integer, parameter, public :: uniform_unique_min_cap = 250
    real(dp), parameter, public :: uniform_cap_shrink = 0.75_dp

    ! Default reference constants for energy comparison and reporting
    real(dp), parameter, public :: conv_ev_to_kjmol = 96.48533212331_dp
    real(dp), parameter, public :: reference_y_energy = -121.812200998519_dp
    real(dp), parameter, public :: reference_x_energy = -128.799143746482_dp

    ! Default output file names
    character(len=*), parameter, public :: default_summary_filename = 'sod_ensemble_summary.csv'
    character(len=*), parameter, public :: default_summary_txt_filename = 'sod_ensemble_summary.txt'
    logical, allocatable, save :: blend_override_valid(:)
    real(dp), allocatable, save :: blend_override_high(:)
    integer, save :: blend_override_total_sites = -1

    public :: mixing_weight, mixing_weights, blend_low_high_energy, blend_low_high_energy_level, reference_relative
    public :: reset_blend_overrides, clear_blend_overrides, set_blend_override

contains

    pure real(dp) function mixing_weight(y_fraction) result(weight)
        real(dp), intent(in) :: y_fraction
        real(dp) :: arg, numerator, denom_base, denominator
        real(dp), parameter :: eps = 1.0e-8_dp

        arg = (y_fraction - mix_d0) / mix_x0
        denom_base = 1.0_dp - arg
        if (abs(denom_base) < eps) then
            denom_base = merge(eps, -eps, denom_base >= 0.0_dp)
        end if
        denominator = denom_base**mix_m
        numerator = 1.0_dp - arg**mix_n
        if (abs(denominator) < eps) then
            weight = merge(1.0_dp, 0.0_dp, numerator >= 0.0_dp)
        else
            weight = numerator / denominator
        end if
        weight = max(0.0_dp, min(1.0_dp, weight))
    end function mixing_weight

    pure subroutine mixing_weights(y_fraction, weight_low, weight_high)
        real(dp), intent(in) :: y_fraction
        real(dp), intent(out) :: weight_low, weight_high

        weight_low = mixing_weight(y_fraction)
        weight_low = max(0.0_dp, min(1.0_dp, weight_low))
        weight_high = 1.0_dp - weight_low
    end subroutine mixing_weights

    pure real(dp) function blend_low_high_energy(y_fraction, low_val, high_val) result(energy_val)
        real(dp), intent(in) :: y_fraction, low_val, high_val
        real(dp), parameter :: limit = huge(1.0_dp) * 0.5_dp
        real(dp) :: weight_low, weight_high
        logical :: has_low, has_high

        has_low = (abs(low_val) < limit) .and. (low_val == low_val)
        has_high = (abs(high_val) < limit) .and. (high_val == high_val)

        if (has_low .and. has_high) then
            call mixing_weights(y_fraction, weight_low, weight_high)
            energy_val = weight_low * low_val + weight_high * high_val
        else if (has_low) then
            energy_val = low_val
        else if (has_high) then
            energy_val = high_val
        else
            energy_val = limit
        end if
    end function blend_low_high_energy

    subroutine clear_blend_overrides()
        if (allocated(blend_override_valid)) deallocate(blend_override_valid)
        if (allocated(blend_override_high)) deallocate(blend_override_high)
        blend_override_total_sites = -1
    end subroutine clear_blend_overrides

    subroutine reset_blend_overrides(total_sites)
        integer, intent(in) :: total_sites

        call clear_blend_overrides()
        if (total_sites < 0) return
        allocate(blend_override_valid(0:total_sites))
        allocate(blend_override_high(0:total_sites))
        blend_override_valid = .false.
        blend_override_high = 0.0_dp
        blend_override_total_sites = total_sites
    end subroutine reset_blend_overrides

    subroutine set_blend_override(level, total_sites, weight_high)
        integer, intent(in) :: level, total_sites
        real(dp), intent(in) :: weight_high

        if (total_sites < 0) return
        if (.not. allocated(blend_override_valid) .or. blend_override_total_sites /= total_sites) then
            call reset_blend_overrides(total_sites)
        end if
        if (level < 0 .or. level > total_sites) return

        blend_override_valid(level) = .true.
        blend_override_high(level) = max(0.0_dp, min(1.0_dp, weight_high))
    end subroutine set_blend_override

    subroutine level_blend_weights(level, total_sites, weight_low, weight_high)
        integer, intent(in) :: level, total_sites
        real(dp), intent(out) :: weight_low, weight_high
        real(dp) :: y_fraction

        if (allocated(blend_override_valid)) then
            if (blend_override_total_sites == total_sites) then
                if (level >= lbound(blend_override_valid, 1) .and. level <= ubound(blend_override_valid, 1)) then
                    if (blend_override_valid(level)) then
                        weight_high = blend_override_high(level)
                        weight_low = 1.0_dp - weight_high
                        return
                    end if
                end if
            end if
        end if

        if (total_sites > 0) then
            y_fraction = real(level, dp) / real(total_sites, dp)
        else
            y_fraction = 0.0_dp
        end if
        call mixing_weights(y_fraction, weight_low, weight_high)
    end subroutine level_blend_weights

    real(dp) function blend_low_high_energy_level(level, total_sites, low_val, high_val) result(energy_val)
        integer, intent(in) :: level, total_sites
        real(dp), intent(in) :: low_val, high_val
        real(dp), parameter :: limit = huge(1.0_dp) * 0.5_dp
        real(dp) :: weight_low, weight_high
        logical :: has_low, has_high

        has_low = (abs(low_val) < limit) .and. (low_val == low_val)
        has_high = (abs(high_val) < limit) .and. (high_val == high_val)

        if (has_low .and. has_high) then
            call level_blend_weights(level, total_sites, weight_low, weight_high)
            energy_val = weight_low * low_val + weight_high * high_val
        else if (has_low) then
            energy_val = low_val
        else if (has_high) then
            energy_val = high_val
        else
            energy_val = limit
        end if
    end function blend_low_high_energy_level

    pure real(dp) function reference_relative(level, total_sites, energy) result(rel_val)
        integer, intent(in) :: level, total_sites
        real(dp), intent(in) :: energy
        real(dp) :: y_atoms, x_atoms, denom

        if (total_sites <= 0) then
            rel_val = 0.0_dp
            return
        end if

        y_atoms = real(level, dp)
        x_atoms = real(total_sites - level, dp)
        denom = real(total_sites, dp)
        rel_val = conv_ev_to_kjmol * (energy - y_atoms * reference_y_energy - x_atoms * reference_x_energy) / denom
    end function reference_relative

end module settings
