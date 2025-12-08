!*******************************************************************************
!  Shared configuration constants and helper routines for ensemble drivers.
!*******************************************************************************

module sod_ensemble_config
    use sod_ensemble_consts
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

    ! Physical constants for energy comparison and reporting
    real(dp), parameter, public :: conv_ev_to_kjmol = 96.48533212331_dp
    real(dp), parameter, public :: quartz_ge_energy = -121.812200998519_dp
    real(dp), parameter, public :: quartz_si_energy = -128.799143746482_dp

    ! Default output file names
    character(len=*), parameter, public :: default_summary_filename = 'sod_ensemble_summary.csv'
    character(len=*), parameter, public :: default_summary_txt_filename = 'sod_ensemble_summary.txt'
    public :: mixing_weight, quartz_relative

contains

    pure real(dp) function mixing_weight(ge_fraction) result(weight)
        real(dp), intent(in) :: ge_fraction
        real(dp) :: arg, numerator, denom_base, denominator
        real(dp), parameter :: eps = 1.0e-8_dp

        arg = (ge_fraction - mix_d0) / mix_x0
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

    pure real(dp) function quartz_relative(level, total_sites, energy) result(rel_val)
        integer, intent(in) :: level, total_sites
        real(dp), intent(in) :: energy
        real(dp) :: ge_atoms, si_atoms, denom

        if (total_sites <= 0) then
            rel_val = 0.0_dp
            return
        end if

        ge_atoms = real(level, dp)
        si_atoms = real(total_sites - level, dp)
        denom = real(total_sites, dp)
        rel_val = conv_ev_to_kjmol * (energy - ge_atoms * quartz_ge_energy - si_atoms * quartz_si_energy) / denom
    end function quartz_relative

end module sod_ensemble_config
