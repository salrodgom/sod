!*******************************************************************************
!    Copyright (c) 2025 
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!******************************************************************************

PROGRAM sod_mc
    USE energy_calc
    USE mc_sampler
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER, PARAMETER :: max_sites = 1000

    ! Variables
    INTEGER :: n_sites, n_steps, ge_count, ios
    REAL(dp) :: temperature, r
    INTEGER, ALLOCATABLE :: initial_config(:)
    REAL(dp), ALLOCATABLE :: properties(:)
    INTEGER :: i, j
    LOGICAL :: file_exists

    ! Check MC_INPUT existence
    INQUIRE(FILE='MC_INPUT', EXIST=file_exists)
    IF (.NOT. file_exists) THEN
        WRITE(*,*) 'Error: MC_INPUT file not found'
        WRITE(*,*) 'Expected format (one value per line):'
        WRITE(*,*) ' n_sites'
        WRITE(*,*) ' temperature (K)'
        WRITE(*,*) ' n_steps'
        WRITE(*,*) ' ge_count'
        STOP
    END IF

    ! Read input parameters
    OPEN(UNIT=10, FILE='MC_INPUT', STATUS='OLD', IOSTAT=ios)
    IF (ios /= 0) THEN
        WRITE(*,*) 'Error opening MC_INPUT'
        STOP
    END IF
    READ(10,*,IOSTAT=ios) n_sites
    IF (ios /= 0 .OR. n_sites <= 0 .OR. n_sites > max_sites) THEN
        WRITE(*,*) 'Error: Invalid number of sites. Must be 1..', max_sites
        STOP
    END IF
    READ(10,*,IOSTAT=ios) temperature
    IF (ios /= 0 .OR. temperature <= 0.0_dp) THEN
        WRITE(*,*) 'Error: Invalid temperature'
        STOP
    END IF
    READ(10,*,IOSTAT=ios) n_steps
    IF (ios /= 0 .OR. n_steps <= 0) THEN
        WRITE(*,*) 'Error: Invalid number of steps'
        STOP
    END IF
    READ(10,*,IOSTAT=ios) ge_count
    IF (ios /= 0 .OR. ge_count < 0 .OR. ge_count > n_sites) THEN
        WRITE(*,*) 'Error: Invalid ge_count (0..', n_sites, ')'
        STOP
    END IF
    CLOSE(10)

    ! Initialize configuration (1 = Si, 2 = Ge)
    ALLOCATE(initial_config(n_sites))
    initial_config = 1

    ! Place ge_count Ge atoms at random positions
    DO i = 1, ge_count
        DO
            CALL RANDOM_NUMBER(r)
            j = INT(r * n_sites) + 1
            IF (initial_config(j) == 1) THEN
                initial_config(j) = 2
                EXIT
            END IF
        END DO
    END DO

    ! Allocate properties (energy history)
    ALLOCATE(properties(n_steps))

    ! Initialize modules
    CALL init_energy_calc()
    CALL init_mc_sampler(initial_config, n_sites, temperature, n_steps)

    ! Run MC sampling
    CALL run_mc_sampling(n_sites, properties)

    ! Write results and cleanup
    CALL write_results(properties, n_steps)
    CALL cleanup_energy_calc()
    CALL cleanup_mc_sampler()

    DEALLOCATE(initial_config)
    DEALLOCATE(properties)

CONTAINS

    SUBROUTINE write_results(props, steps)
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: props(:)
        INTEGER, INTENT(IN) :: steps
        INTEGER :: i

        ! Write energy evolution
        OPEN(UNIT=20, FILE='energy.dat', STATUS='REPLACE')
        DO i = 1, steps
            WRITE(20,*) props(i)
        END DO
        CLOSE(20)

        ! Basic statistics
        OPEN(UNIT=21, FILE='mc_stats.dat', STATUS='REPLACE')
        WRITE(21,*) 'Average energy:', SUM(props)/steps
        WRITE(21,*) 'Minimum energy:', MINVAL(props)
        WRITE(21,*) 'Maximum energy:', MAXVAL(props)
        CLOSE(21)
    END SUBROUTINE write_results

END PROGRAM sod_mc