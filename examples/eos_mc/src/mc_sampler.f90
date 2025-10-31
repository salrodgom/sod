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

MODULE mc_sampler
    USE energy_calc, ONLY: calculate_structure_energy, write_vasp_file
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: run_mc_sampling, init_mc_sampler, cleanup_mc_sampler

    ! Parameters
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    REAL(dp), PARAMETER :: kB = 8.61734E-5_dp ! Boltzmann constant in eV/K
    
    ! MC parameters
    INTEGER, PRIVATE :: n_steps
    REAL(dp), PRIVATE :: temperature
    INTEGER, ALLOCATABLE, PRIVATE :: current_config(:)
    REAL(dp), PRIVATE :: current_energy
    
    ! Acceptance tracking
    INTEGER, PRIVATE :: attempted_moves = 0
    INTEGER, PRIVATE :: accepted_moves = 0

    ! Best (lowest-energy) structure seen
    REAL(dp), PRIVATE :: best_energy
    INTEGER, ALLOCATABLE, PRIVATE :: best_config(:)

    CONTAINS

    SUBROUTINE init_mc_sampler(initial_config, n_sites, temp, steps)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n_sites, steps
        INTEGER, INTENT(IN) :: initial_config(n_sites)
        REAL(dp), INTENT(IN) :: temp
        
        n_steps = steps
        temperature = temp
        
        ALLOCATE(current_config(n_sites))
        current_config = initial_config
        
        ! Calculate initial energy
        CALL calculate_structure_energy(current_config, n_sites, current_energy)
        ! Initialize acceptance counters and best structure
        attempted_moves = 0
        accepted_moves = 0
        best_energy = current_energy
        IF (ALLOCATED(best_config)) DEALLOCATE(best_config)
        ALLOCATE(best_config(n_sites))
        best_config = current_config
    END SUBROUTINE init_mc_sampler

    SUBROUTINE run_mc_sampling(n_sites, properties)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n_sites
        REAL(dp), INTENT(OUT) :: properties(:) ! Array to store calculated properties
        
        INTEGER :: i, site1, site2, struct_count
        REAL(dp) :: trial_energy, delta_e, rand_num
        INTEGER, ALLOCATABLE :: trial_config(:)
        CHARACTER(len=20) :: filename
        LOGICAL :: accepted
        
        ALLOCATE(trial_config(n_sites))
        properties = 0.0_dp ! Initialize properties array
        
        struct_count = 0
        DO i = 1, n_steps
            ! Generate trial configuration by swapping two sites
            trial_config = current_config
            CALL random_number_pair(n_sites, site1, site2)
            
            ! Only swap if exactly one site is Si and one is Ge
            IF ((trial_config(site1) == 1 .AND. trial_config(site2) == 2) .OR. &
                (trial_config(site1) == 2 .AND. trial_config(site2) == 1)) THEN
                ! Count this attempted swap
                attempted_moves = attempted_moves + 1
                
                CALL swap_sites(trial_config, site1, site2)
                
                ! Calculate energy
                CALL calculate_structure_energy(trial_config, n_sites, trial_energy)
                
                ! Metropolis acceptance criterion
                delta_e = trial_energy - current_energy
                accepted = .FALSE.
                
                IF (delta_e < 0.0_dp) THEN
                    accepted = .TRUE.
                ELSE
                    CALL random_number(rand_num)
                    IF (rand_num < EXP(-delta_e/(kB*temperature))) THEN
                        accepted = .TRUE.
                    END IF
                END IF
                
                IF (accepted) THEN
                    ! Debug output
                    WRITE(*,*) 'Step', i, ': Move accepted'
                    WRITE(*,*) '  Sites swapped:', site1, site2
                    WRITE(*,*) '  Old energy:', current_energy
                    WRITE(*,*) '  New energy:', trial_energy
                    WRITE(*,*) '  Delta E:', delta_e
                    
                    CALL accept_move(trial_config, trial_energy)
                    ! Update best structure if improved
                    IF (current_energy < best_energy) THEN
                        best_energy = current_energy
                        best_config = current_config
                    END IF

                    ! Save structure periodically
                    IF (MOD(i, MAX(1,n_steps/100)) == 0) THEN
                        struct_count = struct_count + 1
                        WRITE(filename,'(A,I5.5,A)') 'mc_', struct_count, '.vasp'
                        CALL write_vasp_file(current_config, n_sites, filename)
                    END IF
                END IF
            END IF
            
            ! Store energy in properties array
            properties(i) = current_energy
        END DO
        
        DEALLOCATE(trial_config)
        
        ! Print acceptance statistics and write predicted best structure
        IF (attempted_moves > 0) THEN
            WRITE(*,*) 'Acceptance rate:', REAL(accepted_moves,dp)/REAL(attempted_moves,dp)
        ELSE
            WRITE(*,*) 'Acceptance rate: 0 (no attempted swaps)'
        END IF

        ! Write predicted best structure (VASP)
        WRITE(*,*) 'Writing predicted best structure to predicted_best.vasp'
        CALL write_vasp_file(best_config, n_sites, 'predicted_best.vasp')
    END SUBROUTINE run_mc_sampling

    SUBROUTINE cleanup_mc_sampler()
        IMPLICIT NONE
        IF (ALLOCATED(current_config)) DEALLOCATE(current_config)
        IF (ALLOCATED(best_config)) DEALLOCATE(best_config)
    END SUBROUTINE cleanup_mc_sampler

    ! Helper subroutines
    SUBROUTINE random_number_pair(n, i, j)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: i, j
        REAL(dp) :: r1, r2
        INTEGER :: attempts
        LOGICAL :: valid_pair
        
        attempts = 0
        valid_pair = .FALSE.
        
        DO WHILE (.NOT. valid_pair .AND. attempts < 100)
            CALL RANDOM_NUMBER(r1)
            CALL RANDOM_NUMBER(r2)
            i = INT(r1 * n) + 1
            j = INT(r2 * n) + 1
            
            ! Check if sites are different and not symmetry-equivalent
            IF (i /= j) THEN
                valid_pair = .TRUE.
                ! Add symmetry check here if needed
            END IF
            
            attempts = attempts + 1
        END DO
        
        IF (.NOT. valid_pair) THEN
            ! If we couldn't find a valid pair, just ensure they're different
            i = 1
            j = 2
        END IF
    END SUBROUTINE random_number_pair

    SUBROUTINE swap_sites(config, i, j)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: config(:)
        INTEGER, INTENT(IN) :: i, j
        INTEGER :: temp
        
        temp = config(i)
        config(i) = config(j)
        config(j) = temp
    END SUBROUTINE swap_sites

    SUBROUTINE accept_move(new_config, new_energy)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: new_config(:)
        REAL(dp), INTENT(IN) :: new_energy
        
        current_config = new_config
        current_energy = new_energy
        accepted_moves = accepted_moves + 1
    END SUBROUTINE accept_move

END MODULE mc_sampler