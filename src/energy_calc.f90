!*******************************************************************************
!    Copyright (c) 2025, Salvador R.G. Balestra
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!******************************************************************************

MODULE energy_calc
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: calculate_structure_energy, init_energy_calc, cleanup_energy_calc, write_vasp_file, write_eqmatrix_file, get_eqmatrix, get_base_energy, get_high_base_energy, get_max_low_order, get_max_high_order

    ! Parameters from spbesod
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp
    !
    ! Energy parameters
    REAL(dp), PRIVATE :: E0              ! Reference energy
    REAL(dp), ALLOCATABLE, PRIVATE :: dE1(:)     ! One-body terms
    REAL(dp), ALLOCATABLE, PRIVATE :: dE2(:,:)   ! Two-body terms
    ! Three- and four-body terms stored sparsely as ordered tuples
    INTEGER, ALLOCATABLE, PRIVATE :: dE3_i(:), dE3_j(:), dE3_k(:)
    REAL(dp), ALLOCATABLE, PRIVATE :: dE3_val(:)
    INTEGER, PRIVATE :: n_dE3 = 0
    INTEGER, ALLOCATABLE, PRIVATE :: dE4_i(:), dE4_j(:), dE4_k(:), dE4_l(:)
    REAL(dp), ALLOCATABLE, PRIVATE :: dE4_val(:)
    INTEGER, PRIVATE :: n_dE4 = 0
    ! Complementary expansion (starting from all-Ge reference)
    REAL(dp), PRIVATE :: E0_high = 0.0_dp
    LOGICAL, PRIVATE :: high_base_loaded = .false.
    REAL(dp), ALLOCATABLE, PRIVATE :: hE1(:), hE2(:,:)
    INTEGER, ALLOCATABLE, PRIVATE :: hE3_i(:), hE3_j(:), hE3_k(:)
    REAL(dp), ALLOCATABLE, PRIVATE :: hE3_val(:)
    INTEGER, PRIVATE :: n_h_dE3 = 0
    INTEGER, ALLOCATABLE, PRIVATE :: hE4_i(:), hE4_j(:), hE4_k(:), hE4_l(:)
    REAL(dp), ALLOCATABLE, PRIVATE :: hE4_val(:)
    INTEGER, PRIVATE :: n_h_dE4 = 0
    ! Book-keeping for expansion orders available
    INTEGER, PRIVATE :: max_low_order = 0
    INTEGER, PRIVATE :: max_high_order = 0
    INTEGER, PRIVATE :: npos              ! Number of positions
    INTEGER, PRIVATE :: nop               ! Number of operators
    INTEGER, ALLOCATABLE, PRIVATE :: eqmatrix(:,:) ! Symmetry operations
    INTEGER, ALLOCATABLE, PRIVATE :: subpos(:) ! Mapping: MC site index -> position index (EQMATRIX column)
    
    ! Structure parameters
    REAL(dp), PRIVATE :: cell_params(6)  ! a, b, c, alpha, beta, gamma
    REAL(dp), ALLOCATABLE, PRIVATE :: coords0(:,:)  ! Asymmetric unit coordinates
    REAL(dp), ALLOCATABLE, PRIVATE :: coords1(:,:)  ! Unit cell coordinates after symmetry
    REAL(dp), ALLOCATABLE, PRIVATE :: coords(:,:)   ! Final structure coordinates
    INTEGER, PRIVATE :: n_ge, n_si, n_o   ! Number of atoms per type
    INTEGER, PRIVATE :: n_si_sites        ! Number of substitutable Si sites (from OUTSOD1)
    INTEGER, PRIVATE :: nat0              ! Total atoms in asymmetric unit
    INTEGER, PRIVATE :: nat1              ! Total atoms in unit cell
    INTEGER, PRIVATE :: natsp0(2)         ! Number of atoms per species in asymm unit
    CHARACTER(len=2), ALLOCATABLE, PRIVATE :: atom_types(:)  ! Atom symbols
    
    ! Space group symmetry 
    INTEGER, PRIVATE :: nop1              ! Number of symmetry operators
    REAL(dp), ALLOCATABLE, PRIVATE :: mgroup1(:,:,:) ! Rotation matrices
    REAL(dp), ALLOCATABLE, PRIVATE :: vgroup1(:,:)   ! Translation vectors
    INTEGER, ALLOCATABLE, PRIVATE :: spat0(:)        ! Species index per atom
    INTEGER, ALLOCATABLE, PRIVATE :: spat1(:)        ! Species index per atom in unit cell (deduped)
    INTEGER, ALLOCATABLE, PRIVATE :: pos2coord(:)    ! Mapping: position index (1:npos) -> coords index
    INTEGER, PRIVATE :: natsp1(2)                    ! Number of atoms per species in unit cell
    
    CONTAINS

    SUBROUTINE init_energy_calc()
        IMPLICIT NONE
    INTEGER :: op, m, m1, m2, aux, i, j, k, io_stat
    INTEGER :: Mm1, Mm2, Mm3, Mm4
    ! Variables for EQMATRIX generation
    INTEGER :: na, nb, nc, sptarget, t, at, pos, attmp, op1
    INTEGER, ALLOCATABLE :: fulleqmatrix(:,:)
    INTEGER, ALLOCATABLE :: coord2pos(:)
    REAL(dp), ALLOCATABLE :: mgroup(:,:,:), vgroup(:,:)
    REAL(dp), ALLOCATABLE :: vt(:,:)
    REAL(dp), ALLOCATABLE :: energies1(:), energies2(:)
    INTEGER, ALLOCATABLE :: conf1(:), conf2(:,:), omega1(:), omega2(:)
    INTEGER, ALLOCATABLE :: conf3(:,:), omega3(:)
    INTEGER, ALLOCATABLE :: conf4(:,:), omega4(:)
    REAL(dp), ALLOCATABLE :: energies3(:), energies4(:)
    INTEGER :: idx, found, ii, jj, kk, p, q
    INTEGER, ALLOCATABLE :: tmpi(:), tmpj(:), tmpk(:)
    REAL(dp), ALLOCATABLE :: tmpv(:)
    INTEGER, ALLOCATABLE :: tmp4_i(:), tmp4_j(:), tmp4_k(:), tmp4_l(:)
    REAL(dp), ALLOCATABLE :: tmp4_v(:)
    INTEGER :: arr4(4)
    REAL(dp) :: triple_sum
    CHARACTER(len=80) :: line
    INTEGER :: cnt1, cnt2
    ! Variables for unit-cell generation and deduplication
    INTEGER :: at0, at1i, nat1r, at1r
    REAL(dp) :: prod, tol0
    REAL(dp), ALLOCATABLE :: coords1r(:,:)
    INTEGER, ALLOCATABLE :: spat1r(:)
    REAL(dp) :: coordstemp(3)
        
    max_low_order = 0
    max_high_order = 0
    E0_high = 0.0_dp
    high_base_loaded = .false.
    n_h_dE3 = 0
    n_h_dE4 = 0
    IF (ALLOCATED(hE1)) DEALLOCATE(hE1)
    IF (ALLOCATED(hE2)) DEALLOCATE(hE2)
    IF (ALLOCATED(hE3_i)) DEALLOCATE(hE3_i)
    IF (ALLOCATED(hE3_j)) DEALLOCATE(hE3_j)
    IF (ALLOCATED(hE3_k)) DEALLOCATE(hE3_k)
    IF (ALLOCATED(hE3_val)) DEALLOCATE(hE3_val)
    IF (ALLOCATED(hE4_i)) DEALLOCATE(hE4_i)
    IF (ALLOCATED(hE4_j)) DEALLOCATE(hE4_j)
    IF (ALLOCATED(hE4_k)) DEALLOCATE(hE4_k)
    IF (ALLOCATED(hE4_l)) DEALLOCATE(hE4_l)
    IF (ALLOCATED(hE4_val)) DEALLOCATE(hE4_val)

    WRITE(*,*) 'init_energy_calc: start'
        ! Read SGO file first to get symmetry operations
        OPEN(UNIT=21, FILE='SGO', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open SGO'
            STOP
        END IF
        WRITE(*,*) 'init_energy_calc: opened SGO'
        
        ! Skip space group number and read first operator index
        READ(21,*)
        READ(21,*) i
        
        ! Count operators and store them (following combsod.f90 format)
        nop1 = 0
        ALLOCATE(mgroup1(48,3,3))  ! Max 48 operators
        ALLOCATE(vgroup1(48,3))
        
        DO WHILE (i > 0)
            nop1 = i
            DO j = 1, 3
                READ(21,*) mgroup1(i,j,1:3), vgroup1(i,j)
            END DO
            READ(21,*) i
        END DO
        CLOSE(21)
        
        WRITE(*,*) 'init_energy_calc: finished reading SGO; reading INSOD'
        ! Read INSOD for structure information (robust parsing)
        OPEN(UNIT=20, FILE='INSOD', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open INSOD'
            STOP
        END IF

        ! Find line with six real cell parameters
        DO
            READ(20,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: unexpected end of INSOD while searching for cell params'
                STOP
            END IF
            READ(line,*,IOSTAT=io_stat) cell_params
            IF (io_stat == 0) EXIT
        END DO

        ! Find next integer: number of species (nsp)
        DO
            READ(20,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: unexpected end of INSOD while searching for nsp'
                STOP
            END IF
            READ(line,*,IOSTAT=io_stat) i
            IF (io_stat == 0) EXIT
        END DO

        ! Read atom symbols line
        DO
            READ(20,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: unexpected end of INSOD while searching for atom symbols'
                STOP
            END IF
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            ALLOCATE(atom_types(2))
            READ(line,*,IOSTAT=io_stat) atom_types
            IF (io_stat == 0) EXIT
        END DO

        ! Read natsp0 line (numbers of atoms per species: initially all Si + O)
        DO
            READ(20,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: unexpected end of INSOD while searching for natsp0'
                STOP
            END IF
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            READ(line,*,IOSTAT=io_stat) natsp0(1), natsp0(2)
            IF (io_stat == 0) EXIT
        END DO
        
        ! Initialize counts - will be set properly when reading OUTSOD1
        n_si_sites = 0  ! Number of substitutable Si sites
        n_ge = 0        ! Initially no Ge
        n_si = natsp0(1) ! Initially all Si in asymmetric unit

    ! natsp0 already set from INSOD line: natsp0(1)=Si, natsp0(2)=O
        
        ! Calculate total atoms in asymmetric unit
        nat0 = natsp0(1) + natsp0(2)
        
        ! Allocate and read asymmetric unit coordinates
        ALLOCATE(coords0(nat0,3))
        ALLOCATE(spat0(nat0))
        
        ! Read coordinates and set species indices
        i = 0
        DO WHILE (i < nat0)
            READ(20,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: unexpected end of INSOD while reading coordinates'
                STOP
            END IF
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            i = i + 1
            READ(line,*,IOSTAT=io_stat) coords0(i,1:3)
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: bad coordinate line in INSOD at index', i
                STOP
            END IF
            ! Set species index (1=Si for i≤natsp0(1), 2=O for i>natsp0(1))
            spat0(i) = MERGE(1, 2, i <= natsp0(1))
        END DO
        CLOSE(20)
    WRITE(*,*) 'init_energy_calc: finished INSOD parsing; nat0=', nat0, 'nop1=', nop1
        
    ! Generate unit cell by applying symmetry operations (redundant list)
    ! Then deduplicate to obtain the same unit-cell coords as combsod
    tol0 = 0.001_dp
        ! Create redundant unit-cell list coords1r
        nat1r = nat0 * nop1
        ALLOCATE(coords1r(nat1r,3))
        ALLOCATE(spat1r(nat1r))
        at1r = 0
        DO at0 = 1, nat0
            DO j = 1, nop1
                at1r = at1r + 1
                coords1r(at1r,1:3) = MATMUL(mgroup1(j,1:3,1:3), coords0(at0,1:3)) + vgroup1(j,1:3)
                coords1r(at1r,1:3) = coords1r(at1r,1:3) - FLOOR(coords1r(at1r,1:3))
                spat1r(at1r) = spat0(at0)
            END DO
        END DO

        ! Deduplicate coords1r -> coords1
        IF (ALLOCATED(coords1)) DEALLOCATE(coords1)
        ALLOCATE(coords1(nat1r,3))
        IF (ALLOCATED(spat1)) DEALLOCATE(spat1)
        ALLOCATE(spat1(nat1r))

        coords1(1,1:3) = coords1r(1,1:3)
        spat1(1) = spat1r(1)
        at1r = 1
        DO at0 = 2, nat1r
            prod = 0.0_dp
            DO at1i = 1, at1r
                prod = DOT_PRODUCT(coords1r(at0,1:3)-coords1(at1i,1:3), coords1r(at0,1:3)-coords1(at1i,1:3))
                IF (prod .LE. tol0) EXIT
            END DO
            IF (prod .GT. tol0) THEN
                at1r = at1r + 1
                coords1(at1r,1:3) = coords1r(at0,1:3)
                spat1(at1r) = spat1r(at0)
            END IF
        END DO
        nat1 = at1r

        ! Compute species counts in unit cell
        natsp1(1) = 0
        natsp1(2) = 0
        DO at1r = 1, nat1
            IF (spat1(at1r) == 1) THEN
                natsp1(1) = natsp1(1) + 1
            ELSE
                natsp1(2) = natsp1(2) + 1
            END IF
        END DO

        ! Set n_si and n_o from deduplicated unit cell (use species 1 as substitutable)
        n_si = natsp1(1)
        n_o = natsp1(2)

        ! Prepare final coords array (copy deduped coords)
        IF (ALLOCATED(coords)) DEALLOCATE(coords)
        ALLOCATE(coords(nat1,3))
        coords(1:nat1,1:3) = coords1(1:nat1,1:3)

        ! Build pos2coord mapping: positions correspond to species==1 entries in coords1
        IF (ALLOCATED(pos2coord)) DEALLOCATE(pos2coord)
        ALLOCATE(pos2coord(n_si))
        at1r = 0
        DO at0 = 1, nat1
            IF (spat1(at0) == 1) THEN
                at1r = at1r + 1
                pos2coord(at1r) = at0
            END IF
        END DO
        ! Free redundant arrays
        IF (ALLOCATED(coords1r)) DEALLOCATE(coords1r)
        IF (ALLOCATED(spat1r)) DEALLOCATE(spat1r)
        
        ! Try to read EQMATRIX from file; if absent, generate from SGO+INSOD
        WRITE(*,*) 'init_energy_calc: attempting to read EQMATRIX file'
        OPEN(UNIT=13, FILE='EQMATRIX', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat == 0) THEN
            READ(13,*) nop, npos
            WRITE(*,*) 'init_energy_calc: read EQMATRIX nop=', nop, 'npos=', npos
            ALLOCATE(eqmatrix(nop,npos))
            DO op = 1, nop
                READ(13,*) eqmatrix(op,1:npos)
            END DO
            CLOSE(13)
        ELSE
            WRITE(*,*) 'init_energy_calc: EQMATRIX not found — generating from SGO+INSOD'
            ! Re-open INSOD and parse na,nb,nc and sptarget
            OPEN(UNIT=22, FILE='INSOD', STATUS='OLD', IOSTAT=io_stat)
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Error: cannot open INSOD to read supercell params'
                STOP
            END IF
            ! Scan for na, nb, nc (first line with three integers)
            na = 1; nb = 1; nc = 1
            DO
                READ(22,'(A)',IOSTAT=io_stat) line
                IF (io_stat /= 0) THEN
                    EXIT
                END IF
                IF (TRIM(line) == '') CYCLE
                IF (line(1:1) == '#') CYCLE
                READ(line,*,IOSTAT=io_stat) na, nb, nc
                IF (io_stat == 0) EXIT
            END DO
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Warning: could not find na nb nc in INSOD, defaulting to 1 1 1'
                na = 1; nb = 1; nc = 1
            END IF
            ! Find sptarget (next integer in file)
            DO
                READ(22,'(A)',IOSTAT=io_stat) line
                IF (io_stat /= 0) THEN
                    EXIT
                END IF
                IF (TRIM(line) == '') CYCLE
                IF (line(1:1) == '#') CYCLE
                READ(line,*,IOSTAT=io_stat) sptarget
                IF (io_stat == 0) EXIT
            END DO
            CLOSE(22)

            ! Build translation vectors vt (na*nb*nc)
            ALLOCATE(vt(na*nb*nc,3))
            t = 0
            DO i = 0, na-1
                DO j = 0, nb-1
                    DO op1 = 0, nc-1
                        t = t + 1
                        vt(t,1) = REAL(i,dp) / REAL(na,dp)
                        vt(t,2) = REAL(j,dp) / REAL(nb,dp)
                        vt(t,3) = REAL(op1,dp) / REAL(nc,dp)
                    END DO
                END DO
            END DO

            ! Build supercell operators mgroup and vgroup
            nop = nop1 * na * nb * nc
            ALLOCATE(mgroup(nop,3,3))
            ALLOCATE(vgroup(nop,3))
            op = 0
            DO op1 = 1, nop1
                DO t = 1, na*nb*nc
                    op = op + 1
                    mgroup(op,:,: ) = mgroup1(op1,:,:)
                    vgroup(op,1) = vgroup1(op1,1)/REAL(na,dp) + vt(t,1)
                    vgroup(op,2) = vgroup1(op1,2)/REAL(nb,dp) + vt(t,2)
                    vgroup(op,3) = vgroup1(op1,3)/REAL(nc,dp) + vt(t,3)
                END DO
            END DO

            ! Allocate fulleqmatrix and build mapping for each operator and each deduped atom
            ALLOCATE(fulleqmatrix(nop,nat1))
            ALLOCATE(coord2pos(nat1))
            coord2pos = 0
            ! Build coord2pos inverse mapping for target species positions
            DO pos = 1, SIZE(pos2coord)
                coord2pos(pos2coord(pos)) = pos
            END DO

            DO op = 1, nop
                DO at = 1, nat1
                    ! Apply operator to coords(at)
                    coordstemp = MATMUL(mgroup(op,:,:), coords(at,1:3)) + vgroup(op,1:3)
                    coordstemp = coordstemp - FLOOR(coordstemp)
                    ! Find matching index in coords (deduped)
                    attmp = 0
                    DO i = 1, nat1
                        prod = DOT_PRODUCT(coordstemp - coords(i,1:3), coordstemp - coords(i,1:3))
                        IF (prod .LE. tol0) THEN
                            attmp = i
                            EXIT
                        END IF
                    END DO
                    IF (attmp == 0) THEN
                        WRITE(*,*) 'Error: operator mapping failed for op=', op, 'at=', at
                        STOP
                    END IF
                    fulleqmatrix(op,at) = attmp
                END DO
            END DO

            ! Extract eqmatrix for target positions (npos = number of positions = SIZE(pos2coord))
            npos = SIZE(pos2coord)
            ALLOCATE(eqmatrix(nop,npos))
            DO op = 1, nop
                DO pos = 1, npos
                    eqmatrix(op,pos) = coord2pos(fulleqmatrix(op,pos2coord(pos)))
                END DO
            END DO

            ! Free temporary arrays
            IF (ALLOCATED(mgroup)) DEALLOCATE(mgroup)
            IF (ALLOCATED(vgroup)) DEALLOCATE(vgroup)
            IF (ALLOCATED(vt)) DEALLOCATE(vt)
            IF (ALLOCATED(fulleqmatrix)) DEALLOCATE(fulleqmatrix)
            IF (ALLOCATED(coord2pos)) DEALLOCATE(coord2pos)
        END IF
        
        ! Allocate energy terms
        ALLOCATE(dE1(npos))
        ALLOCATE(dE2(npos,npos))
        
        WRITE(*,*) 'init_energy_calc: reading n00/ENERGIES'
        ! Read reference energy from n00/ENERGIES
        OPEN(UNIT=10, FILE='n00/ENERGIES', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open n00/ENERGIES'
            STOP
        END IF
        READ(10,*) E0
        CLOSE(10)
        WRITE(*,*) 'init_energy_calc: read E0=', E0
        
        ! Read one-body terms
    WRITE(*,*) 'init_energy_calc: reading n01/OUTSOD'
    ! Read single-substitution OUTSOD from n01/OUTSOD
        OPEN(UNIT=15, FILE='n01/OUTSOD', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open n01/OUTSOD'
            STOP
        END IF
        ! Parse file: find the line containing the number of configurations
        Mm1 = 0
        DO
            READ(15,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                EXIT
            END IF
            IF (INDEX(line,'configuration') /= 0) THEN
                ! line like: "   4  configurations"
                READ(line,*,IOSTAT=io_stat) Mm1
                IF (io_stat == 0) EXIT
            END IF
        END DO
        IF (Mm1 <= 0) THEN
            WRITE(*,*) 'Error: could not determine Mm1 from n01/OUTSOD'
            STOP
        END IF
        WRITE(*,*) 'init_energy_calc: Mm1=', Mm1
        ALLOCATE(conf1(Mm1), omega1(Mm1), energies1(Mm1))
        ! Now rewind and read numeric mapping lines (skip non-numeric lines)
        REWIND(15)
        m = 0
        DO
            READ(15,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) EXIT
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            ! Try to read three integers from the line: aux, omega, conf
            READ(line,*,IOSTAT=io_stat) aux, omega1(m+1), conf1(m+1)
            IF (io_stat == 0) THEN
                m = m + 1
                IF (m == Mm1) EXIT
            END IF
        END DO
    WRITE(*,*) 'init_energy_calc: finished reading OUTSOD1 mappings, read m=', m
    CLOSE(15)

    ! Store substitution site -> position mapping for MC using the full set of Si positions
    IF (ALLOCATED(subpos)) DEALLOCATE(subpos)
    ALLOCATE(subpos(npos))
    subpos = [(i, i=1, npos)]
        
    WRITE(*,*) 'init_energy_calc: reading n01/ENERGIES'
    ! Read ENERGIES1 from n01/ENERGIES
        OPEN(UNIT=11, FILE='n01/ENERGIES', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open n01/ENERGIES'
            STOP
        END IF
        DO m = 1, Mm1
            READ(11,*) energies1(m)
        END DO
    WRITE(*,*) 'init_energy_calc: finished reading ENERGIES1'
    CLOSE(11)
    max_low_order = MAX(max_low_order, 1)
        
        ! Calculate dE1 - single substitution energies
        dE1 = 0.0_dp
        DO m1 = 1, Mm1
            DO op = 1, nop
                dE1(eqmatrix(op,conf1(m1))) = dE1(eqmatrix(op,conf1(m1))) + &
                    (energies1(m1) - E0) / COUNT(eqmatrix(:,conf1(m1)) == eqmatrix(op,conf1(m1)))
            END DO
        END DO
        
        ! Read two-body terms
    WRITE(*,*) 'init_energy_calc: reading n02/OUTSOD'
    ! Read pair-substitution OUTSOD from n02/OUTSOD
        OPEN(UNIT=16, FILE='n02/OUTSOD', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open n02/OUTSOD'
            STOP
        END IF
        Mm2 = 0
        DO
            READ(16,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                EXIT
            END IF
            IF (INDEX(line,'configuration') /= 0) THEN
                READ(line,*,IOSTAT=io_stat) Mm2
                IF (io_stat == 0) EXIT
            END IF
        END DO
        IF (Mm2 < 0) THEN
            WRITE(*,*) 'Error: could not determine Mm2 from n02/OUTSOD'
            STOP
        END IF
        ALLOCATE(conf2(Mm2,2), omega2(Mm2), energies2(Mm2))
        REWIND(16)
        m = 0
        DO
            READ(16,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) EXIT
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            READ(line,*,IOSTAT=io_stat) aux, omega2(m+1), conf2(m+1,1), conf2(m+1,2)
            IF (io_stat == 0) THEN
                m = m + 1
                IF (m == Mm2) EXIT
            END IF
        END DO
    WRITE(*,*) 'init_energy_calc: finished reading OUTSOD2 mappings, read m=', m
    CLOSE(16)
        
    WRITE(*,*) 'init_energy_calc: reading n02/ENERGIES'
    ! Read ENERGIES2 from n02/ENERGIES
        OPEN(UNIT=12, FILE='n02/ENERGIES', STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'Error: cannot open n02/ENERGIES'
            STOP
        END IF
        DO m = 1, Mm2
            READ(12,*) energies2(m)
        END DO
    WRITE(*,*) 'init_energy_calc: finished reading ENERGIES2'
    CLOSE(12)
    max_low_order = MAX(max_low_order, 2)
        
        ! Calculate dE2 - pair interaction energies
        dE2 = 0.0_dp
        DO m2 = 1, Mm2
            DO op = 1, nop
                IF (eqmatrix(op,conf2(m2,1)) < eqmatrix(op,conf2(m2,2))) THEN
                    dE2(eqmatrix(op,conf2(m2,1)), eqmatrix(op,conf2(m2,2))) = &
                        dE2(eqmatrix(op,conf2(m2,1)), eqmatrix(op,conf2(m2,2))) + &
                        (energies2(m2) - E0 - dE1(eqmatrix(op,conf2(m2,1))) - &
                         dE1(eqmatrix(op,conf2(m2,2)))) / &
                        COUNT(eqmatrix(:,conf2(m2,1)) == eqmatrix(op,conf2(m2,1)) .AND. &
                             eqmatrix(:,conf2(m2,2)) == eqmatrix(op,conf2(m2,2)))
                    
                    ! Apply symmetry for (j,i) pairs
                    dE2(eqmatrix(op,conf2(m2,2)), eqmatrix(op,conf2(m2,1))) = &
                        dE2(eqmatrix(op,conf2(m2,1)), eqmatrix(op,conf2(m2,2)))
                END IF
            END DO
        END DO
        
    ! --- Read three-body terms (n03) ------------------------------------
    WRITE(*,*) 'init_energy_calc: reading n03/OUTSOD'
    OPEN(UNIT=17, FILE='n03/OUTSOD', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'Warning: n03/OUTSOD not found; skipping three-body terms'
        Mm3 = 0
    ELSE
        Mm3 = 0
        DO
            READ(17,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                EXIT
            END IF
            IF (INDEX(line,'configuration') /= 0) THEN
                READ(line,*,IOSTAT=io_stat) Mm3
                IF (io_stat == 0) EXIT
            END IF
        END DO
        IF (Mm3 <= 0) THEN
            WRITE(*,*) 'Warning: could not determine Mm3 from n03/OUTSOD; skipping three-body'
            CLOSE(17)
            Mm3 = 0
        ELSE
            WRITE(*,*) 'init_energy_calc: Mm3=', Mm3
            ALLOCATE(conf3(Mm3,3), omega3(Mm3), energies3(Mm3))
            REWIND(17)
            m = 0
            DO
                READ(17,'(A)',IOSTAT=io_stat) line
                IF (io_stat /= 0) EXIT
                IF (TRIM(line) == '') CYCLE
                IF (line(1:1) == '#') CYCLE
                READ(line,*,IOSTAT=io_stat) aux, omega3(m+1), conf3(m+1,1), conf3(m+1,2), conf3(m+1,3)
                IF (io_stat == 0) THEN
                    m = m + 1
                    IF (m == Mm3) EXIT
                END IF
            END DO
            WRITE(*,*) 'init_energy_calc: finished reading OUTSOD3 mappings, read m=', m
            CLOSE(17)

            ! Read ENERGIES3
            WRITE(*,*) 'init_energy_calc: reading n03/ENERGIES'
            OPEN(UNIT=18, FILE='n03/ENERGIES', STATUS='OLD', IOSTAT=io_stat)
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Warning: cannot open n03/ENERGIES; skipping three-body terms'
                DEALLOCATE(conf3, omega3)
                Mm3 = 0
            ELSE
                DO m = 1, Mm3
                    READ(18,*) energies3(m)
                END DO
                CLOSE(18)
                max_low_order = MAX(max_low_order, 3)

                ! Build sparse dE3: lists of ordered triples (i<j<k) and values
                IF (ALLOCATED(dE3_i)) THEN
                    DEALLOCATE(dE3_i, dE3_j, dE3_k, dE3_val)
                END IF
                n_dE3 = 0
                DO m = 1, Mm3
                    DO op = 1, nop
                        ! Map the three representative positions through this op
                        i = eqmatrix(op, conf3(m,1))
                        j = eqmatrix(op, conf3(m,2))
                        k = eqmatrix(op, conf3(m,3))
                        ! Order indices ascending
                        IF (i == j .OR. i == k .OR. j == k) CYCLE
                        IF (i < j .AND. j < k) THEN
                            idx = 1
                        ELSEIF (i < k .AND. k < j) THEN
                            idx = 2
                        ELSEIF (j < i .AND. i < k) THEN
                            idx = 3
                        ELSEIF (j < k .AND. k < i) THEN
                            idx = 4
                        ELSEIF (k < i .AND. i < j) THEN
                            idx = 5
                        ELSE
                            idx = 6
                        END IF
                        SELECT CASE (idx)
                        CASE (1)
                            aux = i; i = i; j = j; k = k
                        CASE (2)
                            aux = i; i = i; j = k; k = j
                        CASE (3)
                            aux = i; i = j; j = i; k = k
                        CASE (4)
                            aux = i; i = j; j = k; k = aux
                        CASE (5)
                            aux = i; i = k; j = i; k = j
                        CASE DEFAULT
                            aux = i; i = k; j = j; k = i
                        END SELECT
                        ! At this point we must ensure i<j<k; compute canonical triple
                        ! Simpler: compute ii,jj,kk as ascending values
                        IF (i < j) THEN
                            IF (j < k) THEN
                                idx = 0
                            END IF
                        END IF
                        ! Get ordered triple
                        ii = MIN(i, MIN(j,k))
                        kk = MAX(i, MAX(j,k))
                        jj = i + j + k - ii - kk

                        ! Compute denominator: number of operators mapping to same ordered triple
                        aux = COUNT(eqmatrix(:,conf3(m,1)) == eqmatrix(op,conf3(m,1)) .AND. &
                                    eqmatrix(:,conf3(m,2)) == eqmatrix(op,conf3(m,2)) .AND. &
                                    eqmatrix(:,conf3(m,3)) == eqmatrix(op,conf3(m,3)))
                        IF (aux == 0) CYCLE

                        ! Subtract reference, one-body and pair contributions
                        prod = energies3(m) - E0 - dE1(ii) - dE1(jj) - dE1(kk) - &
                               (dE2(MIN(ii,jj),MAX(ii,jj)) + dE2(MIN(ii,kk),MAX(ii,kk)) + dE2(MIN(jj,kk),MAX(jj,kk)))
                        prod = prod / REAL(aux,dp)

                        ! Append or accumulate into sparse arrays
                        found = 0
                        IF (n_dE3 > 0 .AND. ALLOCATED(dE3_i)) THEN
                            DO idx = 1, n_dE3
                                IF (dE3_i(idx) == ii .AND. dE3_j(idx) == jj .AND. dE3_k(idx) == kk) THEN
                                    dE3_val(idx) = dE3_val(idx) + prod
                                    found = 1
                                    EXIT
                                END IF
                            END DO
                        END IF
                        IF (found == 0) THEN
                            ! append new entry
                            IF (.NOT. ALLOCATED(dE3_i)) THEN
                                ALLOCATE(dE3_i(1), dE3_j(1), dE3_k(1))
                                ALLOCATE(dE3_val(1))
                                n_dE3 = 1
                                dE3_i(1) = ii
                                dE3_j(1) = jj
                                dE3_k(1) = kk
                                dE3_val(1) = prod
                            ELSE
                                ALLOCATE(tmpi(n_dE3+1), tmpj(n_dE3+1), tmpk(n_dE3+1), tmpv(n_dE3+1))
                                tmpi(1:n_dE3) = dE3_i(1:n_dE3)
                                tmpj(1:n_dE3) = dE3_j(1:n_dE3)
                                tmpk(1:n_dE3) = dE3_k(1:n_dE3)
                                tmpv(1:n_dE3) = dE3_val(1:n_dE3)
                                tmpi(n_dE3+1) = ii
                                tmpj(n_dE3+1) = jj
                                tmpk(n_dE3+1) = kk
                                tmpv(n_dE3+1) = prod
                                DEALLOCATE(dE3_i, dE3_j, dE3_k, dE3_val)
                                ALLOCATE(dE3_i(n_dE3+1), dE3_j(n_dE3+1), dE3_k(n_dE3+1))
                                ALLOCATE(dE3_val(n_dE3+1))
                                dE3_i = tmpi
                                dE3_j = tmpj
                                dE3_k = tmpk
                                dE3_val = tmpv
                                DEALLOCATE(tmpi, tmpj, tmpk, tmpv)
                                n_dE3 = n_dE3 + 1
                            END IF
                        END IF
                    END DO
                END DO
            END IF
        END IF
    END IF

    ! --- Read four-body terms (n04) ------------------------------------
    IF (ALLOCATED(dE4_i)) DEALLOCATE(dE4_i)
    IF (ALLOCATED(dE4_j)) DEALLOCATE(dE4_j)
    IF (ALLOCATED(dE4_k)) DEALLOCATE(dE4_k)
    IF (ALLOCATED(dE4_l)) DEALLOCATE(dE4_l)
    IF (ALLOCATED(dE4_val)) DEALLOCATE(dE4_val)
    n_dE4 = 0

    WRITE(*,*) 'init_energy_calc: reading n04/OUTSOD'
    OPEN(UNIT=19, FILE='n04/OUTSOD', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'Warning: n04/OUTSOD not found; skipping four-body terms'
        Mm4 = 0
    ELSE
        Mm4 = 0
        DO
            READ(19,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) EXIT
            IF (INDEX(line,'configuration') /= 0) THEN
                READ(line,*,IOSTAT=io_stat) Mm4
                IF (io_stat == 0) EXIT
            END IF
        END DO
        IF (Mm4 <= 0) THEN
            WRITE(*,*) 'Warning: could not determine Mm4 from n04/OUTSOD; skipping four-body'
            CLOSE(19)
            Mm4 = 0
        ELSE
            WRITE(*,*) 'init_energy_calc: Mm4=', Mm4
            ALLOCATE(conf4(Mm4,4), omega4(Mm4), energies4(Mm4))
            REWIND(19)
            m = 0
            DO
                READ(19,'(A)',IOSTAT=io_stat) line
                IF (io_stat /= 0) EXIT
                IF (TRIM(line) == '') CYCLE
                IF (line(1:1) == '#') CYCLE
                READ(line,*,IOSTAT=io_stat) aux, omega4(m+1), conf4(m+1,1), conf4(m+1,2), conf4(m+1,3), conf4(m+1,4)
                IF (io_stat == 0) THEN
                    m = m + 1
                    IF (m == Mm4) EXIT
                END IF
            END DO
            WRITE(*,*) 'init_energy_calc: finished reading OUTSOD4 mappings, read m=', m
            CLOSE(19)

            WRITE(*,*) 'init_energy_calc: reading n04/ENERGIES'
            OPEN(UNIT=28, FILE='n04/ENERGIES', STATUS='OLD', IOSTAT=io_stat)
            IF (io_stat /= 0) THEN
                WRITE(*,*) 'Warning: cannot open n04/ENERGIES; skipping four-body terms'
                DEALLOCATE(conf4, omega4, energies4)
                Mm4 = 0
            ELSE
                DO m = 1, Mm4
                    READ(28,*,IOSTAT=io_stat) energies4(m)
                    IF (io_stat /= 0) EXIT
                END DO
                CLOSE(28)
                IF (io_stat /= 0) THEN
                    WRITE(*,*) 'Warning: incomplete n04/ENERGIES data; skipping four-body terms'
                    DEALLOCATE(conf4, omega4, energies4)
                    Mm4 = 0
                ELSE
                    n_dE4 = 0
                    max_low_order = MAX(max_low_order, 4)
                    DO m = 1, Mm4
                        quad_operator: DO op = 1, nop
                            arr4 = (/ eqmatrix(op, conf4(m,1)), eqmatrix(op, conf4(m,2)), &
                                      eqmatrix(op, conf4(m,3)), eqmatrix(op, conf4(m,4)) /)
                            DO p = 1, 3
                                DO q = p+1, 4
                                    IF (arr4(p) == arr4(q)) CYCLE quad_operator
                                END DO
                            END DO

                            aux = COUNT(eqmatrix(:, conf4(m,1)) == eqmatrix(op, conf4(m,1)) .AND. &
                                        eqmatrix(:, conf4(m,2)) == eqmatrix(op, conf4(m,2)) .AND. &
                                        eqmatrix(:, conf4(m,3)) == eqmatrix(op, conf4(m,3)) .AND. &
                                        eqmatrix(:, conf4(m,4)) == eqmatrix(op, conf4(m,4)))
                            IF (aux == 0) CYCLE

                            CALL sort_int_small(arr4, 4)

                            prod = energies4(m) - E0
                            DO idx = 1, 4
                                prod = prod - dE1(arr4(idx))
                            END DO
                            DO p = 1, 3
                                DO q = p+1, 4
                                    prod = prod - dE2(arr4(p), arr4(q))
                                END DO
                            END DO

                            triple_sum = 0.0_dp
                            IF (n_dE3 > 0 .AND. ALLOCATED(dE3_i)) THEN
                                triple_sum = triple_sum + get_low_triplet(arr4(1), arr4(2), arr4(3))
                                triple_sum = triple_sum + get_low_triplet(arr4(1), arr4(2), arr4(4))
                                triple_sum = triple_sum + get_low_triplet(arr4(1), arr4(3), arr4(4))
                                triple_sum = triple_sum + get_low_triplet(arr4(2), arr4(3), arr4(4))
                            END IF

                            prod = (prod - triple_sum) / REAL(aux, dp)

                            found = 0
                            IF (n_dE4 > 0 .AND. ALLOCATED(dE4_i)) THEN
                                DO idx = 1, n_dE4
                                    IF (dE4_i(idx) == arr4(1) .AND. dE4_j(idx) == arr4(2) .AND. &
                                        dE4_k(idx) == arr4(3) .AND. dE4_l(idx) == arr4(4)) THEN
                                        dE4_val(idx) = dE4_val(idx) + prod
                                        found = 1
                                        EXIT
                                    END IF
                                END DO
                            END IF
                            IF (found == 0) THEN
                                IF (.NOT. ALLOCATED(dE4_i)) THEN
                                    ALLOCATE(dE4_i(1), dE4_j(1), dE4_k(1), dE4_l(1), dE4_val(1))
                                    dE4_i(1) = arr4(1)
                                    dE4_j(1) = arr4(2)
                                    dE4_k(1) = arr4(3)
                                    dE4_l(1) = arr4(4)
                                    dE4_val(1) = prod
                                    n_dE4 = 1
                                ELSE
                                    ALLOCATE(tmp4_i(n_dE4+1), tmp4_j(n_dE4+1), tmp4_k(n_dE4+1), tmp4_l(n_dE4+1))
                                    ALLOCATE(tmp4_v(n_dE4+1))
                                    tmp4_i(1:n_dE4) = dE4_i(1:n_dE4)
                                    tmp4_j(1:n_dE4) = dE4_j(1:n_dE4)
                                    tmp4_k(1:n_dE4) = dE4_k(1:n_dE4)
                                    tmp4_l(1:n_dE4) = dE4_l(1:n_dE4)
                                    tmp4_v(1:n_dE4) = dE4_val(1:n_dE4)
                                    tmp4_i(n_dE4+1) = arr4(1)
                                    tmp4_j(n_dE4+1) = arr4(2)
                                    tmp4_k(n_dE4+1) = arr4(3)
                                    tmp4_l(n_dE4+1) = arr4(4)
                                    tmp4_v(n_dE4+1) = prod
                                    DEALLOCATE(dE4_i, dE4_j, dE4_k, dE4_l, dE4_val)
                                    ALLOCATE(dE4_i(n_dE4+1), dE4_j(n_dE4+1), dE4_k(n_dE4+1), dE4_l(n_dE4+1), dE4_val(n_dE4+1))
                                    dE4_i = tmp4_i
                                    dE4_j = tmp4_j
                                    dE4_k = tmp4_k
                                    dE4_l = tmp4_l
                                    dE4_val = tmp4_v
                                    DEALLOCATE(tmp4_i, tmp4_j, tmp4_k, tmp4_l, tmp4_v)
                                    n_dE4 = n_dE4 + 1
                                END IF
                            END IF
                        END DO quad_operator
                    END DO
                END IF
            END IF
        END IF
    END IF

    ! Clean up temporary arrays
    IF (Mm3 == 0) THEN
        ! ensure dE3 is allocated to avoid unallocated references later; allocate zero-sized or skip
        ! If no three-body data, do nothing (dE3 remains unallocated)
    END IF

    CALL load_high_order_expansion()

    ! Log summary of computed terms
    cnt1 = COUNT(ABS(dE1) > 1.0E-12_dp)
    cnt2 = 0
    DO i = 1, npos
        DO j = i+1, npos
            IF (ABS(dE2(i,j)) > 1.0E-12_dp) cnt2 = cnt2 + 1
        END DO
    END DO
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') 'init_energy_calc: summary: npos=', npos, &
        ', nop=', nop, ', nonzero dE1=', cnt1, ', nonzero dE2_pairs=', cnt2, &
        ', dE3_triples=', n_dE3, ', dE4_quads=', n_dE4
    WRITE(*,'(A,I0)') 'init_energy_calc: low-side max order = ', max_low_order
    IF (high_base_loaded) THEN
        WRITE(*,'(A,F12.6,A,I0,A,I0,A,I0)') 'init_energy_calc: high-side E0 = ', E0_high, &
            ', max hole order = ', max_high_order, ', hole triples=', n_h_dE3, &
            ', hole quads=', n_h_dE4
    ELSE
        WRITE(*,'(A)') 'init_energy_calc: high-side cluster data not available'
    END IF
    DEALLOCATE(conf1, conf2, omega1, omega2, energies1, energies2)
    IF (ALLOCATED(conf3)) DEALLOCATE(conf3)
    IF (ALLOCATED(omega3)) DEALLOCATE(omega3)
    IF (ALLOCATED(energies3)) DEALLOCATE(energies3)
    IF (ALLOCATED(conf4)) DEALLOCATE(conf4)
    IF (ALLOCATED(omega4)) DEALLOCATE(omega4)
    IF (ALLOCATED(energies4)) DEALLOCATE(energies4)
        
    END SUBROUTINE init_energy_calc

    SUBROUTINE load_high_order_expansion()
        IMPLICIT NONE
        CHARACTER(len=32) :: dirname
        CHARACTER(len=256) :: filename
        INTEGER :: hole_count, ge_count, io_stat

    IF (npos <= 0) RETURN
    IF (.NOT. ALLOCATED(eqmatrix)) RETURN
        IF (.NOT. find_cluster_dir(npos, dirname)) RETURN

        filename = TRIM(dirname)//'/ENERGIES'
        OPEN(UNIT=24, FILE=TRIM(filename), STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'load_high_order_expansion: cannot open ', TRIM(filename)
            RETURN
        END IF
        READ(24,*,IOSTAT=io_stat) E0_high
        CLOSE(24)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'load_high_order_expansion: failed reading ', TRIM(filename)
            E0_high = 0.0_dp
            RETURN
        END IF
        high_base_loaded = .true.

    IF (ALLOCATED(hE1)) DEALLOCATE(hE1)
    IF (ALLOCATED(hE2)) DEALLOCATE(hE2)
    IF (ALLOCATED(hE3_i)) DEALLOCATE(hE3_i)
    IF (ALLOCATED(hE3_j)) DEALLOCATE(hE3_j)
    IF (ALLOCATED(hE3_k)) DEALLOCATE(hE3_k)
    IF (ALLOCATED(hE3_val)) DEALLOCATE(hE3_val)
    IF (ALLOCATED(hE4_i)) DEALLOCATE(hE4_i)
    IF (ALLOCATED(hE4_j)) DEALLOCATE(hE4_j)
    IF (ALLOCATED(hE4_k)) DEALLOCATE(hE4_k)
    IF (ALLOCATED(hE4_l)) DEALLOCATE(hE4_l)
    IF (ALLOCATED(hE4_val)) DEALLOCATE(hE4_val)
    n_h_dE3 = 0
    n_h_dE4 = 0
        max_high_order = 0

    DO hole_count = 1, MIN(4, npos)
            ge_count = npos - hole_count
            IF (ge_count < 0) EXIT
            IF (.NOT. find_cluster_dir(ge_count, dirname)) CYCLE
            CALL process_high_level(TRIM(dirname), ge_count, hole_count)
            max_high_order = MAX(max_high_order, hole_count)
        END DO
    END SUBROUTINE load_high_order_expansion

    SUBROUTINE process_high_level(dirname, ge_count, hole_count)
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(IN) :: dirname
        INTEGER, INTENT(IN) :: ge_count, hole_count
        CHARACTER(len=256) :: filename, line
        INTEGER :: io_stat, Mm, m, aux, omega_tmp
        INTEGER, ALLOCATABLE :: hole_conf(:,:)
        REAL(dp), ALLOCATABLE :: energies(:)
        INTEGER, ALLOCATABLE :: ge_positions(:)
        LOGICAL, ALLOCATABLE :: is_ge(:)
        INTEGER, ALLOCATABLE :: tmpi(:), tmpj(:), tmpk(:)
        REAL(dp), ALLOCATABLE :: tmpv(:)
        INTEGER, ALLOCATABLE :: tmp4_i(:), tmp4_j(:), tmp4_k(:), tmp4_l(:)
        REAL(dp), ALLOCATABLE :: tmp4_v(:)
        INTEGER :: j, hole_idx, p, q
        INTEGER :: denom, op
        INTEGER :: mapped_i, mapped_j, mapped_k
        INTEGER :: ii, jj, kk, idx, found_idx
        REAL(dp) :: contrib
        INTEGER :: arr4(4)
        REAL(dp) :: triple_sum
        LOGICAL :: exists

        IF (hole_count <= 0) RETURN

        filename = TRIM(dirname)//'/OUTSOD'
        INQUIRE(FILE=TRIM(filename), EXIST=exists)
        IF (.NOT. exists) RETURN
        OPEN(UNIT=25, FILE=TRIM(filename), STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'process_high_level: cannot open ', TRIM(filename)
            RETURN
        END IF

        Mm = 0
        DO
            READ(25,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) EXIT
            IF (INDEX(line,'configuration') /= 0) THEN
                READ(line,*,IOSTAT=io_stat) Mm
                IF (io_stat == 0) EXIT
            END IF
        END DO
        IF (Mm <= 0) THEN
            CLOSE(25)
            RETURN
        END IF

        ALLOCATE(hole_conf(Mm, hole_count))
        ALLOCATE(energies(Mm))
        ALLOCATE(ge_positions(MAX(ge_count,1)))
        ALLOCATE(is_ge(npos))

        REWIND(25)
        m = 0
        DO
            READ(25,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) EXIT
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            IF (INDEX(line,'substitutions') /= 0) CYCLE
            IF (INDEX(line,'configuration') /= 0) CYCLE
            BACKSPACE(25)
            ge_positions = 0
            omega_tmp = 0
            READ(25,*,IOSTAT=io_stat) aux, omega_tmp, (ge_positions(j), j=1, ge_count)
            IF (io_stat /= 0) EXIT
            IF (omega_tmp <= 0) omega_tmp = 1
            m = m + 1
            hole_conf(m, :) = 0
            is_ge = .FALSE.
            DO j = 1, ge_count
                IF (ge_positions(j) >= 1 .AND. ge_positions(j) <= npos) THEN
                    is_ge(ge_positions(j)) = .TRUE.
                END IF
            END DO
            hole_idx = 0
            DO j = 1, npos
                IF (.NOT. is_ge(j)) THEN
                    hole_idx = hole_idx + 1
                    IF (hole_idx <= hole_count) hole_conf(m, hole_idx) = j
                END IF
            END DO
            IF (hole_idx /= hole_count) THEN
                WRITE(*,*) 'Warning: inconsistent high-level configuration in ', TRIM(dirname)
            END IF
            IF (m == Mm) EXIT
        END DO
        CLOSE(25)
        DEALLOCATE(is_ge)
        DEALLOCATE(ge_positions)

        filename = TRIM(dirname)//'/ENERGIES'
        OPEN(UNIT=26, FILE=TRIM(filename), STATUS='OLD', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'process_high_level: cannot open ', TRIM(filename)
            DEALLOCATE(hole_conf, energies)
            RETURN
        END IF
        DO m = 1, Mm
            READ(26,*,IOSTAT=io_stat) energies(m)
            IF (io_stat /= 0) EXIT
        END DO
        CLOSE(26)
        IF (io_stat /= 0) THEN
            DEALLOCATE(hole_conf, energies)
            RETURN
        END IF

        SELECT CASE (hole_count)
        CASE (1)
            IF (.NOT. ALLOCATED(hE1)) THEN
                ALLOCATE(hE1(npos))
            END IF
            hE1 = 0.0_dp
            DO m = 1, Mm
                DO op = 1, nop
                    mapped_i = eqmatrix(op, hole_conf(m,1))
                    denom = COUNT(eqmatrix(:, hole_conf(m,1)) == mapped_i)
                    IF (denom <= 0) CYCLE
                    contrib = (energies(m) - E0_high) / REAL(denom, dp)
                    hE1(mapped_i) = hE1(mapped_i) + contrib
                END DO
            END DO
        CASE (2)
            IF (.NOT. ALLOCATED(hE1)) THEN
                DEALLOCATE(hole_conf, energies)
                RETURN
            END IF
            IF (.NOT. ALLOCATED(hE2)) THEN
                ALLOCATE(hE2(npos,npos))
            END IF
            hE2 = 0.0_dp
            DO m = 1, Mm
                DO op = 1, nop
                    mapped_i = eqmatrix(op, hole_conf(m,1))
                    mapped_j = eqmatrix(op, hole_conf(m,2))
                    IF (mapped_i < mapped_j) THEN
                        denom = COUNT(eqmatrix(:, hole_conf(m,1)) == mapped_i .AND. &
                                      eqmatrix(:, hole_conf(m,2)) == mapped_j)
                        IF (denom <= 0) CYCLE
                        contrib = energies(m) - E0_high - hE1(mapped_i) - hE1(mapped_j)
                        contrib = contrib / REAL(denom, dp)
                        hE2(mapped_i, mapped_j) = hE2(mapped_i, mapped_j) + contrib
                        hE2(mapped_j, mapped_i) = hE2(mapped_i, mapped_j)
                    END IF
                END DO
            END DO
        CASE (3)
            IF (.NOT. ALLOCATED(hE1) .OR. .NOT. ALLOCATED(hE2)) THEN
                DEALLOCATE(hole_conf, energies)
                RETURN
            END IF
            IF (ALLOCATED(hE3_i)) THEN
                DEALLOCATE(hE3_i, hE3_j, hE3_k, hE3_val)
            END IF
            n_h_dE3 = 0
            DO m = 1, Mm
                DO op = 1, nop
                    mapped_i = eqmatrix(op, hole_conf(m,1))
                    mapped_j = eqmatrix(op, hole_conf(m,2))
                    mapped_k = eqmatrix(op, hole_conf(m,3))
                    IF (mapped_i == mapped_j .OR. mapped_i == mapped_k .OR. mapped_j == mapped_k) CYCLE
                    ii = MIN(mapped_i, MIN(mapped_j, mapped_k))
                    kk = MAX(mapped_i, MAX(mapped_j, mapped_k))
                    jj = mapped_i + mapped_j + mapped_k - ii - kk
                    denom = COUNT(eqmatrix(:, hole_conf(m,1)) == eqmatrix(op, hole_conf(m,1)) .AND. &
                                  eqmatrix(:, hole_conf(m,2)) == eqmatrix(op, hole_conf(m,2)) .AND. &
                                  eqmatrix(:, hole_conf(m,3)) == eqmatrix(op, hole_conf(m,3)))
                    IF (denom <= 0) CYCLE
                    contrib = energies(m) - E0_high - hE1(ii) - hE1(jj) - hE1(kk) - &
                              (hE2(MIN(ii,jj),MAX(ii,jj)) + hE2(MIN(ii,kk),MAX(ii,kk)) + hE2(MIN(jj,kk),MAX(jj,kk)))
                    contrib = contrib / REAL(denom, dp)
                    found_idx = 0
                    IF (n_h_dE3 > 0 .AND. ALLOCATED(hE3_i)) THEN
                        DO idx = 1, n_h_dE3
                            IF (hE3_i(idx) == ii .AND. hE3_j(idx) == jj .AND. hE3_k(idx) == kk) THEN
                                hE3_val(idx) = hE3_val(idx) + contrib
                                found_idx = 1
                                EXIT
                            END IF
                        END DO
                    END IF
                    IF (found_idx == 0) THEN
                        IF (.NOT. ALLOCATED(hE3_i)) THEN
                            ALLOCATE(hE3_i(1), hE3_j(1), hE3_k(1), hE3_val(1))
                            n_h_dE3 = 1
                            hE3_i(1) = ii
                            hE3_j(1) = jj
                            hE3_k(1) = kk
                            hE3_val(1) = contrib
                        ELSE
                            ALLOCATE(tmpi(n_h_dE3+1), tmpj(n_h_dE3+1), tmpk(n_h_dE3+1), tmpv(n_h_dE3+1))
                            tmpi(1:n_h_dE3) = hE3_i(1:n_h_dE3)
                            tmpj(1:n_h_dE3) = hE3_j(1:n_h_dE3)
                            tmpk(1:n_h_dE3) = hE3_k(1:n_h_dE3)
                            tmpv(1:n_h_dE3) = hE3_val(1:n_h_dE3)
                            tmpi(n_h_dE3+1) = ii
                            tmpj(n_h_dE3+1) = jj
                            tmpk(n_h_dE3+1) = kk
                            tmpv(n_h_dE3+1) = contrib
                            n_h_dE3 = n_h_dE3 + 1
                            DEALLOCATE(hE3_i, hE3_j, hE3_k, hE3_val)
                            ALLOCATE(hE3_i(n_h_dE3), hE3_j(n_h_dE3), hE3_k(n_h_dE3), hE3_val(n_h_dE3))
                            hE3_i = tmpi
                            hE3_j = tmpj
                            hE3_k = tmpk
                            hE3_val = tmpv
                            DEALLOCATE(tmpi, tmpj, tmpk, tmpv)
                        END IF
                    END IF
                END DO
            END DO
        CASE (4)
            IF (.NOT. ALLOCATED(hE1) .OR. .NOT. ALLOCATED(hE2)) THEN
                DEALLOCATE(hole_conf, energies)
                RETURN
            END IF
            IF (ALLOCATED(hE4_i)) THEN
                DEALLOCATE(hE4_i, hE4_j, hE4_k, hE4_l, hE4_val)
            END IF
            n_h_dE4 = 0
            DO m = 1, Mm
                quad_op: DO op = 1, nop
                    arr4 = (/ eqmatrix(op, hole_conf(m,1)), eqmatrix(op, hole_conf(m,2)), &
                              eqmatrix(op, hole_conf(m,3)), eqmatrix(op, hole_conf(m,4)) /)
                    DO j = 1, 3
                        DO p = j+1, 4
                            IF (arr4(j) == arr4(p)) CYCLE quad_op
                        END DO
                    END DO
                    denom = COUNT(eqmatrix(:, hole_conf(m,1)) == eqmatrix(op, hole_conf(m,1)) .AND. &
                                   eqmatrix(:, hole_conf(m,2)) == eqmatrix(op, hole_conf(m,2)) .AND. &
                                   eqmatrix(:, hole_conf(m,3)) == eqmatrix(op, hole_conf(m,3)) .AND. &
                                   eqmatrix(:, hole_conf(m,4)) == eqmatrix(op, hole_conf(m,4)))
                    IF (denom <= 0) CYCLE
                    CALL sort_int_small(arr4, 4)
                    contrib = energies(m) - E0_high
                    DO j = 1, 4
                        contrib = contrib - hE1(arr4(j))
                    END DO
                    DO j = 1, 3
                        DO p = j+1, 4
                            contrib = contrib - hE2(arr4(j), arr4(p))
                        END DO
                    END DO
                    triple_sum = 0.0_dp
                    IF (n_h_dE3 > 0 .AND. ALLOCATED(hE3_i)) THEN
                        triple_sum = triple_sum + get_high_triplet(arr4(1), arr4(2), arr4(3))
                        triple_sum = triple_sum + get_high_triplet(arr4(1), arr4(2), arr4(4))
                        triple_sum = triple_sum + get_high_triplet(arr4(1), arr4(3), arr4(4))
                        triple_sum = triple_sum + get_high_triplet(arr4(2), arr4(3), arr4(4))
                    END IF
                    contrib = (contrib - triple_sum) / REAL(denom, dp)
                    found_idx = 0
                    IF (n_h_dE4 > 0 .AND. ALLOCATED(hE4_i)) THEN
                        DO idx = 1, n_h_dE4
                            IF (hE4_i(idx) == arr4(1) .AND. hE4_j(idx) == arr4(2) .AND. &
                                hE4_k(idx) == arr4(3) .AND. hE4_l(idx) == arr4(4)) THEN
                                hE4_val(idx) = hE4_val(idx) + contrib
                                found_idx = 1
                                EXIT
                            END IF
                        END DO
                    END IF
                    IF (found_idx == 0) THEN
                        IF (.NOT. ALLOCATED(hE4_i)) THEN
                            ALLOCATE(hE4_i(1), hE4_j(1), hE4_k(1), hE4_l(1), hE4_val(1))
                            hE4_i(1) = arr4(1)
                            hE4_j(1) = arr4(2)
                            hE4_k(1) = arr4(3)
                            hE4_l(1) = arr4(4)
                            hE4_val(1) = contrib
                            n_h_dE4 = 1
                        ELSE
                            ALLOCATE(tmp4_i(n_h_dE4+1), tmp4_j(n_h_dE4+1), tmp4_k(n_h_dE4+1), tmp4_l(n_h_dE4+1))
                            ALLOCATE(tmp4_v(n_h_dE4+1))
                            tmp4_i(1:n_h_dE4) = hE4_i(1:n_h_dE4)
                            tmp4_j(1:n_h_dE4) = hE4_j(1:n_h_dE4)
                            tmp4_k(1:n_h_dE4) = hE4_k(1:n_h_dE4)
                            tmp4_l(1:n_h_dE4) = hE4_l(1:n_h_dE4)
                            tmp4_v(1:n_h_dE4) = hE4_val(1:n_h_dE4)
                            tmp4_i(n_h_dE4+1) = arr4(1)
                            tmp4_j(n_h_dE4+1) = arr4(2)
                            tmp4_k(n_h_dE4+1) = arr4(3)
                            tmp4_l(n_h_dE4+1) = arr4(4)
                            tmp4_v(n_h_dE4+1) = contrib
                            DEALLOCATE(hE4_i, hE4_j, hE4_k, hE4_l, hE4_val)
                            ALLOCATE(hE4_i(n_h_dE4+1), hE4_j(n_h_dE4+1), hE4_k(n_h_dE4+1), hE4_l(n_h_dE4+1), hE4_val(n_h_dE4+1))
                            hE4_i = tmp4_i
                            hE4_j = tmp4_j
                            hE4_k = tmp4_k
                            hE4_l = tmp4_l
                            hE4_val = tmp4_v
                            DEALLOCATE(tmp4_i, tmp4_j, tmp4_k, tmp4_l, tmp4_v)
                            n_h_dE4 = n_h_dE4 + 1
                        END IF
                    END IF
                END DO quad_op
            END DO
        END SELECT

        WRITE(*,*) 'load_high_order_expansion: processed holes=', hole_count, ' from ', TRIM(dirname)

        DEALLOCATE(hole_conf, energies)
    END SUBROUTINE process_high_level

    LOGICAL FUNCTION find_cluster_dir(level, dirname)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: level
        CHARACTER(len=*), INTENT(OUT) :: dirname
        CHARACTER(len=32) :: candidate
        LOGICAL :: exists

        dirname = ''

        WRITE(candidate,'("n",I2.2)') level
        candidate = TRIM(ADJUSTL(candidate))
        INQUIRE(FILE=TRIM(candidate)//'/OUTSOD', EXIST=exists)
        IF (.NOT. exists) INQUIRE(FILE=TRIM(candidate)//'/ENERGIES', EXIST=exists)
        IF (exists) THEN
            dirname = TRIM(candidate)
            find_cluster_dir = .TRUE.
            RETURN
        END IF

        WRITE(candidate,'("n",I3.3)') level
        candidate = TRIM(ADJUSTL(candidate))
        INQUIRE(FILE=TRIM(candidate)//'/OUTSOD', EXIST=exists)
        IF (.NOT. exists) INQUIRE(FILE=TRIM(candidate)//'/ENERGIES', EXIST=exists)
        IF (exists) THEN
            dirname = TRIM(candidate)
            find_cluster_dir = .TRUE.
            RETURN
        END IF

        WRITE(candidate,'(A,I0)') 'n', level
        candidate = TRIM(ADJUSTL(candidate))
        INQUIRE(FILE=TRIM(candidate)//'/OUTSOD', EXIST=exists)
        IF (.NOT. exists) INQUIRE(FILE=TRIM(candidate)//'/ENERGIES', EXIST=exists)
        IF (exists) THEN
            dirname = TRIM(candidate)
            find_cluster_dir = .TRUE.
        ELSE
            find_cluster_dir = .FALSE.
        END IF
    END FUNCTION find_cluster_dir

        SUBROUTINE sort_int_small(arr, n)
            IMPLICIT NONE
            INTEGER, INTENT(INOUT) :: arr(:)
            INTEGER, INTENT(IN) :: n
            INTEGER :: i, j, tmp, limit

            limit = MIN(n, SIZE(arr))
            IF (limit <= 1) RETURN

            DO i = 1, limit-1
                DO j = i+1, limit
                    IF (arr(j) < arr(i)) THEN
                        tmp = arr(i)
                        arr(i) = arr(j)
                        arr(j) = tmp
                    END IF
                END DO
            END DO
        END SUBROUTINE sort_int_small

        REAL(dp) FUNCTION get_low_triplet(i, j, k) RESULT(val)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i, j, k
            INTEGER :: ii, jj, kk, idx_local

            val = 0.0_dp
            IF (n_dE3 <= 0) RETURN
            IF (.NOT. ALLOCATED(dE3_i)) RETURN

            ii = MIN(i, MIN(j, k))
            kk = MAX(i, MAX(j, k))
            jj = i + j + k - ii - kk

            DO idx_local = 1, n_dE3
                IF (dE3_i(idx_local) == ii .AND. dE3_j(idx_local) == jj .AND. dE3_k(idx_local) == kk) THEN
                    val = dE3_val(idx_local)
                    RETURN
                END IF
            END DO
        END FUNCTION get_low_triplet

        REAL(dp) FUNCTION get_high_triplet(i, j, k) RESULT(val)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i, j, k
            INTEGER :: ii, jj, kk, idx_local

            val = 0.0_dp
            IF (n_h_dE3 <= 0) RETURN
            IF (.NOT. ALLOCATED(hE3_i)) RETURN

            ii = MIN(i, MIN(j, k))
            kk = MAX(i, MAX(j, k))
            jj = i + j + k - ii - kk

            DO idx_local = 1, n_h_dE3
                IF (hE3_i(idx_local) == ii .AND. hE3_j(idx_local) == jj .AND. hE3_k(idx_local) == kk) THEN
                    val = hE3_val(idx_local)
                    RETURN
                END IF
            END DO
        END FUNCTION get_high_triplet

    SUBROUTINE calculate_structure_energy(config, n_sites, energy, energy_low_side, energy_high_side, low_contrib, high_contrib)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n_sites
        INTEGER, INTENT(IN) :: config(n_sites)
        REAL(dp), INTENT(OUT) :: energy
        REAL(dp), INTENT(OUT), OPTIONAL :: energy_low_side
        REAL(dp), INTENT(OUT), OPTIONAL :: energy_high_side
        REAL(dp), INTENT(OUT), OPTIONAL :: low_contrib(:)
        REAL(dp), INTENT(OUT), OPTIONAL :: high_contrib(:)
        INTEGER :: i, j, k, l, op
        INTEGER :: mapped_i, mapped_j, mapped_k, mapped_l
        INTEGER :: ii, jj, kk, idx
        INTEGER :: arr4(4)
        REAL(dp) :: min_energy_low, min_energy_high, energy_tmp, energy_high_tmp
        REAL(dp) :: huge_val
        INTEGER :: nge, nsi
        INTEGER :: nge_counter, nsi_counter
        LOGICAL :: can_use_high
        REAL(dp) :: x_ge, w_low, w_high, mix_beta
        REAL(dp) :: best_low_contrib(4), best_high_contrib(4)
        REAL(dp) :: low_contrib_local(4), high_contrib_local(4)
        REAL(dp) :: min_energy_low_thread, min_energy_high_thread
        REAL(dp) :: best_low_contrib_thread(4), best_high_contrib_thread(4)
        INTEGER, ALLOCATABLE :: ge_pos_buf(:), si_pos_buf(:)
        INTEGER, ALLOCATABLE :: ge_map_buf_local(:), si_map_buf_local(:)

        huge_val = HUGE(1.0_dp)
        min_energy_low = huge_val
        min_energy_high = huge_val
        best_low_contrib = 0.0_dp
        best_high_contrib = 0.0_dp

        ! Ensure we have a site->position mapping
        IF (.NOT. ALLOCATED(subpos)) THEN
            WRITE(*,*) 'Error: substitution position mapping (subpos) not initialized.'
            STOP
        END IF
        IF (SIZE(subpos) < n_sites) THEN
            WRITE(*,*) 'Error: subpos size', SIZE(subpos), 'is smaller than n_sites', n_sites
            STOP
        END IF

        nge = COUNT(config == 2)
        nsi = n_sites - nge

        IF (nge > 0) ALLOCATE(ge_pos_buf(nge))
        IF (nsi > 0) ALLOCATE(si_pos_buf(nsi))

        IF (nge > 0 .OR. nsi > 0) THEN
            nge_counter = 0
            nsi_counter = 0
            DO i = 1, n_sites
                IF (config(i) == 2) THEN
                    nge_counter = nge_counter + 1
                    ge_pos_buf(nge_counter) = subpos(i)
                ELSE
                    nsi_counter = nsi_counter + 1
                    si_pos_buf(nsi_counter) = subpos(i)
                END IF
            END DO
        END IF

        can_use_high = high_base_loaded
        IF (nsi > 0 .AND. .NOT. ALLOCATED(hE1)) can_use_high = .FALSE.

        ! Try all symmetry operations; OpenMP distributes operations when available
        !$omp parallel default(shared) private(op, energy_tmp, energy_high_tmp, i, j, k, l, mapped_i, mapped_j, mapped_k, mapped_l, ii, jj, kk, idx, arr4, &
        !$omp&     min_energy_low_thread, min_energy_high_thread, best_low_contrib_thread, best_high_contrib_thread, low_contrib_local, high_contrib_local, &
        !$omp&     ge_map_buf_local, si_map_buf_local)

            min_energy_low_thread = huge_val
            min_energy_high_thread = huge_val
            best_low_contrib_thread = 0.0_dp
            best_high_contrib_thread = 0.0_dp

            IF (nge > 0) ALLOCATE(ge_map_buf_local(nge))
            IF (can_use_high .AND. nsi > 0) ALLOCATE(si_map_buf_local(nsi))

            !$omp do schedule(static)
            DO op = 1, nop
                energy_tmp = E0
                low_contrib_local = 0.0_dp

                IF (nge > 0) THEN
                    DO i = 1, nge
                        ge_map_buf_local(i) = eqmatrix(op, ge_pos_buf(i))
                        low_contrib_local(1) = low_contrib_local(1) + dE1(ge_map_buf_local(i))
                        energy_tmp = energy_tmp + dE1(ge_map_buf_local(i))
                    END DO
                END IF

                IF (nge > 1) THEN
                    DO i = 1, nge - 1
                        mapped_i = ge_map_buf_local(i)
                        DO j = i + 1, nge
                            mapped_j = ge_map_buf_local(j)
                            IF (mapped_i < mapped_j) THEN
                                low_contrib_local(2) = low_contrib_local(2) + dE2(mapped_i, mapped_j)
                                energy_tmp = energy_tmp + dE2(mapped_i, mapped_j)
                            ELSE
                                low_contrib_local(2) = low_contrib_local(2) + dE2(mapped_j, mapped_i)
                                energy_tmp = energy_tmp + dE2(mapped_j, mapped_i)
                            END IF
                        END DO
                    END DO
                END IF

                IF (nge > 2 .AND. ALLOCATED(dE3_i) .AND. n_dE3 > 0) THEN
                    DO i = 1, nge - 2
                        mapped_i = ge_map_buf_local(i)
                        DO j = i + 1, nge - 1
                            mapped_j = ge_map_buf_local(j)
                            DO k = j + 1, nge
                                mapped_k = ge_map_buf_local(k)
                                ii = MIN(mapped_i, MIN(mapped_j, mapped_k))
                                kk = MAX(mapped_i, MAX(mapped_j, mapped_k))
                                jj = mapped_i + mapped_j + mapped_k - ii - kk
                                DO idx = 1, n_dE3
                                    IF (dE3_i(idx) == ii .AND. dE3_j(idx) == jj .AND. dE3_k(idx) == kk) THEN
                                        low_contrib_local(3) = low_contrib_local(3) + dE3_val(idx)
                                        energy_tmp = energy_tmp + dE3_val(idx)
                                        EXIT
                                    END IF
                                END DO
                            END DO
                        END DO
                    END DO
                END IF

                IF (nge > 3 .AND. ALLOCATED(dE4_i) .AND. n_dE4 > 0) THEN
                    DO i = 1, nge - 3
                        mapped_i = ge_map_buf_local(i)
                        DO j = i + 1, nge - 2
                            mapped_j = ge_map_buf_local(j)
                            DO k = j + 1, nge - 1
                                mapped_k = ge_map_buf_local(k)
                                DO l = k + 1, nge
                                    mapped_l = ge_map_buf_local(l)
                                    arr4 = (/ mapped_i, mapped_j, mapped_k, mapped_l /)
                                    CALL sort_int_small(arr4, 4)
                                    DO idx = 1, n_dE4
                                        IF (dE4_i(idx) == arr4(1) .AND. dE4_j(idx) == arr4(2) .AND. &
                                            dE4_k(idx) == arr4(3) .AND. dE4_l(idx) == arr4(4)) THEN
                                            low_contrib_local(4) = low_contrib_local(4) + dE4_val(idx)
                                            energy_tmp = energy_tmp + dE4_val(idx)
                                            EXIT
                                        END IF
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END IF

                IF (energy_tmp < min_energy_low_thread) THEN
                    min_energy_low_thread = energy_tmp
                    best_low_contrib_thread = low_contrib_local
                END IF

                IF (can_use_high) THEN
                    energy_high_tmp = E0_high
                    high_contrib_local = 0.0_dp

                    IF (nsi > 0) THEN
                        DO i = 1, nsi
                            si_map_buf_local(i) = eqmatrix(op, si_pos_buf(i))
                        END DO
                    END IF

                    IF (nsi > 0 .AND. ALLOCATED(hE1)) THEN
                        DO i = 1, nsi
                            high_contrib_local(1) = high_contrib_local(1) + hE1(si_map_buf_local(i))
                            energy_high_tmp = energy_high_tmp + hE1(si_map_buf_local(i))
                        END DO
                    END IF

                    IF (nsi > 1 .AND. ALLOCATED(hE2)) THEN
                        DO i = 1, nsi - 1
                            mapped_i = si_map_buf_local(i)
                            DO j = i + 1, nsi
                                mapped_j = si_map_buf_local(j)
                                IF (mapped_i < mapped_j) THEN
                                    high_contrib_local(2) = high_contrib_local(2) + hE2(mapped_i, mapped_j)
                                    energy_high_tmp = energy_high_tmp + hE2(mapped_i, mapped_j)
                                ELSE
                                    high_contrib_local(2) = high_contrib_local(2) + hE2(mapped_j, mapped_i)
                                    energy_high_tmp = energy_high_tmp + hE2(mapped_j, mapped_i)
                                END IF
                            END DO
                        END DO
                    END IF

                    IF (nsi > 2 .AND. ALLOCATED(hE3_i) .AND. n_h_dE3 > 0) THEN
                        DO i = 1, nsi - 2
                            mapped_i = si_map_buf_local(i)
                            DO j = i + 1, nsi - 1
                                mapped_j = si_map_buf_local(j)
                                DO k = j + 1, nsi
                                    mapped_k = si_map_buf_local(k)
                                    ii = MIN(mapped_i, MIN(mapped_j, mapped_k))
                                    kk = MAX(mapped_i, MAX(mapped_j, mapped_k))
                                    jj = mapped_i + mapped_j + mapped_k - ii - kk
                                    DO idx = 1, n_h_dE3
                                        IF (hE3_i(idx) == ii .AND. hE3_j(idx) == jj .AND. hE3_k(idx) == kk) THEN
                                            high_contrib_local(3) = high_contrib_local(3) + hE3_val(idx)
                                            energy_high_tmp = energy_high_tmp + hE3_val(idx)
                                            EXIT
                                        END IF
                                    END DO
                                END DO
                            END DO
                        END DO
                    END IF

                    IF (nsi > 3 .AND. ALLOCATED(hE4_i) .AND. n_h_dE4 > 0) THEN
                        DO i = 1, nsi - 3
                            mapped_i = si_map_buf_local(i)
                            DO j = i + 1, nsi - 2
                                mapped_j = si_map_buf_local(j)
                                DO k = j + 1, nsi - 1
                                    mapped_k = si_map_buf_local(k)
                                    DO l = k + 1, nsi
                                        mapped_l = si_map_buf_local(l)
                                        arr4 = (/ mapped_i, mapped_j, mapped_k, mapped_l /)
                                        CALL sort_int_small(arr4, 4)
                                        DO idx = 1, n_h_dE4
                                            IF (hE4_i(idx) == arr4(1) .AND. hE4_j(idx) == arr4(2) .AND. &
                                                hE4_k(idx) == arr4(3) .AND. hE4_l(idx) == arr4(4)) THEN
                                                high_contrib_local(4) = high_contrib_local(4) + hE4_val(idx)
                                                energy_high_tmp = energy_high_tmp + hE4_val(idx)
                                                EXIT
                                            END IF
                                        END DO
                                    END DO
                                END DO
                            END DO
                        END DO
                    END IF

                    IF (energy_high_tmp < min_energy_high_thread) THEN
                        min_energy_high_thread = energy_high_tmp
                        best_high_contrib_thread = high_contrib_local
                    END IF
                END IF
            END DO
            !$omp end do

            IF (ALLOCATED(ge_map_buf_local)) DEALLOCATE(ge_map_buf_local)
            IF (ALLOCATED(si_map_buf_local)) DEALLOCATE(si_map_buf_local)

            !$omp critical
                IF (min_energy_low_thread < min_energy_low) THEN
                    min_energy_low = min_energy_low_thread
                    best_low_contrib = best_low_contrib_thread
                END IF
                IF (can_use_high .AND. min_energy_high_thread < min_energy_high) THEN
                    min_energy_high = min_energy_high_thread
                    best_high_contrib = best_high_contrib_thread
                END IF
            !$omp end critical
        !$omp end parallel

        IF (ALLOCATED(ge_pos_buf)) DEALLOCATE(ge_pos_buf)
        IF (ALLOCATED(si_pos_buf)) DEALLOCATE(si_pos_buf)

        IF (.NOT. can_use_high .OR. min_energy_high >= huge_val) THEN
            energy = min_energy_low
        ELSE
            ! Blend low/high predictions smoothly; bias Si-side when Ge fraction is small
            x_ge = REAL(nge, dp) / REAL(n_sites, dp)
            mix_beta = 8.0_dp
            w_high = 0.5_dp * (1.0_dp + TANH(mix_beta * (x_ge - 0.5_dp)))
            w_low = 1.0_dp - w_high
            energy = w_low * min_energy_low + w_high * min_energy_high
        END IF

        IF (PRESENT(energy_low_side)) energy_low_side = min_energy_low
        IF (PRESENT(energy_high_side)) THEN
            IF (.NOT. can_use_high) THEN
                energy_high_side = huge_val
            ELSE
                energy_high_side = min_energy_high
            END IF
        END IF

        IF (PRESENT(low_contrib)) THEN
            IF (SIZE(low_contrib) >= 4) THEN
                low_contrib(1:4) = best_low_contrib
            ELSE
                low_contrib(1:SIZE(low_contrib)) = best_low_contrib(1:SIZE(low_contrib))
            END IF
        END IF

        IF (PRESENT(high_contrib)) THEN
            IF (SIZE(high_contrib) >= 4) THEN
                high_contrib(1:4) = best_high_contrib
            ELSE
                high_contrib(1:SIZE(high_contrib)) = best_high_contrib(1:SIZE(high_contrib))
            END IF
        END IF

    END SUBROUTINE calculate_structure_energy

    FUNCTION get_base_energy() RESULT(base_energy)
        REAL(dp) :: base_energy
        base_energy = E0
    END FUNCTION get_base_energy

    FUNCTION get_high_base_energy() RESULT(base_energy)
        REAL(dp) :: base_energy
        base_energy = E0_high
    END FUNCTION get_high_base_energy

    INTEGER FUNCTION get_max_low_order() RESULT(order)
        order = max_low_order
    END FUNCTION get_max_low_order

    INTEGER FUNCTION get_max_high_order() RESULT(order)
        order = max_high_order
    END FUNCTION get_max_high_order

    SUBROUTINE cleanup_energy_calc()
        IMPLICIT NONE
        IF (ALLOCATED(dE1)) DEALLOCATE(dE1)
        IF (ALLOCATED(dE2)) DEALLOCATE(dE2)
        IF (ALLOCATED(dE3_i)) DEALLOCATE(dE3_i)
        IF (ALLOCATED(dE3_j)) DEALLOCATE(dE3_j)
        IF (ALLOCATED(dE3_k)) DEALLOCATE(dE3_k)
        IF (ALLOCATED(dE3_val)) DEALLOCATE(dE3_val)
        n_dE3 = 0
    IF (ALLOCATED(dE4_i)) DEALLOCATE(dE4_i)
    IF (ALLOCATED(dE4_j)) DEALLOCATE(dE4_j)
    IF (ALLOCATED(dE4_k)) DEALLOCATE(dE4_k)
    IF (ALLOCATED(dE4_l)) DEALLOCATE(dE4_l)
    IF (ALLOCATED(dE4_val)) DEALLOCATE(dE4_val)
    n_dE4 = 0
        IF (ALLOCATED(eqmatrix)) DEALLOCATE(eqmatrix)
        IF (ALLOCATED(coords)) DEALLOCATE(coords)
        IF (ALLOCATED(coords0)) DEALLOCATE(coords0)
        IF (ALLOCATED(coords1)) DEALLOCATE(coords1)
        IF (ALLOCATED(mgroup1)) DEALLOCATE(mgroup1)
        IF (ALLOCATED(vgroup1)) DEALLOCATE(vgroup1)
        IF (ALLOCATED(spat0)) DEALLOCATE(spat0)
        IF (ALLOCATED(atom_types)) DEALLOCATE(atom_types)
        IF (ALLOCATED(subpos)) DEALLOCATE(subpos)
        IF (ALLOCATED(spat1)) DEALLOCATE(spat1)
        IF (ALLOCATED(pos2coord)) DEALLOCATE(pos2coord)
        IF (ALLOCATED(hE1)) DEALLOCATE(hE1)
        IF (ALLOCATED(hE2)) DEALLOCATE(hE2)
        IF (ALLOCATED(hE3_i)) DEALLOCATE(hE3_i)
        IF (ALLOCATED(hE3_j)) DEALLOCATE(hE3_j)
        IF (ALLOCATED(hE3_k)) DEALLOCATE(hE3_k)
        IF (ALLOCATED(hE3_val)) DEALLOCATE(hE3_val)
        n_h_dE3 = 0
    IF (ALLOCATED(hE4_i)) DEALLOCATE(hE4_i)
    IF (ALLOCATED(hE4_j)) DEALLOCATE(hE4_j)
    IF (ALLOCATED(hE4_k)) DEALLOCATE(hE4_k)
    IF (ALLOCATED(hE4_l)) DEALLOCATE(hE4_l)
    IF (ALLOCATED(hE4_val)) DEALLOCATE(hE4_val)
    n_h_dE4 = 0
        E0_high = 0.0_dp
        high_base_loaded = .false.
        max_low_order = 0
        max_high_order = 0
    ! conf3/omega3/energies3 are local to init_energy_calc and were deallocated there if needed
    END SUBROUTINE cleanup_energy_calc
    
    SUBROUTINE write_vasp_file(config, n_sites, filename)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n_sites
        INTEGER, INTENT(IN) :: config(n_sites)
        CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER :: i, j
    INTEGER :: n_ge_this ! Number of Ge in this configuration
    LOGICAL, ALLOCATABLE :: is_sub(:)
        
        ! Count number of Ge substitutions in this configuration
        n_ge_this = COUNT(config == 2)

        ! Prepare substitution boolean map for Si positions (1..n_si)
        IF (ALLOCATED(is_sub)) DEALLOCATE(is_sub)
        ALLOCATE(is_sub(n_si))
        is_sub = .FALSE.
        DO j = 1, n_sites
            IF (config(j) == 2) THEN
                IF (subpos(j) >= 1 .AND. subpos(j) <= n_si) THEN
                    is_sub(subpos(j)) = .TRUE.
                ELSE
                    WRITE(*,*) 'Warning: subpos(', j,') =', subpos(j), 'is out of range [1,', n_si,']'
                END IF
            END IF
        END DO
        
        ! Write VASP file
        OPEN(UNIT=30, FILE=filename, STATUS='REPLACE')
        
        ! Write header
        WRITE(30,'(A)') 'MC generated'
        WRITE(30,'(A)') '1.0'
        
        ! Write cell parameters with correct angle transforms
        WRITE(30,'(3F12.6)') cell_params(1), 0.0, 0.0
        WRITE(30,'(3F12.6)') 0.0, cell_params(2), 0.0
        WRITE(30,'(3F12.6)') cell_params(3)*COS(cell_params(5)*PI/180.0), &
                            0.0, cell_params(3)*SIN(cell_params(5)*PI/180.0)
        
        ! Write atom types and counts - note n_si is the total Si sites in full cell
        WRITE(30,'(A)') 'Ge Si O'
        WRITE(30,'(3I6)') n_ge_this, n_si-n_ge_this, n_o
        
        ! Write 'Direct' for fractional coordinates
        WRITE(30,'(A)') 'Direct'
        
        ! First write Ge positions (iterate over position indices 1..n_si, map to coords via pos2coord)
        DO i = 1, n_si
            IF (is_sub(i)) THEN
                WRITE(30,'(3F12.6)') coords(pos2coord(i),:)
            END IF
        END DO

        ! Then write remaining Si positions
        DO i = 1, n_si
            IF (.NOT. is_sub(i)) THEN
                WRITE(30,'(3F12.6)') coords(pos2coord(i),:)
            END IF
        END DO

        ! Finally write O positions (scan deduped coords and write those with spat1==2)
        DO i = 1, nat1
            IF (spat1(i) == 2) THEN
                WRITE(30,'(3F12.6)') coords(i,:)
            END IF
        END DO
        
        CLOSE(30)
        IF (ALLOCATED(is_sub)) DEALLOCATE(is_sub)
    END SUBROUTINE write_vasp_file

    SUBROUTINE write_eqmatrix_file(fname)
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: fname
        CHARACTER(len=80) :: filename
        INTEGER :: op
        IF (.NOT. ALLOCATED(eqmatrix)) THEN
            WRITE(*,*) 'write_eqmatrix_file: eqmatrix not allocated'
            RETURN
        END IF
        IF (PRESENT(fname)) THEN
            filename = fname
        ELSE
            filename = 'EQMATRIX'
        END IF
        OPEN(UNIT=99, FILE=TRIM(filename), STATUS='REPLACE', IOSTAT=op)
        IF (op /= 0) THEN
            WRITE(*,*) 'write_eqmatrix_file: cannot open file', TRIM(filename)
            RETURN
        END IF
        WRITE(99,*) nop, npos
        DO op = 1, nop
            WRITE(99,*) eqmatrix(op,1:npos)
        END DO
        CLOSE(99)
        WRITE(*,*) 'write_eqmatrix_file: wrote', TRIM(filename)
    END SUBROUTINE write_eqmatrix_file

    SUBROUTINE get_eqmatrix(out_mat, out_nop, out_npos)
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: out_nop, out_npos
        INTEGER, ALLOCATABLE, INTENT(OUT) :: out_mat(:,:)
        INTEGER :: i
        IF (.NOT. ALLOCATED(eqmatrix)) THEN
            out_nop = 0
            out_npos = 0
            RETURN
        END IF
        out_nop = nop
        out_npos = npos
        ALLOCATE(out_mat(out_nop, out_npos))
        DO i = 1, out_nop
            out_mat(i,1:out_npos) = eqmatrix(i,1:out_npos)
        END DO
    END SUBROUTINE get_eqmatrix

END MODULE energy_calc