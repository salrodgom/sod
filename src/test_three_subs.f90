PROGRAM test_three_subs
    USE energy_calc
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER :: io_stat
    CHARACTER(len=80) :: line

    ! Local variables for parsing SGO/INSOD and building eqmatrix
    INTEGER :: i,j,k,at0,at1i,at1r,nat0,nat1,nat1r
    REAL(dp) :: tol0, prod
    INTEGER :: nop1, na, nb, nc, op, op1, t, sgo_i, pos, attmp
    REAL(dp), ALLOCATABLE :: mgroup1(:,:,:), vgroup1(:,:), coords0(:,:), coords1r(:,:), coords1(:,:), coords(:,:), vt(:,:)
    REAL(dp), ALLOCATABLE :: mgroup(:,:,:), vgroup(:,:)
    INTEGER, ALLOCATABLE :: spat0(:), spat1r(:), spat1(:), pos2coord(:), coord2pos(:)
    INTEGER, ALLOCATABLE :: fulleqmatrix(:,:)
    INTEGER :: natsp0(2), natsp1(2)
    REAL(dp) :: cell_params(6)

    ! eqmatrix and sizes
    INTEGER, ALLOCATABLE :: eqmatrix(:,:)
    INTEGER :: nop_total, npos

    ! n01/n03 structures
    INTEGER :: Mm1, Mm3, m, aux
    INTEGER, ALLOCATABLE :: conf1(:), omega1(:), omega3(:)
    INTEGER, ALLOCATABLE :: conf3(:,:)
    REAL(dp), ALLOCATABLE :: energies3(:)
    INTEGER, ALLOCATABLE :: subpos(:)
    INTEGER :: iii, jjj, kkk
    INTEGER :: ti, tj, tk, ri, rj, rk
    REAL(dp) :: coord(3)
    INTEGER, ALLOCATABLE :: config(:)

    ! For comparison
    REAL(dp) :: e_calc, e_ref
    INTEGER :: s1,s2,s3
    INTEGER :: found, opcheck, mcheck
    INTEGER :: p1, p2, p3
    LOGICAL :: stop_now
    INTEGER :: max_checks, checked_count

    ! Initialize energy model (builds internal eqmatrix/dE arrays)
    CALL init_energy_calc()
    ! Ensure EQMATRIX file exists by asking energy_calc to write it (module has generated it)
    CALL write_eqmatrix_file()

    ! Now independently build eqmatrix and subpos by parsing local files so we can find reps
    tol0 = 0.001_dp

    ! Read SGO
    OPEN(UNIT=21, FILE='SGO', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'test_three_subs: cannot open SGO'
        STOP
    END IF
    READ(21,*)   ! skip space group num
    READ(21,*) sgo_i
    nop1 = 0
    ALLOCATE(mgroup1(48,3,3))
    ALLOCATE(vgroup1(48,3))
    DO WHILE (sgo_i > 0)
        nop1 = sgo_i
        DO j = 1, 3
            READ(21,*) mgroup1(sgo_i,j,1:3), vgroup1(sgo_i,j)
        END DO
        READ(21,*) sgo_i
    END DO
    CLOSE(21)

    ! Read INSOD
    OPEN(UNIT=20, FILE='INSOD', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'test_three_subs: cannot open INSOD'
        STOP
    END IF
    ! find cell params
    DO
        READ(20,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'test_three_subs: unexpected end of INSOD'
            STOP
        END IF
        READ(line,*,IOSTAT=io_stat) cell_params
        IF (io_stat == 0) EXIT
    END DO
    ! find nsp (ignore)
    DO
        READ(20,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'test_three_subs: unexpected end of INSOD'
            STOP
        END IF
        READ(line,*,IOSTAT=io_stat) aux
        IF (io_stat == 0) EXIT
    END DO
    ! atom symbols
    DO
        READ(20,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'test_three_subs: unexpected end of INSOD'
            STOP
        END IF
        IF (TRIM(line) == '') CYCLE
        IF (line(1:1) == '#') CYCLE
        ! skip symbols
        EXIT
    END DO
    ! read natsp0
    DO
        READ(20,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'test_three_subs: unexpected end of INSOD'
            STOP
        END IF
        IF (TRIM(line) == '') CYCLE
        IF (line(1:1) == '#') CYCLE
        READ(line,*,IOSTAT=io_stat) natsp0(1), natsp0(2)
        IF (io_stat == 0) EXIT
    END DO
    nat0 = natsp0(1) + natsp0(2)
    ALLOCATE(coords0(nat0,3))
    ALLOCATE(spat0(nat0))
    i = 0
    DO WHILE (i < nat0)
        READ(20,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) THEN
            WRITE(*,*) 'test_three_subs: unexpected end reading coords'
            STOP
        END IF
        IF (TRIM(line) == '') CYCLE
        IF (line(1:1) == '#') CYCLE
        i = i + 1
        READ(line,*,IOSTAT=io_stat) coords0(i,1:3)
        spat0(i) = MERGE(1,2, i <= natsp0(1))
    END DO
    CLOSE(20)

    ! Build redundant unit cell coords1r
    nat1r = nat0 * nop1
    ALLOCATE(coords1r(nat1r,3))
    ALLOCATE(spat1r(nat1r))
    at1r = 0
    DO at0 = 1, nat0
        DO op = 1, nop1
            at1r = at1r + 1
            coords1r(at1r,1:3) = MATMUL(mgroup1(op,1:3,1:3), coords0(at0,1:3)) + vgroup1(op,1:3)
            coords1r(at1r,1:3) = coords1r(at1r,1:3) - FLOOR(coords1r(at1r,1:3))
            spat1r(at1r) = spat0(at0)
        END DO
    END DO

    ! Deduplicate
    ALLOCATE(coords1(nat1r,3))
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

    ! Count species in unit cell
    natsp1(1)=0; natsp1(2)=0
    DO at1r = 1, nat1
        IF (spat1(at1r) == 1) THEN
            natsp1(1) = natsp1(1) + 1
        ELSE
            natsp1(2) = natsp1(2) + 1
        END IF
    END DO

    ! Build pos2coord for species 1
    ALLOCATE(pos2coord(natsp1(1)))
    at1r = 0
    DO at0 = 1, nat1
        IF (spat1(at0) == 1) THEN
            at1r = at1r + 1
            pos2coord(at1r) = at0
        END IF
    END DO
    npos = SIZE(pos2coord)

    ! Find na,nb,nc and build vt (supercell) from INSOD (scan file)
    na = 1; nb = 1; nc = 1
    OPEN(UNIT=22, FILE='INSOD', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat == 0) THEN
        DO
            READ(22,'(A)',IOSTAT=io_stat) line
            IF (io_stat /= 0) EXIT
            IF (TRIM(line) == '') CYCLE
            IF (line(1:1) == '#') CYCLE
            READ(line,*,IOSTAT=io_stat) na, nb, nc
            IF (io_stat == 0) EXIT
        END DO
        CLOSE(22)
    END IF

    ALLOCATE(vt(na*nb*nc,3))
    t = 0
    DO i = 0, na-1
        DO j = 0, nb-1
            DO k = 0, nc-1
                t = t + 1
                vt(t,1) = REAL(i,dp)/REAL(na,dp)
                vt(t,2) = REAL(j,dp)/REAL(nb,dp)
                vt(t,3) = REAL(k,dp)/REAL(nc,dp)
            END DO
        END DO
    END DO

    ! Build supercell operators mgroup and vgroup
    nop_total = nop1 * na * nb * nc
    ALLOCATE(mgroup(nop_total,3,3))
    ALLOCATE(vgroup(nop_total,3))
    op = 0
    DO op1 = 1, nop1
        DO t = 1, na*nb*nc
            op = op + 1
            mgroup(op,:,:) = mgroup1(op1,:,:)
            vgroup(op,1) = vgroup1(op1,1)/REAL(na,dp) + vt(t,1)
            vgroup(op,2) = vgroup1(op1,2)/REAL(nb,dp) + vt(t,2)
            vgroup(op,3) = vgroup1(op1,3)/REAL(nc,dp) + vt(t,3)
        END DO
    END DO

    ALLOCATE(fulleqmatrix(nop_total,nat1))
    ALLOCATE(coord2pos(nat1))
    coord2pos = 0
    ! Build inverse mapping coord2pos for target positions
    DO pos = 1, SIZE(pos2coord)
        coord2pos(pos2coord(pos)) = pos
    END DO

    DO op = 1, nop_total
        DO at0 = 1, nat1
            coord = MATMUL(mgroup(op,:,:), coords1(at0,1:3)) + vgroup(op,1:3)
            coord = coord - FLOOR(coord)
            attmp = 0
            DO i = 1, nat1
                prod = DOT_PRODUCT(coord - coords1(i,1:3), coord - coords1(i,1:3))
                IF (prod .LE. tol0) THEN
                    attmp = i
                    EXIT
                END IF
            END DO
            IF (attmp == 0) THEN
                WRITE(*,*) 'test_three_subs: operator mapping failed for op=', op, 'at=', at0
                STOP
            END IF
            fulleqmatrix(op,at0) = attmp
        END DO
    END DO

    ! Extract eqmatrix for target positions
    ALLOCATE(eqmatrix(nop_total, SIZE(pos2coord)))
    DO op = 1, nop_total
        DO pos = 1, SIZE(pos2coord)
            eqmatrix(op,pos) = coord2pos(fulleqmatrix(op,pos2coord(pos)))
        END DO
    END DO

    ! eqmatrix constructed from INSOD+SGO above (nop_total, npos set)

    ! Read n01/OUTSOD for subpos mapping
    OPEN(UNIT=15, FILE='n01/OUTSOD', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'test_three_subs: cannot open n01/OUTSOD'
        STOP
    END IF
    Mm1 = 0
    DO
        READ(15,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) EXIT
        IF (INDEX(line,'configuration') /= 0) THEN
            READ(line,*,IOSTAT=io_stat) Mm1
            IF (io_stat == 0) EXIT
        END IF
    END DO
    IF (Mm1 <= 0) THEN
        WRITE(*,*) 'test_three_subs: could not determine Mm1'
        STOP
    END IF
    ALLOCATE(conf1(Mm1), omega1(Mm1))
    REWIND(15)
    m = 0
    DO
        READ(15,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) EXIT
        IF (TRIM(line) == '') CYCLE
        IF (line(1:1) == '#') CYCLE
        READ(line,*,IOSTAT=io_stat) aux, omega1(m+1), conf1(m+1)
        IF (io_stat == 0) THEN
            m = m + 1
            IF (m == Mm1) EXIT
        END IF
    END DO
    CLOSE(15)
    ALLOCATE(subpos(Mm1))
    subpos = conf1

    ! Read n03/OUTSOD
    OPEN(UNIT=17, FILE='n03/OUTSOD', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'test_three_subs: cannot open n03/OUTSOD'
        STOP
    END IF
    Mm3 = 0
    DO
        READ(17,'(A)',IOSTAT=io_stat) line
        IF (io_stat /= 0) EXIT
        IF (INDEX(line,'configuration') /= 0) THEN
            READ(line,*,IOSTAT=io_stat) Mm3
            IF (io_stat == 0) EXIT
        END IF
    END DO
    IF (Mm3 <= 0) THEN
        WRITE(*,*) 'test_three_subs: could not determine Mm3'
        STOP
    END IF
    ALLOCATE(conf3(Mm3,3))
    ALLOCATE(omega3(Mm3))
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
    CLOSE(17)

    ! Read n03/ENERGIES
    OPEN(UNIT=18, FILE='n03/ENERGIES', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,*) 'test_three_subs: cannot open n03/ENERGIES'
        STOP
    END IF
    ALLOCATE(energies3(Mm3))
    DO m = 1, Mm3
        READ(18,*) energies3(m)
    END DO
    CLOSE(18)

    ! Iterate combinations of site indices (1..Mm1 choose 3)
    max_checks = 25
    checked_count = 0
    stop_now = .FALSE.
    WRITE(*,*) 'test_three_subs: iterating combinations and matching to n03 representatives...'
    ! Enumerate triples over unit-cell positions (npos) and map to MC sites
    DO p1 = 1, npos-2
        DO p2 = p1+1, npos-1
            DO p3 = p2+1, npos
                ! positions
                i = p1
                j = p2
                k = p3
                ! canonical triple (target)
                ti = MIN(i, MIN(j,k))
                tk = MAX(i, MAX(j,k))
                tj = i + j + k - ti - tk

                found = 0
                mcheck = -1
                DO m = 1, Mm3
                    DO op = 1, SIZE(eqmatrix,1)
                        ! get mapped triple from rep m under op
                        aux = eqmatrix(op, conf3(m,1))
                        opcheck = eqmatrix(op, conf3(m,2))
                        t = eqmatrix(op, conf3(m,3))
                        ! canonicalize
                        IF (aux == opcheck .OR. aux == t .OR. opcheck == t) CYCLE
                        ri = MIN(aux, MIN(opcheck,t))
                        rk = MAX(aux, MAX(opcheck,t))
                        rj = aux + opcheck + t - ri - rk
                        IF (ri == ti .AND. rj == tj .AND. rk == tk) THEN
                            found = 1
                            mcheck = m
                            EXIT
                        END IF
                    END DO
                    IF (found == 1) EXIT
                END DO

                ! Build config array for calculate_structure_energy by mapping positions
                IF (.NOT. ALLOCATED(config)) ALLOCATE(config(Mm1))
                IF (SIZE(config) /= Mm1) THEN
                    IF (ALLOCATED(config)) DEALLOCATE(config)
                    ALLOCATE(config(Mm1))
                END IF
                config = 1
                ! For each MC site, set to 2 if its subpos equals one of the chosen positions
                DO s1 = 1, Mm1
                    IF (subpos(s1) == i .OR. subpos(s1) == j .OR. subpos(s1) == k) THEN
                        config(s1) = 2
                    END IF
                END DO

                CALL calculate_structure_energy(config, Mm1, e_calc)

                IF (found == 1) THEN
                    e_ref = energies3(mcheck)
                    WRITE(*,'(A,I3,A,I3,A,I3,A,I4,A,F15.8,A,F15.8,A,F12.8)') 'Positions:', i, ',', j, ',', k, '  rep=', mcheck, '  E_calc=', e_calc, '  E_ref=', e_ref, '  diff=', e_calc - e_ref
                ELSE
                    WRITE(*,'(A,I3,A,I3,A,I3,A,F15.8)') 'Positions:', i, ',', j, ',', k, '  NO_MATCH   E_calc=', e_calc
                END IF
                ! keep config allocated for reuse
                checked_count = checked_count + 1
                IF (checked_count >= max_checks) THEN
                    WRITE(*,'(A,I4)') 'test_three_subs: reached max checks=', checked_count
                    stop_now = .TRUE.
                    EXIT
                END IF
            END DO
            IF (stop_now) EXIT
        END DO
        IF (stop_now) EXIT
    END DO

    WRITE(*,'(A,I6)') 'test_three_subs: total checked=', checked_count

    WRITE(*,*) 'test_three_subs: done.'

    STOP
END PROGRAM test_three_subs
