!   
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\ Force - for local cell interactions  ////////////////////
! Asynchronously done when MPI_send/receive positions
SUBROUTINE Force_local_1(NS, q, cell, param, sys, ion_data, table)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CL):: cell(NS%Nlocal)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(IC):: ion_data
TYPE(VF):: table
!
!
INTEGER :: i, j, k, n, id, jd, i_id, j_id, id_table, id_radial
REAL(KIND=DP):: r, r2, x, dx(3), dv(3), df, E, Eself(Nparam)
!
sys%Epot = 0.0D0
q(:)%ff(1) = 0.0D0
q(:)%ff(2) = 0.0D0
q(:)%ff(3) = 0.0D0
!Eself(:) = - eps* param%z(:)**2*param%a/sqrtpi
!DO i=1, NS%Npt
!   id = q(i)%id
!   q(i)%pot = Eself(id)
!END DO
q(:)%pot = 0.0D0
ion_data%r_min = 100.D0
!
DO n=1, NS%Nlocal/2
   !
   ! Particle interaction for inside of each cell
   DO i=1, cell(n)%Npt-1
      id = cell(n)%node(i)      
      i_id = q(id)%id
      DO j=i+1, cell(n)%Npt
         jd = cell(n)%node(j)
         j_id = q(jd)%id
         
         r2 = 0.0D0
         DO k=1,3
            x = q(id)%xx(k) - q(jd)%xx(k)
            dx(k) = x - param%box(k)*DNINT(x/param%box(k))
            r2 = r2 + dx(k)**2
         END DO
         IF (r2 < param%Rc2) THEN
            r = DSQRT(r2)
            x = r/table%dr
            id_radial = INT(x)
            IF (id_radial < Ntabular) THEN
               !
               ! Pair interaction id
               IF (i_id == 1 .AND. j_id == 1) THEN ! U-U
                  id_table = 1 
               ELSE IF (i_id == 2 .AND. j_id == 2) THEN ! O-O
                  id_table = 2
               ELSE IF ((i_id +j_id) == 3) THEN ! U-O
                  id_table = 3
               ELSE IF (i_id == 1 .OR. j_id == 1) THEN
                  STOP "Table index error 4" !id_table = 4
               ELSE IF (i_id == 2 .OR. j_id == 2) THEN
                  STOP "Table index error 5" !id_table = 5
               ELSE
                  STOP "Table index error N"
               END IF
               !
               x = x - DBLE(id_radial)
               df = table%F(id_table, id_radial)  *(1.D0 - x) + &
                  & table%F(id_table, id_radial+1)*x
               E  = table%V(id_table, id_radial)  *(1.D0 - x) + &
                  & table%V(id_table, id_radial+1)*x
               DO k = 1,3
                  q(id)%ff(k) = q(id)%ff(k) + df*dx(k)
                  q(jd)%ff(k) = q(jd)%ff(k) - df*dx(k)
               END DO
               q(id)%pot = q(id)%pot + E/2.D0
               q(jd)%pot = q(jd)%pot + E/2.D0
               sys%Epot = sys%Epot + E
               ion_data%r_min = MIN(r, ion_data%r_min)
               !
               IF (r < Rcut_firsov .AND. DRAG) THEN
                  df = table%FIRSV(id_table, id_radial)  *(1.D0 - x) + &
                     & table%FIRSV(id_table, id_radial+1)*x
                  DO k=1,3
                     dv(k) = q(id)%xv(k) - q(jd)%xv(k)
                     q(id)%ff(k) = q(id)%ff(k) - dv(k)*df
                     q(jd)%ff(k) = q(jd)%ff(k) + dv(k)*df
                  END DO
               END IF
               IF (id_table == 4 .OR. id_table == 5) THEN 
                  IF (j_id /= 3) THEN 
                     ion_data%id = id
                     j_id = i_id
                  ELSE
                     ion_data%id = jd
                  END IF
                  CALL E_RHO(j_id, r, ion_data)
                  ion_data%tag = .TRUE.
               END IF
            END IF
         END IF
      END DO
   END DO
END DO
RETURN
END SUBROUTINE Force_local_1
!   
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\ Force - for local cell interactions  ////////////////////
! Asynchronously done when MPI_send/receive forces
SUBROUTINE Force_local_2(NS, q, cell, param, sys, ion_data, table)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CL):: cell(NS%Nlocal)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(IC):: ion_data
TYPE(VF):: table
!
!
INTEGER :: i, j, k, n, id, jd, i_id, j_id, id_table, id_radial
REAL(KIND=DP):: r, r2, x, dx(3), df, dv(3), E, v
!
!
DO n=NS%Nlocal/2+1, NS%Nlocal
   !
   ! Particle interaction for inside of each cell
   DO i=1, cell(n)%Npt-1
      id = cell(n)%node(i)   
      i_id = q(id)%id
      DO j=i+1, cell(n)%Npt
         jd = cell(n)%node(j)
         j_id = q(jd)%id
         !
         r2 = 0.0D0
         DO k=1,3
            x = q(id)%xx(k) - q(jd)%xx(k)
            dx(k) = x - param%box(k)*DNINT(x/param%box(k))
            r2 = r2 + dx(k)**2
         END DO
         IF (r2 < param%Rc2) THEN
            r = DSQRT(r2)
            x = r/table%dr
            id_radial = INT(x)
            IF (id_radial < Ntabular) THEN
               !
               ! Pair interaction id
               IF (i_id == 1 .AND. j_id == 1) THEN
                  id_table = 1
               ELSE IF (i_id == 2 .AND. j_id == 2) THEN
                  id_table = 2
               ELSE IF ((i_id +j_id) == 3) THEN
                  id_table = 3
               ELSE IF (i_id == 1 .OR. j_id == 1) THEN
                  STOP "Table index error 4" !id_table = 4
               ELSE IF (i_id == 2 .OR. j_id == 2) THEN
                  STOP "Table index error 5" !id_table = 5
               ELSE
                  STOP "Table index error N"
               END IF
               !
               x = x - DBLE(id_radial)
               df = table%F(id_table, id_radial)  *(1.D0 - x) + &
                  & table%F(id_table, id_radial+1)*x
               E  = table%V(id_table, id_radial)  *(1.D0 - x) + &
                  & table%V(id_table, id_radial+1)*x
               DO k = 1,3
                  q(id)%ff(k) = q(id)%ff(k) + df*dx(k)
                  q(jd)%ff(k) = q(jd)%ff(k) - df*dx(k)
               END DO
               q(id)%pot = q(id)%pot + E/2.D0
               q(jd)%pot = q(jd)%pot + E/2.D0
               sys%Epot = sys%Epot + E
               ion_data%r_min = MIN(r, ion_data%r_min)
               !
               IF (r < Rcut_firsov .AND. DRAG) THEN
                  df = table%FIRSV(id_table, id_radial)  *(1.D0 - x) + &
                     & table%FIRSV(id_table, id_radial+1)*x
                  DO k=1,3
                     dv(k) = q(id)%xv(k) - q(jd)%xv(k)
                     q(id)%ff(k) = q(id)%ff(k) - dv(k)*df
                     q(jd)%ff(k) = q(jd)%ff(k) + dv(k)*df
                  END DO
               END IF
               
               IF (id_table == 4 .OR. id_table == 5) THEN
                  IF (j_id /= 3) THEN 
                     ion_data%id = id
                     j_id = i_id
                  ELSE
                     ion_data%id = jd
                  END IF
                  CALL E_RHO(j_id, r, ion_data)
                  ion_data%tag = .TRUE.
               END IF
            END IF
         END IF
      END DO
   END DO
END DO
!
! Addition of electron stopping
IF (DRAG) THEN
   DO i=1, NS%Npt
      id = q(i)%id
      IF (id /=3 )THEN
         v = 0.0
         DO k = 1, 3
            v = v + q(i)%xv(k)**2
         END DO
         v = DSQRT(v)
         IF (v < table%v_max) THEN
            x = v/table%dv
            id_radial = INT(x)
            x = x - DBLE(id_radial)
            df = table%ESTOP(id, id_radial)    *(1.D0 - x) + &
                 & table%ESTOP(id, id_radial+1)*x
         ELSE
            CALL ESTOP(param%q(q(i)%id), v, df)
         END IF
         DO k=1,3
            q(i)%ff(k) = q(i)%ff(k) - df*q(i)%xv(k)
         END DO
      END IF
   END DO
END IF
RETURN
END SUBROUTINE Force_local_2
!
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\ Force - for local cell interactions  ////////////////////
SUBROUTINE Force_neighbor(NS, q, grec, fsend, cell, gcell, param, sys, &
     & ion_data, table)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(GP):: grec(NS%Gpt)
TYPE(GF):: fsend(NS%Gpt) 
TYPE(CL):: cell(NS%Nlocal), gcell(NS%Nghost)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(IC):: ion_data
TYPE(VF):: table
!
!
INTEGER :: i, j, k, l, m, n, id, jd, i_id, j_id, id_table, id_radial
REAL(KIND=DP):: r, r2, x, dx(3), dv(3), E, df
!
DO k=1, 3
   fsend(:)%ff(k) = 0.0D0
END DO
fsend(:)%pot = 0.0D0
!
DO n=1, NS%Nlocal
   !
   ! Particle interactions with neighboring cells
   DO i=1, cell(n)%Npt
      id = cell(n)%node(i)
      i_id = q(id)%id
      DO l=1, 19 ! NOT 13 !!!
         m = cell(n)%ngbr(l)
         IF (m <= NS%Nlocal .AND. m /=0) THEN
            !
            ! Local cell - including PBC
            DO j=1, cell(m)%Npt
               jd = cell(m)%node(j)
               j_id = q(jd)%id
               !
               r2 = 0.0D0
               DO k=1,3
                  x = q(id)%xx(k) - q(jd)%xx(k)
                  dx(k) = x - param%box(k)*DNINT(x/param%box(k))
                  r2 = r2 + dx(k)**2
               END DO
               IF (r2 < param%Rc2) THEN
                  r = DSQRT(r2)
                  x = r/table%dr
                  id_radial = INT(x)
                  IF (id_radial < Ntabular) THEN
                     !
                     ! Pair interaction id
                     IF (i_id == 1 .AND. j_id == 1) THEN
                        id_table = 1
                     ELSE IF (i_id == 2 .AND. j_id == 2) THEN
                        id_table = 2
                     ELSE IF ((i_id +j_id) == 3) THEN
                        id_table = 3
                     ELSE IF (i_id == 1 .OR. j_id == 1) THEN
                        STOP "Table index error 4" !id_table = 4
                     ELSE IF (i_id == 2 .OR. j_id == 2) THEN
                        STOP "Table index error 5" !id_table = 5
                     ELSE
                        STOP "Table index error N"
                     END IF
                     !
                     x = x - DBLE(id_radial)
                     df = table%F(id_table, id_radial)  *(1.D0 - x) + &
                        & table%F(id_table, id_radial+1)*x
                     E  = table%V(id_table, id_radial)  *(1.D0 - x) + &
                        & table%V(id_table, id_radial+1)*x
                     DO k = 1,3
                        q(id)%ff(k) = q(id)%ff(k) + df*dx(k)
                        q(jd)%ff(k) = q(jd)%ff(k) - df*dx(k)
                     END DO
                     q(id)%pot = q(id)%pot + E/2.D0
                     q(jd)%pot = q(jd)%pot + E/2.D0
                     sys%Epot = sys%Epot + E
                     ion_data%r_min = MIN(r, ion_data%r_min)
                     !
                     IF (r < Rcut_firsov .AND. DRAG) THEN
                        df = table%FIRSV(id_table, id_radial)  *(1.D0 - x) + &
                           & table%FIRSV(id_table, id_radial+1)*x
                        DO k=1,3
                           dv(k) = q(id)%xv(k) - q(jd)%xv(k)
                           q(id)%ff(k) = q(id)%ff(k) - dv(k)*df
                           q(jd)%ff(k) = q(jd)%ff(k) + dv(k)*df
                        END DO
                     END IF
                     
                     IF (id_table == 4 .OR. id_table == 5) THEN
                        IF (j_id /= 3) THEN 
                           ion_data%id = id
                           j_id = i_id
                        ELSE
                           ion_data%id = jd
                        END IF
                        CALL E_RHO(j_id, r, ion_data)
                        ion_data%tag = .TRUE.
                     END IF
                  END IF
               END IF
            END DO
         ELSE IF ( m/=0 ) THEN
            !
            ! Ghost cell
            m = m - NS%Nlocal
            DO j=1, gcell(m)%Npt
               jd = gcell(m)%node(j)
               j_id = grec(jd)%id
               
               r2 = 0.0D0
               DO k=1,3
                  x = q(id)%xx(k) - grec(jd)%xx(k)
                  dx(k) = x - param%box(k)*DNINT(x/param%box(k))
                  r2 = r2 + dx(k)**2
               END DO
               IF (r2 < param%Rc2) THEN
                  r = DSQRT(r2)
                  x = r/table%dr
                  id_radial = INT(x)
                  IF (id_radial < Ntabular) THEN
                     !
                     ! Pair interaction id
                     IF (i_id == 1 .AND. j_id == 1) THEN
                        id_table = 1
                     ELSE IF (i_id == 2 .AND. j_id == 2) THEN
                        id_table = 2
                     ELSE IF ((i_id +j_id) == 3) THEN
                        id_table = 3
                     ELSE IF (i_id == 1 .OR. j_id == 1) THEN
                        STOP "Table index error 4" !id_table = 4
                     ELSE IF (i_id == 2 .OR. j_id == 2) THEN
                        STOP "Table index error 5" !id_table = 5
                     ELSE
                        STOP "Table index error"
                     END IF
                     !
                     x = x - DBLE(id_radial)
                     df = table%F(id_table, id_radial)  *(1.D0 - x) + &
                        & table%F(id_table, id_radial+1)*x
                     E  = table%V(id_table, id_radial)  *(1.D0 - x) + &
                        & table%V(id_table, id_radial+1)*x
                     DO k = 1,3
                        q(    id)%ff(k) = q(    id)%ff(k) + df*dx(k)
                        fsend(jd)%ff(k) = fsend(jd)%ff(k) - df*dx(k)
                     END DO
                     q(id)%pot = q(id)%pot + E/2.D0
                     fsend(jd)%pot = fsend(jd)%pot + E/2.D0
                     sys%Epot = sys%Epot + E
                     ion_data%r_min = MIN(r, ion_data%r_min)
                     !
                     IF (r < Rcut_firsov .AND. DRAG) THEN
                        df = table%FIRSV(id_table, id_radial)  *(1.D0 - x) + &
                           & table%FIRSV(id_table, id_radial+1)*x
                        DO k=1,3
                           dv(k) = q(id)%xv(k) - grec(jd)%xv(k)
                           q(    id)%ff(k) = q(    id)%ff(k) - dv(k)*df
                           fsend(jd)%ff(k) = fsend(jd)%ff(k) + dv(k)*df
                        END DO
                     END IF
                     !
                     IF (id_table == 4 .OR. id_table == 5) THEN
                        IF (i_id == 3) THEN
                           ion_data%id = id
                           ion_data%tag = .TRUE.
                           CALL E_RHO(j_id, r, ion_data)
                        ELSE
                           CALL E_RHO(i_id, r, ion_data)
                        END IF
                     END IF
                  END IF
               END IF
            END DO
         END IF
         !
      END DO
   END DO
   !
END DO
END SUBROUTINE Force_neighbor

