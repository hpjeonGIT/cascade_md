!
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\\\\\\ Cell soring of particles  //////////////////////////
SUBROUTINE CELL_dist(NS, q, cell)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CL):: cell(NS%Nlocal)
!
!
INTEGER :: i, ii, jj, kk, id
REAL(KIND=DP):: x, y, z

cell(:)%Npt = 0

DO i=1, NS%Npt
   x = q(i)%xx(1) - NS%ix(1)
   y = q(i)%xx(2) - NS%ix(2)
   z = q(i)%xx(3) - NS%ix(3)

   ii = INT(x/NS%lc(1))
   IF (ii < 0) ii = 0
   IF (ii >= NS%Nx) ii = ii - 1
   jj = INT(y/NS%lc(2))
   IF (jj < 0) jj = 0
   IF (jj >= NS%Ny) jj = jj - 1
   kk = INT(z/NS%lc(3))
   IF (kk < 0) kk = 0
   IF (kk >= NS%Nz) kk = kk - 1

   id = kk*NS%Nx*NS%Ny + jj*NS%Nx + ii + 1

   cell(id)%Npt = cell(id)%Npt + 1
   cell(id)%node(cell(id)%Npt) = i

END DO
!
RETURN
END SUBROUTINE CELL_dist
!
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\\\ Copy boundary cell into buffer ///////////////////////
SUBROUTINE Copy_buffer(NS, COMM, q, gsend, cell, Nsend)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(MP):: COMM
TYPE(PT):: q(NS%Npt)
TYPE(GP):: gsend(NS%Npt_ghost)
TYPE(CL):: cell(NS%Nlocal)
INTEGER :: Nsend(NS%Nghost)
!
!
INTEGER :: i, n, id, x, y, z, i_npt, i_cell
!
i_npt  = 0
i_cell = 0
IF (COMM%tag(1)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + NS%Nx
   i_cell = i_cell + 1
   Nsend(i_cell) = cell(id)%Npt
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      gsend(i_npt)%xx(:) = q(i)%xx(:)
      gsend(i_npt)%xv(:) = q(i)%xv(:)
      gsend(i_npt)%id    = q(i)%id
   END DO
END IF
IF (COMM%tag(2)) THEN
   DO x = 1,NS%Nx
      id = (NS%Nz-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + x
      i_cell = i_cell + 1
      Nsend(i_cell) = cell(id)%Npt
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         gsend(i_npt)%xx(:) = q(i)%xx(:)
         gsend(i_npt)%xv(:) = q(i)%xv(:)
         gsend(i_npt)%id    = q(i)%id
      END DO
   END DO
END IF
IF (COMM%tag(3)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + 1
   i_cell = i_cell + 1
   Nsend(i_cell) = cell(id)%Npt
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      gsend(i_npt)%xx(:) = q(i)%xx(:)
      gsend(i_npt)%xv(:) = q(i)%xv(:)
      gsend(i_npt)%id    = q(i)%id
   END DO
END IF
IF (COMM%tag(4)) THEN
   DO y = 1,NS%Ny
      id = (NS%Nz-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + NS%Nx
      i_cell = i_cell + 1
      Nsend(i_cell) = cell(id)%Npt
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         gsend(i_npt)%xx(:) = q(i)%xx(:)
         gsend(i_npt)%xv(:) = q(i)%xv(:)
         gsend(i_npt)%id    = q(i)%id
      END DO
   END DO
END IF
IF (COMM%tag(5)) THEN
   DO y = 1,NS%Ny
      DO x= 1, Ns%Nx
         id = (NS%Nz-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + x
         i_cell = i_cell + 1
         Nsend(i_cell) = cell(id)%Npt
         DO n=1, cell(id)%Npt
            i = cell(id)%node(n)
            i_npt = i_npt + 1
            gsend(i_npt)%xx(:) = q(i)%xx(:)
            gsend(i_npt)%xv(:) = q(i)%xv(:)
            gsend(i_npt)%id    = q(i)%id
         END DO
      END DO
   END DO
END IF
IF (COMM%tag(6)) THEN
   DO y = 1,NS%Ny
      id = (NS%Nz-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + 1
      i_cell = i_cell + 1
      Nsend(i_cell) = cell(id)%Npt
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         gsend(i_npt)%xx(:) = q(i)%xx(:)
         gsend(i_npt)%xv(:) = q(i)%xv(:)
         gsend(i_npt)%id    = q(i)%id
      END DO
   END DO
END IF
IF (COMM%tag(7)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + NS%Nx
   i_cell = i_cell + 1
   Nsend(i_cell) = cell(id)%Npt
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      gsend(i_npt)%xx(:) = q(i)%xx(:)
      gsend(i_npt)%xv(:) = q(i)%xv(:)
      gsend(i_npt)%id    = q(i)%id
   END DO
END IF
IF (COMM%tag(8)) THEN
   DO x = 1,NS%Nx
      id = (NS%Nz-1)*NS%Nx*NS%Ny + x
      i_cell = i_cell + 1
      Nsend(i_cell) = cell(id)%Npt
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         gsend(i_npt)%xx(:) = q(i)%xx(:)
         gsend(i_npt)%xv(:) = q(i)%xv(:)
         gsend(i_npt)%id    = q(i)%id
      END DO
   END DO
END IF
IF (COMM%tag(9)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + 1
   i_cell = i_cell + 1
   Nsend(i_cell) = cell(id)%Npt
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      gsend(i_npt)%xx(:) = q(i)%xx(:)
      gsend(i_npt)%xv(:) = q(i)%xv(:)
      gsend(i_npt)%id    = q(i)%id
   END DO
END IF
IF (COMM%tag(10)) THEN
   DO z = 1,NS%Nz
      id = (z-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + NS%Nx
      i_cell = i_cell + 1
      Nsend(i_cell) = cell(id)%Npt
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         gsend(i_npt)%xx(:) = q(i)%xx(:)
         gsend(i_npt)%xv(:) = q(i)%xv(:)
         gsend(i_npt)%id    = q(i)%id
      END DO
   END DO
END IF
IF (COMM%tag(11)) THEN
   DO z = 1,NS%Nz
      DO x= 1, NS%Nx
         id = (z-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + x
         i_cell = i_cell + 1
         Nsend(i_cell) = cell(id)%Npt
         DO n=1, cell(id)%Npt
            i = cell(id)%node(n)
            i_npt = i_npt + 1
            gsend(i_npt)%xx(:) = q(i)%xx(:)
            gsend(i_npt)%xv(:) = q(i)%xv(:)
            gsend(i_npt)%id    = q(i)%id
         END DO
      END DO
   END DO
END IF
IF (COMM%tag(12)) THEN
   DO z = 1,NS%Nz
      id = (z-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + 1
      i_cell = i_cell + 1
      Nsend(i_cell) = cell(id)%Npt
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         gsend(i_npt)%xx(:) = q(i)%xx(:)
         gsend(i_npt)%xv(:) = q(i)%xv(:)
         gsend(i_npt)%id    = q(i)%id
      END DO
   END DO
END IF
IF (COMM%tag(13)) THEN
   DO z = 1,NS%Nz
      DO y= 1, Ns%Ny
         id = (z-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + NS%Nx
         i_cell = i_cell + 1
         Nsend(i_cell) = cell(id)%Npt
         DO n=1, cell(id)%Npt
            i = cell(id)%node(n)
            i_npt = i_npt + 1
            gsend(i_npt)%xx(:) = q(i)%xx(:)
            gsend(i_npt)%xv(:) = q(i)%xv(:)
            gsend(i_npt)%id    = q(i)%id
         END DO
      END DO
   END DO
END IF
!
RETURN
END SUBROUTINE Copy_buffer
!
!##############################################################################
!\\\\\\\\\\\\\ Estimate number of particles to send and receive  //////////////
SUBROUTINE Buffer_estimate(NS, Nsend, Nrec, SBUFF, RBUFF, COMM)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
INTEGER :: Nsend(NS%Nghost), Nrec(NS%Nghost), SBUFF(13), RBUFF(13)
TYPE(MP):: COMM
!
INTEGER:: i_pair, i_cell, x, y, z
! 
i_pair = 0
i_cell = 0
SBUFF = 0
RBUFF = 0
IF (COMM%tag(1)) THEN
   i_pair = i_pair + 1
   i_cell = i_cell + 1
   SBUFF(i_pair) = Nsend(i_cell)
   RBUFF(i_pair) = Nrec(i_cell)
END IF
IF (COMM%tag(2)) THEN
   i_pair = i_pair + 1
   DO x = 1,NS%Nx
      i_cell = i_cell + 1
      SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
      RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
   END DO
END IF
IF (COMM%tag(3)) THEN
   i_pair = i_pair + 1
   i_cell = i_cell + 1
   SBUFF(i_pair) = Nsend(i_cell)
   RBUFF(i_pair) = Nrec(i_cell)
END IF
IF (COMM%tag(4)) THEN
   i_pair = i_pair + 1
   DO y = 1,NS%Ny
      i_cell = i_cell + 1
      SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
      RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
   END DO
END IF
IF (COMM%tag(5)) THEN
   i_pair = i_pair + 1
   DO y = 1,NS%Ny
      DO x= 1, Ns%Nx
         i_cell = i_cell + 1
         SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
         RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
      END DO
   END DO
END IF
IF (COMM%tag(6)) THEN
   i_pair = i_pair + 1
   DO y = 1,NS%Ny
      i_cell = i_cell + 1
      SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
      RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
   END DO
END IF
IF (COMM%tag(7)) THEN
   i_pair = i_pair + 1
   i_cell = i_cell + 1
   SBUFF(i_pair) = Nsend(i_cell)
   RBUFF(i_pair) = Nrec(i_cell)
END IF
IF (COMM%tag(8)) THEN
   i_pair = i_pair + 1
   DO x = 1,NS%Nx
      i_cell = i_cell + 1
      SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
      RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
   END DO
END IF
IF (COMM%tag(9)) THEN
   i_pair = i_pair + 1
   i_cell = i_cell + 1
   SBUFF(i_pair) = Nsend(i_cell)
   RBUFF(i_pair) = Nrec(i_cell)
END IF
IF (COMM%tag(10)) THEN
   i_pair = i_pair + 1
   DO z = 1,NS%Nz
      i_cell = i_cell + 1
      SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
      RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
   END DO
END IF
IF (COMM%tag(11)) THEN
   i_pair = i_pair + 1
   DO z = 1,NS%Nz
      DO x= 1, Ns%Nx
         i_cell = i_cell + 1
         SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
         RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
      END DO
   END DO
END IF
IF (COMM%tag(12)) THEN
   i_pair = i_pair + 1
   DO z = 1,NS%Nz
      i_cell = i_cell + 1
      SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
      RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)      
   END DO
END IF
IF (COMM%tag(13)) THEN
   i_pair = i_pair + 1
   DO z = 1,NS%Nz
      DO y= 1, Ns%Ny
         i_cell = i_cell + 1
         SBUFF(i_pair) = SBUFF(i_pair) + Nsend(i_cell)
         RBUFF(i_pair) = RBUFF(i_pair) + Nrec(i_cell)
      END DO
   END DO
END IF
!
RETURN
END SUBROUTINE Buffer_estimate
!
!##############################################################################
!\\\\\\\\\\\\\\\\ Decode imported buffer and assign forces  /////////////////
SUBROUTINE Decode_buffer(NS, COMM, q, frec, cell)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(MP):: COMM
TYPE(PT):: q(NS%Npt)
TYPE(GF):: frec(NS%Npt_ghost)
TYPE(CL):: cell(NS%Nlocal)
!
!
INTEGER :: i, n, i_npt, id, x, y, z

i_npt = 0
IF (COMM%tag(1)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + NS%Nx
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
      q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
   END DO
END IF
IF (COMM%tag(2)) THEN
   DO x = 1,NS%Nx
      id = (NS%Nz-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + x
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
      q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
      q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
      END DO
   END DO
END IF
IF (COMM%tag(3)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + 1
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
      q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
   END DO
END IF
IF (COMM%tag(4)) THEN
   DO y = 1,NS%Ny
      id = (NS%Nz-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + NS%Nx
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
         q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
      END DO
   END DO
END IF
IF (COMM%tag(5)) THEN
   DO y = 1,NS%Ny
      DO x= 1, Ns%Nx
         id = (NS%Nz-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + x
         DO n=1, cell(id)%Npt
            i = cell(id)%node(n)
            i_npt = i_npt + 1
            q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
            q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
         END DO
      END DO
   END DO
END IF
IF (COMM%tag(6)) THEN
   DO y = 1,NS%Ny
      id = (NS%Nz-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + 1
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
         q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
      END DO
   END DO
END IF
IF (COMM%tag(7)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + NS%Nx
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
      q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
   END DO
END IF
IF (COMM%tag(8)) THEN
   DO x = 1,NS%Nx
      id = (NS%Nz-1)*NS%Nx*NS%Ny + x
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
         q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
      END DO
   END DO
END IF
IF (COMM%tag(9)) THEN
   id = (NS%Nz-1)*NS%Nx*NS%Ny + 1
   DO n=1, cell(id)%Npt
      i = cell(id)%node(n)
      i_npt = i_npt + 1
      q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
      q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
   END DO
END IF
IF (COMM%tag(10)) THEN
   DO z = 1,NS%Nz
      id = (z-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + NS%Nx
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
         q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
      END DO
   END DO
END IF
IF (COMM%tag(11)) THEN
   DO z = 1,NS%Nz
      DO x= 1, Ns%Nx
         id = (z-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + x
         DO n=1, cell(id)%Npt
            i = cell(id)%node(n)
            i_npt = i_npt + 1
            q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
            q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
         END DO
      END DO
   END DO
END IF
IF (COMM%tag(12)) THEN
   DO z = 1,NS%Nz
      id = (z-1)*NS%Nx*NS%Ny + (NS%Ny-1)*NS%Nx + 1
      DO n=1, cell(id)%Npt
         i = cell(id)%node(n)
         i_npt = i_npt + 1
         q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
         q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
      END DO
   END DO
END IF
IF (COMM%tag(13)) THEN
   DO z = 1,NS%Nz
      DO y= 1, Ns%Ny
         id = (z-1)*NS%Nx*NS%Ny + (y-1)*NS%Nx + NS%Nx
         DO n=1, cell(id)%Npt
            i = cell(id)%node(n)
            i_npt = i_npt + 1
            q(i)%ff(:) = q(i)%ff(:) + frec(i_npt)%ff(:)
            q(i)%pot   = q(i)%pot   + frec(i_npt)%pot
         END DO
      END DO
   END DO
END IF
!
RETURN
END SUBROUTINE Decode_buffer
!
!##############################################################################
!\\\\\\\\\\\\\\\\\\\ Ghost ceell soring of ghost particles  ///////////////////
SUBROUTINE gCELL_dist(NS, gcell, Nrec, COMM)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(CL):: gcell(NS%Nghost)
TYPE(MP):: COMM
INTEGER :: Nrec(NS%Nghost)
!
!
INTEGER :: n, x, y, z, i_cell, i_npt
!
i_cell = 0
i_npt = 0

IF (COMM%tag(1)) THEN
   i_cell = i_cell + 1
   gcell(i_cell)%Npt = Nrec(i_cell)
   DO n=1, gcell(i_cell)%Npt
      i_npt = i_npt + 1
      gcell(i_cell)%node(n) = i_npt
   END DO
END IF
IF (COMM%tag(2)) THEN
   DO x = 1,NS%Nx
      i_cell = i_cell + 1
      gcell(i_cell)%Npt = Nrec(i_cell)
      DO n=1, gcell(i_cell)%Npt
         i_npt = i_npt + 1
         gcell(i_cell)%node(n) = i_npt      
      END DO
   END DO
END IF
IF (COMM%tag(3)) THEN
   i_cell = i_cell + 1
   gcell(i_cell)%Npt = Nrec(i_cell)
   DO n=1, gcell(i_cell)%Npt
      i_npt = i_npt + 1
      gcell(i_cell)%node(n) = i_npt
   END DO
END IF
IF (COMM%tag(4)) THEN
   DO y = 1,NS%Ny
      i_cell = i_cell + 1
      gcell(i_cell)%Npt = Nrec(i_cell)
      DO n=1, gcell(i_cell)%Npt
         i_npt = i_npt + 1
         gcell(i_cell)%node(n) = i_npt      
      END DO
   END DO
END IF
IF (COMM%tag(5)) THEN
   DO y = 1,NS%Ny
      DO x= 1, Ns%Nx
         i_cell = i_cell + 1
         gcell(i_cell)%Npt = Nrec(i_cell)
         DO n=1, gcell(i_cell)%Npt
            i_npt = i_npt + 1
            gcell(i_cell)%node(n) = i_npt      
         END DO
      END DO
   END DO
END IF
IF (COMM%tag(6)) THEN
   DO y = 1,NS%Ny
      i_cell = i_cell + 1
      gcell(i_cell)%Npt = Nrec(i_cell)
      DO n=1, gcell(i_cell)%Npt
         i_npt = i_npt + 1
         gcell(i_cell)%node(n) = i_npt      
      END DO
   END DO
END IF
IF (COMM%tag(7)) THEN
   i_cell = i_cell + 1
   gcell(i_cell)%Npt = Nrec(i_cell)
   DO n=1, gcell(i_cell)%Npt
      i_npt = i_npt + 1
      gcell(i_cell)%node(n) = i_npt
   END DO
END IF
IF (COMM%tag(8)) THEN
   DO x = 1,NS%Nx
      i_cell = i_cell + 1
      gcell(i_cell)%Npt = Nrec(i_cell)
      DO n=1, gcell(i_cell)%Npt
         i_npt = i_npt + 1
         gcell(i_cell)%node(n) = i_npt      
      END DO
   END DO
END IF
IF (COMM%tag(9)) THEN
   i_cell = i_cell + 1
   gcell(i_cell)%Npt = Nrec(i_cell)
   DO n=1, gcell(i_cell)%Npt
      i_npt = i_npt + 1
      gcell(i_cell)%node(n) = i_npt
   END DO
END IF
IF (COMM%tag(10)) THEN
   DO z = 1,NS%Nz
      i_cell = i_cell + 1
      gcell(i_cell)%Npt = Nrec(i_cell)
      DO n=1, gcell(i_cell)%Npt
         i_npt = i_npt + 1
         gcell(i_cell)%node(n) = i_npt      
      END DO
   END DO
END IF
IF (COMM%tag(11)) THEN
   DO z = 1,NS%Nz
      DO x= 1, Ns%Nx
         i_cell = i_cell + 1
         gcell(i_cell)%Npt = Nrec(i_cell)
         DO n=1, gcell(i_cell)%Npt
            i_npt = i_npt + 1
            gcell(i_cell)%node(n) = i_npt      
         END DO
      END DO
   END DO
END IF
IF (COMM%tag(12)) THEN
   DO z = 1,NS%Nz
      i_cell = i_cell + 1
      gcell(i_cell)%Npt = Nrec(i_cell)
      DO n=1, gcell(i_cell)%Npt
         i_npt = i_npt + 1
         gcell(i_cell)%node(n) = i_npt      
      END DO
   END DO
END IF
IF (COMM%tag(13)) THEN
   DO z = 1,NS%Nz
      DO y= 1, Ns%Ny
         i_cell = i_cell + 1
         gcell(i_cell)%Npt = Nrec(i_cell)
         DO n=1, gcell(i_cell)%Npt
            i_npt = i_npt + 1
            gcell(i_cell)%node(n) = i_npt      
         END DO
      END DO
   END DO
END IF
NS%Gpt = i_npt
RETURN
END SUBROUTINE gCELL_dist
