!
! Redistribute PME force onto each particle ###################################
!
SUBROUTINE SPME_SUM(NS, sys, COMM, rec, gff, Ngpt)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
TYPE(NM):: NS
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(GF):: rec(NS%Npme), gff(NS%Npt_local,4)
INTEGER :: Ngpt(COMM%lncpu)
!
!
INTEGER :: i, IERR, Nacc(COMM%lncpu), nreq(COMM%lncpu), N_cell, &
     & nstat(MPI_STATUS_SIZE, COMM%lncpu)
N_cell = COMM%lncpu - 4
!
! Reciprocal energy
CALL MPI_BCAST(sys%Um, 1, MPI_REAL8, COMM%lncpu-1, COMM%local, IERR)
! Calculated forces
IF (COMM%TAP) THEN   
   DO i=1,4      
      CALL MPI_Irecv(gff(1,i), NS%Ngpt(i), COMM%GFC, N_cell+i-1, 2, &
           & COMM%local, nreq(i), IERR)
   END DO
   CALL MPI_Waitall(4, nreq, nstat, IERR)
ELSE
   Nacc(1) = 1
   DO i=2,N_cell
      Nacc(i) = Nacc(i-1) + Ngpt(i-1)
   END DO
   DO i=1,N_cell
      CALL MPI_Send(rec(Nacc(i)), Ngpt(i), COMM%GFC, i-1, 2, COMM%local, &
           & IERR)
   END DO
END IF
!
RETURN
END SUBROUTINE SPME_SUM
!
! Redistribute PME force onto each particle ###################################
!
SUBROUTINE Decode_SPME(NS, q, COMM, gff, Narray)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(MP):: COMM
TYPE(GF):: gff(NS%Npt_local,4)
INTEGER :: Narray(NS%Npt_local)
!
!
INTEGER :: i, IERR, n, npt(4)
REAL(KIND=DP):: Fx(3), ifx(3)
!
Fx = 0.0D0
npt(:) = 0
DO i=1, NS%Npt
   n = Narray(i)
   npt(n) = npt(n) + 1
   !Fx(:) = Fx(:) + gff(npt(n),n)%ff(:)
   q(i)%ff(:) = q(i)%ff(:) + gff(npt(n),n)%ff(:)
   Fx(:) = Fx(:) + q(i)%ff(:)
   q(i)%pot   = q(i)%pot   + gff(npt(n),n)%pot
END DO
!
! Remove rigid motion
CALL MPI_ALLREDUCE(Fx, ifx, 3, MPI_REAL8, MPI_SUM, COMM%comm, IERR)
!
! Find average rigid force
ifx(:) = ifx(:)/DBLE(NS%NTpt)
!
! Subtraction for all of the particles
DO i=1,NS%Npt
   q(i)%ff(:) = q(i)%ff(:) - ifx(:)
END DO
!
RETURN
END SUBROUTINE Decode_SPME
!
! ################ Routine for SPME HPC version ##############################
! ################ Using FFTE library ########################################
!
SUBROUTINE Force_SPME_HPC(NS, param, sys, COMM, PSI, gpt, pme, rec, Ngpt)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! Charge unit: e = 1.6022E-19 C
! Length unit: A = 1.E-10 m
! Mass unit: amu = 1.6605E-27 kg
! time unit: 1 = 10fs
! Below C constant is 1/4pi e_0. C times 1/r whose has Angstrom unit will be
! energy in eV
!
TYPE(NM):: NS
TYPE(PM):: param
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(GC):: gpt(NS%Npt_local,4), pme(NS%Npme)
TYPE(GF):: rec(NS%Npme)
REAL(KIND=DP):: PSI(param%K(1), param%K(2), param%K(3))
INTEGER :: Ngpt(COMM%lncpu)
!
!
INTEGER :: i, j, k, m, n, ix, iy, iz, IERR, Ny, Nz, fy, fz, SRC, DEST, &
     & ISTATUS(MPI_STATUS_SIZE, Ndim_FFT), stag, rtag, ireq(Ndim_FFT), &
     & ini_y, ini_z, fin_y, fin_z, AllocateStatus, nreq(COMM%lncpu), N_cell, &
     & ini_y_mdiff, ini_z_mdiff, fin_y_pdiff, fin_z_pdiff, & 
     & move_y1, move_y2, move_y3, move_z1, move_z2, move_z3, Nacc(COMM%lncpu),&
     & nstat(MPI_STATUS_SIZE, COMM%lncpu), ISTAT(MPI_STATUS_SIZE)
REAL(KIND=DP):: V, ux, uy, uz, Epme, Esum, z1, z(3), Kx(3), &
     Qlocal, dQdx(3), coeff
REAL(KIND=DP),  ALLOCATABLE::  Qsum(:,:,:), Qrec(:,:), AX(:,:,:), AXX(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE::Q(:,:,:), Qinv(:,:,:), Qfin(:,:,:), Qalg(:,:,:),&
     & Qacc(:,:)
!
! Number of particles to send
N_cell = COMM%lncpu - 4
IF (COMM%TAP) THEN
   DO i=1,4
      CALL MPI_SEND(NS%Ngpt(i), 1, MPI_INTEGER, N_cell+i-1 , 0, COMM%local, &
           & IERR)
   END DO
ELSE
   DO i=1,N_cell
      CALL MPI_RECV(Ngpt(i), 1, MPI_INTEGER, i-1, 0, COMM%local, ISTAT, IERR)
   END DO
END IF
! Send/receive of actual particles
IF (COMM%TAP) THEN
   DO i=1,4
      CALL MPI_Send(gpt(1,i), NS%Ngpt(i), COMM%GCG, N_cell+i-1 , 1, &
           & COMM%local, IERR)
   END DO
ELSE
   Nacc(1) = 1
   DO i=2,4
      Nacc(i) = Nacc(i-1) + Ngpt(i-1)
   END DO
   DO i=1,N_cell
      CALL MPI_Irecv(pme(Nacc(i)), Ngpt(i), COMM%GCG, i-1, 1, COMM%local, &
           & nreq(i), IERR)
   END DO
   CALL MPI_Waitall(N_cell, nreq, nstat, IERR)
END IF
!
IF (.NOT. COMM%TAP) THEN
   !
   ! Initialization and build Q matrix
   NS%Nlpt = 0
   DO i=1, N_cell
      NS%Nlpt = NS%Nlpt + Ngpt(i)
   END DO
   Ny = param%K(2)/COMM%Ncpu_y
   Nz = param%K(3)/COMM%Ncpu_z
   ini_y = COMM%iy; fin_y = COMM%fy;
   ini_z = COMM%iz; fin_z = COMM%fz;
   ini_y_mdiff = ini_y - param%Ndiff
   ini_z_mdiff = ini_z - param%Ndiff
   fin_y_pdiff = fin_y + param%Ndiff
   fin_z_pdiff = fin_z + param%Ndiff
   move_y1 = Ny + param%Ndiff 
   move_y2 = param%K(2) - param%Ndiff 
   move_y3 = ini_y - param%Ndiff - 1
   move_z1 = Nz + param%Ndiff
   move_z2 = param%K(3) - param%Ndiff
   move_z3 = ini_z - param%Ndiff - 1
        !dQdx(4, 4, 4, NS%Nlpt, 3), Qlocal(4,4,4,NS%Nlpt), &
   ALLOCATE(AX(4,3, NS%Nlpt), AXX(4,3, NS%Nlpt), &
        & Qsum(param%K(1), param%Ny_buffer, param%Nz_buffer), &
        & Qrec(param%Nfft_buffer, Ndim_FFT), &
        STAT = AllocateStatus)
   IF (AllocateStatus /=0 ) STOP "=== Not enough memory for FFT ==="
   V = param%box(1)*param%box(2)*param%box(3)
   Qsum = 0.0D0
   DO m=1, NS%Nlpt
      z1 = param%z(pme(m)%id)
      ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5D0)/param%box(1)
      uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5D0)/param%box(2)
      uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5D0)/param%box(3)
      z(1) = DINT(ux) - ux + 1.D0
      z(2) = DINT(uy) - uy + 1.D0
      z(3) = DINT(uz) - uz + 1.D0
      !
      DO i=1,3
         AX(1,i,m)  = z(i)**3/6.D0
         AXX(1,i,m) = -3.D0*z(i)**2/6.D0
         z(i)       = z(i) + 1.D0
         AX(2,i,m)  = (  4.D0 + z(i)*(-12.D0 + z(i)*( 12.D0 - 3.D0*z(i))))/6.D0
         AXX(2,i,m) = ( 12.D0 + z(i)*(-24.D0 + z(i)*9.D0))/6.D0
         z(i)       = z(i) + 1.D0
         AX(3,i,m)  = (-44.D0 + z(i)*( 60.D0 + z(i)*(-24.D0 + 3.D0*z(i))))/6.D0
         AXX(3,i,m) = (-60.D0 + z(i)*( 48.D0 - z(i)*9.D0))/6.D0
         z(i)       = z(i) + 1.D0
         AX(4,i,m)  = ( 64.D0 + z(i)*(-48.D0 + z(i)*( 12.D0 -      z(i))))/6.D0
         AXX(4,i,m) = ( 48.D0 + z(i)*(-24.D0 + z(i)*3.D0))/6.D0
      END DO
      !
      Kx(:) = z1*DBLE(param%K(:))/param%box(:)
      DO i=1,4
         DO j=1,4
            DO k=1,4
               ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
               iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
               iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
               IF (iy < ini_y_mdiff) THEN
                  iy = iy + move_y1
               ELSE IF (iy > fin_y_pdiff) THEN
                  iy = iy - move_y2
               ELSE
                  iy = iy - move_y3
               END IF
               IF (iz < ini_z_mdiff) THEN
                  iz = iz + move_z1
               ELSE IF (iz > fin_z_pdiff) THEN
                  iz = iz - move_z2
               ELSE
                  iz = iz - move_z3
               END IF
               Qsum(ix,iy,iz) = Qsum(ix,iy,iz) + z1*AX(i,1,m)*AX(j,2,m)*AX(k,3,m)
            END DO
         END DO
      END DO
   END DO
   !
   stag = 0
   rtag = 0
   DO i=1,COMM%Npair
      SRC = COMM%FFT_SRC(i)
      DEST = COMM%FFT_DEST(i)
      iy = COMM%OUTgrid_i(i,1); iz = COMM%OUTgrid_i(i,2)
      fy = COMM%OUTgrid_f(i,1); fz = COMM%OUTgrid_f(i,2)
      CALL MPI_Irecv(Qrec(:,i), COMM%FFT_Nset(i), MPI_REAL8, SRC, rtag, &
           & COMM%comm, ireq(i), IERR)
      CALL MPI_Send(Qsum(:,iy:fy,iz:fz), COMM%FFT_Nset(i), MPI_REAL8, &
           & DEST, stag, COMM%comm, IERR)
   END DO
   CALL MPI_Waitall(COMM%Npair, ireq, ISTATUS, IERR)
   !
   DO m=1, COMM%Npair
      iy = COMM%INNgrid_i(m,1); iz = COMM%INNgrid_i(m,2)
      fy = COMM%INNgrid_f(m,1); fz = COMM%INNgrid_f(m,2)
      n = 0
      DO k=iz, fz
         DO j=iy, fy
            DO i=1, param%K(1)
               n = n + 1
               Qsum(i,j,k) = Qsum(i,j,k) + Qrec(n,m)
            END DO
         END DO
      END DO
   END DO
   !
   DEALLOCATE(Qrec)
   ALLOCATE(Q(param%K(1), Ny, Nz), Qinv(param%K(1), Ny, Nz))
   !
   iy = COMM%INNgrid_i(8,1); iz = COMM%INNgrid_i(8,2)
   fy = COMM%INNgrid_f(1,1); fz = COMM%INNgrid_f(1,2)
   Q(:,:,:) = 0.0D0
   Q(:,1:Ny,1:Nz) = DCMPLX(Qsum(:,iy:fy,iz:fz))
   !
   CALL PZFFT3DV(Q, Qinv, param%K(1), param%K(2), param%K(3), &
        & COMM%comy, COMM%comz, COMM%Ncpu_y, COMM%Ncpu_z, 0)
   CALL PZFFT3DV(Q, Qinv, param%K(1), param%K(2), param%K(3), &
        & COMM%comy, COMM%comz, COMM%Ncpu_y, COMM%Ncpu_z, 1)
   !
   ! Convolution
   Kx(:) = DBLE(param%K(:))
   DO i=1, param%K(1)
      z(1) = (4.D0 + 2.D0*DCOS(2.D0*pi*DBLE(i-1)/Kx(1)))**2
      DO j=1, Ny
         z(2) = (4.D0 + 2.D0*DCOS(2.D0*pi*DBLE(j-1+Ny*COMM%id_y)/Kx(2)))**2
         DO k=1, Nz
            z(3) = (4.D0 + 2.D0*DCOS(2.D0*pi*DBLE(k-1+Nz*COMM%id_z)/Kx(3)))**2
            Qinv(i,j,k) = Qinv(i,j,k)*PSI(i,j+Ny*COMM%id_y,k+Nz*COMM%id_z) &
            & *36.D0**3/z(1)/z(2)/z(3)
         END DO
      END DO
   END DO
   !
   DEALLOCATE(Q, Qsum)
   ALLOCATE(Qfin(param%K(1), Ny, Nz))
   !
   CALL PZFFT3DV(Qinv, Qfin, param%K(1), param%K(2), param%K(3), &
        & COMM%comy, COMM%comz, COMM%Ncpu_y, COMM%Ncpu_z,  0)
   CALL PZFFT3DV(Qinv, Qfin, param%K(1), param%K(2), param%K(3), &
        & COMM%comy, COMM%comz, COMM%Ncpu_y, COMM%Ncpu_z, -1)
   !
   DEALLOCATE(Qinv)
   ALLOCATE(Qacc(param%Nfft_buffer, Ndim_FFT), &
        & Qalg(param%K(1), param%Ny_buffer, param%Nz_buffer))
   !
   DO i=1,COMM%Npair
      SRC = COMM%FFT_DEST(i)
      DEST = COMM%FFT_SRC(i)
      iy = COMM%INNgrid_i(i,1) - param%Ndiff
      fy = COMM%INNgrid_f(i,1) - param%Ndiff
      iz = COMM%INNgrid_i(i,2) - param%Ndiff
      fz = COMM%INNgrid_f(i,2) - param%Ndiff
      CALL MPI_Irecv(Qacc(:,i), COMM%FFT_Nset(i), MPI_DOUBLE_COMPLEX, SRC, &
           & rtag, COMM%comm, ireq(i), IERR)
      CALL MPI_Send(Qfin(:,iy:fy,iz:fz), COMM%FFT_Nset(i), &
           & MPI_DOUBLE_COMPLEX, DEST, stag, COMM%comm, IERR)
   END DO
   !   
   Qalg(:,:,:) = 0.D0
   iy = COMM%INNgrid_i(8,1); iz = COMM%INNgrid_i(8,2)
   fy = COMM%INNgrid_f(1,1); fz = COMM%INNgrid_f(1,2)
   Qalg(:,iy:fy,iz:fz) = Qfin(:,1:Ny,1:Nz)
   CALL MPI_Waitall(COMM%Npair, ireq, ISTATUS, IERR)   
   DEALLOCATE(Qfin)
   !
   DO m=1, COMM%Npair
      iy = COMM%OUTgrid_i(m,1); iz = COMM%OUTgrid_i(m,2)
      fy = COMM%OUTgrid_f(m,1); fz = COMM%OUTgrid_f(m,2)
      n = 0
      DO k=iz, fz
         DO j=iy, fy
            DO i=1,param%K(1)
               n = n + 1
               Qalg(i,j,k) = Qalg(i,j,k) + Qacc(n,m)
            END DO
         END DO
      END DO
   END DO
   !
   Epme = 0.0D0
   coeff = eps*DBLE(param%K(1)*param%K(2)*param%K(3))/V/pi
   Kx(:) = DBLE(param%K(:))/param%box(:)
   !
   ! Force estimation
   DO m=1, NS%Nlpt
      ux = Kx(1)*(pme(m)%xx(1)+param%box(1)*.5D0)
      uy = Kx(2)*(pme(m)%xx(2)+param%box(2)*.5D0)
      uz = Kx(3)*(pme(m)%xx(3)+param%box(3)*.5D0)
      rec(m)%ff(:) = 0.0D0
      rec(m)%pot = 0.0D0
      z1 = param%z(pme(m)%id)
      DO i=1,4
         DO j=1,4
            DO k=1,4
               Qlocal  = z1* AX(i,1,m)* AX(j,2,m)* AX(k,3,m)/2.D0
               dQdx(1) = z1*AXX(i,1,m)* AX(j,2,m)* AX(k,3,m)*Kx(1)
               dQdx(2) = z1* AX(i,1,m)*AXX(j,2,m)* AX(k,3,m)*Kx(2)
               dQdx(3) = z1* AX(i,1,m) *AX(j,2,m)*AXX(k,3,m)*Kx(3)
               
               ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
               iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
               iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
               
               IF (iy < ini_y_mdiff) THEN
                  iy = iy + move_y1
               ELSE IF (iy > fin_y_pdiff) THEN
                  iy = iy - move_y2
               ELSE
                  iy = iy  - move_y3
               END IF
               IF (iz < ini_z_mdiff) THEN
                  iz = iz + move_z1
               ELSE IF (iz > fin_z_pdiff) THEN
                  iz = iz - move_z2
               ELSE
                  iz = iz  - move_z3
               END IF
               rec(m)%pot   = rec(m)%pot   + coeff*Qlocal *Qalg(ix,iy,iz)
               rec(m)%ff(1) = rec(m)%ff(1) - coeff*dQdx(1)*Qalg(ix,iy,iz)
               rec(m)%ff(2) = rec(m)%ff(2) - coeff*dQdx(2)*Qalg(ix,iy,iz)
               rec(m)%ff(3) = rec(m)%ff(3) - coeff*dQdx(3)*Qalg(ix,iy,iz)
            END DO
         END DO
      END DO
      Epme = Epme + rec(m)%pot
   END DO
   CALL MPI_ALLREDUCE(Epme, Esum, 1, MPI_REAL8, MPI_SUM, COMM%comm, IERR)
   sys%Um = eps*Esum/V/pi/2.D0
   DEALLOCATE(Qalg, Qacc, AX, AXX)
   !
END IF
!
RETURN
!
END SUBROUTINE FORCE_SPME_HPC
!
! ################## Self-interaction term of EWALD sum #######################
!
SUBROUTINE EWALD_SELF(NS, q, param, sys, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(MP):: COMM
!
! Internal variables
INTEGER:: IERR, i, id
REAL(KIND=DP):: z1, Uo, Uosum
!
!
param%a = 3.1D0/param%Rc
!
!
Uo = 0.0D0
DO i=1, NS%Npt
   id = q(i)%id
   z1 = param%z(id)
   Uo = Uo - z1**2
END DO
CALL MPI_ALLREDUCE(Uo, Uosum, 1, MPI_REAL8, MPI_SUM, COMM%comm,IERR)
sys%Uo = Uosum * param%a/sqrtpi * eps
!
RETURN
END SUBROUTINE EWALD_SELF
!
! ################## Self-interaction term of EWALD sum #######################
!
SUBROUTINE EWALD_PSI(param, PSI)
USE DATASTR
IMPLICIT NONE
TYPE(PM):: param
REAL(KIND=DP):: PSI(param%K(1), param%K(2), param%K(3))
!
! Internal variables
INTEGER:: i, j, k, ix, iy, iz, m
REAL(KIND=DP):: m2
!
!
param%a = 3.1D0/param%Rc
!
!
DO i=1, param%K(1)
   ix = MOD(i-1+param%K(1)/2, param%K(1)) - param%K(1)/2
   DO j=1, param%K(2)
      iy = MOD(j-1+param%K(2)/2, param%K(2)) - param%K(2)/2
      DO k=1, param%K(3)
         iz = MOD(k-1+param%K(3)/2, param%K(3)) - param%K(3)/2
         m = ix**2 + iy**2 + iz**2
         IF( m /= 0) THEN
            m2 = DBLE(ix)**2/param%box(1)**2 + DBLE(iy)**2/param%box(2)**2+ &
                 & DBLE(iz)**2/param%box(3)**2
            PSI(i,j,k) = EXP(-m2*pi**2/param%a**2)/m2
         END IF
      END DO
   END DO
END DO
!
RETURN
END SUBROUTINE EWALD_PSI
!
! ######## Communication setup for volumetric parallel FFT  ##################
!
SUBROUTINE FFT_GRID(param, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(PM):: param
TYPE(MP):: COMM
!
INTEGER:: Ncpu_y, Ncpu_z, id_y, id_z, Ny, Nz, Ky, Kz, iy(8), iz(8)
INTEGER, PARAMETER:: Ndiff = 4
INTEGER:: Ntmp, ICOLOR, N, IERR
!
Ky = param%K(2)
Kz = param%K(3)
Ncpu_y = COMM%icpu(2)*2
Ncpu_z = COMM%icpu(3)*2
IF (Ndiff >= Ky/Ncpu_y .OR. Ndiff >= Kz/Ncpu_z) &
& PRINT *, "==== Size of FFT grids is too small for the given MPI mapping ===="
id_y = MOD(COMM%id, Ncpu_y) ! 0 <= id_y < Ncpu_y
id_z = COMM%id/Ncpu_y       ! 0 <= id_z < Ncpu_z
Ntmp = MAX(param%K(2)/Ncpu_y, param%K(3)/Ncpu_z)
param%Nfft_buffer = param%K(1)*Ntmp*Ndiff
param%Ny_buffer = Ky/Ncpu_y + Ndiff*2
param%Nz_buffer = Kz/Ncpu_z + Ndiff*2
param%Ndiff = Ndiff
!
! | 6 | 7 | 8 |
! | 4 | y | 5 | ~~~f_z
! | 1 | 2 | 3 | ~~~i_z
!     iy  fy
! iy(:) = [1,Ky], 
! for 64 grid with 4 cpus, iy(:) = 1,4, 5,6, 17,20, 21,24
iy(1) = 1
iy(2) = Ndiff 
iy(3) = Ndiff + 1
iy(4) = Ndiff*2
iy(5) = Ky/Ncpu_y + 1
iy(6) = Ky/Ncpu_y + Ndiff
iy(7) = Ky/Ncpu_y + Ndiff + 1
iy(8) = Ky/Ncpu_y + Ndiff*2
iz(1) = 1
iz(2) = Ndiff 
iz(3) = Ndiff + 1
iz(4) = Ndiff*2
iz(5) = Kz/Ncpu_z + 1
iz(6) = Kz/Ncpu_z + Ndiff
iz(7) = Kz/Ncpu_z + Ndiff + 1
iz(8) = Kz/Ncpu_z + Ndiff*2
!
N = 1
Nz = MOD(id_z - 1 + Ncpu_z, Ncpu_z)
Ny = MOD(id_y - 1 + Ncpu_y, Ncpu_y)
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Ndiff
COMM%OUTgrid_i(N,1) = iy(1)
COMM%OUTgrid_f(N,1) = iy(2)
COMM%OUTgrid_i(N,2) = iz(1)
COMM%OUTgrid_f(N,2) = iz(2)
COMM%INNgrid_i(N,1) = iy(5)
COMM%INNgrid_f(N,1) = iy(6)
COMM%INNgrid_i(N,2) = iz(5)
COMM%INNgrid_f(N,2) = iz(6) 
N = 2
Ny = id_y
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Ky/Ncpu_y
COMM%OUTgrid_i(N,1) = iy(3)
COMM%OUTgrid_f(N,1) = iy(6)
COMM%OUTgrid_i(N,2) = iz(1)
COMM%OUTgrid_f(N,2) = iz(2)
COMM%INNgrid_i(N,1) = iy(3)
COMM%INNgrid_f(N,1) = iy(6)
COMM%INNgrid_i(N,2) = iz(5)
COMM%INNgrid_f(N,2) = iz(6) 
N = 3
Ny = MOD(id_y + 1 + Ncpu_y, Ncpu_y)
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Ndiff
COMM%OUTgrid_i(N,1) = iy(7)
COMM%OUTgrid_f(N,1) = iy(8)
COMM%OUTgrid_i(N,2) = iz(1)
COMM%OUTgrid_f(N,2) = iz(2)
COMM%INNgrid_i(N,1) = iy(3)
COMM%INNgrid_f(N,1) = iy(4)
COMM%INNgrid_i(N,2) = iz(5)
COMM%INNgrid_f(N,2) = iz(6) 
N = 4
Nz = id_z
Ny = MOD(id_y - 1 + Ncpu_y, Ncpu_y)
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Kz/Ncpu_z
COMM%OUTgrid_i(N,1) = iy(1)
COMM%OUTgrid_f(N,1) = iy(2)
COMM%OUTgrid_i(N,2) = iz(3)
COMM%OUTgrid_f(N,2) = iz(6)
COMM%INNgrid_i(N,1) = iy(5)
COMM%INNgrid_f(N,1) = iy(6)
COMM%INNgrid_i(N,2) = iz(3)
COMM%INNgrid_f(N,2) = iz(6) 
N = 5
Ny = MOD(id_y + 1 + Ncpu_y, Ncpu_y)
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Kz/Ncpu_z
COMM%OUTgrid_i(N,1) = iy(7)
COMM%OUTgrid_f(N,1) = iy(8)
COMM%OUTgrid_i(N,2) = iz(3)
COMM%OUTgrid_f(N,2) = iz(6)
COMM%INNgrid_i(N,1) = iy(3)
COMM%INNgrid_f(N,1) = iy(4)
COMM%INNgrid_i(N,2) = iz(3)
COMM%INNgrid_f(N,2) = iz(6) 
!
N = 6
Nz = MOD(id_z + 1 + Ncpu_z, Ncpu_z)
Ny = MOD(id_y - 1 + Ncpu_y, Ncpu_y)
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Ndiff
COMM%OUTgrid_i(N,1) = iy(1)
COMM%OUTgrid_f(N,1) = iy(2)
COMM%OUTgrid_i(N,2) = iz(7)
COMM%OUTgrid_f(N,2) = iz(8)
COMM%INNgrid_i(N,1) = iy(5)
COMM%INNgrid_f(N,1) = iy(6)
COMM%INNgrid_i(N,2) = iz(3)
COMM%INNgrid_f(N,2) = iz(4) 
N = 7
Ny = id_y
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Ky/Ncpu_y
COMM%OUTgrid_i(N,1) = iy(3)
COMM%OUTgrid_f(N,1) = iy(6)
COMM%OUTgrid_i(N,2) = iz(7)
COMM%OUTgrid_f(N,2) = iz(8)
COMM%INNgrid_i(N,1) = iy(3)
COMM%INNgrid_f(N,1) = iy(6)
COMM%INNgrid_i(N,2) = iz(3)
COMM%INNgrid_f(N,2) = iz(4) 
N = 8
Ny = MOD(id_y + 1 + Ncpu_y, Ncpu_y)
COMM%FFT_DEST(N) = Nz*Ncpu_y + Ny
COMM%FFT_Nset(N) = param%K(1)*Ndiff*Ndiff
COMM%OUTgrid_i(N,1) = iy(7)
COMM%OUTgrid_f(N,1) = iy(8)
COMM%OUTgrid_i(N,2) = iz(7)
COMM%OUTgrid_f(N,2) = iz(8)
COMM%INNgrid_i(N,1) = iy(3)
COMM%INNgrid_f(N,1) = iy(4)
COMM%INNgrid_i(N,2) = iz(3)
COMM%INNgrid_f(N,2) = iz(4) 
!
COMM%FFT_SRC(1) = COMM%FFT_DEST(8);COMM%FFT_SRC(5) = COMM%FFT_DEST(4);
COMM%FFT_SRC(2) = COMM%FFT_DEST(7);COMM%FFT_SRC(6) = COMM%FFT_DEST(3);
COMM%FFT_SRC(3) = COMM%FFT_DEST(6);COMM%FFT_SRC(7) = COMM%FFT_DEST(2);
COMM%FFT_SRC(4) = COMM%FFT_DEST(5);COMM%FFT_SRC(8) = COMM%FFT_DEST(1);
!
COMM%Npair = Ndim_FFT
COMM%Ncpu_y = Ncpu_y
COMM%Ncpu_z = Ncpu_z
COMM%iy = id_y*Ky/Ncpu_y + 1; COMM%fy = COMM%iy + Ky/Ncpu_y - 1
COMM%iz = id_z*Kz/Ncpu_z + 1; COMM%fz = COMM%iz + Kz/Ncpu_z - 1
ICOLOR = COMM%id/Ncpu_y
CALL MPI_COMM_SPLIT(COMM%comm, ICOLOR, COMM%id, COMM%comy, IERR)
CALL MPI_COMM_RANK(COMM%comy, COMM%id_y, IERR)
ICOLOR = MOD(COMM%id,Ncpu_y)
CALL MPI_COMM_SPLIT(COMM%comm, ICOLOR, COMM%id, COMM%comz, IERR)
CALL MPI_COMM_RANK(COMM%comz, COMM%id_z, IERR)
!
RETURN
END SUBROUTINE FFT_GRID
