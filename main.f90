!
! CASCADE simulation code
! Base code is the TCP analysis code from NAGATO code
! UO2 potential and ZBL components are implemented
!
PROGRAM CASCADE_SIMULATOR
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(PT), POINTER:: q(:)
TYPE(GP), POINTER:: gsend(:), grec(:)
TYPE(GC), POINTER:: gpt(:,:), pme(:)
TYPE(GF), POINTER:: frec(:), fsend(:), rec(:), gff(:,:)
TYPE(CL), POINTER:: cell(:), gcell(:)
REAL(KIND=DP), POINTER:: PSI(:,:,:)
TYPE(TM):: ttime
TYPE(PM):: param
TYPE(ST):: sys
TYPE(NM):: NS
TYPE(MP):: COMM
TYPE(IC):: ion_data
TYPE(VF):: table
REAL(KIND=DP):: dt
!
INTEGER, ALLOCATABLE:: Nsend(:), Nrec(:), Narray(:), Ngpt(:)
INTEGER :: i, j, BUFF, SBUFF(13), RBUFF(13), Nloop_max, Ninit, Nend, &
     & SRC, DEST, stag, rtag, IERR, ISTATUS(MPI_STATUS_SIZE,26), ireq(26), &
     & Sinit, Send, Rinit, Rend, ISTAT(MPI_STATUS_SIZE), id_array(3,3,3)
LOGICAL :: TAG, MD_RUN
REAL   :: time0, time1, time2, secnds
REAL(KIND=DP):: z_time, tmp
!
time1 = secnds(time0)
z_time = MPI_Wtime()
!
!
CALL Parameter_init(NS, param, ttime, COMM, dt, ion_data)
CALL MPI_CELL_init (NS, COMM, cell, gcell, id_array)
!
! Random number initializer
DO i=1, COMM%GID+10
   CALL RANDOM_NUMBER(tmp)
END DO
!
IF (COMM%TAP) THEN
   CALL Position_init (NS, q, param, COMM, ttime)
ELSE
   ALLOCATE(q(1))
END IF
IF (NS%Nghost == 0) THEN
   TAG = .TRUE.
   NS%Nghost = 1
ELSE
   TAG = .FALSE.
END IF
!
CALL MPI_BCAST(ttime%tnow, 1, MPI_REAL8, 0, COMM%GLOBAL, IERR)
!   
NS%Npt_local = NS%Nptmax*NS%Nlocal
Ns%Npt_ghost = NS%Nptmax*NS%Nghost
NS%Nslnt = 0
ion_data%rho_acc = 0.0D0
ion_data%r_min = 2.D0
ion_data%tag = .FALSE.
!
ALLOCATE(Nsend(NS%Nghost), Nrec(NS%Nghost), gsend(NS%Npt_ghost), &
     & grec(NS%Npt_ghost), gpt(NS%Npt_local,4), pme(NS%Npme), &
     & frec(NS%Npt_ghost), fsend(NS%Npt_ghost), rec(NS%Npme), &
     & gff(NS%Npt_local,4), PSI(param%K(1), param%K(2), param%K(3)), &
     & Narray(NS%Npt_local), Ngpt(COMM%lncpu))
IF (TAG) NS%Nghost = 0
IF (COMM%TAP) THEN
   CALL CELL_dist(NS, q, cell)
   CALL EWALD_SELF(NS, q, param, sys, COMM)
   CALL VF_TABULAR(param, table)
ELSE
   CALL EWALD_PSI(param, PSI)
   CALL FFT_GRID(param, COMM)
END IF
!
!
IF (param%rest == 'OFF') THEN
   IF (COMM%TAP) THEN
      CALL Vinit(NS, q, param, COMM)
   END IF
   CALL Finit(NS, q, param, sys, COMM, gsend, grec, gpt, fsend, frec, &
        & pme, rec, gff, cell, gcell, PSI, dt, ion_data, table, Narray, Ngpt)
END IF
!
! Parameter initialization
ttime%Nloop = 0
ttime%tdump = ttime%tnow + ttime%trest(1)
Nloop_max = NINT(ttime%tmax/dt)
sys%v2_max = param%T(1)/param%xm(1)/3.0D0
ion_data%rho_acc = 0.0D0
ion_data%tag = .FALSE.
stag = 0
rtag = 0
MD_RUN = .TRUE.
!
IF (COMM%GID == 0) THEN
   OPEN(UNIT=55, FILE = "energy000.dat")
   WRITE(55, 5) 
END IF
5 FORMAT ("#  time (fs)      T_U (K)       T_O (K)       Epot (eV)   &
       &    Etotal (eV)")
!
!
DO WHILE (MD_RUN)
   !
   ! loop numbering
   ttime%Nloop = ttime%Nloop +  1
   ttime%tnow  = ttime%tnow  + dt
   !
   ! Verlet update 
   ! Migration/cell position check
   ! Collect particle for PME
   IF (COMM%TAP) THEN
      CALL TIME_STEP(ttime%tnow, dt, sys, ion_data, COMM)
      CALL MPI_BCAST(dt, 1, MPI_REAL8, 0, COMM%comm, IERR)
   END IF
   IF (COMM%TAP) THEN
      IF (param%thermo == 'STOCHA')  THEN
         CALL VVerletStoch1(NS,  q, param, dt)
      ELSE
         CALL VVerletNotemp1(NS, q, param, dt)
      END IF
      IF (MOD(ttime%Nloop,ttime%Nsort) == 0) THEN
         CALL MIGRATION(NS, q, COMM, param, id_array)
         CALL CELL_dist(NS, q, cell)
      END IF
      CALL PME_collect(NS, q, gpt, Narray)
   END IF
   !
   CALL Force_SPME_HPC(NS, param, sys, COMM, PSI, gpt, pme, rec, Ngpt)
   !
   ! Force routine
   ! Dump particles onto PME routine
   ! Neighboring cells for direct sums
   ! Sum up direct/reciprocal
   ! Remove rigid motions
   IF (COMM%TAP) THEN      
      CALL Copy_buffer(NS, COMM, q, gsend, cell, Nsend)
      Ninit = 1
      DO i=1, COMM%Npair
         SRC = COMM%SRC(i)
         DEST= COMM%DEST(i)
         BUFF= COMM%Nbuff(i)
         Nend = Ninit + BUFF - 1
         CALL MPI_SENDRECV(Nsend(Ninit:Nend), BUFF, MPI_INTEGER, DEST, stag, &
              & Nrec(Ninit:Nend), BUFF, MPI_INTEGER, SRC, rtag, COMM%comm, &
              & ISTAT, IERR)
         Ninit = Nend + 1
      END DO
      CALL Buffer_estimate(NS, Nsend, Nrec, SBUFF, RBUFF, COMM)
      Sinit = 1
      Rinit = 1
      DO i=1, COMM%Npair
         SRC = COMM%SRC(i)
         DEST= COMM%DEST(i)
         Send = Sinit + SBUFF(i) - 1
         Rend = Rinit + RBUFF(i) - 1
         CALL MPI_Isend(gsend(Sinit:Send), SBUFF(i), COMM%GPT, DEST, stag, &
              & COMM%comm, ireq(2*i-1), IERR)
         CALL MPI_Irecv(grec(Rinit:Rend), RBUFF(i), COMM%GPT, SRC, rtag, &
              & COMM%comm, ireq(2*i), IERR)
         Sinit = Send + 1
         Rinit = Rend + 1
      END DO
      CALL Force_local_1(NS, q, cell, param, sys, ion_data, table)
      CALL MPI_Waitall(COMM%Npair*2, ireq, ISTATUS, IERR)
      CALL gCELL_dist(NS, gcell, Nrec, COMM)
      CALL Force_neighbor(NS, q, grec, fsend, cell, gcell, param, sys, &
           & ion_data,table)
      Sinit = 1
      Rinit = 1
      DO i=1, COMM%Npair
         SRC = COMM%DEST(i)
         DEST= COMM%SRC(i)
         Send = Sinit + RBUFF(i) - 1
         Rend = Rinit + SBUFF(i) - 1
         CALL MPI_Isend(fsend(Sinit:Send), RBUFF(i), COMM%GFC, DEST, stag, &
              & COMM%comm, ireq(i*2-1), IERR)
         CALL MPI_Irecv(frec(Rinit:Rend), SBUFF(i), COMM%GFC, SRC, rtag, &
              & COMM%comm, ireq(i*2), IERR)
         Sinit = Send + 1
         Rinit = Rend + 1
      END DO
      CALL Force_local_2(NS, q, cell, param, sys, ion_data, table)
      CALL MPI_Waitall(COMM%Npair*2, ireq, ISTATUS, IERR)
      CALL Decode_buffer(NS, COMM, q, frec, cell)
   END IF
   CALL SPME_SUM(NS, sys, COMM, rec, gff, Ngpt)
   !
   ! Verlet update
   IF (COMM%TAP) THEN
      CALL Decode_SPME(NS, q, COMM, gff, Narray)
      IF (DRAG) CALL ION_ESTOP(NS, ion_data, param%q(3), q, COMM)
      IF (param%thermo == 'ISOKIN')  THEN
         CALL VVerletIsokin2(NS, q, param, sys, dt, COMM)
      ELSE IF (param%thermo == 'NOTEMP')  THEN      
         CALL VVerletNotemp2(NS, q, param, sys, dt)
      ELSE IF (param%thermo == 'STOCHA')  THEN
         CALL VVerletStoch2(NS, q, param, sys, dt)
      ELSE IF (param%thermo == 'SILENT')  THEN
         CALL VVerletSilnt2(NS, q, param, sys, dt, COMM)
      ELSE
         CALL MPI_FINALIZE(IERR)
         STOP "Verlet Error"
      END IF
      !
      ! Post-processing
      IF (MOD(ttime%Nloop,ttime%Nsamp) == 0) THEN
         CALL SAMPLING(ttime, sys, param, COMM)
      END IF
      IF (ttime%tnow >= ttime%tdump) THEN
         IF (ttime%tnow > ttime%trnge(ttime%Nstage) .AND. &
              & ttime%Nstage < Nstage_max) ttime%Nstage = ttime%Nstage + 1
         CALL Restart_serial(NS, q, ttime, COMM)
         ttime%tdump = ttime%tdump + ttime%trest(ttime%Nstage)
      END IF
   END IF
   !
   IF (COMM%GID == 0 .AND. (ttime%tnow > ttime%tmax)) MD_RUN = .FALSE.
   CALL MPI_BCAST(MD_RUN, 1, MPI_LOGICAL, 0, COMM%GLOBAL, IERR)
END DO
!
! Write final snapshot
!IF (COMM%TAP) CALL Restart_serial(NS, q, ttime, COMM)
!
! Finish
DEALLOCATE(q, gsend, grec, cell, gcell, Nsend, Nrec, gpt, pme, rec, &
     & fsend, frec, gff, PSI, Narray, Ngpt)
time2 = secnds(time1)
IF (COMM%GID == 0) THEN
   CLOSE(55)
   WRITE(*,200) time2
   WRITE(*,220) MPI_Wtime() - z_time
END IF
200 FORMAT ("Wall time is", ES14.4, "sec in secnds function")
220 FORMAT ("Wall time is", ES14.4, "sec in MPI_Wtime function")
CALL MPI_FINALIZE(IERR)
STOP

CONTAINS
!
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\\ Simulation parameter parsing ///////////////////////////
SUBROUTINE Parameter_init(NS, param, ttime, COMM, dt, ion_data)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PM):: param
TYPE(TM):: ttime
TYPE(MP):: COMM
TYPE(IC):: ion_data
REAL(KIND=DP):: dt
!
INTEGER:: i, IERR, OpenStatus
CHARACTER*100:: dummy
OPEN(UNIT=11, FILE="config.prm", status="old", IOSTAT = OpenStatus)
IF (OpenStatus > 0) STOP "==== Cannot open config.prm ===="
READ(11,*)
!
! Time parameter
! 1.0 = 10.18fs
READ(11,*) dummy
READ(11,*) ttime%tmax, dt
ttime%tmax  = ttime%tmax / tps
dt = dt / tps
!
READ(11,*) dummy
DO i=1, Nstage_max
   READ(11,*) ttime%trnge(i), ttime%trest(i)
   ttime%trnge(i) = ttime%trnge(i) / tps
   ttime%trest(i) = ttime%trest(i) / tps
END DO
ttime%Nstage = 1
!
READ(11,*) dummy
READ(11,*) ttime%Nsort, ttime%Nsamp
ttime%Nfreq = 1
!IF (MOD(ttime%Nsort,Nstep_rec) /= 0) STOP "=== Resorting frequency error ==="
!
! Box size
READ(11,*) dummy
READ(11,*) (param%box(j), j=1,3)
!
! Initial temperature
READ(11,*) dummy
DO i=1, Nparam
   READ(11,*) param%T(i)
   param%T(i) = param%T(i)*kB
END DO
!
! Number of all particles
READ(11,*) dummy
NS%NTpt = 0
DO i=1, Nparam
   READ(11,*) param%NTpt(i)
   NS%NTpt = NS%NTpt + param%NTpt(i)
END DO
!
! Mass for each particle kind
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) param%xm(i)
END DO
!
! charge for long range interactions
READ(11,*) dummy
DO i=1, Nparam
   READ(11,*) param%z(i)
END DO
!
! nuclear charge for ZBL interaction
READ(11,*) dummy
DO i=1, Nparam
   READ(11,*) param%q(i)
END DO
!
! Cutoff radius
READ(11,*) dummy
READ(11,*) param%Rc
param%Rc2 = param%Rc**2
!
! Berendsen thermostat/damping constant
READ(11,*) dummy
READ(11,*) param%thermo
READ(11,*) param%alpha
!
! Restart option
READ(11,*) dummy
READ(11,*) param%rest
!
! FFT parameter
READ(11,*) dummy
READ(11,*) param%K(1), param%K(2), param%K(3)
!
! CPU parameter
READ(11,*) dummy
READ(11,*) COMM%icpu(1), COMM%icpu(2), COMM%icpu(3), COMM%icpu(4)
IF (COMM%icpu(2)*COMM%icpu(3)*4 /= COMM%icpu(4)) THEN
   CALL MPI_FINALIZE(IERR)
   CLOSE(11)
   STOP "====== CPU mapping is not working for real/rec. sum ====="
END IF
!
! File close
CLOSE(11)
!
! Electron stopping data parsing
OPEN(UNIT=11, FILE="8_O.rho", status="old", IOSTAT = OpenStatus)
IF (OpenStatus > 0) STOP "==== Cannot open 8_O.rho ===="
READ(11,*) dummy
READ(11,*) dummy
DO i=1, Nline
   READ(11,*) ion_data%rho(1,i), ion_data%rho(2,i)
END DO
CLOSE(11)
OPEN(UNIT=11, FILE="92_U.rho", status="old", IOSTAT = OpenStatus)
IF (OpenStatus > 0) STOP "==== Cannot open 92_U.rho ===="
READ(11,*) dummy
READ(11,*) dummy
DO i=1, Nline
   READ(11,*) ion_data%rho(3,i), ion_data%rho(4,i)
END DO
CLOSE(11)
!
RETURN
END SUBROUTINE Parameter_init
!
!##############################################################################
!\\\\\\\\\\\\\\\\\\\\\\\ Initial position parsing  ////////////////////////////
SUBROUTINE Position_init(NS, q, param, COMM, ttime)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(PT), POINTER:: q(:)
TYPE(NM):: NS
TYPE(PM):: param
TYPE(MP):: COMM
TYPE(TM):: ttime
!
INTEGER :: i, Ntemp, ix, iy, iz, id, OpenStatus, IERR
REAL(KIND=DP):: x, y, z, vx, vy, vz, fx, fy, fz, density, tinit
CHARACTER*256:: dummy, tmp(7)
!
Ntemp = 10*NS%NTpt/COMM%icpu(1)/COMM%icpu(2)/COMM%icpu(3)
ALLOCATE(q(Ntemp))
NS%Npt = 0
IF (param%rest == "OFF") THEN
   OPEN(UNIT=15, FILE="input.xyz", STATUS = "OLD", IOSTAT= OpenStatus)
   IF (OpenStatus > 0) THEN
      CALL MPI_Finalize(IERR)
      STOP " *** CANNOT OPEN input.xyz *** "
   END IF
   READ(15,*) dummy
   READ(15,*) dummy

   DO i=1,NS%NTpt
      READ(15,*) dummy,  x, y, z

      ix = INT((x + 0.5D0*param%box(1))/NS%lx(1))
      IF (ix >= COMM%icpu(1)) ix = ix - 1
      IF (ix < 0) ix = ix + 1
      iy = INT((y + 0.5D0*param%box(2))/NS%lx(2))
      IF (iy >= COMM%icpu(2)) iy = iy - 1
      IF (iy < 0) iy = iy + 1
      iz = INT((z + 0.5D0*param%box(3))/NS%lx(3))
      IF (iz >= COMM%icpu(3)) iz = iz - 1
      IF (iz < 0) iz = iz + 1

      id = iz*COMM%icpu(1)*COMM%icpu(2) + iy*COMM%icpu(1) + ix

      IF (id == COMM%id) THEN
         NS%Npt = NS%Npt + 1
         q(NS%Npt)%xx(1) = x
         q(NS%Npt)%xx(2) = y
         q(NS%Npt)%xx(3) = z
         q(NS%Npt)%nd    = i
         IF (dummy == 'U') THEN 
            q(NS%Npt)%id = 1
         ELSE IF (dummy =='O' ) THEN
            q(NS%Npt)%id = 2
         ELSE
            q(NS%Npt)%id = 3
         END IF
      END IF
   END DO
   CLOSE(15)
   ttime%tnow = 0.0D0
ELSE
   OPEN(UNIT=15, FILE="restart.xyz", STATUS = "OLD", IOSTAT= OpenStatus)
   IF (OpenStatus > 0) THEN
      CALL MPI_Finalize(IERR)
      STOP " *** CANNOT OPEN input.xyz *** "
   END IF
   READ(15,*) dummy
   READ(15,*) tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tinit

   DO i=1,NS%NTpt
      READ(15,*) dummy, x, y, z, vx, vy, vz, fx, fy, fz

      ix = INT((x + 0.5D0*param%box(1))/NS%lx(1))
      IF (ix >= COMM%icpu(1)) ix = ix - 1
      IF (ix < 0) ix = ix + 1
      iy = INT((y + 0.5D0*param%box(2))/NS%lx(2))
      IF (iy >= COMM%icpu(2)) iy = iy - 1
      IF (iy < 0) iy = iy + 1
      iz = INT((z + 0.5D0*param%box(3))/NS%lx(3))
      IF (iz >= COMM%icpu(3)) iz = iz - 1
      IF (iz < 0) iz = iz + 1

      id = iz*COMM%icpu(1)*COMM%icpu(2) + iy*COMM%icpu(1) + ix

      IF (id == COMM%id) THEN
         NS%Npt = NS%Npt + 1
         q(NS%Npt)%xx(1) = x
         q(NS%Npt)%xx(2) = y
         q(NS%Npt)%xx(3) = z
         q(NS%Npt)%xv(1) = vx
         q(NS%Npt)%xv(2) = vy
         q(NS%Npt)%xv(3) = vz
         q(NS%Npt)%ff(1) = fx
         q(NS%Npt)%ff(2) = fy
         q(NS%Npt)%ff(3) = fz
         q(NS%Npt)%nd    = i
         IF (dummy == 'U') THEN 
            q(NS%Npt)%id = 1
         ELSE IF (dummy =='O' ) THEN
            q(NS%Npt)%id = 2
         ELSE
            q(NS%Npt)%id = 3
         END IF
      END IF
   END DO
   CLOSE(15)
   ttime%tnow = tinit/tps
END IF
q => REALLOC_PT(q, NS%Npt)
!   
!
CALL MPI_ALLREDUCE(NS%Npt, Ntemp, 1, MPI_INTEGER, MPI_SUM, COMM%comm, IERR)
IF (NS%NTpt /= Ntemp .AND. COMM%TAP) &
     STOP "======= particle position parsing error ======="
!
!
density = DBLE(NS%NTpt)/param%box(1)/param%box(2)/param%box(3)
IF (COMM%GID == 0) WRITE(*,100) density*param%box(1)*param%box(2)*param%box(3)&
     /DBLE(param%K(1)*param%K(2)*param%K(3))
100 FORMAT("For a single FFT grid box, there are ", ES14.4, "particles")
!
RETURN
END SUBROUTINE Position_init
!
!#############################################################################
!\\\\\\\\\\\\\\ MPI variables and CELL neighbor initialization ///////////////
SUBROUTINE MPI_CELL_init(NS, COMM, cell, gcell, id_array)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(MP):: COMM
TYPE(CL), POINTER:: cell(:), gcell(:)
INTEGER :: id_array(3,3,3)
!
CHARACTER(len=256)  :: hostname
INTEGER, ALLOCATABLE:: cell_id(:,:,:)
LOGICAL:: PERIOD, REORDER
INTEGER:: IERR, MYID, NUMPROC, offset(2), Oldtype(2), blockcnt(2), extent, &
     i, j, k, Nx, Ny, Nz, Nleft, Nright, Ndown, Nup, Nbottom, Ntop, &
     neighbor(26), n, ii, jj, kk, id, i_cell, Ncpu_total, Ncpu_cell, &
     Ncpu_fft, COLORR, KEYY, Cpu_rate, Nfft_y, Nfft_z
REAL(KIND=DP):: density
!
!\\\\\\\\\\\\\\\\\\\\\ MPI parameter configurations ////////////////////////
!
CALL MPI_INIT(IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC, IERR)
Ncpu_total = COMM%icpu(1)*COMM%icpu(2)*COMM%icpu(3) + COMM%icpu(4)
Ncpu_cell  = COMM%icpu(1)*COMM%icpu(2)*COMM%icpu(3)
Ncpu_fft   = COMM%icpu(4)
Cpu_rate = Ncpu_cell / Ncpu_fft + 1
IF (NUMPROC /= Ncpu_total .AND. MYID==0) STOP "====== CPU allocation error"
CALL MPI_GET_PROCESSOR_NAME(hostname, i, IERR)
PERIOD  = .TRUE.
REORDER = .TRUE.
!
! Global communicator
!
CALL MPI_CART_CREATE(MPI_COMM_WORLD, 1, NUMPROC, PERIOD, REORDER, &
     COMM%GLOBAL, IERR)
CALL MPI_COMM_RANK(COMM%GLOBAL, COMM%GID, IERR)
!
!
KEYY = COMM%GID
IF (COMM%GID < Ncpu_cell) THEN
!IF (MOD(COMM%GID,Cpu_rate) < (Cpu_rate-1)) THEN
   COLORR = 0
   COMM%TAP = .TRUE.
ELSE 
   COLORR = 1
   COMM%TAP = .FALSE.
END IF
!
! Cell sorting/FFT communicator
!
CALL MPI_COMM_SPLIT(COMM%GLOBAL, COLORR, KEYY, COMM%comm, IERR)
CALL MPI_COMM_RANK(COMM%comm, COMM%id, IERR)
CALL MPI_COMM_SIZE(COMM%comm, COMM%ncpu, IERR)
!
! Local communicator for collecting FFT particles
!
Ny = COMM%icpu(2)
Nz = COMM%icpu(3)
Nfft_y = COMM%icpu(2)*2
Nfft_z = COMM%icpu(3)*2
DO i=1, 4
   IF (COMM%TAP) THEN
      COLORR = COMM%id/COMM%icpu(1)
   ELSE
      j = MOD(COMM%id, Nfft_y)/2
      k = COMM%id/(Nfft_y*2)
      COLORR = k*Ny + j
   END IF
   CALL MPI_COMM_SPLIT(COMM%GLOBAL, COLORR, KEYY, COMM%local, IERR)
   CALL MPI_COMM_RANK(COMM%local, COMM%lid, IERR)
   CALL MPI_COMM_SIZE(COMM%local, COMM%lncpu, IERR)
END DO
!
! COMM%id = [0, NUMPROC-1]
! i = [0, icpu(:)-1]
Nx = COMM%icpu(1)
Ny = COMM%icpu(2)
Nz = COMM%icpu(3)
i = MOD(COMM%id,Nx)
j = MOD(INT(COMM%id/Nx), Ny)
k = INT(COMM%id/(Nx*Ny))

NS%lx(1:3) = param%box(1:3)/DBLE(COMM%icpu(1:3))
NS%ix(1) = -0.5D0*param%box(1) + DBLE(i)*NS%lx(1)
NS%ix(2) = -0.5D0*param%box(2) + DBLE(j)*NS%lx(2)
NS%ix(3) = -0.5D0*param%box(3) + DBLE(k)*NS%lx(3)
NS%cx(:) = NS%ix(:) + NS%lx(:)*0.5D0

Nleft  = MOD(i-1+Nx,Nx)
Nright = MOD(i+1+Nx,Nx)
Nup    = MOD(j+1+Ny,Ny)
Ndown  = MOD(j-1+Ny,Ny)
Ntop   = MOD(k+1+Nz,Nz)
Nbottom= MOD(k-1+Nz,Nz)

neighbor( 1) = Nleft  + Ndown*Nx + Nbottom*Nx*Ny
neighbor( 2) =     i  + Ndown*Nx + Nbottom*Nx*Ny
neighbor( 3) = Nright + Ndown*Nx + Nbottom*Nx*Ny
neighbor( 4) = Nleft  +     j*Nx + Nbottom*Nx*Ny
neighbor( 5) =     i  +     j*Nx + Nbottom*Nx*Ny
neighbor( 6) = Nright +     j*Nx + Nbottom*Nx*Ny
neighbor( 7) = Nleft  +   Nup*Nx + Nbottom*Nx*Ny
neighbor( 8) =     i  +   Nup*Nx + Nbottom*Nx*Ny
neighbor( 9) = Nright +   Nup*Nx + Nbottom*Nx*Ny
neighbor(10) = Nleft  + Ndown*Nx +       k*Nx*Ny
neighbor(11) =     i  + Ndown*Nx +       k*Nx*Ny
neighbor(12) = Nright + Ndown*Nx +       k*Nx*Ny
neighbor(13) = Nleft  +     j*Nx +       k*Nx*Ny
neighbor(14) = Nright +     j*Nx +       k*Nx*Ny
neighbor(15) = Nleft  +   Nup*Nx +       k*Nx*Ny
neighbor(16) =     i  +   Nup*Nx +       k*Nx*Ny
neighbor(17) = Nright +   Nup*Nx +       k*Nx*Ny
neighbor(18) = Nleft  + Ndown*Nx +    Ntop*Nx*Ny
neighbor(19) =     i  + Ndown*Nx +    Ntop*Nx*Ny
neighbor(20) = Nright + Ndown*Nx +    Ntop*Nx*Ny
neighbor(21) = Nleft  +     j*Nx +    Ntop*Nx*Ny
neighbor(22) =     i  +     j*Nx +    Ntop*Nx*Ny
neighbor(23) = Nright +     j*Nx +    Ntop*Nx*Ny
neighbor(24) = Nleft  +   Nup*Nx +    Ntop*Nx*Ny
neighbor(25) =     i  +   Nup*Nx +    Ntop*Nx*Ny
neighbor(26) = Nright +   Nup*Nx +    Ntop*Nx*Ny

NS%Nx = INT(NS%lx(1)/param%Rc)
NS%Ny = INT(NS%lx(2)/param%Rc)
NS%Nz = INT(NS%lx(3)/param%Rc)
!
! Maximum number of particles in the cell - 5 times margin
density = DBLE(NS%NTpt)/param%box(1)/param%box(2)/param%box(3)
NS%Nptmax = NINT(density*NS%lx(1)*NS%lx(2)*NS%lx(3)/DBLE(NS%Nx*NS%Ny*NS%Nz))*5
IF (NS%Nptmax > NMAX) STOP "=== maximum number of array exceeds the limit ===="

COMM%tag(:) = .FALSE.
COMM%Npair = 0
DO i=1,26
   IF (neighbor(i) /= COMM%id) THEN
      COMM%tag(i)  = .TRUE.
      COMM%Npair   = COMM%Npair + 1
      COMM%SRC(COMM%Npair)  = neighbor(i)
      COMM%DEST(COMM%Npair) = neighbor(27-i)
   END IF
END DO
IF (INT(COMM%Npair/2)*2 /= COMM%Npair) STOP "=== Neighboring CPUs are odd ==="
COMM%Npair = COMM%Npair/2
j = 0
DO i=1, 13
   IF (COMM%tag(i)) THEN
      j = j + 1
      IF (i==1 .OR. i==3 .OR. i==7 .OR. i==9) THEN
         COMM%Nbuff(j) = 1
      ELSE IF (i==2 .OR. i==8) THEN
         COMM%Nbuff(j) = NS%Nx
      ELSE IF (i==4 .OR. i==6) THEN
         COMM%Nbuff(j) = NS%Ny
      ELSE IF (i==5) THEN         
         COMM%Nbuff(j) = NS%Nx*NS%Ny
      ELSE IF (i==10 .OR. i==12) THEN
         COMM%Nbuff(j) = NS%Nz
      ELSE IF (i==11) THEN                  
         COMM%Nbuff(j) = NS%Nx*NS%Nz
      ELSE ! i == 13
         COMM%Nbuff(j) = NS%Ny*NS%Nz
      END IF
   END IF
END DO        
!
!
offset(1) = 0
Oldtype(1) = MPI_INTEGER
blockcnt(1) = 1
CALL MPI_TYPE_EXTENT(MPI_REAL8, extent, IERR)
offset(2) = extent
Oldtype(2) = MPI_REAL8
blockcnt(2) = 6
!
! New data type for ghost particles
CALL MPI_TYPE_STRUCT(2, blockcnt, offset, Oldtype, COMM%GPT, IERR)
CALL MPI_TYPE_COMMIT(COMM%GPT, IERR)
!
! New data type for ghost charges
blockcnt(2) = 3
CALL MPI_TYPE_STRUCT(2, blockcnt, offset, Oldtype, COMM%GCG, IERR)
CALL MPI_TYPE_COMMIT(COMM%GCG, IERR)
!
! New data type for moving particles
blockcnt(1) = 2
blockcnt(2) = 10
CALL MPI_TYPE_STRUCT(2, blockcnt, offset, Oldtype, COMM%NEW, IERR)
CALL MPI_TYPE_COMMIT(COMM%NEW, IERR)
!
! New data type for ghost force only
blockcnt(2) = 4
CALL MPI_TYPE_STRUCT(1, blockcnt(2), offset(1), Oldtype(2), COMM%GFC, IERR)
CALL MPI_TYPE_COMMIT(COMM%GFC, IERR)
!
!\\\\\\\\\\\\\\\\\\\\\ Neighboring cell configurations ////////////////////////
!
NS%lc(1) = NS%lx(1)/DBLE(NS%Nx)
NS%lc(2) = NS%lx(2)/DBLE(NS%Ny)
NS%lc(3) = NS%lx(3)/DBLE(NS%Nz)
!
!
ALLOCATE(cell_id(NS%Nx+2, NS%Ny+2, NS%Nz+2))
cell_id = 0
!
! Assign local cells into whole groups
DO i=1, NS%Nx
   DO j=1, NS%Ny
      DO k=1, NS%Nz
         id =  i + (j-1)*NS%Nx + (k-1)*NS%Nx*NS%Ny
         cell_id(i+1,j+1,k+1) = id
      END DO
   END DO
END DO
!
! Assign ghost cells and PBC cells depending on neighboring CPUs
NS%Nghost = 0
i_cell = NS%Nx*NS%Ny*NS%Nz
i = 0
IF (neighbor(1) /= COMM%id) THEN
   i_cell = i_cell + 1
   cell_id(1,1,1) = i_cell
   NS%Nghost = NS%Nghost + 1
   i = i + 1
ELSE
   cell_id(1,1,1) = cell_id(NS%Nx+1,NS%Ny+1,NS%Nz+1)
   cell_id(NS%Nx+2,NS%Ny+2,NS%Nz+2)  = cell_id(2,2,2)
END IF
IF (neighbor(2) /= COMM%id) THEN
   DO i=1,NS%Nx
      i_cell = i_cell + 1
      cell_id(1+i,1,1) = i_cell
      cell_id(1+i,NS%Ny+2,NS%Nz+2) = cell_id(1+i,2,2)
   END DO
   NS%Nghost = NS%Nghost + NS%Nx
   i = i + 1
ELSE
   DO i=1, NS%Nx
      cell_id(1+i,1,1) = cell_id(1+i,NS%Ny+1,NS%Nz+1)
   END DO
END IF
IF (neighbor(3) /= COMM%id) THEN
   i_cell = i_cell + 1
   cell_id(NS%Nx+2,1,1) = i_cell
   NS%Nghost = NS%Nghost + 1
   i = i + 1
ELSE
   cell_id(NS%Nx+2,1,1) = cell_id(2,NS%Ny+1,NS%Nz+1)
   cell_id(1,NS%Ny+2,NS%Nz+2)  = cell_id(NS%Nx+1,2,2)
END IF
IF (neighbor(4) /= COMM%id) THEN
   DO j=1,NS%Ny
      i_cell = i_cell + 1
      cell_id(1,1+j,1) = i_cell
   END DO
   NS%Nghost = NS%Nghost + NS%Ny
   i = i + 1
ELSE
   DO j=1, NS%Ny
      cell_id(1,1+j,1) = cell_id(NS%Nx+1,1+j,NS%Nz+1)
      cell_id(NS%Nx+2,1+j,NS%Nz+2) = cell_id(2,1+j,2)
   END DO
END IF
IF (neighbor(5) /= COMM%id) THEN
   DO j=1,NS%Ny
      DO i=1,NS%Nx
         i_cell = i_cell + 1
         cell_id(1+i,1+j,1) = i_cell
      END DO
   END DO
   NS%Nghost = NS%Nghost + NS%Nx*NS%Ny
   i = i + 1
ELSE
   DO j=1, NS%Ny
      DO i=1, NS%Nx
         cell_id(1+i,1+j,1) = cell_id(1+i,1+j,NS%Nz+1)
         cell_id(1+i,1+j,NS%Nz+2) = cell_id(1+i,1+j,2)
      END DO
   END DO
END IF
IF (neighbor(6) /= COMM%id) THEN
   DO j=1,NS%Ny
      i_cell = i_cell + 1
      cell_id(NS%Nx+2,1+j,1) = i_cell
   END DO
   NS%Nghost = NS%Nghost + NS%Ny
   i = i + 1
ELSE
   DO j=1, NS%Ny
      cell_id(NS%Nx+2,1+j,1) = cell_id(2,1+j,NS%Nz+1)
      cell_id(1,1+j,NS%Nz+2) = cell_id(NS%Nx+1,1+j,2)
   END DO
END IF
IF (neighbor(7) /= COMM%id) THEN
   i_cell = i_cell + 1
   cell_id(1,NS%Ny+2,1) = i_cell
   NS%Nghost = NS%Nghost + 1
   i = i + 1
ELSE
   cell_id(1,NS%Ny+2,1) = cell_id(NS%Nx+1,2,NS%Nz+1)
   cell_id(NS%Nx+2,1,NS%Nz+2) = cell_id(2,NS%Ny+1,2)
END IF
IF (neighbor(8) /= COMM%id) THEN
   DO i=1,NS%Nx
      i_cell = i_cell + 1
      cell_id(1+i,NS%Ny+2,1) = i_cell
   END DO
   NS%Nghost = NS%Nghost + NS%Nx
   i = i + 1
ELSE
   DO i=1, NS%Nx
      cell_id(1+i,NS%Ny+2,1) = cell_id(1+i,2,NS%Nz+1)
      cell_id(1+i,1,NS%Nz+2) = cell_id(1+i,NS%Ny+1,2)
   END DO
END IF
IF (neighbor(9) /= COMM%id) THEN
   i_cell = i_cell + 1
   cell_id(NS%Nx+2,NS%Ny+2,1) = i_cell
   NS%Nghost = NS%Nghost + 1
   i = i + 1
ELSE
   cell_id(NS%Nx+2,NS%Ny+2,1) = cell_id(2,2,NS%Nz+1)
   cell_id(1,1,NS%Nz+2) = cell_id(NS%Nx+1,NS%Ny+1,2)
END IF
IF (neighbor(10) /= COMM%id) THEN
   DO k=1,NS%Nz
      i_cell = i_cell + 1
      cell_id(1,1,1+k) = i_cell
   END DO
   NS%Nghost = NS%Nghost + NS%Nz
   i = i + 1
ELSE
   DO k=1, NS%Nz
      cell_id(1,1,1+k) = cell_id(NS%Nx+1,NS%Ny+1,1+k)
      cell_id(NS%Nx+2,NS%Ny+2,1+k) = cell_id(2,2,1+k)
   END DO
END IF
IF (neighbor(11) /= COMM%id) THEN
   DO k=1,NS%Nz
      DO i=1,NS%Nx
         i_cell = i_cell + 1
         cell_id(1+i,1,1+k) = i_cell
      END DO
   END DO
   NS%Nghost = NS%Nghost + NS%Nx*NS%Nz
   i = i + 1
ELSE
   DO k=1, NS%Nz
      DO i=1, NS%Nx
         cell_id(1+i,1,1+k) = cell_id(1+i,NS%Ny+1,1+k)
         cell_id(1+i,NS%Ny+2,1+k) = cell_id(1+i,2,1+k)
      END DO
   END DO
END IF
IF (neighbor(12) /= COMM%id) THEN
   DO k=1,NS%Nz
      i_cell = i_cell + 1
      cell_id(NS%Nx+2,1,1+k) = i_cell
   END DO
   NS%Nghost = NS%Nghost + NS%Nz
   i = i + 1
ELSE
   DO k=1, NS%Nz
      cell_id(NS%Nx+2,1,1+k) = cell_id(2,NS%Ny+1,1+k)
      cell_id(1,NS%Ny+2,1+k) = cell_id(NS%Nx+1,2,1+k)
   END DO
END IF
IF (neighbor(13) /= COMM%id) THEN
   DO k=1,NS%Nz
      DO j=1,NS%Ny
         i_cell = i_cell + 1
         cell_id(1,1+j,1+k) = i_cell
      END DO
   END DO
   NS%Nghost = NS%Nghost + NS%Ny*NS%Nz
   i = i + 1
ELSE
   DO k=1, NS%Nz
      DO j=1, NS%Ny
         cell_id(1,1+j,1+k) = cell_id(NS%Nx+1,1+j,1+k)
         cell_id(NS%Nx+2,1+j,1+k) = cell_id(2,1+j,1+k)
      END DO
   END DO
END IF

!
!
NS%Nlocal = NS%Nx*NS%Ny*NS%Nz
NS%Ntotal = NS%Nlocal + NS%Nghost
ALLOCATE(cell(NS%Nlocal), gcell(NS%Nghost))
!
! Allocate neighboring cell id
DO i=1, NS%Nx
   DO j=1, NS%Ny
      DO k=1, NS%Nz
         id = i + (j-1)*NS%Nx + (k-1)*NS%Nx*NS%Ny
         DO n=1,19
            
            SELECT CASE(n)
            CASE(1)
               ii = i -1
               jj = j -1
               kk = k -1
            CASE(2)
               ii = i 
               jj = j -1
               kk = k -1
            CASE(3)
               ii = i +1
               jj = j -1
               kk = k -1
            CASE(4)
               ii = i -1
               jj = j 
               kk = k -1
            CASE(5)
               ii = i 
               jj = j 
               kk = k -1
            CASE(6)
               ii = i +1
               jj = j 
               kk = k -1
            CASE(7)
               ii = i -1
               jj = j +1
               kk = k -1
            CASE(8)
               ii = i 
               jj = j +1
               kk = k -1
            CASE(9)
               ii = i +1
               jj = j +1
               kk = k -1
            CASE(10)
               ii = i -1
               jj = j -1
               kk = k
            CASE(11)
               ii = i
               jj = j -1
               kk = k
            CASE(12)
               ii = i +1
               jj = j -1
               kk = k
            CASE(13)
               ii = i -1
               jj = j 
               kk = k 

            CASE(14)
               ii = i -1
               jj = j +1
               kk = k


            CASE(15)
               ii = i -1
               jj = j -1
               kk = k +1 
            CASE(16)
               ii = i -1
               jj = j
               kk = k +1 
            CASE(17)
               ii = i -1
               jj = j +1 
               kk = k +1
 
            CASE(18)
               ii = i
               jj = j -1 
               kk = k +1 
            CASE(19)
               ii = i +1
               jj = j -1
               kk = k +1

            END SELECT

            i_cell = cell_id(ii+1,jj+1,kk+1)
            IF (n < 14) THEN
               cell(id)%ngbr(n) = i_cell
            ELSE IF ( n > 13 .AND. i_cell > NS%Nlocal) THEN
               cell(id)%ngbr(n) = i_cell
            ELSE
               cell(id)%ngbr(n) = 0
            END IF
            
         END DO
      END DO
   END DO
END DO
!
!
DEALLOCATE(cell_id)
!

IF (COMM%GID == 0) WRITE(*,100) NS%Nx, NS%Ny, NS%Nz, &
     density*NS%lc(1)*NS%lc(2)*NS%lc(3)
100 FORMAT("Each processor will have ",I5,"x",I5,"x",I5 "cells with", ES14.4, &
         & "particles per a cell")

NS%Npme = NS%Nptmax*NS%Nlocal*COMM%icpu(1)*COMM%icpu(2)*COMM%icpu(3)/COMM%icpu(4)
IF (COMM%GID >= Ncpu_cell) THEN
   NS%Nghost = 1
   NS%Npt = 1
END IF
!
! Indexing neighboring CPUs
id_array = 0
id = 0
DO n=1,13
   IF (COMM%tag(n)) THEN

      i = MOD(n-1,3) + 1
      j = MOD((n-1)/3,3) + 1
      k = (n-1)/9 + 1

      id = id + 1
      id_array(i,j,k) = id
   END IF
END DO
DO n=14,26
   IF (COMM%tag(n)) THEN

      i = MOD(n,3) + 1
      j = MOD(n/3,3) + 1
      k = n/9 + 1

      id = id + 1
      id_array(i,j,k) = id
   END IF
END DO
!
RETURN
END SUBROUTINE MPI_CELL_init
!
! Random number seeding
SUBROUTINE INIT_RANDOM_SEED(NID)
  INTEGER :: NID, i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  clock = NID
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!
END PROGRAM CASCADE_SIMULATOR
