!
! ###################### Veolcity Initialization ##############################
!
! With given temperature, initial velocity is given to each particle.
! But distrubtion follows Gaussian(Normalized) random distribution
!
SUBROUTINE Vinit(NS, q, param, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(PM):: param
TYPE(MP):: COMM
!
INTEGER:: i, j, IERR, id
REAL(KIND=DP):: xm, lambda(Nparam), T0(Nparam), xv(NS%Npt,3), &
     & mv2(Nparam), mv2sum(Nparam), sumv(3), rand
!
!
sumv(:) = 0.0D0
DO i = 1, NS%Npt
   DO j=1, 3
      CALL RANDOM_NUMBER(rand)
      xv(i,j) = rand - 0.5D0
      sumv(j) = sumv(j) + xv(i,j)
     END DO
END DO
!
sumv(:) = sumv(:)/DBLE(NS%Npt)
mv2(:) = 0.0D0
DO i=1, NS%Npt
   xv(i,:)= xv(i,:) - sumv(:)
   id = q(i)%id
   xm = param%xm(q(i)%id)   
  mv2(id) = mv2(id) + xm*(xv(i,1)**2 + xv(i,2)**2 + xv(i,3)**2)
END DO
CALL MPI_ALLREDUCE(mv2, mv2sum, Nparam, MPI_REAL8, MPI_SUM, COMM%comm, IERR)
DO i=1, Nparam
   IF (param%NTpt(i) > 0) THEN
      T0(i) = mv2sum(i)/DBLE(param%NTpt(i)*3)
   ELSE
      T0(i) = param%T(i)
   END IF
END DO
lambda(:) = DSQRT(param%T(:)/T0(:))
mv2 = 0.0D0
DO i = 1, NS%Npt
   DO j = 1,3
      id = q(i)%id
      q(i)%xv(j) = xv(i,j)*lambda(id)
      q(i)%ff(j) = 0.0D0
   END DO
END DO
!
RETURN
END SUBROUTINE Vinit
!
!
! ################### Initializing routine for scratch run ####################
!
SUBROUTINE Finit(NS, q, param, sys, COMM, gsend, grec, gpt, fsend, frec, &
     & pme, rec, gff, cell, gcell, PSI, dt, ion_data, table, Narray, Ngpt)

USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(GP):: gsend(NS%Npt_ghost), grec(NS%Npt_ghost)
TYPE(GC):: gpt(NS%Npt_local,4), pme(NS%Npme)
TYPE(GF):: fsend(NS%Npt_ghost), frec(NS%Npt_ghost), rec(NS%Npme), &
     & gff(NS%Npt_local,4)
TYPE(CL):: cell(NS%Nlocal), gcell(NS%Nghost)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(IC):: ion_data
TYPE(VF):: table
REAL(KIND=DP):: PSI(param%K(1), param%K(2), param%K(3))
REAL(KIND=DP):: dt
INTEGER:: Narray(NS%Npt_local), Ngpt(COMM%lncpu)
!
INTEGER:: i, k, BUFF, SBUFF(13), RBUFF(13), SRC, DEST, Send, Sinit, Rend, Rinit
INTEGER:: Ninit, Nend, ireq(26)
INTEGER:: rtag, stag, IERR, ISTAT(MPI_STATUS_SIZE), ISTATUS(MPI_STATUS_SIZE,26)
INTEGER:: Nsend(NS%Nghost), Nrec(NS%Nghost)
REAL(KIND=DP):: x, xm, v2
!
rtag = 0
stag = 0
IF (COMM%TAP) THEN
   DO i = 1, NS%Npt
      xm = param%xm(q(i)%id)
      DO k = 1, 3
         x = q(i)%xx(k) - dt*q(i)%xv(k)
         q(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
      END DO
   END DO
   CALL PME_collect(NS, q, gpt, Narray)
END IF
CALL MPI_BARRIER(COMM%GLOBAL, IERR)
!
CALL Force_SPME_HPC(NS, param, sys, COMM, PSI, gpt, pme, rec, Ngpt)
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
   CALL Force_neighbor(NS, q, grec, fsend, cell, gcell, param, sys, ion_data, table)
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
!
IF (COMM%TAP) THEN
   CALL Decode_SPME(NS, q, COMM, gff, Narray)
   sys%v2_max = 0.0D0
   DO i = 1, NS%Npt
      xm = param%xm(q(i)%id)
      v2 = 0.0D0
      DO k = 1, 3
         x = q(i)%xx(k) + dt*q(i)%xv(k)
         q(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
         v2 = v2 + q(i)%xv(k)**2
      END DO
      sys%v2_max = MAX(v2, sys%v2_max)
   END DO
END IF
!
RETURN
END SUBROUTINE Finit
!
! Verlet routine without thermostat - for microcanonical ensemble (NVE) ######
! ############################################################################
! Also can be used for isokinetic/berendsen thermostat as 1st verlet routine
SUBROUTINE VVerletNotemp1(NS, q, param, dt)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL(KIND=DP):: x, xm
!
DO i = 1, NS%Npt
   xm = param%xm(q(i)%id)
   DO k = 1, 3
      q(i)%xv(k) = q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm
      x = q(i)%xx(k) + dt*q(i)%xv(k)
      q(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
   END DO
END DO
!
RETURN
END SUBROUTINE VVerletNotemp1
!
! ############################################################################
! 2nd routine for Verlet routine - without thermostat #########################
SUBROUTINE VVerletNotemp2(NS, q, param, sys, dt)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, id
REAL(KIND=DP):: xm, v2
!
!
sys%mv2  = 0.0D0
sys%v2_max = 0.0D0
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   v2 = 0.0D0
   DO k = 1, 3
      q(i)%xv(k) = q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm
      v2 = v2 + q(i)%xv(k)**2
   END DO
   sys%mv2(id) = sys%mv2(id) + v2*xm
   sys%v2_max = MAX(v2, sys%v2_max)
END DO
!
RETURN
END SUBROUTINE VVerletNotemp2
!
! ############################################################################
! Isokinetic thermostat - linear velocity scaling thermostat ##############
SUBROUTINE VVerletIsokin2(NS, q, param, sys, dt, COMM)
USE DataStr
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, IERR, id
REAL(KIND=DP):: xm, lambda(Nparam), T0(Nparam), v2, mv2(Nparam), mv2sum(Nparam)
!
!
mv2 = 0.0D0
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)   
   v2 = 0.0D0
   DO k = 1, 3
      q(i)%xv(k) = q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm
      v2 = v2 + q(i)%xv(k)**2
   END DO
   mv2(id) = mv2(id) + xm*v2
END DO
!
CALL MPI_ALLREDUCE(mv2, mv2sum, Nparam, MPI_REAL8, MPI_SUM, COMM%comm, IERR)
DO i=1, Nparam
   IF (param%NTpt(i) > 0) THEN
      T0(i) = mv2sum(i)/DBLE(param%NTpt(i)*3)
   ELSE
      T0(i) = 1.0D0
   END IF
END DO
lambda(:) = DSQRT(param%T(:)/T0(:))
sys%mv2  = 0.0D0
sys%v2_max = 0.0D0
!
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   v2 = 0.0D0
   DO k = 1, 3
      q(i)%xv(k) = q(i)%xv(k)*lambda(id)
      v2 = v2 + q(i)%xv(k)**2
   END DO
   sys%mv2(id) = sys%mv2(id) + v2*xm
   sys%v2_max = MAX(v2, sys%v2_max)
END DO
!
RETURN
END SUBROUTINE VVerletIsokin2
!
! ############################################################################
! Stochastic thermostat - from Green-Kubo dissipation fluctuation theorem ####
SUBROUTINE VVerletStoch1(NS,  q, param, dt)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL(KIND=DP):: x, xm
!
! 
DO i = 1, NS%Npt
   xm = param%xm(q(i)%id)
   DO k = 1, 3
      q(i)%xv(k) = q(i)%xv(k)*(1.D0-0.5D0*param%alpha*dt/xm) + &
           & 0.5D0*dt*q(i)%ff(k)/xm
      x = q(i)%xx(k) + dt*q(i)%xv(k)
      q(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
   END DO
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch1
!
! ############################################################################
! Stochastic thermostat - Second part ########################################
SUBROUTINE VVerletStoch2(NS, q, param, sys, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     USE DATASTR
     IMPLICIT NONE
     REAL(KIND=DP):: fluct, x, r, v1, v2
     REAL(KIND=DP):: rand1, rand2, ran2
     REAL:: ranf
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, id
REAL(KIND=DP):: xm, sigma, beta(Nparam), eta, v2
sigma = 1.D0       ! Deviation for Gaussian distribution
beta(:) = DSQRT(2.D0*param%alpha*param%T(:)/dt) ! Temperature unit is eV
sys%mv2(:) = 0.0D0
sys%v2_max = 0.0D0
!
!
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   v2 = 0.0D0
   DO k = 1, 3
      eta = fluct(sigma)
      q(i)%ff(k) = q(i)%ff(k) + eta*beta(id)
      q(i)%xv(k) = (q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm) / &
           & (1.D0 + 0.5D0*param%alpha*dt/xm)
      v2 = v2 + q(i)%xv(k)**2
   END DO
   sys%mv2(id) = sys%mv2(id) + xm*v2
   sys%v2_max = MAX(v2, sys%v2_max)
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch2
!
! ############################################################################
! Stochastic thermostat - Second part ########################################
SUBROUTINE VVerletSilnt2(NS, q, param, sys, dt, COMM)
USE DataStr
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, IERR, id, Npt(Nparam), Nptsum(Nparam)
REAL(KIND=DP):: xm, lambda(Nparam), T0(Nparam), mv2(Nparam), mv2sum(Nparam)
REAL(KIND=DP):: dx, tv2(Nparam), v2
LOGICAL:: EDGE(NS%Npt)
REAL(KIND=DP), PARAMETER:: r_cut = 3.D0
!
!
tv2(:) = 0.0D0
mv2(:) = 0.0D0
Npt(:) = 0
sys%v2_max = 0.0D0
NS%Nslnt = NS%Nslnt + 1
!
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)   
   v2 = 0.0D0
   EDGE(i) = .FALSE.
   DO k = 1, 3
      dx = 0.5D0*param%box(k) - DABS(q(i)%xx(k))
      IF (DABS(dx) < r_cut) EDGE(i) = .TRUE.
      q(i)%xv(k) = q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm
      v2 = v2 + q(i)%xv(k)**2
   END DO
   IF (EDGE(i)) THEN
      tv2(id) = tv2(id) + v2*xm
      Npt(id) = Npt(id) + 1
   ELSE
      mv2(id) = mv2(id) + v2*xm
      sys%v2_max = MAX(v2, sys%v2_max)
   END IF
END DO
!
IF (MOD(NS%Nslnt,10) == 9) THEN 
   CALL MPI_ALLREDUCE(tv2, mv2sum, Nparam, MPI_REAL8,   MPI_SUM, COMM%comm, IERR)
   CALL MPI_ALLREDUCE(Npt, Nptsum, Nparam, MPI_INTEGER, MPI_SUM, COMM%comm, IERR)
   DO i=1, Nparam
      IF (Nptsum(i) > 0) THEN
         T0(i) = mv2sum(i)/DBLE(Nptsum(i)*3)
      ELSE
         T0(i) = param%T(i)
      END IF
   END DO
   lambda(:) = DSQRT(param%T(:)/T0(:))
   !lambda(:) = DSQRT(1.D0 + 0.1D0*(param%T(:)/T0(:)-1.0D0))
   IF(COMM%GID == 0) PRINT '(3(F5.3, 1X))', lambda(:)
ELSE
   lambda(:) = 1.0D0
END IF
sys%mv2  = 0.0D0
sys%mv2(:) = mv2(:)
!
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   v2 = 0.0D0
   IF (EDGE(i)) THEN
      DO k = 1, 3
         q(i)%xv(k) = q(i)%xv(k)*lambda(id)
         v2 = v2 + q(i)%xv(k)**2
      END DO
      sys%mv2(id) = sys%mv2(id) + v2*xm
      sys%v2_max = MAX(v2, sys%v2_max)
   END IF
END DO
!
!
END SUBROUTINE VVerletSilnt2
!
! ############################################################################
! Stochastic thermostat - Second part ########################################
SUBROUTINE VVerletSilnt3(NS, q, param, sys, dt, COMM)
USE DataStr
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, IERR, id, Npt(Nparam), Nptsum(Nparam)
REAL(KIND=DP):: xm, lambda, T0(Nparam), mv2(Nparam), mv2sum(Nparam)
REAL(KIND=DP):: dx, tv2(Nparam), v2, dx_min(NS%Npt)
LOGICAL:: EDGE(NS%Npt)
REAL(KIND=DP), PARAMETER:: r_cut = 10.D0
!
!
tv2(:) = 0.0D0
mv2(:) = 0.0D0
Npt(:) = 0
sys%v2_max = 0.0D0
!
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)   
   v2 = 0.0D0
   EDGE(i) = .FALSE.
   dx_min(i) = r_cut
   DO k = 1, 3
      dx = 0.5D0*param%box(k) - DABS(q(i)%xx(k))
      dx = DABS(dx)
      IF (dx < r_cut) THEN
         EDGE(i) = .TRUE.
         dx_min(i) = MIN(dx_min(i), dx)
      END IF
      q(i)%xv(k) = q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm
      v2 = v2 + q(i)%xv(k)**2
   END DO
   IF (EDGE(i)) THEN
      tv2(id) = tv2(id) + v2*xm
      Npt(id) = Npt(id) + 1
   ELSE
      mv2(id) = mv2(id) + v2*xm
      sys%v2_max = MAX(v2, sys%v2_max)
   END IF
END DO
!
CALL MPI_ALLREDUCE(tv2, mv2sum, Nparam, MPI_REAL8,   MPI_SUM, COMM%comm, IERR)
CALL MPI_ALLREDUCE(Npt, Nptsum, Nparam, MPI_INTEGER, MPI_SUM, COMM%comm, IERR)
DO i=1, Nparam
   IF (Nptsum(i) > 0) THEN
      T0(i) = mv2sum(i)/DBLE(Nptsum(i)*3)
   ELSE
      T0(i) = param%T(i)
   END IF
END DO
sys%mv2  = 0.0D0
sys%mv2(:) = mv2(:)
!
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   v2 = 0.0D0
   IF (EDGE(i)) THEN
      lambda = DSQRT(1.D0 + (1.D0-dx_min(i)/r_cut)*(param%T(id)/T0(id)-1.D0))
      DO k = 1, 3
         q(i)%xv(k) = q(i)%xv(k)*lambda
         v2 = v2 + q(i)%xv(k)**2
      END DO
      sys%mv2(id) = sys%mv2(id) + v2*xm
      sys%v2_max = MAX(v2, sys%v2_max)
   END IF
END DO
!
!
END SUBROUTINE VVerletSilnt3
!
!
! ############## Random Gaussian(Normal) Distribution Function ################
! 
! For stochastic thermostat, fluctuation dissipation theorem is implemented.
! Basically, random number generator which follows Gaussian distribution is
! needed to implement this thermal noise force.
! Random number is generated using FORTRAN intrinsic fucntion - RANDOM_SEED
! and RANDOM_NUMBER. But during the implementation, it is found that those
! intrinsic functions may not work well under false seeding - to use this
! routine on new machine or new compiler, please make sure that the 
! distribution follows zero-mean with suitable deviation.
!
! This function provides a random number with zero-mean and deviation x 
! along Gaussian distribution.
! <p> = 0.0
! <p**2> = x**2
! 
FUNCTION fluct(x)
  USE DATASTR
  IMPLICIT NONE
  REAL(KIND=DP):: fluct, x, r, v1, v2
  REAL(KIND=DP):: rand1, rand2
  !
  ! Initialization
  r=1.D0
  DO WHILE (r.ge.1.D0)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.D0*rand1 - 1.D0
     v2 = 2.D0*rand2 - 1.D0
     r = v1*v1+v2*v2
  END DO
  fluct = v1*DSQRT(-2.D0*log(r)/r)*x
  RETURN
END FUNCTION fluct
!
! Collect position data for PME ###############################################
SUBROUTINE PME_collect(NS, q, gpt, Narray)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(GC)::gpt(NS%Npt_local,4)
INTEGER ::Narray(NS%Npt_local)
!
! INTERNAL VARIABLES
INTEGER:: i, npt(4), n
!
npt(:) = 0
DO i= 1, NS%Npt
   IF (q(i)%xx(2) < NS%cx(2)) THEN
      IF (q(i)%xx(3) < NS%cx(3)) THEN
         n = 1
         npt(n) = npt(n) + 1
         gpt(npt(n),n)%id = q(i)%id
         gpt(npt(n),n)%xx(:) = q(i)%xx(:)
         Narray(i) = n
      ELSE
         n = 3
         npt(n) = npt(n) + 1
         gpt(npt(n),n)%id = q(i)%id
         gpt(npt(n),n)%xx(:) = q(i)%xx(:)
         Narray(i) = n
      END IF
   ELSE 
      IF (q(i)%xx(3) < NS%cx(3)) THEN
         n = 2
         npt(n) = npt(n) + 1
         gpt(npt(n),n)%id = q(i)%id
         gpt(npt(n),n)%xx(:) = q(i)%xx(:)
         Narray(i) = n
      ELSE
         n = 4
         npt(n) = npt(n) + 1
         gpt(npt(n),n)%id = q(i)%id
         gpt(npt(n),n)%xx(:) = q(i)%xx(:)
         Narray(i) = n
      END IF
   END IF
END DO
NS%Ngpt(:) = npt(:)
!
RETURN
END SUBROUTINE PME_collect

