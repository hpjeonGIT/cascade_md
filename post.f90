!
! ################ Intermediate temperature and energy print #################
SUBROUTINE SAMPLING(ttime, sys, param, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(TM):: ttime
TYPE(ST):: sys
TYPE(PM):: param
TYPE(MP):: COMM
!
!
INTEGER:: IERR, i
REAL(KIND=DP):: Ndof(Nparam), mv2(Nparam), Epot, Ekin
!
Ndof(:) = DBLE(param%NTpt(:)*3 - Nref)
!
! Collect kinetic energy
CALL MPI_REDUCE(sys%mv2, mv2, Nparam, MPI_REAL8, MPI_SUM, 0, COMM%comm, IERR)
CALL MPI_REDUCE(sys%Epot,Epot, 1,     MPI_REAL8, MPI_SUM, 0, COMM%comm,IERR)
!
! Print data
Ekin = 0.0
DO i=1,Nparam
   Ekin = Ekin + mv2(i)
END DO
sys%Epot = Epot + sys%Uo + sys%Um
sys%Etot = sys%Epot + Ekin/2.D0
IF (COMM%GID == 0) THEN
!   PRINT *,  sys%Uo , sys%Um , Ur
   WRITE(55, 200) ttime%tnow*tps, mv2(1)/Ndof(1)/kB, mv2(2)/Ndof(2)/kB, &
        & sys%Epot, sys%Etot
END IF
200 FORMAT(5(ES11.4, 1X))
!
RETURN
END SUBROUTINE SAMPLING
!
! ############ restart file generation ###################
SUBROUTINE Restart_serial(NS, q, ttime, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(TM):: ttime
TYPE(MP):: COMM
!
!
TYPE(PT):: local(NS%Npt)
TYPE(PT), ALLOCATABLE:: global(:)
INTEGER:: i, Nfreq, CNT(COMM%ncpu), DPL(COMM%ncpu), IERR
CHARACTER(LEN=20):: FILENAME, FILE2, dummy
!
! File naming
Nfreq = ttime%Nfreq
ttime%Nfreq = ttime%Nfreq + 1
WRITE(FILENAME, 10) Nfreq
WRITE(FILE2   , 20) Nfreq
10 FORMAT("snap",  I3.3, ".xyz") 
20 FORMAT("energy", I3.3, ".dat")
!
! Dump particle information
DO i=1, NS%Npt
   local(i)%xx(:) = q(i)%xx(:)
   local(i)%xv(:) = q(i)%xv(:)
   local(i)%ff(:) = q(i)%ff(:)
   local(i)%pot   = q(i)%pot
   local(i)%id    = q(i)%id
   local(i)%nd    = q(i)%nd
END DO
!
! Number of particles for each processor
CALL MPI_ALLGATHER(NS%Npt, 1, MPI_INTEGER, CNT, 1, MPI_INTEGER, COMM%comm,IERR)
DPL(1) = 0
DO i=1, COMM%ncpu-1
   DPL(i+1) = CNT(i) + DPL(i)
END DO
!
IF(COMM%id == 0) ALLOCATE(global(NS%NTpt))

CALL MPI_GatherV(local, NS%Npt, COMM%NEW, global, CNT, DPL, COMM%NEW, 0, &
     & COMM%comm, IERR)
!
IF (COMM%id == 0) THEN
   !! Sorting
   CALL hpsort(NS%NTpt, global)        

   OPEN (UNIT=30,file=FILENAME)
   WRITE(30,*) NS%NTpt
   WRITE(30,*) "frame = ", Nfreq, " energy= 0 # at ", ttime%tnow*tps, "fs"
   DO i=1, NS%NTpt
      IF (global(i)%id == 1) THEN
         dummy = "U "
      ELSE IF (global(i)%id == 2) THEN
         dummy = "O "
      ELSE 
         dummy = "Kr "
      END IF
      
      WRITE(30,110) dummy, global(i)%xx(1), global(i)%xx(2), global(i)%xx(3), &
           & global(i)%xv(1), global(i)%xv(2), global(i)%xv(3), &
           & global(i)%ff(1), global(i)%ff(2), global(i)%ff(3), global(i)%pot
   
   END DO
   CLOSE(30)
   CLOSE(55)
   OPEN(UNIT=55,file=FILE2)
   WRITE(55,100)
   DEALLOCATE(global)
END IF
100 FORMAT ("#  time (fs)      T_U (K)       T_O (K)       Epot (eV)   &
       &    Etotal (eV)")
110 FORMAT(A3, 3(1X, ES15.8), 7(1X, ES11.4))
RETURN
!
END SUBROUTINE Restart_serial
!!
!! HEAPSORTING routine
SUBROUTINE hpsort(n, ra)
USE DATASTR
IMPLICIT NONE
INTEGER :: n
TYPE(PT):: ra(n), rra
INTEGER :: i, ir, j, l

IF (n < 2) RETURN
l = n/2 + 1
ir = n

10 CONTINUE
IF(l > 1) THEN
   l = l - 1
   rra = ra(l)
ELSE
   rra = ra(ir)
   ra(ir) = ra(1)
   ir = ir-1
   IF (ir == 1) THEN
      ra(1) = rra
      RETURN
   END IF
END IF
i=l
j=l+l
20 IF (j <= ir) THEN
   IF (j < ir) THEN
      IF (ra(j)%nd < ra(j+1)%nd) j = j + 1
   END IF
   IF (rra%nd < ra(j)%nd) THEN
      ra(i) = ra(j)
      i = j
      j =j+j
   ELSE
      j = ir + 1
   END IF
GOTO 20
END IF
ra(i) = rra
GOTO 10
END SUBROUTINE hpsort

