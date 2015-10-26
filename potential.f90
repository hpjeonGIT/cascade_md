!
! /// Electronic stopping on the implanted ion \\\\\\\\\\
! --- Neighboring charged electron density is accumulated in the loop
! --- of pair interaction 
SUBROUTINE ION_ESTOP(NS, ion_data, z1, q, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(IC):: ion_data
TYPE(PT):: q(NS%Npt)
TYPE(MP):: COMM
INTEGER:: IERR, id, k
REAL(KIND=DP):: rho, rs, v1, v2, v4, v_f, v_f2, v_f4, vr, yr, qeff, lambda
REAL(KIND=DP):: gamma, C, z_eff, z1, dEdx, G, x0
REAL(KIND=DP), PARAMETER:: aB = 0.52917720859D0
REAL(KIND=DP), PARAMETER:: rs0 = 1.2D0, v_f_coeff = 1.9191582926775133D0
REAL(KIND=DP), PARAMETER:: vB = 222.71802408586265D0, a00 = 0.24005195D0
REAL(KIND=DP), PARAMETER:: one_3 = 0.3333333333333333D0
REAL(KIND=DP), PARAMETER:: two_3 = 0.6666666666666666D0
REAL(KIND=DP), PARAMETER:: x0_coeff = 6.0292135933516056D0
REAL(KIND=DP), PARAMETER:: C_force_hart2ev = 51.422086

CALL MPI_ALLREDUCE(ion_data%rho_acc, rho, 1, MPI_REAL8, MPI_SUM, COMM%comm, &
     & IERR)
IF (ion_data%tag) THEN
   id = ion_data%id
   v1 = 0.0D0
   DO k=1,3
      v1 = v1 + q(id)%xv(k)**2
   END DO
   v1 = DSQRT(v1)/vB
   rs = rs0/aB
   v_f = v_f_coeff/rs
   IF (v1 >= v_f) THEN
      vr = v1*(1.D0 + v_f*v_f/5./v1/v1)
   ELSE
      v2 = v1**2
      v4 = v2**2
      v_f2 = v_f**2
      v_f4 = v_f2**2
      vr = 0.75D0*v_f*(1.D0 + (2.D0*v2/(3.D0*v_f2)) - (v4/(15.D0*v_f4)))
   END IF
   yr = vr/z1**(one_3)
   IF (yr < 0.1D0) yr = 0.1D0
   qeff = 1.D0 - DEXP(-0.95D0*(yr-0.07D0))
   lambda = 2.D0*a00*(1.D0-qeff)**two_3/z1**one_3/(1.D0-(1.D0-qeff)/7.D0)
   C = 0.5
   gamma = qeff + C*(1.D0-qeff)*DLOG(1.D0+16.D0*lambda**2/rs**2)
   z_eff = z1*gamma
   IF (rho > 1.D-6) rs = (3.D0/4.D0/pi/rho)**one_3/aB
   x0 = x0_coeff/rs
   dEdx = 2.*v1*(DLOG(1.D0+x0) - x0/(1.D0 + x0))/3.D0/pi
   IF (rs < 6.88505) THEN
      G = 1.0D0 + rs*(0.717D0 - rs*(0.125D0 + rs*(0.0124D0 - rs*0.00212D0)));
   ELSE
      G = 0.72791831295D0;
   END IF
   dEdx = dEdx*z_eff**2*G*C_force_hart2ev/v1/vB
   DO k=1,3
      q(id)%ff(k) = q(id)%ff(k) - dEdx*q(id)%xv(k)
   END DO
END IF

ion_data%tag = .FALSE.
RETURN
END SUBROUTINE ION_ESTOP
!
! /// Variable time step determination \\\
! --- Using the variable time step routine of the REED-MD,
! --- an appropriate time step is calculated by the highest kinetic energy
! --- and the closest pair distance
SUBROUTINE TIME_STEP(tnow, dt, sys, ion_data, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(MP):: COMM
TYPE(ST):: sys
TYPE(IC):: ion_data
INTEGER:: IERR
REAL(KIND=DP):: tnow, dt, dt_tmp, dt_new, dt_min, dt_old, v_max, r_min
REAL(KIND=DP), PARAMETER:: C_diff = 0.1D0
!
CALL MPI_ALLREDUCE(sys%v2_max, v_max, 1, MPI_REAL8, MPI_MAX, COMM%comm, IERR)
CALL MPI_ALLREDUCE(ion_data%r_min,r_min,1,MPI_REAL8, MPI_MIN, COMM%comm, IERR)
v_max = DSQRT(v_max)
dt_tmp = C_diff/v_max
dt_old = dt
dt_new = MIN(dt_old*1.05D0, dt_tmp*0.25D0 + 0.75D0*dt_old)
IF (r_min < 0.5D0) THEN
   dt_min = 0.005D0/v_max
ELSE IF (r_min < 1.0D0) THEN
   dt_min = 0.01D0/v_max
ELSE
   dt_min = 0.02D0/v_max
END IF
dt = MIN(dt_new,dt_min,0.05D0) 
IF (COMM%GID==0) PRINT '("At=", F8.2, " vmax= ", F8.4, " r_min= ", F8.4, &
     & " dt=", F10.7)', tnow*tps, v_max, r_min, dt
ion_data%rho_acc = 0.0D0
ion_data%r_min = 100.D0
ion_data%id = 0
ion_data%tag = .FALSE.
RETURN
END SUBROUTINE TIME_STEP
!
SUBROUTINE E_RHO(id, r, ion_data)
USE DATASTR
IMPLICIT NONE
INTEGER::id, i
REAL(KIND=DP):: r
TYPE(IC):: ion_data
LOGICAL:: loop_tag
!
!
loop_tag = .TRUE.
i=Nline
IF (id == 1 .AND. r < 1.45D0) THEN
   DO WHILE(loop_tag .OR. i == 1)
      IF (ion_data%rho(1,i) < r) loop_tag = .FALSE.
      i = i - 1
   END DO
   ion_data%rho_acc = ion_data%rho_acc + ion_data%rho(2,i)
ELSE IF (id == 2 .AND. r < 0.95D0) THEN
   DO WHILE(loop_tag .OR. i == 1)
      IF (ion_data%rho(3,i) < r) loop_tag = .FALSE.
      i = i - 1
   END DO
   ion_data%rho_acc = ion_data%rho_acc + ion_data%rho(4,i)
END IF
!
RETURN
END SUBROUTINE E_RHO
!
! /// Potential and force tabularization \\\
! --- Short range potential and forces are included such as:
! --- direct sum of Coulomb, ZBL, and the potential from P. Tiwary, A. van
! --- de Walle, and N. Gr\onbech-Jensen
! --- Energy is calculated in eV unit.
! --- Force data is stored for the use as df*dx[k] = eV/\AA. df is the value.
SUBROUTINE VF_TABULAR(param, table)
USE DATASTR
IMPLICIT NONE
TYPE(PM):: param
TYPE(VF):: table
INTEGER:: i
REAL(KIND=DP):: q1, q2, z1, z2, dr, r, df_zbl, df_direct, df_mrln, df_firs
REAL(KIND=DP):: E_zbl, E_direct, E_mrln, df_estop, dv, v
!
!
dr = param%Rc/DBLE(Ntabular)
!
! U-U
q1 = 92.D0
q2 = 92.D0
z1 = param%z(1)
z2 = param%z(1)
DO i=1, Ntabular
   r = DBLE(i)*dr
   CALL COULOMB_direct(z1, z2, r, param%a, df_direct, E_direct)
   CALL V_UU(r, df_mrln, E_mrln)
   table%V(1,i) = E_direct  + E_mrln
   table%F(1,i) = df_direct + df_mrln
   CALL FIRSOV(q1, q2, r, df_firs)
   table%FIRSV(1,i) = df_firs
END DO
!
! O-O
q1 = 8.D0
q2 = 8.D0
z1 = param%z(2)
z2 = param%z(2)
DO i=1, Ntabular
   r = DBLE(i)*dr
   CALL COULOMB_direct(z1, z2, r, param%a, df_direct, E_direct)
   CALL V_OO(r, df_mrln, E_mrln)
   table%V(2,i) = E_direct +  E_mrln
   table%F(2,i) = df_direct + df_mrln
   CALL FIRSOV(q1, q2, r, df_firs)
   table%FIRSV(2,i) = df_firs
END DO
!
! U-O
q1 = 92.D0
q2 =  8.D0
z1 = param%z(1)
z2 = param%z(2)
DO i=1, Ntabular
   r = DBLE(i)*dr
   CALL COULOMB_direct(z1, z2, r, param%a, df_direct, E_direct)
   CALL V_UO(r, df_mrln, E_mrln)
   table%V(3,i) = E_direct +  E_mrln
   table%F(3,i) = df_direct + df_mrln
   CALL FIRSOV(q1, q2, r, df_firs)
   table%FIRSV(3,i) = df_firs
END DO
!
! X-U
q1 = param%q(3)
q2 = param%q(1)
DO i=1, Ntabular
   r = DBLE(i)*dr
   CALL ZBL(q1, q2, q1, q2, r, df_zbl, E_zbl)
   table%V(4,i) = E_zbl
   table%F(4,i) = df_zbl
   CALL FIRSOV(q1, q2, r, df_firs)
   table%FIRSV(4,i) = df_firs
END DO
!
! X-O
q1 = param%q(3)
q2 = param%q(2)
DO i=1, Ntabular
   r = DBLE(i)*dr
   CALL ZBL(q1, q2, q1, q2, r, df_zbl, E_zbl)
   table%V(5,i) = E_zbl
   table%F(5,i) = df_zbl
   CALL FIRSOV(q1, q2, r, df_firs)
   table%FIRSV(5,i) = df_firs
END DO
!
! Electron stopping
table%v_max = v_max_estop
dv = table%v_max/DBLE(Ntabular)
q1 = param%q(1)
q2 = param%q(2)
DO i=1, Ntabular
   v = DBLE(i)*dv
   CALL ESTOP(q1, v, df_estop)
   table%ESTOP(1,i) = df_estop
   CALL ESTOP(q2, v, df_estop)
   table%ESTOP(2,i) = df_estop
END DO
!
table%dr = dr
table%dv = dv
RETURN
END SUBROUTINE VF_TABULAR
!
! /// ZBL nuclear-nuclear collision potential \\\
! This is a hickup for Tiwary potential - screening is using different charge.
SUBROUTINE ZBL(q1, q2, z1, z2, r, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: q1, q2, z1, z2, r, df, E
REAL(KIND=DP):: phi(4), x, a
REAL(KIND=DP), PARAMETER:: aB = 0.52917720859D0
!
a = 0.8854D0*aB/(z1**0.23D0 + z2**0.23D0)
x = r/a
phi(1) = 0.18175D0*DEXP(-3.1998D0*x)
phi(2) = 0.50986D0*DEXP(-0.94229D0*x)
phi(3) = 0.28022D0*DEXP(-0.4029D0*x)
phi(4) = 0.028171D0*DEXP(-0.20162D0*x)

df = ( (x*3.1998D0+1.0D0)*phi(1) + (x*0.94229D0+1.0D0)*phi(2)+ &
     & (x*0.4029D0+1.0D0)*phi(3) + (x*0.20162D0+1.0D0)*phi(4))*eps*q1*q2/r**3
E = eps*q1*q2*(phi(1) + phi(2) + phi(3) + phi(4))/r
RETURN
END SUBROUTINE ZBL
!
! /// Firsov inelastic collision \\\
SUBROUTINE FIRSOV(q1, q2, r, df)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: q1, q2, r, df
REAL(KIND=DP):: z1, z2, alpha
REAL(KIND=DP), PARAMETER:: aB = 0.52917720859D0, a = 0.4685024785301D0, &
     & C = 0.051444896719640566D0

z1 = MAX(q1,q2)
z2 = MIN(q1,q2)
alpha = 1.0D0/(1.D0 + (z2/z1)**(1./6.))

df = C*(z1**2/ (1.D0 + 0.8D0*alpha*z1**(1./3.)*r/a)**4 + &
     &  z2**2/ (1.D0 + 0.8D0*(1.D0-alpha)*z2**(1./3.)*r/a)**4)
RETURN
END SUBROUTINE FIRSOV
!
! /// Morelon et al's potential \\\
! --- Splined by N. Gr\onbech-Jensen
SUBROUTINE V_UU(r, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: r, df, E
REAL(KIND=DP):: r1, r2, q1, df_zbl, E_zbl, E_off, x, Rs, z_u, dRdr
!
r1 = 1.6D0
r2 = 2.1D0
q1 = 92.D0
z_u = 3.227252D0
CALL ZBL(q1, q1, q1, q1, r2, df_zbl, E_zbl)
E_off = E_zbl - eps*z_u*z_u/r2
CALL ZBL(q1, q1, q1, q1, r,  df_zbl, E_zbl)
!
IF (r <r1) THEN
   E = E_zbl - E_off - eps*z_u*z_u/r
   df = df_zbl - eps*z_u*z_u/r**3
ELSE IF (r < r2) THEN
   x = (r - r1)/(r2 - r1)
   Rs = 3.D0*x**2 - 2.D0*x**3
   dRdr = (6.D0*x - 6.D0*x**2)/(r2 - r1)
   E = (1.D0 - Rs)*(E_zbl - E_off - eps*z_u*z_u/r)
   df = (1.D0 - Rs)*(df_zbl - eps*z_u*z_u/r**3)+ &
        & dRdr*(E_zbl - E_off - eps*z_u*z_u/r)/r
ELSE
   E = 0.0
   df = 0.0
END IF
!
RETURN
END SUBROUTINE V_UU
SUBROUTINE V_OO(r, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: r, df, E, dRdr
REAL(KIND=DP):: r1, r2, q1, df_zbl, df_oo, E_zbl, E_off, E_oo, x, Rs, z_u, z_o
!
r1 = 1.0D0
r2 = 1.5D0
q1 = 8.D0
z_u = 3.227252D0
z_o = -z_u/2.D0
CALL ZBL(q1, q1, q1, q1, r2, df_zbl, E_zbl)
CALL FOO(r2, df_oo, E_oo)
E_off = E_zbl - E_oo - eps*z_o*z_o/r2
CALL ZBL(q1, q1, q1, q1, r,  df_zbl, E_zbl)
CALL FOO(r, df_oo, E_oo)
!
IF (r <r1) THEN
   E = E_zbl - E_off - eps*z_o*z_o/r
   df = df_zbl - eps*z_o*z_o/r**3
ELSE IF (r < r2) THEN
   x = (r - r1)/(r2 - r1)
   Rs = 3.D0*x**2 - 2.D0*x**3
   dRdr = (6.D0*x - 6.D0*x**2)/(r2 - r1)
   E = (1.D0 - Rs)*(E_zbl - E_off - eps*z_o*z_o/r) + Rs*E_oo
   df = (1.D0 - Rs)*(df_zbl - eps*z_o*z_o/r**3) + Rs*df_oo + &
        & dRdr*(E_zbl - E_off - eps*z_o*z_o/r)/r - dRdr*E_oo/r
ELSE
   E = E_oo
   df = df_oo
END IF
!
RETURN
END SUBROUTINE V_OO
SUBROUTINE V_UO(r, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: r, df, E, dRdr
REAL(KIND=DP):: r1, r2, q1, q2, df_zbl, df_uo, E_zbl, E_off, E_uo, x, Rs, z_u, z_o
!
r1 = 1.1D0
r2 = 1.6D0
q1 = 8.D0
q2 = 92.D0
z_u = 3.227252D0
z_o = -z_u/2.D0
CALL ZBL(q1, q2, q1, q2, r2, df_zbl, E_zbl)
CALL FUO(r2, df_uo, E_uo)
E_off = E_zbl - E_uo - eps*z_u*z_o/r2
CALL ZBL(q1, q2, q1, q2, r,  df_zbl, E_zbl)
CALL FUO(r, df_uo, E_uo)
!
IF (r <r1) THEN
   E = E_zbl - E_off - eps*z_u*z_o/r
   df = df_zbl - eps*z_u*z_o/r**3
ELSE IF (r < r2) THEN
   x = (r - r1)/(r2 - r1)
   Rs = 3.D0*x**2 - 2.D0*x**3
   dRdr = (6.D0*x - 6.D0*x**2)/(r2 - r1)
   E = (1.D0 - Rs)*(E_zbl - E_off - eps*z_u*z_o/r) + Rs*E_uo
   df = (1.D0 - Rs)*(df_zbl - eps*z_u*z_o/r**3) + Rs*df_uo + &
        & dRdr*(E_zbl - E_off - eps*z_u*z_o/r)/r - dRdr*E_uo/r
ELSE
   E = E_uo
   df = df_uo
END IF
!
RETURN
END SUBROUTINE V_UO
!
! /// Morelon potential \\\
SUBROUTINE FOO(r, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: r, df, E
REAL(KIND=DP):: aoo, boo, coo, c_5(6), c_3(4)
aoo = 11272.6D0
boo = 0.1363D0
coo = 134.0D0
c_5 =RESHAPE((/479.955, -1372.53, 1562.223, -881.968, 246.435, -27.245/),(/6/))
c_3 =RESHAPE((/42.71316, -55.2905, 22.9981, -3.1212/), (/4/))
!
IF (r < 1.2D0) THEN
   E = aoo*DEXP(-r/boo)
   df = E/boo/r
ELSE IF (r < 2.1D0) THEN
   E = c_5(1) + r*(c_5(2) + r*(c_5(3) + r*(c_5(4) + r*(c_5(5) + r*c_5(6)))))
   df = c_5(2)/r + 2.D0*c_5(3) + r*(3.D0*c_5(4)+r*(4.D0*c_5(5)+r*5.D0*c_5(6)))
   df= -df 
ELSE IF (r < 2.6D0) THEN
   E = c_3(1) + r*(c_3(2) + r*(c_3(3) + r*c_3(4)))
   df = c_3(2)/r + 2.D0*c_3(3) + r*3.D0*c_3(4)
   df = -df 
ELSE
   E = -coo/r**6
   df = 6.D0*E/r**2;
END IF
!
RETURN
END SUBROUTINE FOO
!
! /// Morelon potential \\\
SUBROUTINE FUO(r, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: r, df, E
REAL(KIND=DP):: aou, bou
aou = 566.498D0
bou = 0.42056D0
!
E = aou*DEXP(-r/bou)
df = E/bou/r
!
RETURN
END SUBROUTINE FUO
!
! /// direct sum part of Coulomb long-range \\\
! Conventional complementary error function is implemented
SUBROUTINE COULOMB_direct(z1, z2, r, a, df, E)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: z1, z2, r, a, df, E
REAL(KIND=DP):: ar, x, r3, DERF
!
ar = a*r
r3 = r**3
x = 1.D0 - DERF(ar)
df = eps*z1*z2*(2.D0*DEXP(-ar*ar)*ar/sqrtpi + x)/r3
E = eps*z1*z2*x/r
!
RETURN
END SUBROUTINE COULOMB_direct
!
! /// Electronic stopping on the target atom \\\\\\\\\\
! --- the effect of accumulated charge is neglected.
! --- only the plane wave electron drag force
SUBROUTINE ESTOP(z1, v, df_stop)
USE DATASTR
IMPLICIT NONE
REAL(KIND=DP):: z1, v, df_stop
REAL(KIND=DP):: rs, v1, v2, v4, v_f, v_f2, v_f4, vr, yr, qeff, lambda
REAL(KIND=DP):: gamma, C, z_eff, G, x0, dEdx
REAL(KIND=DP), PARAMETER:: aB = 0.52917720859D0
REAL(KIND=DP), PARAMETER:: rs0 = 1.2D0, v_f_coeff = 1.9191582926775133D0
REAL(KIND=DP), PARAMETER:: vB = 222.71802408586265D0, a00 = 0.24005195D0
REAL(KIND=DP), PARAMETER:: one_3 = 0.3333333333333333D0
REAL(KIND=DP), PARAMETER:: two_3 = 0.6666666666666666D0
REAL(KIND=DP), PARAMETER:: x0_coeff = 6.0292135933516056D0
REAL(KIND=DP), PARAMETER:: C_force_hart2ev = 51.422086

v1 = v/vB
rs = rs0/aB
v_f = v_f_coeff/rs
IF (v1 >= v_f) THEN
   vr = v1*(1.D0 + v_f*v_f/5./v1/v1)
ELSE
   v2 = v1**2
   v4 = v2**2
   v_f2 = v_f**2
   v_f4 = v_f2**2
   vr = 0.75D0*v_f*(1.D0 + (2.D0*v2/(3.D0*v_f2)) - (v4/(15.D0*v_f4)))
END IF
yr = vr/z1**(one_3)
IF (yr < 0.1D0) yr = 0.1D0
qeff = 1.D0 - DEXP(-0.95D0*(yr-0.07D0))
lambda = 2.D0*a00*(1.D0-qeff)**two_3/z1**one_3/(1.D0-(1.D0-qeff)/7.D0)
C = 0.5
gamma = qeff + C*(1.D0-qeff)*DLOG(1.D0+16.D0*lambda**2/rs**2)
z_eff = z1*gamma

x0 = x0_coeff/rs
dEdx = 2.*v1*(DLOG(1.D0+x0) - x0/(1.D0 + x0))/3.D0/pi
IF (rs < 6.88505) THEN
   G = 1.0D0 + rs*(0.717D0 - rs*(0.125D0 + rs*(0.0124D0 - rs*0.00212D0)));
ELSE
   G = 0.72791831295D0;
END IF
df_stop = dEdx*z_eff**2*G*C_force_hart2ev/v1/vB
RETURN
END SUBROUTINE ESTOP
