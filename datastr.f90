MODULE DATASTR
IMPLICIT NONE
!
! PRECISION
! SI = single precision for integer
! DI = double precision for integer
! SP = single precision for real
! DP = double precision for real
!INTEGER, PARAMETER:: SI = SELECTED_INT_KIND(7),  DI = SELECTED_INT_KIND(15)
!INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(7), DP = SELECTED_REAL_KIND(14)
INTEGER, PARAMETER:: SI = 4, DI = 8
INTEGER, PARAMETER:: SP = 4, DP = 8
INTEGER, PARAMETER:: Nparam=3, Nref=3, NMAX=10000, Ntabular = 5000
INTEGER, PARAMETER:: Nline = 80
INTEGER, PARAMETER:: Ndim_FFT = 8 ! 8 neighboring communication
INTEGER, PARAMETER:: Nstage_max = 3 ! Range for different restart frequency
REAL(KIND=DP), PARAMETER:: eps = 14.399644154020082D0 
REAL(KIND=DP), PARAMETER:: tps = 10.180505306645389D0
REAL(KIND=DP), PARAMETER:: pi  = 3.1415926535897931D0
REAL(KIND=DP), PARAMETER:: sqrtpi = 1.7724538509055159D0
REAL(KIND=DP), PARAMETER:: kB = 8.617343D-5
REAL(KIND=DP), PARAMETER:: Rcut_firsov = 3.0D0 ! cutoff radius for Firsov
REAL(KIND=DP), PARAMETER:: v_max_estop = 1.D3  ! max data for tabular E_stop 
LOGICAL, PARAMETER:: DRAG = .FALSE.
!
! V/F tabular data
TYPE VF
   ! 1 = U-U
   ! 2 = O-O
   ! 3 = U-O
   ! 4 = X-U
   ! 5 = X-O
   REAL(KIND=DP):: V(5,Ntabular), F(5,Ntabular), dr, dv
   REAL(KIND=DP):: FIRSV(5, Ntabular), ESTOP(3, Ntabular), v_max
END TYPE VF
!
! Ion Collision data set
TYPE IC
   REAL(KIND=DP):: rho(4,80), r_min, rho_acc
   INTEGER:: id
   LOGICAL:: tag
END TYPE IC
!
! particle data type
TYPE PT
   !
   ! id = kind of particle
   ! xx = position
   ! xv = velocity
   ! ff = force
   INTEGER:: id, nd
   REAL(KIND=DP):: xx(3), xv(3), ff(3), pot
END TYPE PT
!
! ghost particle data type (for cell-sorting type direct sum)
TYPE GP
   !
   ! id = kind of particle
   ! xx = position (for communication)
   INTEGER:: id
   REAL(KIND=DP):: xx(3), xv(3)
END TYPE GP
!
! ghost charge data type (for reciprocal sum only)
TYPE GC
   !
   ! id = kind of particle
   ! xx = position (for communication)
   INTEGER:: id
   REAL(KIND=DP):: xx(3)
END TYPE GC
!
! ghost force data type (for data feed-back)
TYPE GF
   REAL(KIND=DP):: ff(3), pot
END TYPE GF
!
! cell data type
TYPE CL
   !
   ! Npt  = number of particle in the cell
   ! node = corresponding particle number
   ! ngbr = neighbor cell list
   !      = basically 13 neighbors are needed, but ghost cell treatment 
   !        modifies search engine and demands max 19.
   INTEGER:: Npt, node(NMAX), ngbr(19)
END TYPE CL
!
! time data
TYPE TM
   !
   ! Nloop = current loop number
   ! Nrest = frequency of restart file generation
   ! Ndump = frequency of RDF file generation
   ! Nsamp = frequency of temporal data generation
   ! Nsort = frequency of cell sorting
   ! tmax  = time of maximum time steps
   ! trest = time interval of restart
   ! tdist = time interval of cell distribution
   ! tsamp = time interval of sampling temporal data
   ! tnow  = current time
   INTEGER:: Nloop, Nrest, Ndump, Nsort, Nsamp, Nfreq, Nstage
   REAL(KIND=DP):: tnow,  trest(Nstage_max), trnge(Nstage_max), &
        & tdump, tdist, tsamp, tmax
END TYPE TM
!
! parameter data
TYPE PM
   !
   ! xm = mass of particle
   ! box = simulation PBC box size
   ! T = given temperature
   ! rest = restart option (ON/OFF)
   INTEGER:: NTpt(Nparam), K(3), Nfft_buffer, Ny_buffer, Nz_buffer, Ndiff
   REAL(KIND=DP) :: xm(Nparam), box(3), T(Nparam), q(Nparam), z(Nparam), &
        & Rc, Rc2, alpha, a, rs, dr
   CHARACTER*3:: rest
   CHARACTER*6:: thermo
END TYPE PM
!
! system data
TYPE ST
   !
   ! temp = temperature
   ! mv2  = twice of kinetic energy
   ! Epot = potential energy
   ! Etot = total energy
   REAL(KIND=DP):: temp(Nparam), mv2(Nparam), Epot, Etot, v2_max, Uo, Um
END TYPE ST
!
! number data
TYPE NM
   !
   ! Npt = number of particles per processor
   ! NTpt = total number of particles of system
   ! Nx/y/z = number of cells along x/y/z
   ! Nlocal = Number of local cells = Nx*Ny*Nz
   ! Nghost = number of ghost cells
   ! Gpt    = number of ghost particles
   ! Nlpt   = number of particles for PME
   ! lc = length of cell box
   ! ix = left/down/bottm position of each CPU domain
   ! cx = center of each CPU domain
   ! lx = length of each CPU domain
   INTEGER:: Npt, NTpt, Nx, Ny, Nz, Nlocal, Nghost, Ntotal, Gpt, Npme, &
        & Nptmax, Nlpt, Npt_local, Npt_ghost, Ngpt(4), Nslnt
   REAL(KIND=DP) :: lc(3), ix(3), lx(3), dx_ee, dx_ei, dx_ii, cx(3)
END TYPE NM
!
! MPI data
TYPE MP
   !
   ! id = id of each processor
   ! comm = communicator
   ! icpu = structure of allocated processors
   ! Npair = Number of communication pairs for neighboring cells
   ! SRC/DEST = source/destiny node of communications for neighboring cells
   !          = [0, Npair*2]
   ! tag = Tag for neighboring CPUs = [0, 26]
   INTEGER:: id, comm, icpu(4), SRC(26), DEST(26), Npair, Nbuff(13), &
        & NEW, GPT, GCG, GFC, GID, GLOBAL, lid, local, lncpu, ncpu
   LOGICAL:: tag(26), TAP
   !
   ! FFT_SRC/DEST = Communicating nodes for FFT communication
   ! OUTG_i(i,x) = initial mesh (grid) id to send at x-coordinate 
   ! OUTG_f(i,x) = final mesh (grid) id to send at x-coodrdinate
   ! ING_i(i,x) = initial mesh (grid) id to be updated at x-coordinate
   ! ING_f(i,x) = final mesh (grid) id to be updated at x-coodrdinate
   ! SELF_i/f(:) = Mesh index for self at x/y-coordinate
   INTEGER:: FFT_SRC(8), FFT_DEST(8), FFT_Nset(8), id_y, id_z, comy, comz, &
        & OUTgrid_i(8,2), OUTgrid_f(8,2), INNgrid_i(8,2), INNgrid_f(8,2), &
        & Ncpu_y, Ncpu_z, iy, iz, fy, fz
END TYPE MP
!
CONTAINS
!
! ############### Particle set reconfiguration functions ###################
! ##########################################################################
!
! Resize pointer p into n length
FUNCTION REALLOC_PT(p, n)
  TYPE(PT), POINTER:: p(:), REALLOC_PT(:)
  INTEGER, intent(in):: n
  INTEGER:: nold, ierr
  ALLOCATE(REALLOC_PT(1:n), STAT = ierr)
  IF(ierr /=0) STOP "allocate error"
  IF(.NOT.ASSOCIATED(p)) RETURN
  nold = MIN(size(p),n)
  REALLOC_PT(1:nold) = p(1:nold)
  DEALLOCATE(p)
END FUNCTION REALLOC_PT
!
! Remove some particles from pointer p
FUNCTION REMOVE_PT(p, l, n, m)
  TYPE(PT), POINTER:: p(:), REMOVE_PT(:)
  INTEGER, intent(in):: l, n, m(n)
  INTEGER:: ierr, ii, kk, ni, nf
  ALLOCATE(REMOVE_PT(1:l-n), STAT = ierr)
  IF(ierr /=0) STOP "allocate error"
  ni = 1
  nf = 0
  kk = 1
  DO ii=1, n
     ni = nf + 1
     nf = m(ii) - ii
     REMOVE_PT(ni:nf) = p(kk:m(ii)-1)
     kk = m(ii)+1
  END DO
  IF (m(n) < l) REMOVE_PT(nf+1:l-n ) = p(kk:l)
  DEALLOCATE(p)
END FUNCTION REMOVE_PT
!
! Add some particles to pointer p
FUNCTION ADD_PT(x, ll, nl, q)
  TYPE(PT), POINTER:: x(:), ADD_PT(:)
  INTEGER, intent(in):: ll, nl
  TYPE(PT)           :: q(nl)
  INTEGER:: ierr
  ALLOCATE(ADD_PT(1:ll+nl), STAT = ierr)
  IF(ierr /=0) STOP "allocate error"
  ADD_PT(1:ll) = x(1:ll)
  ADD_PT(ll+1:ll+nl) = q(1:nl)
  DEALLOCATE(x)
END FUNCTION ADD_PT
!
!#############################################################################
!\\\\\\\\\\\\\\\\\ Management of exporting/importing particles ///////////////
SUBROUTINE MIGRATION(NS, q, COMM, param, id_array)
IMPLICIT NONE
INCLUDE 'mpif.h'
!
INTEGER,  PARAMETER    :: Nsize = NMAX
TYPE(NM), INTENT(INOUT):: NS
TYPE(MP), INTENT(IN)   :: COMM
TYPE(PT), POINTER      :: q(:)
TYPE(PM)               :: param
INTEGER                :: id_array(3,3,3)
!
TYPE(PT):: qsend(Nsize,26), qrec(Nsize,26)
INTEGER :: Nsend(26), Nrec(26)
INTEGER :: i, j,  k, n, id, Nprune, Dprune(Nsize), IERR, &
     & ISTATUS(MPI_STATUS_SIZE,52), ireq(52), SRC, DEST, stag, rtag
REAL(KIND=DP):: x(3), rx(3)
!
! Nprune = number of pruned particles = sum of Nsend(:)
Nprune = 0
Nsend = 0
DO n=1, NS%Npt
   i = 0
   j = 0
   k = 0
   x(:) = q(n)%xx(:) - NS%ix(:)
   rx(:) = x(:) - param%box(:)*DNINT(x(:)/param%box(:))
   IF (rx(1) < 0.0) i = -1
   IF (rx(2) < 0.0) j = -1
   IF (rx(3) < 0.0) k = -1
   IF (rx(1) > NS%lx(1)) i = +1
   IF (rx(2) > NS%lx(2)) j = +1
   IF (rx(3) > NS%lx(3)) k = +1
   id = id_array(i+2,j+2,k+2)
   IF ((id /= 0)) THEN
      Nsend(id) = Nsend(id) + 1
      Nprune = Nprune + 1
      qsend(Nsend(id),id)%id    = q(n)%id
      qsend(Nsend(id),id)%nd    = q(n)%nd
      qsend(Nsend(id),id)%xx(:) = q(n)%xx(:)
      qsend(Nsend(id),id)%xv(:) = q(n)%xv(:)
      qsend(Nsend(id),id)%ff(:) = q(n)%ff(:)
      qsend(Nsend(id),id)%pot   = q(n)%pot
      Dprune(Nprune) = n      
   END IF
END DO

IF (Nprune > 0) THEN
   q => REMOVE_PT(q, NS%Npt, Nprune, Dprune(1:Nprune))
   NS%Npt = NS%Npt - Nprune
END IF
!
! export/import
stag = 0
rtag = 0
!
! for DEST/SRC, think about it!
DO i=1, COMM%Npair*2
   DEST = COMM%SRC(i)
   SRC  = COMM%DEST(i)
   CALL MPI_SENDRECV(Nsend(i), 1, MPI_INTEGER, DEST, stag, &
        & Nrec(i), 1, MPI_INTEGER, SRC, rtag, COMM%comm, ISTATUS, IERR)
END DO
DO i=1, COMM%Npair*2
   DEST = COMM%SRC(i)
   SRC  = COMM%DEST(i)
   CALL MPI_Isend(qsend(1:Nsend(i),i), Nsend(i), COMM%NEW, DEST, stag, &
        & COMM%comm, ireq(2*i-1), IERR)
   CALL MPI_Irecv(qrec( 1:Nrec(i), i), Nrec(i),  COMM%NEW, SRC,  rtag, &
        & COMM%comm, ireq(2*i), IERR)
END DO
CALL MPI_Waitall(COMM%Npair*4, ireq, ISTATUS, IERR)
!
! Add data set
DO i=1, COMM%Npair*2
   IF (Nrec(i) > 0) THEN
      q => ADD_PT(q, NS%Npt, Nrec(i), qrec(1:Nrec(i),i))
      NS%Npt = NS%Npt + Nrec(i)
   END IF
END DO
!
RETURN
END SUBROUTINE MIGRATION
!
END MODULE DATASTR
