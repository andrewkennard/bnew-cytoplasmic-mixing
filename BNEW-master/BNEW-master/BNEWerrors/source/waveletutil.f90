MODULE WAVELETUTIL
  ! utilities for dealing with wavelet smoothing 
  IMPLICIT NONE

  ! global arrays for dealing with linearized indexing of wavelet coefficients  

CONTAINS

SUBROUTINE GETALLCOEFF(NMAX,DEG,WAVETYPE,ALLVELC,ALLPOSC,NKMAX,NISTART,AIMAX,JMAX,AVALS,BVALS)
  ! get the full matrix of coefficients for positions and velocities
  ! for all different wavelet half-spans (n) and separations (k)
  ! ---------
  ! input:
  ! --------
  ! NMAX : maximal wavelet span
  ! DEG: degree for Savitzky-Golay wavelet; only use spans N>=DEG
  ! WAVETYPE: wavelet type; can be HAAR, SVG (savitzky-golay),
  ! MEAN (sliding mean), or MSD (ordinary MSD, null wavelet)
  ! NKMAX: maximal K value for each N
  ! NISTART: linearized index for the N,1 coefficient for each N
  ! AIMAX: maximal linearized index
  ! JMAX: kmax+2*nmax (maximal index for the data step)
  ! --------
  ! output:
  ! --------
  ! ALLVELC : c coefficients for each trajectory step
  ! ALLPOSC: c-hat coefficients for each trajectory point
  ! optional output:
  ! AVALS, BVALS: rescaling coefficients (linearized indexing)
  
  INTEGER, INTENT(IN) :: NMAX, DEG, NKMAX(NMAX), NISTART(NMAX),AIMAX,JMAX
  CHARACTER(LEN=*), INTENT(IN) :: WAVETYPE
  DOUBLE PRECISION, INTENT(OUT) :: ALLVELC(JMAX,AIMAX), &
       & ALLPOSC(JMAX+1,AIMAX)
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: AVALS(AIMAX),BVALS(AIMAX)
  DOUBLE PRECISION :: WS(2*NMAX+1), WHAT(2*NMAX),VELC(JMAX),POSC(JMAX+1)
  INTEGER :: N, K, A
  
  ALLVELC = 0; ALLPOSC = 0;
  IF (PRESENT(AVALS)) AVALS=0D0
  IF (PRESENT(BVALS)) BVALS = 0D0
  
  DO N = 1,NMAX
     IF (WAVETYPE.EQ.'SVG') THEN ! Savitzky-Golay wavelets
        IF (NISTART(N)+NKMAX(N).LT.DEG) THEN
           ! for wavelet coefficients to be well defined, need at least DEG+1 number of datapoints
           CYCLE
        END IF
        CALL GETWAVELETCOEFF_SVG(N,DEG,WS(1:2*N+1),WHAT(1:2*N))
     ELSE IF (WAVETYPE.EQ.'HAAR') THEN ! Haar wavelets
        CALL GETWAVELETCOEFF_HAAR(N,WS(1:2*N+1),WHAT(1:2*N))
     ELSE IF (WAVETYPE.EQ.'MEAN') THEN ! sliding mean wavelets
        CALL GETWAVELETCOEFF_MEAN(N,WS(1:2*N+1),WHAT(1:2*N))
     ELSE IF (WAVETYPE.EQ.'MSD') THEN ! ordinary MSD (null wavelet)
        WS = 0D0; WHAT = 0D0
     ELSE
        PRINT*, 'ERROR: BAD WAVETYPE', WAVETYPE
        STOP 3
     ENDIF

     DO K = 1, NKMAX(N)
        A = NISTART(N)+K-1
        CALL GETVELPOSCOEFF(N,K,WHAT(1:2*N),VELC(1:2*N+K-1),POSC(1:2*N+K))

        ALLVELC(-N+NMAX+1:N+K+NMAX-1,A) = VELC(1:2*N+K-1)
        IF (PRESENT(AVALS)) THEN
           AVALS(A) = SUM(VELC(1:2*N+K-1)**2)
        ENDIF

        ALLPOSC(-N+NMAX+1:N+K+NMAX,A) = POSC(1:2*N+K)
        IF (PRESENT(BVALS)) THEN
           BVALS(A) = 0.5D0*SUM(POSC(1:2*N+K)**2)
        ENDIF
     END DO
  END DO

END SUBROUTINE GETALLCOEFF

SUBROUTINE GETVELPOSCOEFF(N,K,WHAT,VELC,POSC)
  ! for a given (xa(k)-xa(0)) obtained using wavelet of halfspan N
  ! VELC: the coefficient of all the velocities from -N to N+K-1
  ! POSC: the coefficient of all the localization errors from -N to N+K
  ! WHAT: velocity coefficients used to define the wavelet shape

  IMPLICIT NONE
INTEGER, INTENT(IN) :: N, K
DOUBLE PRECISION, INTENT(IN) :: WHAT(2*N)
DOUBLE PRECISION, INTENT(OUT) :: VELC(2*N+K-1), POSC(2*N+K)
INTEGER :: J, I

VELC = 0
POSC = 0

! coefficients from wavelet analysis itself
DO J = -N,N+K-2
   DO I = MAX(J-K+1,-N),MIN(J,N-1)
      VELC(J+N+1) = VELC(J+N+1)+WHAT(I+N+1)
   END DO
END DO

! plain velocities minus corrected ones
VELC(N+1:N+K) = VELC(N+1:N+K)-1
VELC = -VELC

! localization error coefficients
POSC(1) = -VELC(1)
POSC(2:2*N+K-1) = VELC(1:2*N+K-2) - VELC(2:2*N+K-1)
POSC(2*N+K) = VELC(2*N+K-1)

END SUBROUTINE GETVELPOSCOEFF


SUBROUTINE GETWAVELETCOEFF_SVG(N,DEG,WS,WHAT)
  ! get wavelet coefficient for a Savitsky-Golay wavelet of degree DEG (degree in terms of positions. should be odd)
  ! N is wavelet span
  ! WS is the coefficients for position values
  ! WHAT is the coefficients for velocities
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, DEG
  DOUBLE PRECISION, INTENT(OUT) :: WS(2*N+1), WHAT(2*N)
  DOUBLE PRECISION :: A(2*N+1,DEG+1), GMAT(DEG+1,DEG+1), M(2*N,2*N)
  DOUBLE PRECISION :: LAM(DEG+1),G(DEG+1)
  INTEGER :: IPIV(DEG+1), INFO, I, IPIV2(2*N)
  DOUBLE PRECISION :: JVALS(2*N+1)

  WHAT = 0D0
  WS = 0D0

  JVALS = (/(DBLE(I), I=-N,N)/)

  A(:,1)=1D0
  DO I = 1,DEG
     A(:,I+1) = JVALS**I
  END DO

  ! get GMAT=A'*A
  CALL DGEMM('T','N',DEG+1,DEG+1,2*N+1,1D0,A,2*N+1,A,2*N+1,0D0,GMAT,DEG+1)

  G = 0
  G(2)=1D0

  ! solve system of linear equations
  LAM = G
  CALL DGESV(DEG+1,1,GMAT,DEG+1,IPIV,LAM,DEG+1,INFO)

  IF (INFO.NE.0) THEN
     PRINT*, 'FAILED TO SOLVE LINEAR SYSTEM FOR WAVELET COEFFICIENTS', INFO
     STOP 2
  END IF

  ! get WS=A*LAM
  CALL DGEMV('N',2*N+1,DEG+1,1D0,A,2*N+1,LAM,1,0D0,WS,1)

  ! convert to step coefficients
  WHAT(1) = -WS(1)
  DO I = 2,2*N
     WHAT(I) = WHAT(I-1)-WS(I)
  ENDDO

END SUBROUTINE GETWAVELETCOEFF_SVG

SUBROUTINE GETWAVELETCOEFF_HAAR(N,WS,WHAT)
  ! get wavelet coefficient for a Haar wavelet 
  ! N is wavelet half-span
  ! WS is the coefficients for position values
  ! WHAT is the coefficients for velocities
  ! translated from matlab: getoptwavelet_sgolay_pos.m

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(OUT) :: WS(2*N+1), WHAT(2*N)
  DOUBLE PRECISION :: M(2*N,2*N)
  INTEGER :: INFO, I, IPIV2(2*N)

  WHAT = 0D0
  WS = 0D0
  WS(1:N) = -1D0;
  WS(N+2:2*N+1) = 1D0;
  WS = WS/(N*(N+1));

  ! convert to step coefficients
  WHAT(1) = -WS(1)
  DO I = 2,2*N
     WHAT(I) = WHAT(I-1)-WS(I)
  ENDDO

END SUBROUTINE GETWAVELETCOEFF_HAAR

SUBROUTINE GETWAVELETCOEFF_MEAN(N,WS,WHAT)
  ! get wavelet coefficient for a wavelet that estimates velocity by just 
  ! subtracting x_n-x(-n)
  ! N is wavelet half-span
  ! WS is the coefficients for position values
  ! WHAT is the coefficients for velocities
  ! translated from matlab: getoptwavelet_sgolay_pos.m

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(OUT) :: WS(2*N+1), WHAT(2*N)
  DOUBLE PRECISION :: M(2*N,2*N)
  INTEGER :: INFO, I, IPIV2(2*N)

  WHAT = 0D0
  WS = 0D0
  WS(1) = -1D0;
  WS(2*N+1) = 1D0;
  WS = WS/(2*N);

  ! convert to step coefficients
  WHAT(1) = -WS(1)
  DO I = 2,2*N
     WHAT(I) = WHAT(I-1)-WS(I)
  ENDDO

END SUBROUTINE GETWAVELETCOEFF_MEAN

END MODULE WAVELETUTIL
