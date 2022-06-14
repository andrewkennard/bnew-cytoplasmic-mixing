MODULE INDEXARRAYS
  ! arrays for dealing with linearized indexing of coefficients
  ! where we're allowed to have a different range of time separations k for each
  ! wavelet span n
  
  IMPLICIT NONE

  ! NKMAX lists the maximal K value for each N
  ! NISTART gives the index at which each N starts in the linearized n,k indices
  ! ANVALS lists the N value associated with each linearized (n,k) index A
  ! AKVALS lists the K value associated with each linearized index 
  ! AIMAX = maximal linearized index
  INTEGER, ALLOCATABLE :: NKMAX(:), NISTART(:), ANVALS(:), AKVALS(:)
  INTEGER :: AIMAX, MAXKMAX

CONTAINS
SUBROUTINE SETUPINDEXARRAYS(KSCL,NMAX)
  ! allocate and initialize arrays for keeping track of n and k indices
  ! assume that for each n value, we have 1<= k <= floor(kscl*n)
  ! and that n goes from 1 to NMAX
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NMAX
  DOUBLE PRECISION, INTENT(IN) :: KSCL
  INTEGER :: N, K

  ALLOCATE(NKMAX(NMAX), NISTART(NMAX))

  NISTART(1) = 1;
  DO N = 1,NMAX
     ! maximal k for each N
     NKMAX(N) = FLOOR(KSCL*N)
     !NKMAX(N) = 2*N+1

     ! starting point for each n in the linearized indices
     IF (N.GT.1) NISTART(N) = NISTART(N-1)+NKMAX(N-1)
  END DO

  ! max linearized index
  AIMAX = NISTART(NMAX)+NKMAX(NMAX)-1

  ! maximum k index over all n
  MAXKMAX = MAXVAL(NKMAX)

  ! n and k associated with each linearized index
  ALLOCATE(ANVALS(AIMAX),AKVALS(AIMAX))

  DO N = 1,NMAX
     DO K = 1,NKMAX(N)
        ANVALS(NISTART(N)+K-1) = N
        AKVALS(NISTART(N)+K-1) = K
     ENDDO
  ENDDO

END SUBROUTINE SETUPINDEXARRAYS

SUBROUTINE CLEANUPINDEXARRAYS
  ! deallocate global indexing arrays
  IMPLICIT NONE

  DEALLOCATE(NKMAX, NISTART, ANVALS, AKVALS)

END SUBROUTINE CLEANUPINDEXARRAYS


END MODULE INDEXARRAYS
