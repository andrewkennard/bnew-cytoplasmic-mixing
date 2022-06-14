PROGRAM MAIN
  USE KEYS
  USE INDEXARRAYS
  USE COVARUTIL, ONLY : GETALLCOVMATS, GETDRIFTCOUPLINGMAT
  IMPLICIT NONE

  ! read in keyword arguments
  CALL READKEY
  
  ! allocate arrays for linearized indices
  CALL SETUPINDEXARRAYS(KSCL,NMAX)

  SELECT CASE(ACTION)
  CASE('GETCOVARMAT')
     ! calculate the different components of the covariance matrices
     ! for estimating error in parameters
     CALL GETCOVMATDRIVER(WAVETYPE,DEG,NMAX,DEL,TRACKLEN,OUTFILE,HMATFROMFILE,SAVEHMATFILE,HMATFILE)

  CASE('GETDRIFTCOUPLEMAT')
     CALL GETDRIFTCOUPLEDRIVER(WAVETYPE,DEG,NMAX,DEL,DEL2,TRACKLEN,OUTFILE,HMATFROMFILE,HMATFILE,HMATFILE2)
  CASE DEFAULT
     PRINT*, 'Action must be: GETCOVARMAT or GETDRIFTCOUPLEMAT'
     STOP 1
  END SELECT

CONTAINS
  SUBROUTINE GETCOVMATDRIVER(WAVETYPE,DEG,NMAX,DEL,TRACKLEN,OUTFILE,HMATFROMFILE,SAVEHMATFILE,HMATFILE)
    ! calculate components of the covariance matrix for estimating errors
    ! INPUTS:
    ! wavetype: wavelet shape; can be SVG, HAAR, MEAN, or MSD
    ! deg: wavelet degree (only used for svg wavelets)
    ! nmax: maximal wavelet span
    ! tracklen: length of each track
    ! outfile: output file for components of covariance matrix
    ! hmatfromfile: load matrices for individual q from file?
    ! savehmatfile: save matrices for individual q to file?
    ! hmatfile: file name to load/save individual q results; # is replaced by q value
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WAVETYPE, OUTFILE,HMATFILE
    INTEGER, INTENT(IN) :: DEG,NMAX,TRACKLEN
    DOUBLE PRECISION, INTENT(IN) :: DEL
    LOGICAL, INTENT(IN) :: HMATFROMFILE,SAVEHMATFILE
    INTEGER :: A1
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX) :: HUMAT,HVMAT,HLMAT,FUVMAT,FULMAT,FVLMAT
    DOUBLE PRECISION, DIMENSION(AIMAX) :: AVALS,BVALS,FVALS
    

    CALL GETALLCOVMATS(NMAX,DEG,WAVETYPE,DEL,TRACKLEN,&
         & HMATFROMFILE,SAVEHMATFILE,HMATFILE,HUMAT,HVMAT,HLMAT,FUVMAT,&
         & FULMAT,FVLMAT,AVALS,BVALS,FVALS)

    PRINT*, 'Writing matrices to file:', TRIM(ADJUSTL(OUTFILE))
    OPEN(UNIT=88,FILE=OUTFILE,STATUS='UNKNOWN')
    WRITE(88,*) NMAX
    WRITE(88,*) NKMAX
    DO A1 = 1,AIMAX
       WRITE(88,*) HUMAT(A1,:)
    ENDDO
    DO A1 = 1,AIMAX
       WRITE(88,*) HVMAT(A1,:)
    ENDDO
    DO A1 = 1,AIMAX
       WRITE(88,*) HLMAT(A1,:)
    ENDDO
    DO A1 = 1,AIMAX
       WRITE(88,*) FUVMAT(A1,:)
    ENDDO
    DO A1 = 1,AIMAX
       WRITE(88,*) FULMAT(A1,:)
    ENDDO
    DO A1 = 1,AIMAX
       WRITE(88,*) FVLMAT(A1,:)
    ENDDO
    WRITE(88,*) AVALS
    WRITE(88,*) BVALS
    WRITE(88,*) FVALS
    CLOSE(88)
  END SUBROUTINE GETCOVMATDRIVER

    SUBROUTINE GETDRIFTCOUPLEDRIVER(WAVETYPE,DEG,NMAX,DEL1,DEL2,TRACKLEN,OUTFILE,HMATFROMFILE,HMATFILE1,HMATFILE2)
    ! calculate components of the covariance matrix for estimating errors
    ! INPUTS:
    ! wavetype: wavelet shape; can be SVG, HAAR, MEAN, or MSD
    ! deg: wavelet degree (only used for svg wavelets)
    ! nmax: maximal wavelet span
    ! tracklen: length of each track
    ! outfile: output file for components of covariance matrix
    ! hmatfromfile: load matrices for individual q from file?
    ! savehmatfile: save matrices for individual q to file?
    ! hmatfile: file name to load/save individual q results; # is replaced by q value
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WAVETYPE, OUTFILE,HMATFILE1,HMATFILE2
    INTEGER, INTENT(IN) :: DEG,NMAX,TRACKLEN
    DOUBLE PRECISION, INTENT(IN) :: DEL1,DEL2
    LOGICAL, INTENT(IN) :: HMATFROMFILE
    INTEGER :: A1
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX) :: FU1U2MAT
    
    PRINT*, 'Calculating coupling matrix between two exponentially decorrelating drifts...', DEL1, DEL2

    IF (HMATFROMFILE) THEN
       CALL GETDRIFTCOUPLINGMAT(NMAX,DEG,WAVETYPE,DEL,DEL2,TRACKLEN,FU1U2MAT,HMATFILE1,HMATFILE2)
    ELSE
       CALL GETDRIFTCOUPLINGMAT(NMAX,DEG,WAVETYPE,DEL,DEL2,TRACKLEN,FU1U2MAT)
    ENDIF

    ! output coupling matrix to file
  PRINT*, 'Writing coupling matrix to file:', TRIM(ADJUSTL(OUTFILE))
  OPEN(UNIT=88,FILE=OUTFILE,STATUS='UNKNOWN')
  DO A1 = 1,AIMAX
     WRITE(88,*) FU1U2MAT(A1,:)
  ENDDO
  CLOSE(88)

  END SUBROUTINE GETDRIFTCOUPLEDRIVER
END PROGRAM MAIN
