MODULE COVARUTIL
  ! utilities for calculating the covariance matrices of datapoints after BNEW analysis
  USE INDEXARRAYS
  IMPLICIT NONE

CONTAINS
  SUBROUTINE GETALLCOVMATS(NMAX,DEG,WAVETYPE,DEL,TRACKLEN,&
       & LOADFILES,SAVEFILES,FILENAME,HUMAT,HVMAT,HLMAT,FUVMAT,FULMAT,FVLMAT,AVALS,BVALS,FVALS)
    ! get all the matrices necessary for calculating covariance between datapoints
    ! sum over Q values
    ! Optionally, can load in drift matrices from previously generated binary files
    ! these files can have a different NMAX but *must* have the same range of 
    ! k values for each N 
    ! Does not do the rescaling by B here
    ! ----------
    ! inputs:
    ! ----------
    ! NMAX: maximal wavelet span
    ! DEG: wavelet degree
    ! WAVETYPE: wavelet type
    ! DEL: 1/tau where tau is the correlation time
    ! TRACKLEN: track length (N+1 using the notation in the manuscript)
    ! LOADFILES: whether to load in drift data from files
    ! SAVEFILES: whether to save the drift data to file
    ! FILENAME: file name to load or save to; should have a single # that will
    ! get replaced by the value of Q
    ! -------------
    ! outputs:
    ! ------------
    ! HUMAT: 4-th order covariance matrix for drift velocities    
    ! HVMAT: 4-th order covariance for diffusion velocities
    ! HLMAT: 4-th order covariance for localization error
    ! FUVMAT: drift-diffusion coupling
    ! FULMAT: drift-localization error coupling
    ! FVLMAT: diffusion-localization coupling
    ! AVALS, BVALS: wavelet rescaling functions
    ! FVALS: bias component from drift 

    USE WAVELETUTIL, ONLY : GETALLCOEFF
    USE GENUTIL, ONLY : REPLACESUBSTR

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NMAX, DEG, TRACKLEN
    CHARACTER(LEN=*), INTENT(IN) :: WAVETYPE
    DOUBLE PRECISION, INTENT(IN) :: DEL
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX), INTENT(OUT) :: HUMAT,HVMAT,HLMAT,FUVMAT,FULMAT,FVLMAT
    DOUBLE PRECISION, INTENT(OUT) :: AVALS(AIMAX),BVALS(AIMAX), FVALS(AIMAX)
    LOGICAL, INTENT(IN) :: LOADFILES, SAVEFILES
    CHARACTER(LEN=100), INTENT(IN) :: FILENAME    
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX) :: HUMATQ, FUMATQ, FUMAT0,&
         & HVMATQ,HLMATQ, FLMATQ, FVMATQ
    INTEGER :: KMAX,JMAX, QMAX, EMAX
    DOUBLE PRECISION, ALLOCATABLE :: EVEC(:), EMAT(:,:), ALLVELC(:,:), ALLPOSC(:,:), CVELMAT(:,:)
    INTEGER :: IC, JC, A1, A2, Q, K1, K2
    LOGICAL :: ISFILE
    CHARACTER*100 :: FNAME
    CHARACTER*10 :: SNUM
    
    DOUBLE PRECISION :: COEF1, COEF2, COEF
    INTEGER :: TOTCOUNTS(AIMAX,AIMAX)
    

    IF (LOADFILES) THEN
       PRINT*, 'WILL LOAD FROM FILES...', FILENAME
    ENDIF
    IF (SAVEFILES) THEN
       PRINT*, 'WILL SAVE TO FILES...', FILENAME
    ENDIF    

    KMAX = MAXKMAX
    JMAX = KMAX+2*NMAX-1   
    QMAX = TRACKLEN-1
    EMAX = KMAX+2*NMAX+1+QMAX

    ALLOCATE(EVEC(EMAX+1),EMAT(EMAX+1,EMAX+1))
    ALLOCATE(CVELMAT(JMAX*JMAX,AIMAX))
    ALLOCATE(ALLVELC(JMAX,AIMAX), ALLPOSC(JMAX+1,AIMAX))  
    
    ! precalculate exponential values
    CALL SETUPEXPARRAYS(EMAX, DEL,EVEC,EMAT)

    ! calculate wavelet coefficients
    CALL GETALLCOEFF(NMAX,DEG,WAVETYPE,ALLVELC,ALLPOSC,NKMAX,NISTART,AIMAX,JMAX,AVALS,BVALS)
    ! linearized array of pairs of coefficients
    DO IC = 1,JMAX
       DO JC = 1,JMAX             
          CVELMAT(JMAX*(IC-1)+JC,:) = ALLVELC(IC,:)*ALLVELC(JC,:)         
       ENDDO
    ENDDO        
    
    HUMAT = 0D0; HVMAT = 0D0; HLMAT = 0D0; 
    FUVMAT = 0D0; FULMAT =0D0; FVLMAT = 0D0
    TOTCOUNTS = 0

    DO Q = 0,QMAX
       print*, 'Q, QMAX: ', Q, QMAX

       ! file name to load and/or save to
       IF (LOADFILES.OR.SAVEFILES) THEN
          FNAME = FILENAME
          WRITE(SNUM,'(I10)') Q
          CALL REPLACESUBSTR(FNAME,'#',SNUM)
       ENDIF

       IF (LOADFILES) THEN
          ! load previously calculated data from file
          CALL READDATAFROMFILE(FNAME,NMAX,AIMAX,ISFILE,HUMATQ,FUMATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ)

          IF (Q.EQ.0) THEN
             DO A1 = 1,AIMAX
                FVALS(A1) = FUMATQ(A1,A1)
             ENDDO
          ENDIF
       ENDIF
       
       IF (.NOT.LOADFILES.OR..NOT.ISFILE) then ! file not read, recalculate
          ! calculate 2-nd order drift correlations
          CALL GETFUMATQ(Q,JMAX,EMAX,ALLVELC,EVEC,FUMATQ)
          IF (Q.EQ.0) THEN
             DO A1 = 1,AIMAX
                FVALS(A1) = FUMATQ(A1,A1)
             ENDDO
          ENDIF

          ! calculate 4-st order drift correlations
          CALL GETHUMATQ(Q,NMAX,JMAX,EMAX,EVEC,EMAT,CVELMAT,HUMATQ)

          ! subtract to get covariance
          DO A1 = 1,AIMAX
             HUMATQ(A1,:) = HUMATQ(A1,:)-FVALS(A1)*FVALS
          ENDDO

          ! calculate correlations involving only diffusion and loc. error
          CALL GETCORRVXI(Q,NMAX,ALLVELC,ALLPOSC,HVMATQ,HLMATQ,FVMATQ,FLMATQ)  
       ENDIF

       IF (SAVEFILES) THEN
          ! save calculated data to a binary file
          CALL SAVEDATATOFILE(FNAME,NMAX,NKMAX,HUMATQ,FUMATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ)
       ENDIF
        
       ! add this q value to the total
       DO A1 = 1,AIMAX
          K1 = AKVALS(A1)
          COEF1 = TRACKLEN-K1-2*NMAX
          DO A2 = 1,AIMAX
             K2 = AKVALS(A2)
             COEF2 = TRACKLEN-K2-2*NMAX

             IF (COEF2.GT.Q) THEN       
                TOTCOUNTS(A1,A2) = TOTCOUNTS(A1,A2) + MIN(COEF1,COEF2-Q)

                COEF = MIN(COEF1,COEF2-Q)/(COEF1*COEF2);
                HUMAT(A1,A2) = HUMAT(A1,A2) + COEF*HUMATQ(A1,A2)
                HVMAT(A1,A2) = HVMAT(A1,A2) + COEF*HVMATQ(A1,A2)
                HLMAT(A1,A2) = HLMAT(A1,A2) + COEF*HLMATQ(A1,A2)
                FUVMAT(A1,A2) = FUVMAT(A1,A2) + COEF*FUMATQ(A1,A2)*FVMATQ(A1,A2)
                FULMAT(A1,A2) = FULMAT(A1,A2) + COEF*FUMATQ(A1,A2)*FLMATQ(A1,A2)
                FVLMAT(A1,A2) = FVLMAT(A1,A2) + COEF*FLMATQ(A1,A2)*FVMATQ(A1,A2)
             ENDIF

             IF (Q.GT.0.AND.COEF1.GT.Q) THEN ! add the negative q case
                TOTCOUNTS(A1,A2) = TOTCOUNTS(A1,A2) + MIN(COEF2,COEF1-Q)
                COEF = MIN(COEF2,COEF1-Q)/(COEF1*COEF2)

                HUMAT(A1,A2) = HUMAT(A1,A2) + COEF*HUMATQ(A2,A1)
                HVMAT(A1,A2) = HVMAT(A1,A2) + COEF*HVMATQ(A2,A1)
                HLMAT(A1,A2) = HLMAT(A1,A2) + COEF*HLMATQ(A2,A1)
                FUVMAT(A1,A2) = FUVMAT(A1,A2) + COEF*FUMATQ(A2,A1)*FVMATQ(A2,A1)
                FULMAT(A1,A2) = FULMAT(A1,A2) + COEF*FUMATQ(A2,A1)*FLMATQ(A2,A1)
                FVLMAT(A1,A2) = FVLMAT(A1,A2) + COEF*FLMATQ(A2,A1)*FVMATQ(A2,A1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(ALLVELC,ALLPOSC,CVELMAT,EVEC,EMAT)

  END SUBROUTINE GETALLCOVMATS

  SUBROUTINE GETDRIFTCOUPLINGMAT(NMAX,DEG,WAVETYPE,DEL1,DEL2,TRACKLEN,FU1U2MAT,&
       & FILENAME1,FILENAME2)
    ! get the covariance matrix coupling two different drift velocity terms
    ! sum over Q values
    ! Optionally, can load in pre-calculated matrices for the individual drifts
    ! these files can have a different NMAX but *must* have the same range of 
    ! k values for each N 
    ! ----------
    ! inputs:
    ! ----------
    ! NMAX: maximal wavelet span
    ! DEG: wavelet degree
    ! WAVETYPE: wavelet type
    ! DEL1,DEL2: 1/tau where tau is the correlation time for each drift
    ! TRACKLEN: track length (N+1 using the notation in the manuscript)
    ! FILENAME1, FILENAME2: optional arguments; if both are supplied then 
    ! read in the precalculated values for individual drifts from these files;
    ! # is replaced by Q value
    ! -------------
    ! outputs:
    ! ------------
    ! FU1U2MAT: coupling matrix, without gamma prefactor or rescaling by B

    USE WAVELETUTIL, ONLY : GETALLCOEFF
    USE GENUTIL, ONLY : REPLACESUBSTR

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NMAX, DEG, TRACKLEN
    CHARACTER(LEN=*), INTENT(IN) :: WAVETYPE
    DOUBLE PRECISION, INTENT(IN) :: DEL1, DEL2
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX), INTENT(OUT) :: FU1U2MAT
    CHARACTER(LEN=100), INTENT(IN),OPTIONAL :: FILENAME1, FILENAME2
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX) :: FU1MATQ, FU2MATQ, HUMATQ,HVMATQ,HLMATQ,FVMATQ,FLMATQ
    DOUBLE PRECISION :: AVALS(AIMAX),BVALS(AIMAX),FVALS(AIMAX)
    INTEGER :: KMAX,JMAX, QMAX, EMAX
    DOUBLE PRECISION, ALLOCATABLE :: EVEC1(:), EVEC2(:), EMAT(:,:), ALLVELC(:,:), ALLPOSC(:,:)
    INTEGER :: IC, JC, A1, A2, Q, K1, K2
    LOGICAL :: ISFILE, LOADFILES
    CHARACTER*100 :: FNAME
    CHARACTER*10 :: SNUM
    
    DOUBLE PRECISION :: COEF1, COEF2, COEF    
    
    LOADFILES = PRESENT(FILENAME1).AND.PRESENT(FILENAME2)
    IF (LOADFILES) THEN
       PRINT*, 'WILL LOAD FROM FILES...', FILENAME1, FILENAME2
    ENDIF

    KMAX = MAXKMAX
    JMAX = KMAX+2*NMAX-1   
    QMAX = TRACKLEN-1
    EMAX = KMAX+2*NMAX+1+QMAX

    ALLOCATE(EVEC1(EMAX+1),EVEC2(EMAX+1),EMAT(EMAX+1,EMAX+1))
    ALLOCATE(ALLVELC(JMAX,AIMAX), ALLPOSC(JMAX+1,AIMAX))  
    
    ! precalculate exponential values
    IF(.NOT.LOADFILES) THEN
       CALL SETUPEXPARRAYS(EMAX, DEL1,EVEC1,EMAT)
       CALL SETUPEXPARRAYS(EMAX, DEL2,EVEC2,EMAT)
    ENDIF

    ! calculate wavelet coefficients
    CALL GETALLCOEFF(NMAX,DEG,WAVETYPE,ALLVELC,ALLPOSC,NKMAX,NISTART,AIMAX,JMAX,AVALS,BVALS)
    
    FU1U2MAT = 0D0
    DO Q = 0,QMAX
       print*, 'Q, QMAX: ', Q, QMAX

       ! --------------------------
       ! calculations for first drift
       IF (LOADFILES) THEN
          FNAME = FILENAME1
          WRITE(SNUM,'(I10)') Q
          CALL REPLACESUBSTR(FNAME,'#',SNUM)

          ! load previously calculated data from file
          CALL READDATAFROMFILE(FNAME,NMAX,AIMAX,ISFILE,HUMATQ,FU1MATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ)
       ENDIF
       
       IF (.NOT.LOADFILES.OR..NOT.ISFILE) then ! file not read, recalculate
          ! calculate 2-nd order drift correlations
          CALL GETFUMATQ(Q,JMAX,EMAX,ALLVELC,EVEC1,FU1MATQ)           
       ENDIF

       ! --------------------------
       ! calculations for second drift
       IF (LOADFILES) THEN
          FNAME = FILENAME2
          WRITE(SNUM,'(I10)') Q
          CALL REPLACESUBSTR(FNAME,'#',SNUM)

          ! load previously calculated data from file
          CALL READDATAFROMFILE(FNAME,NMAX,AIMAX,ISFILE,HUMATQ,FU2MATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ)
       ENDIF
       
       IF (.NOT.LOADFILES.OR..NOT.ISFILE) then ! file not read, recalculate
          ! calculate 2-nd order drift correlations
          CALL GETFUMATQ(Q,JMAX,EMAX,ALLVELC,EVEC2,FU2MATQ)           
       ENDIF
       
       ! add this q value to the total
       DO A1 = 1,AIMAX
          K1 = AKVALS(A1)
          COEF1 = TRACKLEN-K1-2*NMAX
          DO A2 = 1,AIMAX
             K2 = AKVALS(A2)
             COEF2 = TRACKLEN-K2-2*NMAX

             IF (COEF2.GT.Q) THEN                      
                COEF = MIN(COEF1,COEF2-Q)/(COEF1*COEF2);    
                FU1U2MAT(A1,A2) = FU1U2MAT(A1,A2) + COEF*FU1MATQ(A1,A2)*FU2MATQ(A1,A2)
             ENDIF

             IF (Q.GT.0.AND.COEF1.GT.Q) THEN ! add the negative q case
                COEF = MIN(COEF2,COEF1-Q)/(COEF1*COEF2)
                
                FU1U2MAT(A1,A2) = FU1U2MAT(A1,A2) + COEF*FU1MATQ(A2,A1)*FU2MATQ(A2,A1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(ALLVELC,ALLPOSC,EVEC1,EVEC2,EMAT)

  END SUBROUTINE GETDRIFTCOUPLINGMAT

  SUBROUTINE GETCORRVXI(Q,NMAX,ALLVELC,ALLPOSC,HVMAT,HLMAT,FVMAT,FLMAT)
    ! Get the correlation matrices for velocities v and localization errors xi
    ! inputs:
    ! Q: time separation btwn data points for the covariance (always positive)
    ! NMAX: maximal wavelet span
    ! ALLVELC: coefficients for steps c_i
    ! ALLPOSC: coefficients for positions chat_i
    ! outputs:
    ! HVMAT: value of H^v matrix each m1,m2 index (diffusion)
    ! HLMAT: value of H^xi matrix each m1,m2 index (localization errors)
    ! FVMAT: value of F^v matrix (2nd order diffusive correlations)
    ! FLMAT: value of F^xi matrix (2nd order localization error correlations)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Q, NMAX
    DOUBLE PRECISION, INTENT(IN) :: ALLVELC(:,:), ALLPOSC(:,:)
    DOUBLE PRECISION, INTENT(OUT) :: HVMAT(AIMAX,AIMAX), HLMAT(AIMAX,AIMAX), &
         & FVMAT(AIMAX,AIMAX), FLMAT(AIMAX,AIMAX)
    INTEGER :: A1, A2, N1, N2, K1, K2, IMIN, IMAX, IC

    HVMAT = 0D0
    HLMAT = 0D0
    FLMAT = 0D0
    FVMAT = 0D0

    DO A1 = 1,AIMAX
       N1 = ANVALS(A1); K1 = AKVALS(A1)
       DO A2 = 1,AIMAX
          N2 = ANVALS(A2); K2 = AKVALS(A2)

          IMIN = MAX(-N1,-N2+Q)+NMAX+1
          IMAX = MIN(N1+K1-2,N2+K2-2+Q)+NMAX+1

          HVMAT(A1,A2) = 4*SUM((ALLVELC(IMIN:IMAX,A1)*ALLVELC(IMIN-Q:IMAX-Q,A2))**2)
          FVMAT(A1,A2) = 2*SUM(ALLVELC(IMIN:IMAX,A1)*ALLVELC(IMIN-Q:IMAX-Q,A2))

          HLMAT(A1,A2) = 4*SUM((ALLPOSC(IMIN:IMAX+1,A1)*ALLPOSC(IMIN-Q:IMAX-Q+1,A2))**2)
          FLMAT(A1,A2) = 2*SUM(ALLPOSC(IMIN:IMAX+1,A1)*ALLPOSC(IMIN-Q:IMAX-Q+1,A2))

          DO IC = IMIN,IMAX
             HVMAT(A1,A2) = HVMAT(A1,A2)+8*SUM(ALLVELC(IC,A1)*ALLVELC(IC-Q,A2)*&
                  & ALLVELC(IC+1:IMAX,A1)*ALLVELC(IC+1-Q:IMAX-Q,A2))
             HLMAT(A1,A2) = HLMAT(A1,A2)+8*SUM(ALLPOSC(IC,A1)*ALLPOSC(IC-Q,A2)*&
                  & ALLPOSC(IC+1:IMAX+1,A1)*ALLPOSC(IC+1-Q:IMAX-Q+1,A2))
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE GETCORRVXI

  SUBROUTINE GETFUMATQ(Q,JMAX,EMAX,ALLVELC,EVEC,FUMAT)
    ! get 2-point correlations of drift velocity
    ! inputs:
    ! Q: separation between points to be correlated
    ! JMAX: maximal coefficient index for c_i
    ! EMAX: size of exponential stored vector
    ! ALLVELC: c_i indices for all n,k combinations
    ! EVEC: vector of stored exponential values (see SETUPEXPARRAYS)
    ! outputs:
    ! FUMAT(m1,m2): F^u matrix for m1,m2

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Q, JMAX,EMAX
    DOUBLE PRECISION, INTENT(IN) :: ALLVELC(JMAX,AIMAX), EVEC(EMAX+1)
    DOUBLE PRECISION, INTENT(OUT) :: FUMAT(AIMAX,AIMAX)
    INTEGER :: I1, I2
    DOUBLE PRECISION :: EMATQ(JMAX,JMAX), TMPMAT(JMAX,AIMAX)

    ! get the coupling matrix
    DO I1 = 1,JMAX
       DO I2= 1,JMAX
          EMATQ(I1,I2) = EVEC(ABS(I2+Q-I1)+1)
       ENDDO
    ENDDO

    ! multiply by coefficients
    CALL DGEMM('N','N',JMAX,AIMAX,JMAX,1D0,EMATQ,JMAX,ALLVELC,JMAX,0D0,TMPMAT,JMAX)
    CALL DGEMM('T','N',AIMAX,AIMAX,JMAX,1D0,ALLVELC,JMAX,TMPMAT,JMAX,0D0,FUMAT,AIMAX)

    ! For two dimensions
    FUMAT = FUMAT*2

  END SUBROUTINE GETFUMATQ

  SUBROUTINE GETHUMATQ(Q,NMAX,JMAX,EMAX,EVEC,EMAT,CVELMAT,HUMAT,SAVEFILE)
    ! Get H^(u) matrix for the drift velocities
    ! (assumed to be persistent random walk), for a given value of q
    ! where q is the separation between pairs of points in the adjusted MSD
    ! NOTE: this calculates the 4th order correlation, still have to subtract Fu(q=0) term to get covariance
    ! ---------
    ! inputs:
    ! ---------
    ! Q: time separation between points for the covariance
    ! NMAX: maximal wavelet span
    ! JMAX: maximal coefficient c_i 
    ! EMAX: size of saved exponential arrays
    ! EVEC: stored single exponential array (see SETUPEXPARRAYS)
    ! EMAT: stored double exponential array
    ! CVELMAT: matrix of products of c_i coefficients
    ! .... CVELMAT(JMAX*(i-1)+c,m) has c_i*c_j for n and k combination m
    ! SAVEFILE: optionally, save result to this file name
    ! ------------
    ! outputs:
    ! -------------
    ! HUMAT: H^u matrix
        
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Q,NMAX,JMAX,EMAX
    DOUBLE PRECISION, INTENT(IN) :: EVEC(EMAX+1),EMAT(EMAX+1,EMAX+1)
    DOUBLE PRECISION, INTENT(IN) :: CVELMAT(JMAX*JMAX,AIMAX)
    DOUBLE PRECISION, INTENT(OUT) :: HUMAT(AIMAX,AIMAX)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: SAVEFILE
    DOUBLE PRECISION :: COV4U(JMAX*JMAX,JMAX*JMAX)
    DOUBLE PRECISION :: HCMAT(JMAX*JMAX,AIMAX)
    LOGICAL :: DOSAVE
    INTEGER :: A1,A2

    HUMAT = 0D0;

    DOSAVE = PRESENT(SAVEFILE)

    CALL GETCOV4UMAT(Q,NMAX,JMAX,EMAX,EVEC,EMAT,COV4U)
    ! do the matrix multiplication H*C_b
    CALL DGEMM('N','N',JMAX*JMAX,AIMAX,JMAX*JMAX,1D0,COV4U,JMAX*JMAX,CVELMAT,JMAX*JMAX,0D0,HCMAT,JMAX*JMAX)   
        
    ! do the matrix multiplication C_a*H
    CALL DGEMM('T','N',AIMAX,AIMAX,JMAX*JMAX,1D0,CVELMAT,JMAX*JMAX,HCMAT,JMAX*JMAX,0D0,HUMAT,AIMAX)   
    
    IF (DOSAVE) THEN
       OPEN(UNIT=99,FILE=SAVEFILE,FORM='UNFORMATTED')
       WRITE(99) HUMAT
       CLOSE(99)
    ENDIF
    
    
  END SUBROUTINE GETHUMATQ
  
  SUBROUTINE GETCOV4UMAT(Q,NMAX,JMAX,EMAX,EVEC,EMAT,COV4U)
    ! get the 4-point covariance matrix 
    ! cov[u_i1 u_j1, u_(i2+q) u_(j2+q)]
    ! first index is linearized combo of i1,j1 
    ! second index is linearized combo of i2,j2
    ! terms: (x,x,x,x)+(x,x,y,y) (or different orders of x,y depending on indices)
    ! multiply result by 2 to get full correlation of 2D vectors
    ! ---------
    ! inputs:
    ! ----------
    ! Q: time separation between data points used for the covariance
    ! NMAX: maximal span
    ! JMAX: maximal i index for c_i; should be kmax+2*nmax-1
    ! EMAX: size for stored exponential array
    ! EVEC: stored single exponential array (see SETUPEXPARRAYS)
    ! EMAT: stored double exponential array
    ! --------
    ! outputs:
    ! --------
    ! COV4U(jmax*(i1-1)+j1, jmax*(i2-1)+j2): 4-point covariance matrix for i1,j1,i2,j2
    
    
    INTEGER, INTENT(IN) :: Q, EMAX,NMAX,JMAX
    DOUBLE PRECISION, INTENT(IN) :: EVEC(EMAX+1),EMAT(EMAX+1,EMAX+1)
    DOUBLE PRECISION, INTENT(OUT) :: COV4U(JMAX*JMAX,JMAX*JMAX)
    INTEGER :: KMAX
    INTEGER :: IC,JC,LC,II,JJ,LL,MM,JCN,LCN, IJMAX, IJMIN, MSHIFT, LQ, IND1, ICN
    LOGICAL :: LSMALL, LBIG

    KMAX=MAXKMAX    

    DO IC = -NMAX,KMAX+NMAX-2
       ICN = IC+NMAX+1
       DO JC = -NMAX, KMAX+NMAX-2
          JCN = JC+NMAX+1
          IND1 = JMAX*(ICN-1)+JCN ! index in first dimension

          IJMIN = MIN(IC,JC)
          IJMAX = MAX(IC,JC)

          DO LC = -NMAX, KMAX+NMAX-2
             LQ = LC+Q

             LCN = LC+NMAX

             LSMALL = LQ.LE.IJMIN
             LBIG = LQ.GE.IJMAX

             IF (LQ<=IJMIN) THEN
                II=LQ; JJ=IJMIN; LL=IJMAX
             ELSEIF (LQ>IJMIN .AND. LQ<=IJMAX) THEN
                II=IJMIN; JJ=LQ; LL=IJMAX
             ELSE
                II=IJMIN; JJ=IJMAX; LL=LQ
             END IF

             MSHIFT = JMAX*LCN+NMAX+1-Q
             !PRINT*, 'TESTX1:', II,JJ,LL,MSHIFT,Q

             DO MM = -NMAX+Q,KMAX+NMAX-2+Q
                IF (MM<=II) THEN
                   IF (LSMALL) THEN
                      COV4U(IND1,MM+MSHIFT) = EVEC(LL-JJ+II-MM+1)
                   ELSE
                      COV4U(IND1,MM+MSHIFT) = EMAT(LL-JJ+II-MM+1,JJ-II+1)
                   ENDIF
                ELSEIF(MM<=JJ) THEN
                   IF (LSMALL) THEN
                      COV4U(IND1,MM+MSHIFT) = EVEC(LL-JJ+MM-II+1)
                   ELSE
                      COV4U(IND1,MM+MSHIFT) = EMAT(LL-JJ+MM-II+1,JJ-MM+1)
                   ENDIF
                ELSEIF(MM<=LL) THEN
                   IF (LBIG) THEN
                      COV4U(IND1,MM+MSHIFT) =  EVEC(LL-MM+JJ-II+1)
                   ELSE
                      COV4U(IND1,MM+MSHIFT) =  EMAT(LL-MM+JJ-II+1,MM-JJ+1)
                   ENDIF
                ELSE
                   IF (LBIG) THEN
                      COV4U(IND1,MM+MSHIFT) = EVEC(MM-LL+JJ-II+1)
                   ELSE
                      COV4U(IND1,MM+MSHIFT) = EMAT(MM-LL+JJ-II+1,LL-JJ+1)
                   ENDIF
                END IF
             ENDDO
          END DO
       END DO
    ENDDO

    ! times 2 for 2 dimensions
    COV4U = COV4U*2
  END SUBROUTINE GETCOV4UMAT

  SUBROUTINE SETUPEXPARRAYS(EMAX, DEL,EVEC,EMAT)
    ! set up precalculated exponential arrays for covariance calculations    
    ! inputs:
    ! EMAX is maximal time separation
    ! DEL is 1/tau where tau is the correlation time
    ! outputs:
    ! EVEC is a list of two-point correlations <u_0 u_i>, (for one dimension)
    ! also used for 4-point correlations of the form (x,x,x,x) + (x,x,y,y)
    ! where the index is t1+t3
    ! EMAT is a matrix of four point corelations of the form (x,x,x,x)+(x,y,x,y)
    ! first index is t1+t3, second index is t2

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: EMAX
    DOUBLE PRECISION, INTENT(IN) :: DEL
    DOUBLE PRECISION, INTENT(OUT) :: EVEC(EMAX+1),EMAT(EMAX+1,EMAX+1)
    DOUBLE PRECISION :: ELIST(EMAX+1), EVEC4(EMAX+1)
    INTEGER :: IC

    ! list of time separations for correlations
    ELIST = (/(DBLE(IC), IC=0,EMAX)/)

    ! Store 2-point correlations for drift velocity (in one dimension)
    EVEC = 1D0/2D0*EXP(-ELIST*DEL)

    ! store matrix of exponential values
    EVEC4 = EXP(-4*ELIST*DEL)    
    DO IC = 1,EMAX       
       EMAT(IC,:) = 1D0/2D0*EVEC(IC)*(EVEC4+1)       
    END DO

  END SUBROUTINE SETUPEXPARRAYS
  
  SUBROUTINE SAVEDATATOFILE(FNAME,NMAX,NKMAX,HUMATQ,FUMATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ)
    ! save covariance matrices for a particular value of q to a binary file
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    INTEGER, INTENT(IN) :: NMAX, NKMAX(:)
    DOUBLE PRECISION, INTENT(IN) :: HUMATQ(:,:),FUMATQ(:,:),&
         &  HVMATQ(:,:),FVMATQ(:,:),HLMATQ(:,:),FLMATQ(:,:)

    ! save calculated matrices to file
    OPEN(UNIT=88,FILE=FNAME,FORM='UNFORMATTED')
    WRITE(88) NMAX
    WRITE(88) NKMAX(1:NMAX)
    WRITE(88) HUMATQ          
    WRITE(88) FUMATQ
    WRITE(88) HVMATQ
    WRITE(88) FVMATQ         
    WRITE(88) HLMATQ
    WRITE(88) FLMATQ
    CLOSE(88)
  END SUBROUTINE SAVEDATATOFILE

  SUBROUTINE READDATAFROMFILE(FNAME,NMAX,AIMAX,ISFILE,HUMATQ,FUMATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ)
    ! read in precalculated data from a binary file (generated with SAVEDATATOFILE)
    ! NMAX: maximal wavespan to read in
    ! AIMAX: maximal linearized index to read in

    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    INTEGER, INTENT(IN) :: NMAX,AIMAX
    LOGICAL, INTENT(OUT) :: ISFILE
    DOUBLE PRECISION, DIMENSION(AIMAX,AIMAX), INTENT(OUT) :: HUMATQ,FUMATQ,HVMATQ,FVMATQ,HLMATQ,FLMATQ
    INTEGER, ALLOCATABLE :: NKMAXLOAD(:)
    DOUBLE PRECISION, ALLOCATABLE :: DATALOAD(:,:)
    INTEGER :: NMAXORIG, AMAXORIG

    ! check if file exists
    ! if not, return without reading anything
    INQUIRE(FILE=FNAME, EXIST=ISFILE)
    IF (.NOT.ISFILE) RETURN    
    
    ! Read in previously calculated data 
    
    OPEN(UNIT=55,FILE=FNAME,FORM='UNFORMATTED',STATUS='old')
    READ(55) NMAXORIG ! nmax used in calculating data for  file
    ALLOCATE(NKMAXLOAD(NMAXORIG)) ! kmax values used in calculating data
    READ(55) NKMAXLOAD         

    AMAXORIG = SUM(NKMAXLOAD)
    IF (AIMAX.GT.AMAXORIG) THEN
       PRINT*, 'ERROR: ORIGINAL MATRICES FOR LOADING ARE TOO SMALL', AIMAX, AMAXORIG
       STOP 1
    ELSEIF (ANY(NKMAXLOAD(1:NMAX).NE.NKMAX)) THEN
       PRINT*, 'ERROR: ORIGINAL MATRICES HAD DIFFERENT K RANGES FOR SOME N'
       PRINT*, NKMAXLOAD(1:NMAX)
       PRINT*, NKMAX
       STOP 1
    ENDIF

    ALLOCATE(DATALOAD(AMAXORIG,AMAXORIG))
    READ(55) DATALOAD
    HUMATQ = DATALOAD(1:AIMAX,1:AIMAX)
    READ(55) DATALOAD
    FUMATQ = DATALOAD(1:AIMAX,1:AIMAX)
    READ(55) DATALOAD
    HVMATQ = DATALOAD(1:AIMAX,1:AIMAX)
    READ(55) DATALOAD
    FVMATQ = DATALOAD(1:AIMAX,1:AIMAX)
    READ(55) DATALOAD
    HLMATQ = DATALOAD(1:AIMAX,1:AIMAX)
    READ(55) DATALOAD
    FLMATQ = DATALOAD(1:AIMAX,1:AIMAX)
    CLOSE(55)        

    DEALLOCATE(NKMAXLOAD,DATALOAD)
  END SUBROUTINE READDATAFROMFILE
END MODULE COVARUTIL
