SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing stuff
  INTEGER :: TIMEVAL(8), SEED
  !DOUBLE PRECISION :: ROTMAT(3,3)
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM, QRANGESET

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'

  ! maximal wavelet half-span
  NMAX = 5

  ! min and max values of time separation Q
  QMIN = 1
  QMAX = 1
  QRANGESET = .FALSE.
  
  ! output file
  OUTFILE = '*.txt'

  ! file for saving h matrices
  HMATFILE = 'Hmat_Q#.bin'
  HMATFILE2 = 'Hmat_Q#.bin'
  ! read Hmat from files rather than recalculating
  HMATFROMFILE = .FALSE.
  ! save H matrix calculations to file 
  SAVEHMATFILE = .FALSE.

  ! degree of savitsky-golay wavelet
  DEG = 3

  ! dimensionless drift persistence time (as a multiple of time step)
  TAU = 100;

  ! length of track for calculating covariance matrix
  TRACKLEN = 300

  ! Type of wavelet to use
  WAVETYPE = 'SVG'

  ! type of overlap when doing time average
  ! 0: fully overlapping intervals (spacing 1)
  ! 1: intervals spaced at K
  OVERLAPTYPE = 1

  ! kmax as a multiple of n
  KSCL = 1D0


  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO 
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)  
        CASE('DEG')
           CALL READI(DEG)
        CASE('TAU')
           CALL READF(TAU)
           IF (NITEMS.GT.2) THEN
              CALL READF(TAU2)
           ELSE
              TAU2 = TAU
           END IF
        CASE('HMATFROMFILE')
           CALL READO(HMATFROMFILE)
           IF (NITEMS.GT.2) CALL READA(HMATFILE)  
           IF (NITEMS.GT.3) CALL READA(HMATFILE2)
        CASE('KSCL')
           CALL READF(KSCL)
        CASE('SAVEHMATFILE')
           CALL READO(SAVEHMATFILE)
           IF (NITEMS.GT.2) CALL READA(HMATFILE)
        CASE('NMAX')
           CALL READI(NMAX)
        CASE('OUTFILE')
           CALL READA(OUTFILE)        
        CASE('TRACKLEN')
           CALL READI(TRACKLEN)
        CASE('WAVETYPE')
           CALL READA(WAVETYPE,CASESET=1)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO

  ! ----- set some more defaults -----
  IF (.NOT.QRANGESET) THEN
     QMIN = 1; QMAX = 3*NMAX;
  ENDIF

  ! time step relative to correlation time
  DEL = 1D0/TAU;
  DEL2 = 1D0/TAU2;

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  
  IF (WAVETYPE.NE.'SVG'.AND.WAVETYPE.NE.'MSD'&
       & .and.WAVETYPE.NE.'HAAR'.AND.WAVETYPE.NE.'MEAN') THEN
     PRINT*, 'ERROR: WAVETYPE MUST BE SVG, HAAR, MEAN, or MSD FOR NOW: ', WAVETYPE
     STOP 2
  ENDIF  
    
  ! ----------- fill in  file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))  
  ! ---------------------------


  ! print out parameter values
  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)    
  PRINT*, 'Nmax, Deg:', NMAX, DEG
  PRINT*, 'TAU, DEL, TAU2, DEL2:', TAU, DEL, TAU2, DEL2
  print*, 'KSCL:', KSCL
  PRINT*, 'HMATFILE:', trim(adjustl(HMATFILE)), HMATFROMFILE
  PRINT*, 'TRACKLEN:', TRACKLEN
  PRINT*, 'WAVELET TYPE:' ,WAVETYPE  
  print*, '----------------------------------------------------'
  CALL FLUSH(6)

END SUBROUTINE READKEY
