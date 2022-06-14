MODULE KEYS
  ! keyword parameters that are globally used and set from the parameter file
  IMPLICIT NONE

  CHARACTER*100 :: ACTION
  CHARACTER*100 :: OUTFILE, HMATFILE, HMATFILE2
  CHARACTER*20 :: WAVETYPE
  INTEGER :: NMAX, QMIN, QMAX, DEG, TRACKLEN
  DOUBLE PRECISION :: DEL, TAU, DEL2, TAU2, KSCL
  LOGICAL :: HMATFROMFILE, SAVEHMATFILE
  INTEGER :: OVERLAPTYPE
END MODULE KEYS
