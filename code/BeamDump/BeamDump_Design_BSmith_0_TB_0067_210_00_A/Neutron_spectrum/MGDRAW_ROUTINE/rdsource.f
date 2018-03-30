*$ CREATE SOURCE.FOR
*COPY SOURCE
*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2006      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA200x:                                *
*                                                                      *
*     Created on 07 january 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 03-mar-06     by    Alfredo Ferrari               *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*  Read source term from a file                                        *
*     Author: Vasilis.Vlachoudis@cern.ch                               *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST

*      PARAMETER (NMAX=1000000)
      PARAMETER (NMAX=5000000)

*
      SAVE LFIRST
      DATA LFIRST / .TRUE. /

*      CHARACTER*250 LINE
      CHARACTER*500 LINE
      INTEGER    NNN, IJ(NMAX)
      DIMENSION  XXX(NMAX), YYY(NMAX), ZZZ(NMAX)
      DIMENSION  UUU(NMAX), VVV(NMAX), WWW(NMAX)
      DIMENSION  ERG(NMAX), TME(NMAX), WGT(NMAX)

      SAVE XXX, ZZZ, YYY
      SAVE UUU, VVV, WWW
      SAVE ERG, TME, WGT
      SAVE ACUT, RNDDIR

*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***

*  |  WHASOU(1) = UNIT number (File name is linked with OPEN)
         LUNRD  = NINT(WHASOU(1))
*  |  WHASOU(2) = Time cut in (s)
         TCUT   = WHASOU(2)
         IF (WHASOU(3).GT.AZRZRZ) THEN
            ACUT = COS(WHASOU(3) * PIPIPI/180.0)
         ELSE
            ACUT = -ONEONE
         END IF
         RNDDIR = ONEONE - COS(WHASOU(4) * PIPIPI/180.0D0)
         WRITE (LUNOUT,*)
         WRITE (LUNOUT,*) " ** rdsource: TCUT=",TCUT
         WRITE (LUNOUT,*) " ** rdsource: ACUT=",ACUT
         WRITE (LUNOUT,*) " ** rdsource: RNDDIR=",RNDDIR
         IF (TCUT .LE. AZRZRZ) TCUT = AINFNT

         NNN = 0
         TOTWG = ZERZER
 10      CONTINUE
            READ( LUNRD, '(A)', ERR=9999, END=20 ) LINE
            IF (LINE(1:1) .EQ. '*') GO TO 10
            READ (LINE,*,ERR=10) X, Y, Z, U, V, W, E, T, WG, I
* |   Time cut
            IF (T.GE.TCUT) GO TO 10

* |   Angular cut
            IF (U*UBEAM + V*VBEAM + W*WBEAM .LT. ACUT) GO TO 10

* |   increase counter
            NNN = NNN + 1
            IF (NNN.GT.NMAX) CALL FLABRT('SOURCE','Increase NMAX')

            IJ(NNN)  = I
            XXX(NNN) = X
            YYY(NNN) = Y
            ZZZ(NNN) = Z
***ANNES INPUT            
            ZZZ(NNN) = Z + 190.0
****

*  | U = dir[nX], V = dir[nY], W = dir[nZ]
*  |  Normalize direction to 1.0
            UVW = SQRT(U**2 + V**2 + W**2)
            UUU(NNN) =  U / UVW
            VVV(NNN) =  V / UVW
            WWW(NNN) =  W / UVW
* if U,V,W are vz, vy, vz then already cosines: UUU(NNN)=U etc.

            ERG(NNN) = E
            TME(NNN) = T
            WGT(NNN) = WG
            TOTWG = TOTWG + WG
         GOTO 10
 20      CONTINUE
         IF (NNN.EQ.0) CALL FLABRT('SOURCE','Error reading file')
         WRITE (LUNOUT,*) '*** rdsource: ',NNN,' particles loaded'
         WRITE (LUNOUT,*) '*** rdsource: ',TOTWG,' total weight'
         WRITE (LUNOUT,*)
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
* Npflka is the stack counter: of course any time source is called it
* must be =0
      NPFLKA = NPFLKA + 1
*  +-------------------------------------------------------------------*
*  | Choose a random particle
      RNDSIG = FLRNDM (RNDSIG)
      N = INT(NNN*RNDSIG)+1
*      WRITE (LUNOUT,*) '*** N=',N, " ", IJ(N)," ", ERG(N)," ", WGT(N)
*  |
*  +-------------------------------------------------------------------*
* Wt is the weight of the particle
      WTFLK  (NPFLKA) = WGT(N)
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
* Particle type (1=proton.....). Ijbeam is the type set by the BEAM
* card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJ(N)
         ILOFLK (NPFLKA) = IJ(N)
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*

* Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
* User dependent flag:
      LOUSE  (NPFLKA) = 0
* User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE

* User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
* Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
* ... to this point: don't change anything
* Particle age (s)
      AGESTK (NPFLKA) = TME(N)
      AKNSHR (NPFLKA) = -TWOTWO
* Group number for "low" energy neutrons, set to 0 anyway
      IGROUP (NPFLKA) = 0
* Kinetic energy of the particle (GeV)
*      TKEFLK (NPFLKA) = SQRT ( PBEAM**2 + AM (IONID)**2 ) - AM (IONID)
      TKEFLK (NPFLKA) = ERG(N)
* Particle momentum
*      PMOFLK (NPFLKA) = PBEAM
      PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
     &                       + TWOTWO * AM (IONID) ) )

* Cosines (tx,ty,tz)
*      TXFLK  (NPFLKA) = UBEAM
*      TYFLK  (NPFLKA) = VBEAM
*      TZFLK  (NPFLKA) = WBEAM

* Randomize direction based on aperture in WHASOU(4) deg
      IF (RNDDIR.GT.ZERZER) THEN
         PHI    = TWOTWO*PIPIPI*FLRNDM(X)
         CTHETA = ONEONE - RNDDIR*FLRNDM(X)
         STHETA = SQRT(ONEONE-CTHETA**2)
*         PRINT *,PHI,CTHETA
* Random vector
         X = COS(PHI)*STHETA
         Y = SIN(PHI)*STHETA
         Z = CTHETA

* Constuct a temporary axis system
* Find orthogonal to Z
         UX = UUU(N)
         UY = VVV(N)
         UZ = WWW(N)

         IF (ABS(UX) .LT. ABS(UY)) THEN
                IF (ABS(UX) .LT. ABS(UZ)) THEN
                    YX =  ZERZER
                    YY =  UZ
                    YZ = -UY
                ELSE
                    YX =  UY
                    YY = -UX
                    YZ =  ZERZER
                ENDIF
         ELSE
                IF (ABS(UY) .LT. ABS(UZ)) THEN
                    YX = -UZ
                    YY =  ZERZER
                    YZ =  UX
                ELSE
                    YX =  UY
                    YY = -UX
                    YZ =  ZERZER
                ENDIF
         ENDIF

         AL = ONEONE/SQRT(YX**2 + YY**2 + YZ**2)
         YX = YX * AL
         YY = YY * AL
         YZ = YZ * AL

         ! X = Y cross Z
         XX = YY*UZ - YZ*UY
         XY = YZ*UX - YX*UZ
         XZ = YX*UY - YY*UX

*         PRINT *,"V",X,Y,Z
*         PRINT *,"X",XX,XY,XZ
*         PRINT *,"Y",YX,YY,YZ
*         PRINT *,"Z",UX,UY,UZ

         XN = X*XX + Y*YX + Z*UX
         YN = X*XY + Y*YY + Z*UY
         ZN = X*XZ + Y*YZ + Z*UZ
         AL = ONEONE/SQRT(XN**2 + YN**2 + ZN**2)
*         PRINT *,"N",XN,YN,ZN

         TXFLK  (NPFLKA) = XN*AL
         TYFLK  (NPFLKA) = YN*AL
         TZFLK  (NPFLKA) = ZN*AL
      ELSE
         TXFLK  (NPFLKA) = UUU(N)
         TYFLK  (NPFLKA) = VVV(N)
         TZFLK  (NPFLKA) = WWW(N)
      END IF

* Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER
* Particle coordinates
*      XFLK   (NPFLKA) = XXX(N) + XBEAM
*      YFLK   (NPFLKA) = YYY(N) + YBEAM
*      ZFLK   (NPFLKA) =          ZBEAM
      XFLK   (NPFLKA) = XXX(N) + XBEAM
      YFLK   (NPFLKA) = YYY(N) + YBEAM
      ZFLK   (NPFLKA) = ZZZ(N) + ZBEAM

*      WRITE (LUNOUT,*) '*** Running particle # ',N

*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN
*=== End of subroutine Source =========================================*
 9999 CONTINUE
      CALL FLABRT('SOURCE','Error reading source file')
      END
