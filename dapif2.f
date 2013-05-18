C**********************************************************************
C
        SUBROUTINE DALI12 (NBOX,LEN,L,NL,XB,YB,
     *             NPAR,NCHI,L1,L2,LL)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Form (first) interaction list and (second) list of colleagues
C   and the pointers to these lists. The colleagues are the adjacent
C   boxes of the same size. The element of the first list are the
C   separated children of the parent box  and the separated children
C   of the parent's colleagues.
C
C   *** INPUT PARAMETERS:
C
C   NBOX      =  number of boxes
C   L         =  number of levels plus one
C   NL(I)     =  pointer to the first box of ith level
C   XB(I)     =  first coordinate of the center of ith box
C   YB(I)     =  second coordinate of the center of ith box
C   NPAR(I)   = address of ith box's parent (0 for boxes at first level)
C   NCHI(I)   = address of ith box's first child (0 for a childless box)
C
C   *** OUTPUT PARAMETERS:
C
C    L1(I)    = pointer to the first element of the first
C               list for the ith box
C    L2(I)    = pointer to the first element of the scond
C               list for the ith box
C    LL(J)    = address of a box in first list for ith box if
C               J  in [L1(I),L2(I)-1]
C               address of a box in first list for ith box if
C               J  in [L2(I),L1(I+1)-1]
C
C   *** SUBROUTINES CALLED :
C
C       DAIN12
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NL(1),NPAR(1),NCHI(1),L1(1),L2(1),LL(1),
     *          INLI(27),COL(8),ICO,ICOK,ICOK1
        DOUBLE PRECISION  XB(1), YB(1)
        DATA BOXSIZ/128.0/
        D = 0.5 * 1.25 * BOXSIZ
CCCC    CALL PRINF('INSIDE DALI12, NPAR=*',NPAR,NBOX)
CCCC    CALL PRINF('INSIDE DALI12, NCHI=*',NCHI,NBOX)
CCCC    CALL PRINF('INSIDE DALI12, NL=*',NL,L-1)
C
C ----- Initializations
C
        CALL DAIN12 (NBOX,L,NL,XB,YB,NPAR,NCHI,L1,L2,LL)
        IF (L.EQ.2)  THEN
CCC           LEN = L1(NBOX+1)-1
           GOTO 110
        END IF
C
C ----- Loop on the levels
C
        DO 100 LEV = 2,L-1
        D = D/2
           IMIN = NL(LEV)
           IMAX = NL(LEV+1)-1
C
C ----- Loop on the boxes at level lev
C
           DO 90 I=IMIN,IMAX
              NP = NPAR(I)
              LPMIN = L2(NP)
              LPMAX = L1(NP+1)-1
              X = XB(I)
              Y = YB(I)
              NCO = 0
              NLI = 0
              IF (LPMAX .LT. LPMIN) GOTO 65
C
C ----- Loop on parent's colleagues
C
              DO 60 JCO=LPMIN,LPMAX
                 ICO = LL(JCO)
                 ICOK = NCHI(ICO)
                 IF (ICOK .EQ. 0) GOTO 60
                 ICOKPA = NPAR(ICOK)
C
C ----- Loop on children of parent's colleagues
C
                 DO 50 IFTY=1,100000000
                    IF (ICOKPA .NE. ICO) GOTO 51
                    XK = XB(ICOK)
                    YK = YB(ICOK)
                    XD = DABS(X-XK)
                    YD = DABS(Y-YK)
                    IF ((XD.LE.D) .AND. (YD.LE.D)) THEN
C
C ----- They can be colleague => add them to COL
C
                       NCO = NCO + 1
                       COL(NCO) = ICOK
                       ICOK = ICOK+1
                       ICOKPA = NPAR(ICOK)
                    ELSE
C
C ----- Or in (first) interaction list => add them to INLI
C
                       NLI = NLI + 1
                       INLI(NLI) = ICOK
                       ICOK = ICOK + 1
                       ICOKPA = NPAR(ICOK)
                    END IF
 50              CONTINUE
 51              CONTINUE
 60           CONTINUE
C
C ---- Add brothers to colleagues' list  COL
C
 65           CONTINUE
              NBRO = NCHI(NP)
              NPP  = NPAR(NBRO)
              DO 67  IFTY=1,100000000
                 IF(NPP .NE. NP) GOTO 68
                 IF (NBRO .NE. I) THEN
                    NCO = NCO + 1
                    COL(NCO) = NBRO
                 END IF
                    NBRO = NBRO + 1
                    NPP = NPAR(NBRO)
 67           CONTINUE
 68           CONTINUE
C
C ----- Transfer  INLI and COL  in  LL, update L1 and L2,
C       first inli...
C
              LR = L1(I)
              J1 = 0
              DO 70 J=1,NLI
                 J1 = J-1
                 LIND = LR + J1
                 LL(LIND) = INLI(J)
                 J1 = J1+1
 70           CONTINUE
C
C ----- Update L2
C
              L2(I) = LR + J1
C
C ----- Transfer COL
C
              LR = L2(I)
              J1 = 0
              DO 80 J=1,NCO
                 LIND = LR + J1
                 LL(LIND) = COL(J)
                 J1 = J1+1
 80           CONTINUE
C
C ----- Update L1
C
              L1(I+1) = LR + J1
 90        CONTINUE
 100    CONTINUE
        LEN = LIND+1
 110    CONTINUE
        LEN = L1(NBOX+1)-1
CCC     CALL PRINF(' L1 = *',L1,NBOX+1)
CCC     CALL PRINF(' L2 = *',L2,NBOX+1)
CCC     CALL PRINF(' LL = *',LL,LEN)
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DAIN12 (NBOX,L,NL,XB,YB,NPAR,NCHI,L1,L2,LL)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Initialize list of colleagues and
C   (first) interaction list at level 1 (which is empty)
C
C   *** PARAMETERS :
C
C   See subroutine DALI12
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NL(1),NPAR(1),NCHI(1),L1(1),L2(1),LL(1)
        DOUBLE PRECISION XB(1),YB(1)
        IMIN = NL(1)
        IMAX = NL(2)-1
        L2(1) = 1
        L1(1) = 1
        DO 2 I = IMIN,IMAX
C
C ----- (first) interaction list is empty at first level
C       set only  pointers
C
           LC = L2(I)
C
C ----- At the first level every box is
C       the colleague of all the other boxes
C
           DO 1 J=IMIN,IMAX
              IF (I.EQ.J) GOTO 1
              LL(LC) = J
              LC = LC + 1
 1         CONTINUE
           L2(I+1) = LC
           L1(I+1) = LC
 2      CONTINUE
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DALI34 (NAPB,NBOX,LEN3,LEN4,L,NL,
     *                     XB,YB,NBAT,NPAR,NCHI,INDB,L1,L2,JL,L3,L4,
     *                     IDESC,JL3,JL4)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Form two lists (for childless boxes)
C   and the pointers to these lists .
C   LIST 3: adjacent childless descendants of colleagues
C   LIST 4: separated descendents of colleagues only
C           at the coarser possible mesh level
C
C   *** INPUT PARAMETERS :
C
C   NAPB        = maximum number of particles per box
C   NBOX        = number of boxes
C   L         =  number of levels plus one
C   NL(I)     =  pointer to first box of ith level
C   XB(I)     =  first coordinate of the center of ith box
C   YB(I)     =  second coordinate of the center of ith box
C   NBAT(I)   =  number of particles in ith box
C   NPAR(I)   =  address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)   =  address of ith box's first child (0 for a childless box)
C   INDB(I)   = 0 if i is a parent box, 1 if i is a childless box
C   L1(I)     = pointer to the first element of the first
C               list for the ith box
C   L2(I)     = pointer to the first element of the second
C               list for the ith box
C   JL(J)     = address of a box in first list for ith box if
C               J  in [L1(I),L2(I)-1]
C               address of a box in first list for ith box if
C               J  in [L2(I),L1(I+1)-1]
C
C   *** OUTPUT PARAMETERS :
C
C   LEN3     = length of JL3
C   LEN4     = length of JL4
C   L3(I)    = pointer to the first element of third list for ith box
C   L4(I)    = pointer to the first element of fourth list for ith box
C   IDESC(I) = list of elements in descendence of the current
C   JL3(J)   = address of a box in third list for ith box if
C               J  in [L3(I),L3(I+1)-1]
C   JL4(J)   = address of a box in fourth list for ith box if
C               J  in [L4(I),L4(I+1)-1]
C
C   *** SUBROUTINES CALLED :
C
C       DASCAN
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NL(1),NBAT(1),NPAR(1),NCHI(1),INDB(1),L1(1),L2(1),
     *          JL(1),L3(1),L4(1),IDESC(1),JL3(1),JL4(1)
        DOUBLE PRECISION XB(1),YB(1)
C
C ----- Initializations
C
        LEN3 = 0
        LEN4 = 0
        L3(1) = 1
        L4(1) = 1
C
C ----- Loop on the boxes
C
        DO 100 I=1,NBOX
           NB = NBAT(I)
           L2I = L2(I)
           L1I = L1(I+1)-1
C
C ----- If box not childless or childless with no colleague: skip it
C
CCC           IF ((L1I .LT. L2I).OR.(NB .GT. NAPB)) THEN
           IF ((L1I .LT. L2I).OR.(INDB(I) .EQ. 0)) THEN
              L3(I+1) = L3(I)
              L4(I+1) = L4(I)
           ELSE
C
C ----- Initialize IDESC with colleagues
C
              ISTA = 1
              IEND = 0
              NNBAR = 0
              DO 10 J=L2I,L1I
                 IEND = IEND+1
                 IDESC(IEND) = JL(J)
                 NBC = NBAT(JL(J))
CCC                 IF (NBC .GT. NAPB) NNBAR = 1
                 IF (INDB(JL(J)) .EQ. 0) NNBAR = 1
 10           CONTINUE
C
C ----- In case when all colleagues are childless: skip
C
              IF (NNBAR .EQ. 0) THEN
                 L3(I+1) = L3(I)
                 L4(I+1) = L4(I)
                 GOTO 100
              END IF
C
C ----- Compute distance bound between two colleagues
C
              X = XB(I)
              Y = YB(I)
              XC = XB(JL(L2I))
              YC = YB(JL(L2I))
              XD = DABS(X-XC)
              YD = DABS(Y-YC)
              D1 = DMAX1(XD,YD)
C
C ----- Scan the descendence of colleagues and complete
C       lists  JL3 and JL4
C
              ILOOP = 0
              DO  50 IFTY=1,100000000
                 IF(IEND .LT. ISTA) GOTO 51
                 ILOOP = ILOOP+1
                 CALL DASCAN (NAPB,I,X,Y,NB,ILOOP,D1,LEN3,
     *                        LEN4,ISTA,IEND,L,NL,XB,YB,NBAT,
     *                        NPAR,NCHI,INDB,L3,L4,IDESC,JL3,JL4)

 50           CONTINUE
 51           CONTINUE
C
C ----- Update L3 and L4
C
              L3(I+1) = LEN3+1
              L4(I+1) = LEN4+1
           ENDIF
 100    CONTINUE
CCC        CALL PRINF(' L3=*',L3,NBOX+1)
CCC        CALL PRINF(' L4=*',L4,NBOX+1)
CCC        CALL PRINF(' JL3=*',JL3,LEN3+1)
CCC        CALL PRINF(' JL4=*',JL4,LEN4+1)
        RETURN
        END
C
C**********************************************************************
C
           SUBROUTINE DASCAN (NAPB,IA,X,Y,NB,ILOOP,D1,LEN3,LEN4,
     *                      ISTA,IEND,L,NL,XB,YB,NBAT,NPAR,NCHI,
     *                      INDB,L3,L4,IDESC,JL3,JL4)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Scan descendence of colleagues at ILOOP-th level
C   form lists JL3 and JL4
C
C   *** INPUT PARAMETERS :
C
C   NAPB        = maximum number of particles per box
C   IA          = address of curent box
C   X           = first coordinate of the center of the box
C   Y           = second coordinate of the center of the box
C   NB          = number of particles in current box
C   ILOOP       = indicate the level
C                 1: colleagues
C                 2: one level smaller etc...
C   D1        = Distance bound for the first level
C   L         =  number of levels plus one
C   NL(I)     =  pointer to first box of ith level
C   XB(I)     =  first coordinate of the center of ith box
C   YB(I)     =  second coordinate of the center of ith box
C   NBAT(I)   =  number of particles in ith box
C   NPAR(I)   =  address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)   =  address of ith box's first child (0 for a childless box)
C   INDB(I)   = 0 if i is a parent box, 1 if i is a childless box
C
C   *** OUTPUT PARAMETERS :
C
C   LEN3     = length of JL3
C   LEN4     = length of JL4
C   ISTA     = first position of active boxes in IDESC
C   IEND     = last  position in active boxes in IDESC
C   L3(I)    = pointer to the first element of third list for ith box
C   L4(I)    = pointer to the first element of fourth list for ith box
C   IDESC(I) = list of elements in descendence of the current
C   JL3(J)   = address of a box in third list for ith box if
C               J  in [L3(I),L3(I+1)-1]
C   JL4(J)   = address of a box in fourth list for ith box if
C               J  in [L4(I),L4(I+1)-1]
C
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NL(1),NBAT(1),NPAR(1),NCHI(1),INDB(1),L3(1),
     *          L4(1),IDESC(1),JL3(1),JL4(1)
        DOUBLE PRECISION XB(1),YB(1)
C
C ---- Compute distance bound for the current level
C

        D = D1*(0.5+1.25*0.5**ILOOP)
C
C ----- Loop on active elements in idesc
C
        IENDN = IEND
        ISTAN = IEND+1
        DO 100 I=ISTA,IEND
           IBO = IDESC(I)
           NBI = NBAT(IBO)
           INI = INDB(IBO)
C
C ----- Check if adjacent or not.
C
           XI = XB(IBO)
           YI = YB(IBO)
           XD = DABS(X-XI)
           YD = DABS(Y-YI)
           IF ((XD .GT. D) .OR. (YD  .GT. D)) THEN
C
C ---- If not adjacent : add box to list 4
C
              LEN4 = LEN4+1
              JL4(LEN4) = IBO
           ELSE
C
C ----- If adjacent check if childless or not
C
              IF (INI .EQ. 1) THEN
CCC              IF (NBI .LE. NAPB) THEN
C
C ----- If childless : add box to list 3
C
                 LEN3 = LEN3+1
                 JL3(LEN3) = IBO
              ELSE
C
C ----- If not childless : add its children to the list
C       of active elements
C
                 IKID=NCHI(IBO)
                 NPP = NPAR(IKID)
                 DO 50 IFTY=1,100000000
                    IF(NPP.NE.IBO) GOTO 51
                    IENDN = IENDN+1
                    IDESC(IENDN) = IKID
                    IKID = IKID+1
                    NPP = NPAR(IKID)
 50              CONTINUE
 51              CONTINUE
              ENDIF
           ENDIF
 100    CONTINUE
C
C ----- Update delimiters to active elements in idesc
C
        ISTA = ISTAN
        IEND = IENDN
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DAUPWD (N,NAPB,IFLAG,NAT,NBOX,L,NL,
     *               XA,YA,QA,IBOX,XB,YB,NBAT,NPAR,NCHI,
     *               INDB,IAT,JAT,LO,ZI,CHARGI)
C
C**********************************************************************
C
C
C   *** DESCRIPTION :
C
C       Form multipole expansions for childless boxes, shift expansions
C   to all coarser levels
C
C   *** INPUT PARAMETERS:
C
C   N           =  number of terms in the expansions
C   NAPB        =  maximum number of particles per box
C   IFLAG       =  1 for computing fields only or Hilbert matrix
C                  2 for computing potential and electric fields
C   NAT         =  number of particles
C   NBOX        =  number of boxes
C   L           =  number of levels plus one
C   NL(I)       =  pointer to first box of ith level
C   XA(I)       =  first coordinate of ith particle
C   YA(I)       =  second coordinate of ith particle
C   QA(I)       =  charge of ith particle
C   IBOX(I)     =  address of the box containing ith particle
C   XB(I)       =  first coordinate of the center of ith box
C   YB(I)       =  second coordinate of the center of ith box
C   NBAT(I)     =  number of particles in ith box
C   NPAR(I)     =  address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)     =  address of ith box's first child (0 for a childless box)
C   INDB(I)     = 0 if i is a parent box, 1 if i is a childless box
C   IAT(I)      = pointer to the list of particles in ith box
C                 negative for boxes that are not childless
C   JAT(IAT(I)) = is  equal to NBAT(I)
C   JAT(J)      = address of particle contained in ith box
C                 if J in [ IAT(I)+1 , IAT(I) + NBAT(I) ]
C
C   *** OUTPUT PARAMETERS :
C
C   LO(*,I)     = coefficients of i-th box's multipole expansion
C
C   *** LOCAL PARAMETERS :
C
C   ZI       Position of particles
C   CHARGI   Charges of particles
C
C   *** SUBROUTINES CALLED :
C
C       DSTTIN
C       DSEXPI
C       DACRLO
C       DSLOLO
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NL(1),IBOX(1),NBAT(1),NPAR(1),NCHI(1),
     *          INDB(1),IAT(1),JAT(1)
        DOUBLE PRECISION XA(1),YA(1),XB(1),YB(1)
        DOUBLE COMPLEX QA(1),Z,LO(N,1),LO1(61),ZI(1),CHARGI(1),IMAG
        DATA IMAG/(0.0,1.0)/
C
C ----- Call initialization entry points
C
        CALL DSTTIN(IDUM)
        CALL DSEXPI(IDUM)
C
C ----- Initializations of expansions for parent boxes
C
        N1  = N-1
        DO 20 I=1,NBOX
           NB = NBAT(I)
CCC           IF (NB .LE. NAPB) GOTO 20
           IF (INDB(I) .EQ. 1) GOTO 20
           DO 10 J=1,N
              LO(J,I) = 0.0
 10        CONTINUE
 20     CONTINUE
C
C ----- Loop on the levels, fine to coarse
C
        DO 100 LEVINV = 2,L
           LEV = L+1 - LEVINV
           JMIN = NL(LEV)
           JMAX = NL(LEV+1)-1
           DO 90 J = JMIN,JMAX
              NB = NBAT(J)
              NPA = NPAR(J)
C
C ----- If childless create multipole expansion for box
C
CCC              IF (NB .LE. NAPB) THEN
              IF (INDB(J) .EQ. 1) THEN
                 CALL DACRLO (N1,IFLAG,J,XA,YA,QA,
     *                        XB,YB,IAT,JAT,LO(1,J),ZI,CHARGI)
              ENDIF
C
C ----- Compute vector by which the expansion is to be shifted
C
              IF (NPA .EQ. 0) GOTO 90
              XV = XB(J)-XB(NPA)
              YV = YB(J)-YB(NPA)
              Z  = DCMPLX(XV,YV)
C
C ----- Shift multipole expansion to the parent's center
C
              CALL DSLOLO (LO(2,J),Z,LO1(2),N1)
C
C ----- Add to parent's expansion
C
              DO 80 I=1,N
                 LO(I,NPA) = LO(I,NPA)+LO1(I)
 80           CONTINUE
 90        CONTINUE
 100    CONTINUE
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DACRLO (N1,IFLAG,JBOX,XA,YA,QA,
     *                     XB,YB,IAT,JAT,LO,Z,CHARGE)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Create multipole expansions for childless boxes.
C
C   *** INPUT PARAMETERS :
C
C   N1          =  number of terms in the expansions
C   IFLAG       =  1 for computing fields only or Hilbert matrix
C                  2 for computing potential and electric fields
C   JBOX        =  address of the childless box
C   XA(I)       =  first coordinate of ith particle
C   YA(I)       =  second coordinate of ith particle
C   QA(I)       =  charge of ith particle
C   XB(I)       =  first coordinate of the center of ith box
C   YB(I)       =  second coordinate of the center of ith box
C   IAT(I)      = pointer to the list of particles in ith box
C                 negative for boxes that are not childless
C   JAT(IAT(I)) = is  equal to NBAT(I)
C   JAT(J)      = address of particle contained in ith box
C                 if J in [ IAT(I)+1 , IAT(I) + NBAT(I) ]
C
C   *** OUTPUT PARAMETERS :
C
C   LO(*,I)     = coefficients of i-th box's multipole expansion
C
C   *** LOCAL PARAMETERS :
C
C   Z        Position of particles
C   CHARGE   Charges of particles
C
C
C   *** SUBROUTINES CALLED :
C
C       DSDENM
C       DSDEBS
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER IAT(1),JAT(1),N1
        DOUBLE PRECISION XA(1),YA(1),XB(1),YB(1)
        DOUBLE COMPLEX QA(1),X0Y0,LO(1),Z(1),CHARGE(1),IMAG
        DATA IMAG/(0.0,1.0)/
C
C ----- M is the number of particles in the box
C
        IAJ = IAT(JBOX)
        M = JAT(IAJ)
        X = XB(JBOX)
        Y = YB(JBOX)
C
C ----- Form Z and CHARGE
C
        JMIN = IAJ+1
        JMAX = JMIN+M-1
        JJ = 0
        DO 2 J = JMIN,JMAX
           JJ = JJ+1
           JAJ = JAT(J)
           XAT = XA(JAJ)
           YAT = YA(JAJ)
           Z(JJ)= DCMPLX(XAT,YAT)
           CHARGE(JJ) = QA(JAJ)
 2      CONTINUE
        X0Y0 = DCMPLX(XB(JBOX),YB(JBOX))
C
C ----- Create the multipole expansions
C
        IF (IFLAG .EQ. 1)  CALL DSDENM (X0Y0,Z,CHARGE,M,LO(2),N1)
        IF (IFLAG .EQ. 2) CALL DSDEBS (X0Y0,Z,CHARGE,M,LO(2),N1)
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DAFLEX (N,NAPB,IFLAG,NAT,NBOX,L,NL,
     *                     XA,YA,IBOX,XB,YB,NBAT,NPAR,NCHI,
     *                     IAT,JAT,L1,L2,JL,LO,TA,POTEN,FIELD)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       At each level, convert multipole expansion of box to local
C   expansion about the centers of boxes in the first list, and
C   add these to initial local expansions obtained from parent level.
C
C   In order to make use of symmetry considerations, DSLNEW is used
C   to compute four separate expansions b1, b2, b3, b4. By recombining
C   them appropriately, one can generate the full expansions at four
C   (or two) locations obtained by rotations of 90 degrees.
C
C   If a box in the first list has no symetric and has few particles
C   interaction is computed by evaluating the multipole expansion
C   at the locations of all the particles in the box.
C
C   For efficiency, we form the table NT giving the number of terms
C   needed for the same precision when the boxes are farther apart.
C
C   *** INPUT PARAMETERS:
C
C   N           =  number of terms in the expansions
C   NAPB        =  maximum number of particles per box
C   IFLAG       =  1 for computing fields only or Hilbert matrix
C                  2 for computing potential and electric fields
C   NAT         =  number of particles
C   NBOX        =  number of boxes
C   L           =  number of levels plus one
C   NL(I)       =  pointer to first box of ith level
C   XA(I)       =  first coordinate of ith particle
C   YA(I)       =  second coordinate of ith particle
C   IBOX(I)     =  address of the box containing ith particle
C   XB(I)       =  first coordinate of the center of ith box
C   YB(I)       =  second coordinate of the center of ith box
C   NBAT(I)     =  number of particles in ith box
C   NPAR(I)     = address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)     = address of ith box's 1st child (0 for a childless box)
C   IAT(I)      = pointer to the list of particles in ith box
C                 negative for boxes that are not childless
C   JAT(IAT(I)) = is  equal to NBAT(I)
C   JAT(J)      = address of particle contained in ith box
C                 if J in [ IAT(I)+1 , IAT(I) + NBAT(I) ]
C   L1(I)       = pointer to the first element of the first
C                 list for the ith box
C   L2(I)       = pointer to the first element of the second
C                 list for the ith box
C   JL(J)       = address of a box in first list for ith box if
C                 J  in [L1(I),L2(I)-1]
C                 address of a box in first list for ith box if
C                 J  in [L2(I),L1(I+1)-1]
C   LO(*,I)     = coefficients of i-th box's multipole expansion
C
C   *** OUTPUT PARAMETERS :
C
C   TA(*,I)     = local expansion about the center of i-th box
C   POTEN(I)    = potential at the location of i-th particle
C   FIELD(I)    = field at the location of i-th particle
C
C   *** SUBROUTINES CALLED :
C
C       DAINSY
C       DAADDD
C       DASYME
C       DACLLO
C       DSLNEW
C       DAADP1
C       DAADPI
C       DAADM1
C       DAADMI
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NL(1),IBOX(1),NBAT(1),NPAR(1),NCHI(1),IAT(1),JAT(1),
     *          L1(1),L2(1),JL(1),NT(27),ILI(27),ILIA(27)
        DOUBLE PRECISION XB(1),YB(1),XA(1),YA(1),POTEN(1)
        DOUBLE COMPLEX LO(N,1),TA(N,1),FIELD(1),B1(61),B2(61),
     *              B3(61),B4(61),ZBB,ZD,IMAG,ZDA,ZGN
        DOUBLE COMPLEX B1P3(61),B1M3(61),B2P4(61),B2M4(61)
        DATA IMAG/(0.0,1.0)/ , BOXSIZ/128.0/
C
C ----- Initializations
C
        N1 = N-1
        L1MIN = NL(1)
        L2MAX = NL(3)-1
        DO 20 LI = L1MIN,L2MAX
           DO 10 J=1,N
              TA(J,LI)=0.0
 10        CONTINUE
 20     CONTINUE
        D = 0.5 *BOXSIZ
C
C ----- Initialize tables for symetries, form table NT
C
        CALL DAINSY(N,XB,YB,JL,NT,ILI,ILIA)
C
C ----- Initialize tables for additions
C
        CALL DAADDD(N,TA,B1P3,B2P4,B1M3,B2M4)
C
C ----- LOOP ON THE LEVELS
C
        DO 1000 LEV=2,L-1
C
C ----- Distance between two boxes
C
           D2 = D
           D = D/2.
           D3 = D*3.
C
C ----- Loop on the boxes at level LEV
C
           IBMIN = NL(LEV)
           IBMAX = NL(LEV+1)-1
           DO 600 IB = IBMIN,IBMAX
              XBB = XB(IB)
              YBB = YB(IB)
              ZBB = DCMPLX(XBB,YBB)
              NPBB = NPAR(IB)
              XPBB = XB(NPBB)
              YPBB = YB(NPBB)
              L1MIN = L1(IB)
              L1MAX = L2(IB)-1
              IF (L1MAX.LT.L1MIN) GOTO 600
C
C ----- Form tables for symetries
C
              CALL DASYME(IB,XBB,YBB,XPBB,YPBB,D,L1MIN,L1MAX,
     1     ILI,ILIA,XB,YB,JL)
C
C ----- Look for existing boxes and symetries in first list
C       if no symetries  : call ACLLO
C       if any symetries : call DSLNEW
C
C
C ----- Boxes # 1 to 4
C
              ISI = 0
              IF (ILI(1).GT.0) ISI = ISI+1
              IF (ILI(2).GT.0) ISI = ISI+2
              IF (ILI(3).GT.0) ISI = ISI+4
              IF (ILI(4).GT.0) ISI = ISI+8
              IF (ISI.EQ.0) GOTO 41
              IF (ISI.EQ.1) THEN
                 II = ILIA(1)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(1),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.2) THEN
                 II = ILIA(2)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(1),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.4) THEN
                 II = ILIA(3)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(1),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.8) THEN
                 II = ILIA(4)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(1),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE
                 ZD = ZBB + D2
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(1))
                 DO 24 K=1,NT(1)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 24              CONTINUE
                 IF (ILI(1).GT.0)  THEN
                    IIA = ILIA(1)
                    CALL DAADP1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(1))
                 END IF
                 IF (ILI(2).GT.0) THEN
                    IIA = ILIA(2)
                    CALL DAADMI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(1))
                 END IF
                 IF (ILI(3).GT.0) THEN
                    IIA = ILIA(3)
                    CALL DAADM1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(1))
                 END IF
                 IF (ILI(4).GT.0) THEN
                    IIA = ILIA(4)
                    CALL DAADPI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(1))
                 END IF
              END IF
C
C ----- Boxes # 5 to 8
C
 41           CONTINUE
              ISI = 0
              IF (ILI(5).GT.0) ISI = ISI+1
              IF (ILI(6).GT.0) ISI = ISI+2
              IF (ILI(7).GT.0) ISI = ISI+4
              IF (ILI(8).GT.0) ISI = ISI+8
              IF (ISI.EQ.0) GOTO 61
              IF (ISI.EQ.1) THEN
                 II = ILIA(5)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.2) THEN
                 II = ILIA(6)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.4) THEN
                 II = ILIA(7)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.8) THEN
                 II = ILIA(8)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE
                 ZD = ZBB + DCMPLX(D2,D)
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(5))
                 DO 44 K=1,NT(5)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 44              CONTINUE
                 IF (ILI(5).GT.0)  THEN
                    IIA = ILIA(5)
                    CALL DAADP1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
                 IF (ILI(6).GT.0) THEN
                    IIA = ILIA(6)
                    CALL DAADMI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
                 IF (ILI(7).GT.0) THEN
                    IIA = ILIA(7)
                    CALL DAADM1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
                 IF (ILI(8).GT.0) THEN
                    IIA = ILIA(8)
                    CALL DAADPI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
              END IF
C
C ----- Boxes # 9 to 12
C
 61           CONTINUE
              ISI = 0
              IF (ILI(9).GT.0) ISI = ISI+1
              IF (ILI(10).GT.0) ISI = ISI+2
              IF (ILI(11).GT.0) ISI = ISI+4
              IF (ILI(12).GT.0) ISI = ISI+8
              IF (ISI.EQ.0) GOTO 81
              IF (ISI.EQ.1) THEN
                 II = ILIA(9)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(9),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.2) THEN
                 II = ILIA(10)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(9),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.4) THEN
                 II = ILIA(11)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(9),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.8) THEN
                 II = ILIA(12)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(9),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE
                 ZD = ZBB + DCMPLX(D2,D2)
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(9))
                 DO 64 K=1,NT(9)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 64              CONTINUE
                 IF (ILI(9).GT.0)  THEN
                    IIA = ILIA(9)
                    CALL DAADP1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(9))
                 END IF
                 IF (ILI(10).GT.0) THEN
                    IIA = ILIA(10)
                    CALL DAADMI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(9))
                 END IF
                 IF (ILI(11).GT.0) THEN
                    IIA = ILIA(11)
                    CALL DAADM1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(9))
                 END IF
                 IF (ILI(12).GT.0) THEN
                    IIA = ILIA(12)
                    CALL DAADPI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(9))
                 END IF
              END IF
C
C ----- Boxes # 13 to 16
C
 81           CONTINUE
              ISI = 0
              IF (ILI(13).GT.0) ISI = ISI+1
              IF (ILI(14).GT.0) ISI = ISI+2
              IF (ILI(15).GT.0) ISI = ISI+4
              IF (ILI(16).GT.0) ISI = ISI+8
              IF (ISI.EQ.0) GOTO 101
              IF (ISI.EQ.1) THEN
                 II = ILIA(13)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.2) THEN
                 II = ILIA(14)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.4) THEN
                 II = ILIA(15)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE IF (ISI.EQ.8) THEN
                 II = ILIA(16)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(5),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE
                 ZD = ZBB + DCMPLX(D,D2)
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(5))
                 DO 84 K=1,NT(5)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 84              CONTINUE
                 IF (ILI(13).GT.0)  THEN
                    IIA = ILIA(13)
                    CALL DAADP1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
                 IF (ILI(14).GT.0) THEN
                    IIA = ILIA(14)
                    CALL DAADMI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
                 IF (ILI(15).GT.0) THEN
                    IIA = ILIA(15)
                    CALL DAADM1(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
                 IF (ILI(16).GT.0) THEN
                    IIA = ILIA(16)
                    CALL DAADPI(N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(5))
                 END IF
              END IF
C
C ----- Boxes # 17 to 18
C
 101          CONTINUE
              IF ( (ILI(17).GT.0) .AND. (ILI(18).GT.0) ) THEN
                 II = ILIA(17)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(17))
                 DO 102 K=1,NT(17)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 102             CONTINUE
                 IIA = ILIA(18)
                 CALL DAADP1 (N,TA,B1P3,B1M3,B2P4,B2M4,II,NT(17))
                 CALL DAADMI (N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(17))
              ELSE  IF ((ILI(17).GT.0) .AND. (ILI(18).LE.0)) THEN
                 II = ILIA(17)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(17),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE  IF ((ILI(17).LE.0) .AND. (ILI(18).GT.0)) THEN
                 II = ILIA(18)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(17),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              END IF
C
C ----- Boxes # 19 to 20
C
              IF ( (ILI(19).GT.0) .AND. (ILI(20).GT.0) ) THEN
                 II = ILIA(19)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(19))
                 DO 202 K=1,NT(19)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 202             CONTINUE
                 IIA = ILIA(20)
                 CALL DAADP1 (N,TA,B1P3,B1M3,B2P4,B2M4,II,NT(19))
                 CALL DAADMI (N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(19))
              ELSE  IF ((ILI(19).GT.0) .AND. (ILI(20).LE.0)) THEN
                 II = ILIA(19)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(19),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE  IF ((ILI(19).LE.0) .AND. (ILI(20).GT.0)) THEN
                 II = ILIA(20)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(19),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              END IF
C
C ----- Boxes # 21 to 22
C
              IF ( (ILI(21).GT.0) .AND. (ILI(22).GT.0) ) THEN
                 II = ILIA(21)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(21))
                 DO 302 K=1,NT(21)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 302             CONTINUE
                 IIA = ILIA(22)
                 CALL DAADP1 (N,TA,B1P3,B1M3,B2P4,B2M4,II,NT(21))
                 CALL DAADMI (N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(21))
              ELSE  IF ((ILI(21).GT.0) .AND. (ILI(22).LE.0)) THEN
                 II = ILIA(21)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(21),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE  IF ((ILI(21).LE.0) .AND. (ILI(22).GT.0)) THEN
                 II = ILIA(22)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(21),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              END IF
C
C ----- Boxes # 23 to 24
C
              IF ( (ILI(23).GT.0) .AND. (ILI(24).GT.0) ) THEN
                 II = ILIA(23)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(19))
                 DO 402 K=1,NT(19)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 402             CONTINUE
                 IIA = ILIA(24)
                 CALL DAADP1 (N,TA,B1P3,B1M3,B2P4,B2M4,II,NT(19))
                 CALL DAADMI (N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(19))
              ELSE  IF ((ILI(23).GT.0) .AND. (ILI(24).LE.0)) THEN
                 II = ILIA(23)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(19),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE  IF ((ILI(23).LE.0) .AND. (ILI(24).GT.0)) THEN
                 II = ILIA(24)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(19),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              END IF
C
C ----- Boxes # 25 to 26
C
              IF ( (ILI(25).GT.0) .AND. (ILI(26).GT.0) ) THEN
                 II = ILIA(25)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DSLNEW (LO(2,IB),ZBB,ZD,
     *                           B1(2),B2(2),B3(2),B4(2),NT(17))
                 DO 502 K=1,NT(17)+1
                    B1P3(K) =  B1(K) + B3(K)
                    B1M3(K) =  B1(K) - B3(K)
                    B2P4(K) =  B2(K) + B4(K)
                    B2M4(K) = (B2(K) - B4(K))*IMAG
 502             CONTINUE
                 IIA = ILIA(26)
                 CALL DAADP1 (N,TA,B1P3,B1M3,B2P4,B2M4,II,NT(17))
                 CALL DAADMI (N,TA,B1P3,B1M3,B2P4,B2M4,IIA,NT(17))
              ELSE  IF ((ILI(25).GT.0) .AND. (ILI(26).LE.0)) THEN
                 II = ILIA(25)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(17),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              ELSE  IF ((ILI(25).LE.0) .AND. (ILI(26).GT.0)) THEN
                 II = ILIA(26)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(17),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              END IF
C
C ----- Box # 27
C
              IF (ILI(27).GT.0) THEN
                 II = ILIA(27)
                 ZD = DCMPLX(XB(II),YB(II))
                 CALL DACLLO (N,NAPB,IFLAG,IB,II,ZBB,ZD,NT(27),
     *                        XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
              END IF
C
 600        CONTINUE
 1000    CONTINUE
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DAINSY(N,XB,YB,JL,NT,ILI,ILIA)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Initialize tables for symetries. Compute the numbers
C       of terms to be taken in the expansion depending
C       on the distance between boxes.
C
C   *** INPUT PARAMETERS:
C
C   N        = maximum number of terms in expansions
C   XB(I)    = first coordinate of the center of i-th box
C   YB(I)    = second coordinate of the center of i-th box
C   JL(J)    = address of a box in first list for ith box if
C              J  in [L1(I),L2(I)-1]
C              address of a box in first list for ith box if
C              J  in [L2(I),L1(I+1)-1]
C
C   *** OUTPUT PARAMETERS:
C
C   NT(I)     = number of terms in local expansions of i-th box
C               or  symetric boxes
C   ILI(I)    = index on existing boxes in the first list; positive
C               if i-th box exists, i is a local number (see below)
C   ILIA(I)   = give the address of a box having i as local number
C               (see below)
C
C            ***   Local numbers for first list  ***
C
C        10  6 | 2 13 | 9 26         17 10 | 6  2 |13  9
C        14  * | *  * | 5 24         19 14 | *  * | *  5
C        ------|------|-----         ------|------|-----
C         3  * | X  * | 1 22         21  3 | *  X | *  1
C         7  * | *  * |16 20         23  7 | *  * | * 16
C        ------|------|-----         ------|------|-----
C        11 15 | 4  8 |12 18         25 11 |15  4 | 8 12
C        17 19 |21 23 |25 27         27 18 |20 22 |24 26
C
C
C
C        26 24 |22 20 |18 27         27 25 |23 21 |19 17
C        10  6 | 2 13 | 9 25         18 10 | 6  2 |13  9
C        ------|------|-----         ------|------|-----
C        14  * | *  * | 5 23         20 14 | *  * | *  5
C         3  * | X  * | 1 21         22  3 | *  X | *  1
C        ------|------|-----         ------|------|-----
C         7  * | *  * |16 19         24  7 | *  * | * 16
C        11 15 | 4  8 |12 17         26 11 |15  4 | 8 12
C
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION XB(1),YB(1)
        INTEGER JL(1),NT(1),ILI(1),ILIA(1)
        INTEGER TAB1(6,6),TAB2(6,6),TAB3(6,6),TAB4(6,6)
C
C -----  set local numbers : TAB1,TAB2,TAB3,TAB4
C
        DATA TAB1/17,19,21,23,25,27,11,15,4,8,12,18,7,0,0,0,16,20,
     *            3,0,0,0,1,22,14,0,0,0,5,24,10,6,2,13,9,26/
        DATA TAB2/27,18,20,22,24,26,25,11,15,4,8,12,23,7,0,0,0,16,
     *            21,3,0,0,0,1,19,14,0,0,0,5,17,10,6,2,13,9/
        DATA TAB3/11,15,4,8,12,17,7,0,0,0,16,19,3,0,0,0,1,21,
     *            14,0,0,0,5,23,10,6,2,13,9,25,26,24,22,20,18,27/
        DATA TAB4/26,11,15,4,8,12,24,7,0,0,0,16,22,3,0,0,0,1,
     *            20,14,0,0,0,5,18,10,6,2,13,9,27,25,23,21,19,17/
         N1 = N-1
C
C ----- Form table NT
C
        NT(1) = N1
        NN = N1*.935
        NN = NN + 1
        NT(5) = MIN(NN,N1)
        NN = N1*.685
        NN = NN + 1
        NT(9) = NN
        NN = N1*.531
        NN = NN + 1
        NT(17) = NN
        NN = N1*.586
        NN = NN+1
        NT(19) = NN
        NN = N1*.596
        NN = NN +1
        NT(21) = NN
        NN = N1*.468
        NN = NN +1
        NT(27) = NN
        RETURN
C
C**********************************************************************
C
        ENTRY DASYME(I,X,Y,XP,YP,D,L1MIN,L1MAX,
     1     ILI,ILIA,XB,YB,JL)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Form tables ILI and ILIA  (see DASYIN)
C
C   *** INPUT PARAMETERS :
C
C   I        = address of current box
C   X        = first coordinate of center of current box
C   Y        = second coordinate of center of current box
C   XP       = first coordinate of center of current box's parent
C   YP       = second coordinate of center of current box's parent
C   D        = distance between the centers of two colleagues
C              at the current level
C   L1MIN    = pointer to the first position in fist list for box I
C   L1MAX    = pointer to the last  position in fist list for box I
C
C
C**********************************************************************
C
C ----- Initialize ILI and ILIA
C
        DO 5 K=1,27
           ILIA(K)=-1
           ILI(K)=-1
 5     CONTINUE
C
C ----- Loop on the elements in the first list
C       form ILI
C
        IF ((X .LT. XP) .AND. (Y .GT. YP)) THEN
C
C ----- If I is child #1
C
           DO 10 LL = L1MIN,L1MAX
              J = JL(LL)
              XJ = XB(J)
              YJ = YB(J)
              XD = (XJ - X)/D + 3.1
              YD = (YJ - Y)/D + 4.1
              NX = XD
              NY = YD
              NBO = TAB1(NX,NY)
              ILI(NBO) = 1
              ILIA(NBO) = J
 10        CONTINUE
        END IF
        IF ((X .GT. XP) .AND. (Y .GT. YP)) THEN
C
C ----- If I is child #2
C
           DO 20 LL = L1MIN,L1MAX
              J = JL(LL)
              XJ = XB(J)
              YJ = YB(J)
              XD = (XJ - X)/D + 4.1
              YD = (YJ - Y)/D + 4.1
              NX = XD
              NY = YD
              NBO = TAB2(NX,NY)
              ILI(NBO) = 1
              ILIA(NBO) = J
 20        CONTINUE
        END IF
        IF ((X .LT. XP) .AND. (Y .LT. YP)) THEN
C
C ----- If I is child #3
C
           DO 30 LL = L1MIN,L1MAX
              J = JL(LL)
              XJ = XB(J)
              YJ = YB(J)
              XD = (XJ - X)/D + 3.1
              YD = (YJ - Y)/D + 3.1
              NX = XD
              NY = YD
              NBO = TAB3(NX,NY)
              ILI(NBO) = 1
              ILIA(NBO) = J
 30        CONTINUE
        END IF
        IF ((X .GT. XP) .AND. (Y .LT. YP)) THEN
C
C ----- If I is child #4
C
           DO 40 LL = L1MIN,L1MAX
              J = JL(LL)
              XJ = XB(J)
              YJ = YB(J)
              XD = (XJ - X)/D + 4.1
              YD = (YJ - Y)/D + 3.1
              NX = XD
              NY = YD
              NBO = TAB4(NX,NY)
              ILI(NBO) = 1
              ILIA(NBO) = J
 40        CONTINUE
        END IF
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DACLLO (N,NAPB,IFLAG,IB,II,ZB,ZD,NN,
     *                     XA,YA,NBAT,IAT,JAT,LO,TA,POTEN,FIELD)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       If box II has 4 or more particles:
C   convert the multipole expansion of box IB to local expansion
C   about the center of box II.
C       If box II has less than 4 particles:
C   evaluate multipole expansion of IB at every particle in II
C
C   *** INPUT PARAMETERS :
C
C   N           =  maximum number of terms in the expansions
C   NAPB        =  maximum number of particles per box
C   IFLAG       =  1 for computing fields only or Hilbert matrix
C                  2 for computing potential and electric fields
C   IB          =  number of the sending box
C   II          =  number of the receiving box
C   ZB          =  position of the center of box # IB
C   ZD          =  position of the center of box # II
C   NN          =  number of terms in the expansions
C   XA(I)       =  first coordinate of ith particle
C   YA(I)       =  second coordinate of ith particle
C   NBAT(I)     =  number of particles in ith box
C   IAT(I)      = pointer to the list of particles in ith box
C                 negative for boxes that are not childless
C   JAT(IAT(I)) = is  equal to NBAT(I)
C   JAT(J)      = address of particle contained in ith box
C                 if J in [ IAT(I)+1 , IAT(I) + NBAT(I) ]
C                 J  in [L2(I),L1(I+1)-1]
C   LO(*,I)     = coefficients of i-th box's multipole expansion
C
C   *** OUTPUT PARAMETERS :
C
C   TA(*,I)     = local expansion about the center of i-th box
C   POTEN(I)    = potential at the location of i-th particle
C   FIELD(I)    = field at the location of i-th particle
C
C   *** SUBROUTINES CALLED :
C
C       DSLOR1
C       DSLOR0
C       DSLORD
C       DSLOTA
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER NBAT(1),IAT(1),JAT(1)
        DOUBLE PRECISION XA(1),YA(1),POTEN(1)
        DOUBLE COMPLEX LO(N,1),TA(N,1),FIELD(1),B(61),Z,ZB,ZD,
     *               OUT1,OUT2
C
        NBOUN = 3
        IF (NAPB .LT. 3) NBOUN = NAPB -1
        NN1 = NN-1
        NBI = NBAT(II)
C
C ----- If less than 4 particles in box II, evaluate multipole
C       expansion of IB at every particle in II
C
        IF (NBI .LE. NBOUN) THEN
           JMIN = IAT(II)+1
           JMAX = IAT(II)+NBI
           IF (IFLAG .EQ. 1)  THEN
              DO 10 JJ=JMIN,JMAX
                 J = JAT(JJ)
                 X = XA(J)
                 Y = YA(J)
                 Z = DCMPLX(X,Y)
                       CALL DSLOR1(LO(2,IB),ZB,Z,NN1,OUT2)
                       FIELD(J) = FIELD(J) + OUT2
 10           CONTINUE
           ELSE
              DO 20 JJ=JMIN,JMAX
                 J = JAT(JJ)
                 X = XA(J)
                 Y = YA(J)
                 Z = DCMPLX(X,Y)
                       CALL DSLORD(LO(2,IB),ZB,Z,NN1,OUT2)
                       CALL DSLOR0(LO(2,IB),ZB,Z,NN1,OUT1)
                       FIELD(J) = FIELD(J) + OUT2
                       POTEN(J) = POTEN(J) + OUT1
 20           CONTINUE
           END IF
        ELSE
C
C ----- Else convert multipole expansion of IB to a local
C       expansion about the center of of II
C
           CALL DSLOTA (LO(2,IB),ZB,ZD,B(2),NN1)
           DO 100 K=1,NN
              TA(K,II) = TA(K,II) + B(K)
 100       CONTINUE
        END IF
        RETURN
        END
C
C**********************************************************************
C
        SUBROUTINE DAADDD(N,TA,B1P3,B2P4,B1M3,B2M4)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       The for expansion from LOTANEW, B1,B2,B3,B4 have been
C       combined in B1+B3, B2+B4, B1-B3, i*(B2-B4). This subroutine
C       recombine these expansions to generate the full expansions
C       at four (or two) locations obtained by rotations of 90 degrees.
C
C   *** PARAMETERS :
C
C       See DAFLEX
C
C**********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE COMPLEX TA(N,1),B1P3(1),B2P4(1),B1M3(1),B2M4(1)
        DOUBLE COMPLEX IMAG
        DATA IMAG/(0.0,1.0)/
        RETURN
C
C**********************************************************************
C
        ENTRY DAADP1 (N,TA,B1P3,B1M3,B2P4,B2M4,I,NN)
C
C       Rotation angle 0 degree:
C       B = P(1)*(B1+B3 + B2+B4)
C       P(1) = DIAG(1,1,1,...,1,...)
C
C**********************************************************************
C
        NN1 = NN+1
        DO 10 K=1,NN1
           TA(K,I) = TA(K,I)+B1P3(K)+B2P4(K)
 10     CONTINUE
        RETURN
C
C**********************************************************************
C
        ENTRY DAADM1 (N,TA,B1P3,B1M3,B2P4,B2M4,I,NN)
C
C       Rotation angle 180 degrees:
C       B = P(-1)*(B1+B3 - B2+B4)
C       P(-1) = DIAG(1,-1,1,-1,...,1,-1,...)
C
C**********************************************************************
C
        NN1 = NN+1
        DO 20 K=1,NN1,2
           TA(K,I) = TA(K,I)+B1P3(K)-B2P4(K)
 20     CONTINUE
        DO 30 K=2,NN1,2
           TA(K,I) = TA(K,I)-(B1P3(K)-B2P4(K))
 30     CONTINUE
        RETURN
C
C**********************************************************************
C
        ENTRY DAADPI (N,TA,B1P3,B1M3,B2P4,B2M4,I,NN)
C
C       Rotation angle 270 degrees:
C       B = P(i)*((B1-B3) + i*(B2-B4))
C       P(i) = DIAG(1,i,-1,-i,...,1,i,-1,-i,...)
C
C**********************************************************************
C
        NN1 = NN+1
        DO 40 K=1,NN1,4
           TA(K,I) = TA(K,I) + B1M3(K) + B2M4(K)
 40     CONTINUE
        DO 50 K=2,NN1,4
           TA(K,I) = TA(K,I) + (B1M3(K) + B2M4(K))*IMAG
 50     CONTINUE
        DO 60 K=3,NN1,4
           TA(K,I) = TA(K,I) - (B1M3(K) + B2M4(K))
 60     CONTINUE
        DO 70 K=4,NN1,4
           TA(K,I) = TA(K,I) - (B1M3(K) + B2M4(K))*IMAG
 70     CONTINUE
        RETURN
C
C**********************************************************************
C
        ENTRY DAADMI (N,TA,B1P3,B1M3,B2P4,B2M4,I,NN)
C
C       Rotation angle 90 degrees:
C       B = P(-i)*((B1-B3) - i*(B2-B4))
C       P(-i) = DIAG(1,-i,-1,i,...,1,-i,-1,i,...)
C
C**********************************************************************
C
        NN1 = NN+1
        DO 80 K=1,NN1,4
           TA(K,I) = TA(K,I) + B1M3(K) - B2M4(K)
 80     CONTINUE
        DO 90 K=2,NN1,4
           TA(K,I) = TA(K,I) - (B1M3(K) - B2M4(K))*IMAG
 90     CONTINUE
        DO 100 K=3,NN1,4
           TA(K,I) = TA(K,I) - (B1M3(K) - B2M4(K))
 100    CONTINUE
        DO 110 K=4,NN1,4
           TA(K,I) = TA(K,I) + (B1M3(K) - B2M4(K))*IMAG
 110    CONTINUE
        RETURN
        END
C

