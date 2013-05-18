C**********************************************************************
C
	SUBROUTINE DABARR (N,NAPB,IFLAG,NAT,NBOX,L,NL,QTOT,RSCAL,
     *       XAU,YAU,XA,YA,QA,IBOX,XB,YB,NBAT,NPAR,NCHI,
     *       INDB,IAT,JAT,L1,L2,JL,L3,L4,JL3,JL4,
     *       LO,TA,ZI,CHARGI,IBAK,N2,POTEN,FIELD,EPS,CLOSE)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Interactions between childless boxes and their colleagues
C   or the descendence of their colleagues: list #3 and #4.
C   Also compute interactions between particles that are in the
C   same childless box.
C
C   List #3 : childless colleagues and childless adjacent descendents
C             of colleagues.
C   List #4 : Decendents of colleagues that are separated at
C             but whose parent are not .
C   NOTE: If all the colleagues are childless, list #3 is empty,
C         list #2 is then used as if it were list #3.
C
C   Box <--> List #3: direct calculations
C
C   Box  --> List #4: Convert monopole expansion for every particle
C                     in Box to a local expansion about the center
C                     of box in list #4
C
C   List #4 -->  Box: evaluation of multipole expansion of box in
C                     list #4 at every particle in Box.
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
C   QTOT        =  sum of all the charges in the simulation
C   RSCAL       =  scaling parameter
C   XAU(i)      =  same as XA(i) except unscaled
C   YAU(i)      =  same as YA(i) except unscaled
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
C   L3(I)       = pointer to the first element of third list for ith box
C   L4(I)       = pointer to the first element of fourth list for ith box
C   JL3(J)      = address of a box in third list for ith box if
C                  J  in [L3(I),L3(I+1)-1]
C   JL4(J)      = address of a box in fourth list for ith box if
C                  J  in [L4(I),L4(I+1)-1]
C   LO(*,I)     = coefficients of i-th box's multipole expansion
C   EPS         = near approach distance.
C   CLOSE       = is the near approach subroutine to be provided by the
C                 user.It is called when two particles are within EPS of
C                 each other.
C
C   *** OUTPUT PARAMETERS:
C
C   TA(*,I)     = local expansion about the center of i-th box
C   POTEN(I)    = potential at the location of i-th particle
C   FIELD(I)    = field at the location of i-th particle
C
C   *** LOCAL PARAMETERS:
C
C   ZI         = table of position of particles for a given box
C   CHARGI     = table of charges of particles for a given box
C   IBAK       = give the global number of a particle in a given box
C   N2         = number of terms to be actually used in an expansion
C
C   *** SUBROUTINES CALLED :
C
C       DSTABI
C       DSLOT1
C       DSLOT3
C       DSLOR0
C       DSLORD
C       DSLOR1
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER NL(1),IBOX(1),NBAT(1),
     *          NPAR(1),NCHI(1),INDB(1),IAT(1),JAT(1),L1(1),L2(1),
     *          JL(1),L3(1),L4(1),JL3(1),JL4(1),IBAK(1),
     *          N2(1),NTERM(500)
	DOUBLE PRECISION XAU(1),YAU(1),XA(1),YA(1),XB(1),YB(1),POTEN(1)
	DOUBLE COMPLEX QTOT,QA(1),LO(N,1),TA(N,1),FIELD(1),FI1,
     *                 FI2,QP,Q1,Q2,CHARGI(1),ZI(1),Z4,
     *                 A0,B(61),OUT,ZP,IMAG,OUT1,OUT2
	DATA IMAG/(0.0,1.0)/, DHALF/0.5/, COE/180.0/
     *       BOXSIZ/128.0/ ,DONE/1.0/
C
C ----- Call initialization entry point
C
	CALL DSTABI(IFLAG)
	COEF = RSCAL
	COEP = DLOG(RSCAL)
	EPS1 = EPS/RSCAL
	EPS2 =EPS**2
C
C ----- Initializations  form nterm
C
	N1 = N-1
	D = 0.95 * BOXSIZ
	R29=2./9.
	DLR29 = N1*DLOG(R29)
	RD29 = R29/20.
	DD = 0.0
	DO 7 I=1,20
	  DD = DD + RD29
	  DDT = DLR29/DLOG(DD)
	  NTERM(I) = DDT
 7      CONTINUE
C
C ---- Loop on the levels, coarse to fine
C
	DO 400 LEV=1,L-1
	   D = D/2.
	   LMIN = NL(LEV)
	   LMAX = NL(LEV+1)-1
C
C ----- Loop on the boxes of that level
C
	   DO 300 I=LMIN,LMAX
	      NBI = NBAT(I)
C
C ----- If box not childless : skip
C
CCC              IF (NBI .GT. NAPB) GOTO 300
	      IF (INDB(I) .EQ. 0) GOTO 300
	      XI = XB(I)
	      YI = YB(I)
C
C ----- Form tables of positions and charges of particles in i-th box
C
	      IATS = IAT(I)+1
	      IATE = IATS + NBI-1
	      II = 0
	      DO 10 IATI = IATS,IATE
		 II = II+1
		 JATI = JAT(IATI)
		 ZI(II) = DCMPLX(XA(JATI),YA(JATI))
		 CHARGI(II) = QA(JATI)
		 IBAK(II) = JATI
 10           CONTINUE
C
C ---- Fourth list
C
	      L4MIN = L4(I)
	      L4MAX = L4(I+1)-1
C
C ----- If fourth list empty : skip
C
	      IF (L4MAX.LT.L4MIN) GOTO 110
C
C ---- Loop on the elements of fourth list
C
	      DO 100 I4 = L4MIN,L4MAX
		 J4 = JL4(I4)
		 NPJ4 = NPAR(J4)
		 X4 = XB(J4)
		 Y4 = YB(J4)
		 DCP = COE*(X4-XB(NPJ4))**2
		 Z4 = DCMPLX(X4,Y4)
		 NB4 = NBAT(J4)
C
C ----- Interaction I-->J4: convert monopole for every particle in I
C       to a local expansion about the center of J4
C
		 IF (IFLAG .EQ. 1) THEN
		    DO 50 K=1,II
		       IBA=IBAK(K)
		       DCA = (X4-XA(IBA))**2+(Y4-YA(IBA))**2
		       RIND = DCP/DCA
		       IND =  RIND
		       IF (IND .LT. 1)  IND = 1
		       IF (IND .GT. 20) IND = 20
		       N2(K) = NTERM(IND)
		       A0 = CHARGI(K)
		       CALL DSLOT3(A0,ZI(K),Z4,B(2),N2(K))
		       DO 30 J=1,N2(K)+1
			  TA(J,J4) = TA(J,J4) + B(J)
 30                    CONTINUE
 50                 CONTINUE
		 ELSE
		    DO 51 K=1,II
		       IBA=IBAK(K)
		       DCA = (X4-XA(IBA))**2+(Y4-YA(IBA))**2
		       RIND = DCP/DCA
		       IND =  RIND
		       IF (IND .LT. 1)  IND = 1
		       IF (IND .GT. 20) IND = 20
		       N2(K) = NTERM(IND)
		       A0 = CHARGI(K)
		       CALL DSLOT1(A0,ZI(K),Z4,B(2),N2(K))
		       DO 31 J=1,N2(K)+1
			  TA(J,J4) = TA(J,J4) + B(J)
 31                    CONTINUE
 51                 CONTINUE
		 END IF
C
C ----- Interaction J4-->I: evaluation of multipole expansion of J4
C       at every particle in I
C
		 IF (IFLAG .EQ. 1) THEN
		    DO 70 K=1,II
		       CALL DSLOR1(LO(2,J4),Z4,ZI(K),N2(K),OUT2)
		       FIELD(IBAK(K)) = FIELD(IBAK(K)) + OUT2
 70                 CONTINUE
		 ELSE
		    DO 71 K=1,II
		       CALL DSLORD(LO(2,J4),Z4,ZI(K),N2(K),OUT2)
		       CALL DSLOR0(LO(2,J4),Z4,ZI(K),N2(K),OUT1)
		       FIELD(IBAK(K)) = FIELD(IBAK(K)) + OUT2
		       POTEN(IBAK(K)) = POTEN(IBAK(K)) + OUT1
 71                 CONTINUE
		 END IF
 100          CONTINUE
 110          CONTINUE
C
C ----- Adjacent boxes: third list (or second list if all colleagues
C       are childless
C
	      L3MIN = L3(I)
	      L3MAX = L3(I+1)-1
C
C ----- If third list empty: look at second list
C
	      IF (L3MAX.LT.L3MIN) THEN
		 L3MIN = L2(I)
		 L3MAX = L1(I+1)-1
		 IF (L3MAX .LT. L3MIN) GOTO 210
		 LISTNB = 2
	      ELSE
		 LISTNB = 3
	      END IF
C
C ---- Loop on boxes in adjacency (2nd or 3rd) list
C
	      DO 200 I3=L3MIN,L3MAX
		 IF (LISTNB .EQ. 3) J3 = JL3(I3)
		 IF (LISTNB .EQ. 2) J3 = JL(I3)
		 NB3 = NBAT(J3)
CCC                 IF (NB3 .GT. NAPB) GOTO 200
		 IF (INDB(J3) .EQ. 0) GOTO 200
		 X3 = XB(J3)
		 Y3 = YB(J3)
		 XD = DABS(XI-X3)
		 YD = DABS(YI-Y3)
C
C ----- Loop on the particles in J3.
C       Compute their actions on particles in I  directly.
C       Use newton's third law only if boxes have the same size.
C
		 IATS = IAT(J3)+1
		 IATE = IATS + NB3-1
		 IF ((XD .GT. D) .OR. (YD .GT. D)) THEN
		    IF (IFLAG .EQ.1)  THEN
		       DO 140 IP=IATS,IATE
			  JP = JAT(IP)
			  ZP = DCMPLX(XA(JP),YA(JP))
			  QP = QA(JP)
			  DO 120 K=1,II
			     IBA = IBAK(K)
			     XX = XA(IBA)-XA(JP)
			     YY = YA(IBA)-YA(JP)
			     R = XX**2+YY**2
			     IF (R .LE. EPS2) THEN
				I1 = IBA
				I2 = JP
				Q1 = QA(IBA)
				Q2 = QP
				X1 = XAU(IBA)
				Y1 = YAU(IBA)
				X2 = XAU(JP)
				Y2 = YAU(JP)
				CALL CLOSE(EPS1,I1,I2,Q1,Q2,
     *                           X1,Y1,X2,Y2,FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
				FI1 = FI1/COEF
				FIELD(I2)=FIELD(I2)+FI1
			     ELSE
				R = DONE/R
				RX = XX*R
				RY = YY*R
				OUT2 = DCMPLX(RX,-RY)
				FIELD(JP)=FIELD(JP)-OUT2*CHARGI(K)
			     END IF
 120                      CONTINUE
 140                   CONTINUE
		    ELSE
		       DO 141 IP=IATS,IATE
			  JP = JAT(IP)
			  ZP = DCMPLX(XA(JP),YA(JP))
			  QP = QA(JP)
			  DO 121 K=1,II
			     IBA = IBAK(K)
			     XX = XA(IBA)-XA(JP)
			     YY = YA(IBA)-YA(JP)
			     R = XX**2+YY**2
			     IF (R .LE. EPS2) THEN
		                I1 = IBA
				I2 = JP
				Q1 = QA(IBA)
				Q2 = QP
				X1 = XAU(IBA)
				Y1 = YAU(IBA)
				X2 = XAU(JP)
				Y2 = YAU(JP)
				CALL CLOSE(EPS1,I1,I2,Q1,Q2,
     *                          X1,Y1,X2,Y2,FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
CCCCC                           POT1 = POT1 + (QTOT-QA(I2))*COEP
CCCCCC  CALL PRIN2('BEFORE SCALING IN DABARR, QTOT=*',QTOT,1)
CCCC                            POT1 = POT1 +       QA(I1) *COEP
				POT1 = POT1 +       QA(I1) *COEP
				FI1 = FI1/COEF
				POTEN(I2)=POTEN(I2)+POT1
				FIELD(I2)=FIELD(I2)+FI1
			     ELSE
				R = DONE/R
				RX = XX*R
				RY = YY*R
				OUT2 = DCMPLX(RX,-RY)
				OUT3 = -DLOG(R)*DHALF
				OUT3 = OUT3*CHARGI(K)
				POTEN(JP)=POTEN(JP)+OUT3
				FIELD(JP)=FIELD(JP)-OUT2*CHARGI(K)
			     END IF
 121                      CONTINUE
 141                   CONTINUE
		    END IF
		 ELSE
		    IF (IFLAG .EQ. 1)  THEN
		       DO 180 IP=IATS,IATE
			  JP = JAT(IP)
			  ZP = DCMPLX(XA(JP),YA(JP))
			  QP = QA(JP)
			  DO 160 K=1,II
			     IBA = IBAK(K)
			     XX = XA(IBA)-XA(JP)
			     YY = YA(IBA)-YA(JP)
			     R = XX**2+YY**2
			     IF (R .LE. EPS2) THEN
				I1 = IBA
				I2 = JP
				Q1 = QA(IBA)
				Q2 = QP
				X1 = XAU(IBA)
				Y1 = YAU(IBA)
				X2 = XAU(JP)
				Y2 = YAU(JP)
				CALL CLOSE(EPS1,I1,I2,Q1,Q2,
     *                          X1,Y1,X2,Y2,FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
				FI1 = FI1/COEF
				FI2 = FI2/COEF
CCC                             FIELD(I2) = FIELD(I2)-FI1
				FIELD(I2) = FIELD(I2)+FI1
CCC                             FIELD(I1) = FIELD(I1)-FI2
				FIELD(I1) = FIELD(I1)+FI2
			     ELSE
				R = DONE/R
				RX =XX*R
				RY =YY*R
				OUT2 =DCMPLX(RX,-RY)
				FIELD(JP)  = FIELD(JP)-OUT2*CHARGI(K)
				FIELD(IBA) = FIELD(IBA)+OUT2*QP
			     END IF
 160                      CONTINUE
 180                   CONTINUE
		    ELSE
		       DO 181 IP=IATS,IATE
			  JP = JAT(IP)
			  ZP = DCMPLX(XA(JP),YA(JP))
			  QP = QA(JP)
			  DO 161 K=1,II
			     IBA = IBAK(K)
			     XX = XA(IBA)-XA(JP)
			     YY = YA(IBA)-YA(JP)
			     R = XX**2+YY**2
			     IF (R .LE. EPS2) THEN
				I1 = IBA
				I2 = JP
				Q1 = QA(IBA)
				Q2 = QP
				X1 = XAU(IBA)
				Y1 = YAU(IBA)
				X2 = XAU(JP)
				Y2 = YAU(JP)
				CALL CLOSE(EPS1,I1,I2,Q1,Q2,
     *                          X1,Y1,X2,Y2,FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
CCCCC                           POT1 = POT1 + (QTOT-Q2)*COEP
CCCCC                           POT2 = POT2 + (QTOT-Q1)*COEP
				POT1 = POT1 +       Q1 *COEP
				POT2 = POT2 +       Q2 *COEP
				FI1 = FI1/COEF
				FI2 = FI2/COEF
				POTEN(I2) = POTEN(I2)+POT1
				POTEN(I1) = POTEN(I1)+POT2
CCC                             FIELD(I2) = FIELD(I2)-FI1
				FIELD(I2) = FIELD(I2)+FI1
				FIELD(I1) = FIELD(I1)+FI2
			     ELSE
				R = DONE/R
				RX =XX*R
				RY =YY*R
				OUT2 =DCMPLX(RX,-RY)
				OUT3 = -DLOG(R)*DHALF
				OUT3K = OUT3*CHARGI(K)
				OUT3J = OUT3*QP
				POTEN(JP) =POTEN(JP)+OUT3K
				POTEN(IBA)=POTEN(IBA)+OUT3J
				FIELD(JP) =FIELD(JP)-OUT2*CHARGI(K)
				FIELD(IBA)=FIELD(IBA)+OUT2*QP
			     END IF
 161                      CONTINUE
 181                   CONTINUE
		    END IF
		 END IF
 200          CONTINUE
 210          CONTINUE
C
C ----- Compute interactions for the particles
C       that are in the same box using newton's third law
C
	      IF (IFLAG .EQ. 1)  THEN
		 DO 250 K=1,II-1
		    IBA1 = IBAK(K)
		    DO 230 J=K+1,II
		       IBA2 = IBAK(J)
		       XX = XA(IBA1)-XA(IBA2)
		       YY = YA(IBA1)-YA(IBA2)
		       R = XX**2+YY**2
		       IF (R .LE. EPS2) THEN
			  I1 = IBA1
			  I2 = IBA2
			  Q1 = QA(IBA1)
			  Q2 = QA(IBA2)
			  X1 = XAU(IBA1)
			  Y1 = YAU(IBA1)
			  X2 = XAU(IBA2)
			  Y2 = YAU(IBA2)
			  CALL CLOSE(EPS1,I1,I2,Q1,Q2,
     *                    X1,Y1,X2,Y2,FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
			  FI1 = FI1/COEF
			  FI2 = FI2/COEF
CCC                       FIELD(I2) = FIELD(I2)-FI1
			  FIELD(I2) = FIELD(I2)+FI1
			  FIELD(I1) = FIELD(I1)+FI2
		       ELSE
			  R = DONE/R
			  RX = XX*R
			  RY = YY*R
			  OUT2 = DCMPLX(RX,-RY)
			  FIELD(IBA2) = FIELD(IBA2) - OUT2*CHARGI(K)
			  FIELD(IBA1) = FIELD(IBA1) + OUT2*CHARGI(J)
		       END IF
 230                CONTINUE
 250             CONTINUE
	      ELSE
		 DO 251 K=1,II-1
		    IBA1 = IBAK(K)
		    DO 231 J=K+1,II
		       IBA2 = IBAK(J)
		       XX = XA(IBA1)-XA(IBA2)
		       YY = YA(IBA1)-YA(IBA2)
		       R = XX**2+YY**2
		       IF (R .LE. EPS2) THEN
			  I1 = IBA1
			  I2 = IBA2
			  Q1 = QA(IBA1)
			  Q2 = QA(IBA2)
			  X1 = XAU(IBA1)
			  Y1 = YAU(IBA1)
			  X2 = XAU(IBA2)
			  Y2 = YAU(IBA2)
			  CALL CLOSE(EPS1,I1,I2,Q1,Q2,
     *                    X1,Y1,X2,Y2,FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
CCCCC                     POT1 = POT1 + (QTOT-Q2)*COEP
CCCCC                     POT2 = POT2 + (QTOT-Q1)*COEP
			  POT1 = POT1 +       Q1 *COEP
			  POT2 = POT2 +       Q2 *COEP
			  FI1 = FI1/COEF
			  FI2 = FI2/COEF
			  POTEN(I2) = POTEN(I2)+POT1
			  POTEN(I1) = POTEN(I1)+POT2
CCC                       FIELD(I2) = FIELD(I2)-FI1
			  FIELD(I2) = FIELD(I2)+FI1
			  FIELD(I1) = FIELD(I1)+FI2
		       ELSE
			  R = DONE/R
			  RX = XX*R
			  RY = YY*R
			  OUT2 = DCMPLX(RX,-RY)
			  OUT3 = -DLOG(R)*DHALF
			  OUT3K = OUT3*CHARGI(K)
			  OUT3J = OUT3*CHARGI(J)
			  POTEN(IBA2) = POTEN(IBA2) + OUT3K
			  POTEN(IBA1) = POTEN(IBA1) + OUT3J
			  FIELD(IBA2) = FIELD(IBA2) - OUT2*CHARGI(K)
			  FIELD(IBA1) = FIELD(IBA1) + OUT2*CHARGI(J)
		       END IF
 231                CONTINUE
 251             CONTINUE
	      END IF
 300       CONTINUE
 400    CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DADNWD (N,IFLAG,NAT,NBOX,L,NL,XA,YA,QA,IBOX,
     *                     XB,YB,NBAT,NPAR,NCHI,TA,POTEN,FIELD)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Shift all the local expansions to all finest levels.
C   Evaluate the local expansions at the finest possible level
C   (childless boxes).
C
C   *** INPUT PARAMETERS:
C
C   N           =  number of terms in the expansions
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
C   NPAR(I)     = address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)     = address of ith box's 1st child (0 for a childless box)
C
C   *** OUTPUT PARAMETERS:
C
C   TA(*,I)     = local expansion about the center of i-th box
C   POTEN(I)    = potential at the location of i-th particle
C   FIELD(I)    = field at the location of i-th particle
C
C
C   *** SUBROUTINES CALLED :
C
C       DSTATA
C       DSTAYL
C       DSTAYD
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER NL(1),IBOX(1),NBAT(1),
     *          NPAR(1),NCHI(1)
	DOUBLE PRECISION XA(1),YA(1),XB(1),YB(1),POTEN(1)
	DOUBLE COMPLEX QA(1),TA(N,1),FIELD(1),B(61)
	DOUBLE COMPLEX IMAG,Z,ZK,ZV,ZA,ZIB
	DOUBLE COMPLEX Z1,Z2,Z4,Z6,Z9,Z10,OUT1,OUT2
	DATA IMAG/(0.0,1.0)/
	N1 = N-1
CCCC    CALL PRINF('INSIDE DADNWD, NCHI IS*',NCHI,NBOX)
C
C ----- Loop on the boxes
C
	DO 300 I=1,NBOX
CCC        CALL PRINF('I =*',I ,1)
	   NK = NCHI(I)
CCC        CALL PRINF('NK=*',NK,1)
C
C ----- If childless: skip
C
	   IF (NK .EQ. 0) GOTO 300
	   X = XB(I)
	   Y = YB(I)
	   Z = DCMPLX(X,Y)
	   NP0 = I
	   NP  = I
C
C ----- Shift local expansion from parent's center
C       to children's centers, add to children's local expansion
C
	   DO 250 IFTY=1,100000000
	      IF (NP .NE. NP0)  GOTO 251
CCCC       CALL PRINF('NK=*',NK,1)
	      XK = XB(NK)
	      YK = YB(NK)
	      ZK = DCMPLX(XK,YK)
	      ZV = ZK - Z
CCCC          IF( CDABS(ZV) .LT.8 )  GOTO 190
CCCC          CALL PRINF('I=*',I,1)
CCCC          CALL PRINF('NK=*',NK,1)
CCCC          CALL PRINF('NP=*',NP,1)
CCCC          CALL PRIN2('BEFORE DSTATA, ZK=*',ZK,2)
CCCC          CALL PRIN2('BEFORE DSTATA, Z =*',Z ,2)
CCCC          CALL PRIN2('BEFORE DSTATA, ZV=*',ZV,2)
 190       CONTINUE
	      CALL DSTATA(TA(2,I),ZV,B(2),N1)
	      DO 200 K=1,N
		 TA(K,NK) = TA(K,NK) + B(K)
 200          CONTINUE
	      NK = NK+1
	      NP = NPAR(NK)
 250       CONTINUE
 251       CONTINUE
 300    CONTINUE
CCCC        CALL PRINF('FINISHED TATAS, IFLAG=*',IFLAG,1)
C
C ----- Compute the potential for each particle
C
	DO 500 I=1,NAT
	   ZA = DCMPLX(XA(I),YA(I))
	   IB = IBOX(I)
	   ZIB = DCMPLX(XB(IB),YB(IB))
C
C ----- Evaluate local expansion
C
C
C ----- Add what has already been accounted
C       with evaluation of local expansion.
C
	   IF (IFLAG .EQ. 1)  THEN
	      CALL DSTAYL(TA(2,IB),ZIB,ZA,N1,OUT2)
	      FIELD(I) = FIELD(I)+OUT2
	      FIELD(I) = DCONJG(FIELD(I))
	   ELSE
	      CALL DSTAYD(TA(2,IB),ZIB,ZA,N1,OUT2)
	      CALL DSTAYL(TA(2,IB),ZIB,ZA,N1,OUT1)
	      POTEN(I) = POTEN(I)+OUT1
	      FIELD(I) = FIELD(I)+OUT2
	      FIELD(I) = DCONJG(FIELD(I))
	   END IF
 500    CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAAUPO (N,NAPB,IFLAG,NAT,NBOX,L,NL,
     *             XA,YA,XAU,YAU,QA,IBOX,XB,YB,NBAT,NPAR,NCHI,
     *             INDB,IAT,JAT,
     *             LO,XU,YU,X,Y,PO,FI,EPS,RSCAL,QTOT,CLOSEP)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Compute field and/or potential at an arbitrary point
C
C   *** INPUT PARAMETERS:
C
C   N           =  number of terms in the expansions
C   NAPB        =  maximum number of particles per box
C   IFLAG       =  1 for computing fields only or Hilbert Matrix
C                  2 for computing potential and electric fields
C   NAT         =  number of particles
C   NBOX        =  number of boxes
C   L           =  number of levels plus one
C   NL(I)       =  pointer to first box of ith level
C   XA(I)       =  first coordinate of ith particle
C   YA(I)       =  second coordinate of ith particle
C   XAU(I)      =  same as xa(i) except non scaled
C   YAU(I)      =  same as ya(i) except non scaled
C   QA(I)       =  charge of ith particle
C   IBOX(I)     =  address of the box containing ith particle
C   XB(I)       =  first coordinate of the center of ith box
C   YB(I)       =  second coordinate of the center of ith box
C   NBAT(I)     =  number of particles in ith box
C   NPAR(I)     = address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)     = address of ith box's 1st child (0 for a childless box)
C   INDB(I)     = 1 if i is a childless box, 0 if i is a parent box
C   IAT(I)      = pointer to the list of particles in ith box
C                 negative for boxes that are not childless
C   JAT(IAT(I)) = is  equal to NBAT(I)
C   JAT(J)      = address of particle contained in ith box
C                 if J in [ IAT(I)+1 , IAT(I) + NBAT(I) ]
C   LO(*,I)     = coefficients of i-th box's multipole expansion
C   X           = first coordinate of point
C   Y           = second coordinate of point
C   XU          = same as X except unscaled
C   YU          = same as Y except unscaled
C   EPS         = near approach distance.
C   RSCAL       = scaling parameter
C   QTOT        = sum of all charges
C   CLOSEP      = is the near approach subroutine to be provided by the
C                 user.It is called when a particle and the point are
C                 within EPS of each other.
C
C   *** OUTPUT PARAMETERS:
C
C   IERR        =  error code array
C   PO          = potential at the location of point
C   FI          = field at the location of point
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER NL(1),IBOX(1),NBAT(1),NPAR(1),NCHI(1),INDB(1),
     *          IAT(1),JAT(1),KBOX(36),JBOX(36)
	DOUBLE PRECISION XAU(1),YAU(1),XA(1),YA(1),XB(1),YB(1),
     *                   X,Y,PO,POP
	DOUBLE COMPLEX QTOT,QA(1),LO(N,1),FI,FIP,Z,Z1,ZERO,OUT1,
     *                 OUT2,QJ
	DATA ZERO/(0.0,0.0)/ , BOXSIZ/128.0/ ,DONE/1.0/
     *       DHALF/0.5/
C
C ----- Initializations
C
	COEF = RSCAL
	COEP = DLOG(RSCAL)
	EPS2 = EPS**2
	EPS1 = EPS/RSCAL
	Z = DCMPLX(X,Y)
	PO = ZERO
	FI = ZERO
	D  = BOXSIZ*(DONE + DHALF)
	L2MIN = NL(1)
	L2MAX = NL(2)-1
	IEND = 0
	DO 10 K=L2MIN,L2MAX
	   IEN = IEN +1
	   KBOX(IEN) = K
 10     CONTINUE
C
C ----- Loop
C
	DO 250 IFTY=1,100000000
	   IF (IEN .LE. 0) GOTO 251
	   IEND = IEN
	   IEN = 0
	   D  = D/2
	   DO 100 KK = 1,IEND
	      K  = KBOX(KK)
	      NB = NBAT(K)
	      X1 = XB(K)
	      Y1 = YB(K)
	      Z1 = DCMPLX(X1,Y1)
	      XD = DABS(X-X1)
	      YD = DABS(Y-Y1)
C
C ----- If separated, evaluate multipole expansion
C
	      IF ((XD .GE. D) .OR. (YD .GE. D)) THEN
		 IF (IFLAG .EQ. 1)  THEN
		    CALL DSLOR1(LO(2,K),Z1,Z,N,OUT2)
		    FI = FI + OUT2
		 ELSE
		    CALL DSLORD(LO(2,K),Z1,Z,N,OUT2)
		    CALL DSLOR0(LO(2,K),Z1,Z,N,OUT1)
		    FI = FI + OUT2
		    PO = PO + OUT1
		 END IF
C
C ----- Else if box childless: compute direct interactions
C
CCC              ELSE IF (NB .LE. NAPB) THEN
	      ELSE IF (INDB(K) .EQ. 1) THEN
		 JMIN =IAT(K)+1
		 JMAX = IAT(K)+NB
		 IF (IFLAG .EQ. 1)  THEN
		    DO 51 JJ = JMIN,JMAX
		       J = JAT(JJ)
		       XJ = XA(J)
		       YJ = YA(J)
		       QJ = QA(J)
		       XX = X-XJ
		       YY = Y-YJ
		       R = XX**2+YY**2
		       IF (R .LE. EPS2) THEN
			  XJ = XAU(J)
			  YJ = YAU(J)
			  CALL CLOSEP(EPS1,J,QJ,XU,YU,XJ,YJ,FIP,POP)
C
C ----- Scale result from close subroutine
C
			  FIP = FIP/COEF
			  FI = FI + FIP
		       ELSE
			  R = DONE/R
			  RX = XX*R
			  RY = YY*R
			  OUT2 = DCMPLX(RX,-RY)
			  FI = FI + OUT2*QJ
			 END IF
 51                 CONTINUE
		 ELSE
		    DO 52 JJ = JMIN,JMAX
		       J = JAT(JJ)
		       XJ = XA(J)
		       YJ = YA(J)
		       QJ = QA(J)
		       XX = X - XJ
		       YY = Y - YJ
		       R = XX**2+YY**2
		       IF (R .LE. EPS2) THEN
			  XJ = XAU(J)
			  YJ = YAU(J)
			  CALL CLOSEP(EPS1,J,QJ,XU,YU,XJ,YJ,FIP,POP)
C
C ----- Scale result from close subroutine
C
			  FIP = FIP/COEF
CCCCC                     POP = POP + QTOT * COEP
			  POP = POP + QJ   * COEP
			  PO = PO + POP
			  FI = FI + FIP
		       ELSE
			  R = DONE/R
			  RX = XX*R
			  RY = YY*R
			  OUT2 = DCMPLX(RX,-RY)
			  OUT3 = -DLOG(R)*DHALF
			  OUT3 = OUT3*QJ
			  PO = PO + OUT3
			  FI = FI + OUT2*QJ
		       END IF
 52                 CONTINUE
		 END IF
C
C ----- Else if box is not childless: add its children to the
C       non interaction list
C
	      ELSE
		 NK = NCHI(K)
		 NP = K
		 DO 95 IFTY1=1,100000000
		    IF (NP .NE. K) GOTO 96
		    IEN = IEN+1
		    JBOX(IEN)= NK
		    NK = NK+1
		    NP = NPAR(NK)
 95              CONTINUE
 96              CONTINUE
	      END IF
 100       CONTINUE
	   DO 200 J=1,IEN
	      KBOX(J) = JBOX(J)
 200       CONTINUE
 250    CONTINUE
 251    CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAEOUT(IERR)
C
C   *** DESCRIPTION :
C
C       Print out error messages
C
C   *** PARAMETERS :
C
C       IERR  = Error code array
C
C**********************************************************************
C
	INTEGER IERR(1)
C
	IF (IERR(1) .EQ. 1) THEN
	   CALL PRINF(' EPS7 SEEMS TO BE UNREASONABLY LARGE FOR
     *  THE SIMULATION BEING PERFORMED*',I,0)
	   CALL PRINF(' THIS CAUSES AN INCREASE IN THE EXECUTION TIME
     *  FOR THE SUBROUTINE DAPIF2 *',I,0)
	  CALL PRINF(' YOU MIGHT WANT TO SEE WHETHER EPS7 COULD
     *  BE REDUCED*',I,0)
	ELSE IF (IERR(1) .EQ. 4) THEN
	   CALL PRINF(' NUMBER OF TERMS IN EXPANSION IS *',IERR(2),1)
	   CALL PRINF(' WHICH EXCEEDS MAXIMUM CAPACITY OF *',60,1)
	   CALL PRIN2(' AND OBTAINABLE PRECISION.  *',1.0E-6,1)
	   CALL PRINF(' THIS IS A TERMINAL ERROR *',I,0)
	ELSE IF (IERR(1) .EQ. 8) THEN
	   CALL PRINF(' INSUFFICIENT WORKSPACE ALLOTTED *',IERR,0)
	   CALL PRINF(' LENGTH OF WORKSPACE GIVEN IS *',IERR(3),1)
	   CALL PRINF(' AMOUNT OF WORKSPACE NEEDED IS *',IERR(2),1)
	   CALL PRINF(' THIS IS A TERMINAL ERROR *',I,0)
	ELSE IF (IERR(1) .EQ. 16) THEN
	   CALL PRINF
     *  (' TOTAL NUMBER OF PARTICLES IN SIMULATION= *',IERR(2),1)
	   CALL PRINF
     *  (' IS LESS THAN MAXIMUM NUMBER OF PARTICLES
     *  ALLOWED PER BOX:*',IERR(3),1)
	   CALL PRINF(' THIS IS A TERMINAL ERROR *',I,0)
	ELSE IF (IERR(1) .EQ. 32) THEN
	  CALL PRINF (' NUMBER OF LEVELS =*',IERR(2),1)
	  CALL PRINF (' CAN NOT EXCEED =*',IERR(3),1)
	  CALL PRINF(' IN THIS SIMULATION *',I,0)
	  CALL PRINF(' THIS IS A TERMINAL ERROR *',I,0)
	ENDIF
	CALL PRINF('   *',IERR,0)
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DLDDPL(X0,Y0,NORM,X,Y,OUT)
C
C        This entry computes the field "out" at the point
C        (x,y) of a unit dipole located at the point (x0,y0)
C        and oriented in the direction norm.
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION X0,Y0,X,Y,OUT,NORM(2),DX,DY
	DX=X-X0
	DY=Y-Y0
	OUT=DX*NORM(1)+DY*NORM(2)
	OUT=-OUT/(DX**2+DY**2)
	RETURN
C
C**********************************************************************
C
	ENTRY DFLCHG(X0,Y0,X,Y,OUT)
C
C        This entry computes the field  "out"  at the point
C        (x,y) of a unit charge located at the point (x0,y0)
C
	OUT=DLOG((X-X0)**2+(Y-Y0)**2)
	OUT=OUT/2
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DSDEDR(X0Y0,Z,NORM,CHARGE,M,A,N)
C
C
C        This subroutine computes the multiple decomposition of the
C        field of a set of dipoles
C         PARAMETERS
C       X0Y0      - the center of the decomposition
C       Z(I)      - the coordinates of the i-th dipole
C       NORM(I)   - the orientation of the i-th dipole
C       CHARGE(I) - the value of the i-th dipole
C       M         - the number of dipoles in the array z
C       A         - the coefficients of the expansion
C       N         - the order of the expansion
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE COMPLEX X0Y0,Z(1),NORM(1),A(1),A1(61),BUF(2),Z0
	DOUBLE COMPLEX    CHARGE(1),D
	EQUIVALENCE(BUF(2),A1(1))
	DATA  ZERO/0.0/,I0/0/
	A(I0)=ZERO
	DO 1200 I=1,N
	   A(I)=ZERO
 1200   CONTINUE
	DO 1600 I=1,M
	   D=CHARGE(I)
	   Z0=Z(I)-X0Y0
	   CALL DSEXP1(Z0,NORM(I),A1,N)
	   DO 1400 J=1,N
	      A(J)=A(J)+A1(J)*D
 1400      CONTINUE
 1600   CONTINUE
	RETURN
C
C**********************************************************************
C
	ENTRY DSDENM(X0Y0,Z,CHARGE,M,A,N)
C
C       This entry computes the multipole decomposition of the
C       derivative of a set of charges
C       PARAMETERS
C       X0Y0      - the center of the decomposition
C       Z(I)      - the coordinates of the i-th dipole
C       CHARGE(I) - the value of the i-th dipole
C       M         - the number of dipoles in the array z
C       A         - the coefficients of the expansion
C       N         - the order of the expansion
C
C
C ----- Set the array a to zero
C
	A(I0)=ZERO
	DO 2200 I=1,N
	   A(I)=ZERO
 2200   CONTINUE
C
C ----- Create the decomposition of derivatives of charges
C
	DO 2600 I=1,M
	   D=CHARGE(I)
	   Z0=Z(I)-X0Y0
	   CALL DSEXDF(Z0,A1,N)
	   DO 2400 J=1,N
	      A(J)=A(J)+A1(J)*D
 2400      CONTINUE
 2600   CONTINUE
	RETURN
C
C**********************************************************************
C
	ENTRY DSDEBS(X0Y0,Z,CHARGE,M,A,N)
C
C        This entry computes the multipole decomposition of the
C        field of a set of charges
C         PARAMETERS
C       X0Y0      - the center of the decomposition
C       Z(I)      - the coordinates of the i-th dipole
C       CHARGE(I) - the value of the i-th dipole
C       M         - the number of dipoles in the array z
C       A         - the coefficients of the expansion
C       N         - the order of the expansion
C
C ----- Set a to zero
C
	A(I0)=ZERO
	DO 2800 I=1,N
	   A(I)=ZERO
 2800   CONTINUE
C
C ----- Compute the decomposition of the field of charges
C
	DO 3200 I=1,M
	   D=CHARGE(I)
	   Z0=Z(I)-X0Y0
	   CALL DSEXP0(Z0,A1,N)
C          A(I0)=A(I0)+A1(I0)*D
	   A(I0)=A(I0)+       D
	   DO 3000 J=1,N
	      A(J)=A(J)+A1(J)*D
C             A(J)=0
 3000      CONTINUE
 3200   CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DSLOLO(A,Z07,B,N)
C
C        This entry shifts  the origin of a multipole decomposition
C        PARAMETERS
C     A  - the original multipole decomposition
C     Z0 - the vector by which it is to be shifted
C     B  - the shifted decomposition
C     N  - the order of the decompositions a, b
C
C     IMPORTANT : The initialization point DSTTIN must be call
C                  before any other in this subroutine.
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE COMPLEX A(1),B(1),Z0,Z07,Z00,Z0P(61),Z0P2(61),CBUF(2)
	DOUBLE COMPLEX CD,CDD,Z22,Z23,ZS,CBUF2(2)
	DOUBLE COMPLEX A0,ANEW(61),IMAG,B1(1),B2(1),B3(1),B4(1)
	EQUIVALENCE (Z0P(1),CBUF(2))
	EQUIVALENCE (Z0P2(1),CBUF2(2))
CCC        DOUBLE PRECISION FACT(61),BUF(2),C(61,61),RAP(61)
ccc	DOUBLE PRECISION BUF(2),C(81,81),RAP(61)
	DOUBLE PRECISION BUF(2),C(121,121),RAP(61)
	DOUBLE PRECISION SGN ,RZ0(2),RIN(61)
CCC        EQUIVALENCE (FACT(1),BUF(2))
	EQUIVALENCE (Z0,RZ0(1))
	DATA IMAG/(0.0,1.0)/ ,DHALF/0.5/,DONE/1.0/,
     *       ZERO/0.0/, I0/0/,ICALL/0/
C        Z=Z0
ccc      QUICK FIX
	 Z0 = Z07
C
C ----- Create the array of powers of z0
C
	Z0= Z0
	Z0P(I0)=DONE
	Z00=Z0
	CD=DONE/Z0
	CDD=CD
	Z0P2(I0)=DONE

	DO 1200 I=1,N
	   Z0P(I)=Z00
	Z0P2(I) = CDD
	CDD = CDD*CD
	   Z00= Z00*Z0
 1200   CONTINUE
CCC  CALL PRIN2('INSIDE DSLOLO Z0=*',Z0,2)
CCC  CALL PRIN2('AND Z0P=*',Z0P(I0),N+1  )
C
C ----- Set b to zero and create anew
C
	B(I0)=ZERO
	DO 1400 I=1,N
	   B(I)=ZERO
	   ANEW(I) = A(I)*Z0P2(I)
 1400   CONTINUE
C
C ----- Create the array b
C
      B(I0)=A(I0)
      DO 1800 M=1,N
	 DO 1600 K=1,M
	    B(M)=B(M)+ANEW(K)*C(M,K)
 1600    CONTINUE
 1800 CONTINUE
C
C ----- Scale b
C
      DO 1900 M = 1,N
	 B(M) = (B(M) - A(I0)*RIN(M))*Z0P(M)
 1900 CONTINUE
      RETURN
C
C**********************************************************************
C
	ENTRY DSTTIN(IER)
C
C        This is the initialization entry point
C        it precomputes the tables of factorials and polynomial
C        coefficients
C
C
C ----- Create the bynomial coefficients
C
	IF (ICALL.EQ.1) GOTO 2601
	D=DONE
CCC        FACT(I0)=D
	DO 2200 I=1,61
CCC           D=D*I
CCC           FACT(I)=D
	   RIN(I)=DONE/I
	   RAP(I)= I
	   RAP(I) = RAP(I)/(I+1)
 2200   CONTINUE
CCC        DO 2600 I=1,61
CCC           DO 2400 J=1,I
CCC              C(I,J)=FACT(I-1)/(FACT(J-1)*FACT(I-J) )
CCC 2400      CONTINUE
CCC 2600   CONTINUE
	C(1,1) = 1
ccc	DO 2600 I=2,81
	DO 2600 I=2,121
	   I1 = I-1
	   C(I,1) = 1
	   DO 2400 J=2,I
	      J1 = J-1
	      C(I,J) = C(I1,J)+C(I1,J1)
 2400      CONTINUE
 2600   CONTINUE
	ICALL = 1
 2601   CONTINUE
CCC        CALL PRIN2('FACTORIALS ARE *',FACT(I0),60)
ccc        CALL PRIN2('AND C(M,K) ARE*',C,121*121)
	RETURN
C
C**********************************************************************
C
	ENTRY DSLNEW(A,Z22,Z23,B1,B2,B3,B4,N)
C
C        This entry point converts a multipole decomposition with
C        the origin at z22 into 4 expansions.
C        By recombining them appropriately, one can generate
C        the full taylor expansions at four locations:
C        Z23,-Z23,i*Z23,-i*Z23.
C        PARAMETERS
C     Z22 - The origin of the multipole expansion
C     Z23 - The origin of the taylor expansion
C     A   - The input multipole expansion
C     B   - The output taylor expansion
C     N   - The order of the expansions a,b
	Z0=Z23-Z22
C       Z0=-Z0
C
C ----- Create the array of inverse powers of z
C
	CD=DONE/Z0
	CDD=CD
	Z0P(I0)=DONE
	DO 2800 I=1,N
	   Z0P(I)=CDD
	   CDD=CDD*CD
 2800   CONTINUE
C
C -----  Set the arrays b to zero and scale the vector a
C
      DO 3200 I=1,N
	 ANEW(I) = A(I)*Z0P(I)
	 B1(I)=ZERO
	 B2(I)=ZERO
	 B3(I)=ZERO
	 B4(I)=ZERO
 3200 CONTINUE
      B1(I0)=ZERO
      B2(I0)=ZERO
      B3(I0)=ZERO
      B4(I0)=ZERO
C
C -----  Create the array b
C
	DO 3600 M1=1,N
	   M=M1-1
	   NMM = N
	   DO 3410 K = 1,NMM,4
	      KPM=K+M
	      B2(M)=B2(M)+ ANEW(K)*C(KPM,K)
 3410      CONTINUE
	   DO 3420 K = 2,NMM,4
	      KPM=K+M
	      B3(M)=B3(M)+ ANEW(K)*C(KPM,K)
 3420      CONTINUE
	   DO 3430 K = 3,NMM,4
	      KPM=K+M
	      B4(M)=B4(M)+ ANEW(K)*C(KPM,K)
 3430      CONTINUE
	   DO 3440 K = 4,NMM,4
	      KPM=K+M
	      B1(M)=B1(M)+ ANEW(K)*C(KPM,K)
 3440      CONTINUE
 3600   CONTINUE
C
C -----  Scale b
C
      DO 3710 M1=2,N,2
	 M = M1-1
	 ZS =-Z0P(M)
	 B1(M) = B1(M)*ZS
	 B2(M) = B2(M)*ZS
	 B3(M) = B3(M)*ZS
	 B4(M) = B4(M)*ZS
 3710 CONTINUE
      DO 3720 M1=1,N,2
	 M = M1-1
	 ZS = Z0P(M)
	 B1(M) = B1(M)*ZS
	 B2(M) = B2(M)*ZS
	 B3(M) = B3(M)*ZS
	 B4(M) = B4(M)*ZS
 3720 CONTINUE
      B1(I0) = B1(I0) + A(I0)*DLOG(RZ0(1)**2+RZ0(2)**2)*DHALF
CCC      B1(I0) = B1(I0) + A(I0)*DLOG(CDABS(Z0))
      DO 3701 M = 1,N,2
	 B1(M) = B1(M) + A(I0)*Z0P(M)*RIN(M)
 3701 CONTINUE
      DO 3702 M = 2,N,2
	 B1(M) = B1(M) - A(I0)*Z0P(M)*RIN(M)
 3702 CONTINUE
      RETURN
C
C**********************************************************************
C
	ENTRY DSLOTA(A,Z22,Z23,B,N)
C
C        This entry point converts a multipole decomposition with
C        the origin at z22 into a taylor expansion with the origin at
C        z23
C        PARAMETERS
C     Z22 - the origin of the multipole expansion
C     Z23 - the origin of the taylor expansion
C     A   - the input multipole expansion
C     B   - the output taylor expansion
C     N   - the order of the expansions a,b
C
	Z0=Z23-Z22
C       Z0=-Z0
C
C ----- Create the array of inverse powers of z
C
      CD=DONE/Z0
      CDD=CD
      Z0P(I0)=DONE
      DO 3800 I=1,N
	Z0P(I)=CDD
	CDD=CDD*CD
 3800 CONTINUE
C
C ----- Set the array b to zero and create anew
C
      DO 4000 I=1,N
	 ANEW(I) = A(I)*Z0P(I)
	 B(I)=ZERO
 4000 CONTINUE
      B(I0)=ZERO
C
C -----  Create the unscaled array b
C
      DO 4400 M1=1,N
	 M=M1-1
	 NMM = N
	 DO 4200 K=1,NMM
	    KPM=K+M
	    B(M)=B(M)+ ANEW(K)*C(KPM,K)
 4200    CONTINUE
 4400 CONTINUE
C
C -----  Scale b
C
      DO 4599 M1=1,N,2
	 M = M1-1
	 B(M) = B(M)*Z0P(M)
 4599 CONTINUE
      DO 4600 M1=2,N,2
	 M = M1-1
	 B(M) = -B(M)*Z0P(M)
 4600 CONTINUE
      B(I0) = B(I0) + A(I0)*DLOG(RZ0(1)**2+RZ0(2)**2)*DHALF
CCC      B(I0) = B(I0) + A(I0)*DLOG(CDABS(Z0))
      DO 4601 M = 1,N,2
	 B(M) = B(M) + A(I0)*Z0P(M)*RIN(M)
 4601 CONTINUE
      DO 4602 M = 2,N,2
	 B(M) = B(M) - A(I0)*Z0P(M)*RIN(M)
 4602 CONTINUE
      RETURN
C
C**********************************************************************
C
	ENTRY DSLOT1(A0,Z22,Z23,B,N)
C
C        This entry point converts a monopole decomposition with
C        the origin at z22 into a taylor expansion with the origin at
C        Z23 (For fields and potentials only)
C        PARAMETERS
C     Z22 - the origin of the decomposition (position of the particle)
C     Z23 - the origin of the taylor expansion
C     A0  - the input monopole coefficient
C     B   - the output taylor expansion
C     N   - the order of the expansions b
C
       Z0=Z23-Z22
       Z0=-Z0
C
C ----- Create the array of inverse powers of z
C
      CD = DONE/Z0
      B(I0) =  A0*DLOG(RZ0(1)**2+RZ0(2)**2)*DHALF
      B(1)  =  -A0*CD
      DO 3810 I=1,N
	 I1 = I+1
	 B(I1)=B(I)*CD*RAP(I)
 3810 CONTINUE
      RETURN
C
C**********************************************************************
C
	ENTRY DSLOT3(A0,Z22,Z23,B,N)
C
C        This entry point converts a monopole decomposition with
C        The origin at z22 into a taylor expansion with the origin at
C        Z23 (For fields only)
C        PARAMETERS
C     Z22 - the origin of the decomposition (position of the particle)
C     Z23 - the origin of the taylor expansion
C     A0  - the input monopole coefficient
C     B   - the output taylor expansion
C     N   - the order of the expansions b
C
       Z0=Z23-Z22
	Z0 = -Z0
C
C ----- Create the array of inverse powers of z
C
      CD = DONE/Z0
      B(I0) = -A0*CD
      DO 3811 I=1,N
	 I1 = I-1
	 B(I)=B(I1)*CD
 3811 CONTINUE
      RETURN
C
C**********************************************************************
C
	ENTRY DSTATA(A,Z07,B,N)
C
C        This entry shifts the origin of a taylor expansion
C       PARAMETERS
C       A  - the initial taylor expansion
C       Z0 - the vector by which it is to be shifted
C       B  - the shifted expansion
C       N  - the order of both expansions
C
C
CCC     QUICK FIX
	Z0 = Z07
C
C -----  Create the powers of z0
C
	Z00=Z0
	Z0P(I0)=DONE
	CD = DONE/Z0
	CDD = CD
	Z0P2(I0) = DONE
	DO 5200 I=1,N
	   Z0P(I)=Z00
	   Z00= Z00*Z0
	   Z0P2(I) = CDD
	   CDD = CDD*CD
5200    CONTINUE
C
C -----  Set b to zero
C
	DO 5400 I =1,N
	   B(I) = ZERO
	   ANEW(I) = A(I)*Z0P(I)
5400    CONTINUE
	B(I0) = ZERO
	ANEW(I0) = A(I0)
C
C -----  Create b
C
	N1=N+1
	DO 5800 K1=1,N1
	   K=K1-1
	   DO 5600 M1=1,K1
	      M=M1-1
	      B(M)= B(M)+ANEW(K)*C(K1,M1)
 5600      CONTINUE
 5800   CONTINUE
C
C -----  Scale b
C
	DO 5900 M = 1,N
	   B(M) = B(M)*Z0P2(M)
5900    CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DSEXP0(Z07,A,N)
C
C        This entry computes the coefficients of the multipole
C        Decomposition of the field of a unity charge (the origin
C        Of the expansion is located at zero)
C        PARAMETERS
C      Z0 - the location of the charge
C      A  - the coefficients of the expansion
C      N  - the order of the expansion
C
C      IMPORTANT : the initialization entry point DSEXPI
C                  must be called before any other in this
C                  subroutine
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE COMPLEX Z0,Z07,A(1),Z00,BETA
	DOUBLE PRECISION RK1(61)
	DATA I0/0/ ,DONE/1.0/,DHALF/0.5/,ZERO/0.0/
C
C -----  Create the decomposition of the charge at the point z0
C
CCC     QUICK FIX
	Z0 = Z07
	A(I0)=DONE
C       Z00=DONE
	Z00=+Z0
	DO 1200 I=1,N
	A(I)=Z00*RK1(I)
	Z00=Z00*Z0
 1200 CONTINUE
	RETURN
C
C**********************************************************************
C
	ENTRY DSEXPI(IER)
C
C        This is an initialization entry point
C
C -----  Create the array rk1
C
	DO 1400 I=1,61
C          D=(-1)**(I-1)
	   D=-1
	   RK1(I)=D/I
 1400   CONTINUE
	RETURN
C
C**********************************************************************
C
	ENTRY DSEXP1(Z07,BETA,A,N)
C
C
C        this entry computes the multipole decomposition of the
C        field of a dipole (the origin of the decomposition is at zero)
C        PARAMETERS
C     Z0   - the location of the dipole
C     BETA - the orientation of the dipole
C     A    - the coefficients of the decomposition
C     N    - the order of the decomposition
C
C -----  Create the decomposition of the unit dipole at the point
C        z0 oriented in the direction beta
C
CCC   QUICK FIX
	Z0 = Z07
C
	A(I0)=ZERO
	Z00=-   BETA
	DO 1600 I=1,N
	   A(I)=Z00
	   Z00=Z00*Z0
 1600   CONTINUE
	RETURN
C
C**********************************************************************
C
	ENTRY DSEXDF(Z07,A,N)
C
C        This entry computes the coefficients of the multipole
C        decomposition of the derivative of a unity charge (the origin
C        of the expansion is located at zero)
C        PARAMETERS
C      Z0 - the location of the charge
C      A  - the coefficients of the expansion
C      N  - the order of the expansion
C
C
C ----- Create the decomposition of the derivative of the charge
C        located at z0
CCC   QUICK FIX
      Z0 = Z07
C
	A(I0)=ZERO
	Z00=DONE
	DO 1800 I=1,N
	   A(I)=Z00
	   Z00=Z00*Z0
 1800   CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DSTAYL(A,Z07,Z2,N,OUT)
C
C        This entry computes the value of a taylor expansion
C        parameters
C     A   - the coefficients of  the expansion
C     Z0  - the origin of the expansion
C     Z2  - the point at which the value is to be computed
C     N   - the order of the expansion
C     OUT - the computed complex value
C
C     IMPORTANT : The initialization point DSTABI must be
C                 call before any other in this subroutine
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE COMPLEX A(1),Z,CD,Z1,Z11,CD1,OUT,Z0,Z07,Z2
	DOUBLE PRECISION RNS(61)
	DOUBLE PRECISION    RD(2),RD1(2),D
	DATA I0/0/ ,DONE/1.0/,DHALF/0.5/
	EQUIVALENCE (RD(1),CD), (RD1(1),CD1)
C
C -----  Compute the polynomial
C
CCC     QUICK FIX
	Z0 = Z07
C
	Z=Z2-Z0
	Z1=A(N)
	DO 1200 I1=1,N
	   I=N-I1
	   Z1=Z1*Z+A(I)
 1200   CONTINUE
	OUT=Z1
	RETURN
C
C**********************************************************************
C
	ENTRY DSTAYD(A,Z07,Z2,N,OUT)
C
C        This entry computes the value of the derivative of
C         a taylor expansion
C        PARAMETERS
C     A   - the coefficients of  the expansion
C     Z0  - the origin of the expansion
C     Z2  - the point at which the value is to be computed
C     N   - the order of the expansion
C     OUT - the computed complex value
C
C
C -----   Compute the polynomial
C
CCC     QUICK FIX
	Z0 = Z07
C
	Z=Z2-Z0
	Z1=A(N)*RNS(N)
	N1=N-1
	DO 2200 I1=1,N1
	   I=N-I1
	   Z1=Z1*Z+A(I)*RNS(I)
 2200   CONTINUE
	OUT=Z1
	RETURN
C
C**********************************************************************
C
	ENTRY  DSTABI(IER)
C
C  ---- this is an initialization entry point
C
	DO 111 I=1,61
	   RNS(I) = I
 111    CONTINUE
	RETURN
C
C**********************************************************************
C
	 ENTRY DSLOR0(A,Z07,Z2,N,OUT)
C
C        this entry computes the value of a multipole expansion
C        with a zero order term (logarithmic)
C        PARAMETERS
C     A   - the coefficients of  the expansion
C     Z0  - the origin of the expansion
C     Z2  - the point at which the value is to be computed
C     N   - the order of the expansion
C     OUT - the computed complex value
C
C
CCC   QUICK FIX
	Z0 = Z07
C
C -----  Compute the inverse polynomial with the zero term
C
	Z=Z2-Z0
	Z1=A(N)
	Z11=DONE/Z
	N1=N-1
	DO 1600 I1=1,N1
	   I=N-I1
	   Z1=Z1*Z11+A(I)
 1600   CONTINUE
	Z1=Z1*Z11
	CD1=A(I0)
	CD=Z
	D=RD(1)**2+RD(2)**2
	OUT=   DLOG(D)*RD1(1)*DHALF
	OUT=Z1+OUT
CCC         OUT = Z1 + CDLOG(Z)*CD1
	 RETURN
C
C**********************************************************************
C
	 ENTRY DSLORD(A,Z07,Z2,N,OUT)
C
c        This entry computes the value of the derivative of a multipole
C        expansion with a zero order term (logarithmic)
C        PARAMETERS
C     A   - the coefficients of  the expansion
C     Z0  - the origin of the expansion
C     Z2  - the point at which the value is to be computed
C     N   - the order of the expansion
C     OUT - the computed complex value
C
C
CCC   QUICK FIX
	Z0 = Z07
C
C ------ Compute the inverse polynomial with the zero term
C
	Z=Z2-Z0
	Z1=A(N)*RNS(N)
	Z11=DONE/Z
	N1=N-1
	DO 1660 I1=1,N1
	  I=N-I1
	  Z1=Z1*Z11+A(I)*RNS(I)
 1660   CONTINUE
CCC     Z1=Z1*Z11
	Z1=Z1*Z11*Z11
	CD1=A(I0)
	OUT = CD1*Z11
	OUT=OUT-Z1
	RETURN
C
C**********************************************************************
C
	 ENTRY DSLOR1(A,Z07,Z2,N,OUT)
C
C        This entry computes the value of a multipole expansion
C        without a zero order term
C        PARAMETERS
C     A   - the coefficients of  the expansion
C     Z0  - the origin of the expansion
C     Z2  - the point at which the value is to be computed
C     N   - the order of the expansion
C     OUT - the computed complex value
C
CCC   QUICK FIX
      Z0 = Z07
C
C -----  Compute the inverse polynomial without the zero term
C
	Z=Z2-Z0
	Z1=A(N)
	Z11=DONE/Z
	N1=N-1
	DO 1800 I1=1,N1
	   I=N-I1
	   Z1=Z1*Z11+A(I)
 1800   CONTINUE
	OUT=Z1*Z11
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DSLOR2(Z1,A1,Z2,A2,Z3,A3,NEXP)
C
C        This subroutine combines two multipole expansions producing
C        a global one which is equal to the sum of the two
C        PARAMETERS
C     Z1   - the origin of the first input expansion
C     A1   - the coefficients of the first input expansion
C     Z2   - the origin of the second input expansion
C     A2   - the coefficients of the second input expansion
C     Z3   - the origin of the output expansion
C     NEXP - the order of all three expansions
C
C**********************************************************************
C
	DOUBLE COMPLEX Z1,Z2,Z3,Z0,A1(1),A2(1),A3(1),B1(61),
     1             B2(61),CB1(2),CB2(2)
	DATA I0/0/
	EQUIVALENCE (B1(1),CB1(2)), (B2(1),CB2(2))
C
C -----  Recompute the first expansion
C
	Z0=Z1-Z3
	CALL DSLOLO(A1,Z0,B1,NEXP)
C
C -----  Recompute the second expansion
C
	Z0=Z2-Z3
	CALL DSLOLO(A2,Z0,B2,NEXP)
C
C -----  Compute the global expansion
C
	A3(I0)=B1(I0)+B2(I0)
	DO 1400 I=1,NEXP
	   A3(I)=B1(I)+B2(I)
 1400   CONTINUE
	RETURN
C
C**********************************************************************
C
	ENTRY DSTAY2(Z1,A1,Z2,A2,Z3,A3,NEXP)
C
C        This entry from one input taylor expansion produces two
C        output expansions each of which is equal to the initial
C        expansion in its domain of validity
C        PARAMETERS
C     Z1   - the origin of the input expansion
C     A1   - the coefficients of the first expansion
C     Z2   - the origin of the first output expansion
C     A2   - the coefficients of the first output expansion
C     Z3   - the origin og the second output expansion
C     A3   - the coefficients of the second output expansion
C     NEXP - the order of all three expansions
C
C
C -----  Compute the first of new expansions
C
	Z0=Z2-Z1
	CALL DSTATA(A1,Z0,A2,NEXP)
C
C -----  Create the second expansion
C
	 Z0=Z3-Z1
	CALL DSTATA(A1,Z0,A3,NEXP)
	RETURN
	END
C
C**********************************************************************
C
C                     JCG  - 25 Jan 87
C
C**********************************************************************
C
C################   END CODE : FAST MULTIPOLE METHOD ##################
