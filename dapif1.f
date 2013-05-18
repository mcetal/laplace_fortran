C
C####################  BEGIN CODE : NEAR APPROACH  ####################
C
C**********************************************************************
C
	SUBROUTINE CLOSE0(EPS,I1,I2,Q1,Q2,X1,Y1,X2,Y2,
     *                   FI1,FI2,POT1,POT2)
C
C   *** DESCRIPTION :
C
C       Close approch subroutine provided by the user.
C   Called when two particles are closer to each other than EPS.
C
C   *** INPUT PARAMETERS :
C
C   EPS     =  close approach distance
C   I1      =  number of first  particle
C   I2      =  number of second particle
C   Q1      =  charge of first  particle
C   Q2      =  charge of second particle
C   X1      =  x-component of position of first  particle
C   Y1      =  y-component of position of first  particle
C   X2      =  x-component of position of second particle
C   Y2      =  y-component of position of second particle
C
C   *** OUTPUT PARAMETERS :
C
C   FI1     =  field contribution from first particle on second
C   FI2     =  field contribution from second particle on first
C   POT1    =  potential contribution from first particle on second
C   POT2    =  potential contribution from second particle on first
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER II(2)
	DOUBLE PRECISION POS1(2),POS2(2)
	DOUBLE COMPLEX FI1,FI2,OUT2,Q1,Q2
	DATA DONE/1.0/,DHALF/0.5/
	II(1) = I1
	II(2) = I2
	POS1(1) = X1
	POS1(2) = Y1
	POS2(1) = X2
	POS2(2) = Y2
	FI1 = 0
	FI2 = 0
CCC     POT1 = 0
CCC     POT2 = 0
	POT1=DLOG( (X2-X1)**2 + (Y2-Y1)**2 ) /2
	POT2=POT1
ccc	   CALL PRINF('CLOSE PARTICLES ARE:*',II,2)
ccc	   CALL PRIN2(' COORDINATES OF FIRST  PART=*',POS1,2)
ccc	   CALL PRIN2(' COORDINATES OF SECOND PART=*',POS2,2)
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE CLOSE1(EPS,J,QJ,X,Y,XJ,YJ,FIP,POP)
C
C   *** DESCRIPTION :
C
C       Close approch subroutine for arbitrary point provided by
C   the user. Called when a particle is within EPS to arb. point.
C   of a particle.
C
C   *** INPUT PARAMETERS :
C
C   EPS     =  close approach distance
C   J       =  number of particle
C   QJ      =  charge of particle
C   X       =  x-component of position of point
C   Y       =  y-component of position of point
C   XJ      =  x-component of position of particle
C   YJ      =  y-component of position of particle
C
C   *** OUTPUT PARAMETERS :
C
C   FIP     =  field contribution from particle on point
C   POP     =  field contribution from particle on point
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION POSAP(2),POSP(2)
	DOUBLE COMPLEX FIP,QJ
	DATA DONE/1.0/,DHALF/0.5/
	   FIP = 0
	   POP = 0
	   CALL PRINF(' ARBITRARY POINT CLOSE TO PARTICLE :*',J,1)
	   CALL PRIN2(' COORDINATES OF  PARTICLE=*',POSP,2)
	   CALL PRIN2(' COORDINATES OF  ARB. POINT=*',POSAP,2)
	RETURN
	END
C
C####################  END CODE : NEAR APPROACH  ######################
C####################  BEGIN CODE : DIRECT METHOD  ####################
C
C**********************************************************************
C
	SUBROUTINE DADIRE (IFLAG7,NAT,IMAX,XA,YA,QA,POT,
     *                       FIE,EPS7,WKSP,NSP,CLOSE)
C
C   *** DESCRIPTION:
C
C       This subroutine computes the potentials and electric
C   fields by direct calculation in two dimensions for a subset
C   of the NAT charged particles of size IMAX .
C
C
C   *** OTHER ENTRY THAT CAN BE CALLED:
C
C       DADIRP computes the fields and/or potentials at
C       an arbitrary point in the computationnal cell.
C       IMPORTANT: DADIRE must be called before DADIRP.
C
C
C   *** INPUT PARAMETERS:
C
C   IFLAG7     = 1 for computing fields only
C                2 for computing potential and electric fields
C                3 Hilbert Matrix
C   NAT       = number of particles
C   XA(i)     = first coordinate of a particle
C   YA(i)     = second coordinate of a particle
C   QA(i)     = charge of a particle
C   EPS7      = near approach distance.
C               this number is provided by the user to check for very
C               near approaches. If two particles are within EPS7 of
C               each other, the program executes a near approach
C               subroutine CLOSE which the user must declare as
C               EXTERNAL in the calling program.
C   WKSP      = workspace
C   NSP       = length of workspace
C   CLOSE     = is the near approach subroutine to be provided by the
C               user. The calling sequence should be
C               CLOSE(EPS,I1,I2,Q1,Q2,X1,X2,Y1,Y2,FI1,FI2,PO1,PO2)
C               which is called when two particles are within EPS of
C               each other.
C               In many applications, when particles are very close,
C               the nature of the force law changes. For example,
C               the charges may be viewed as continuous distributions
C               about their centers instead of as point charges.
C
C   *** OUTPUT PARAMETERS:
C
C   FIE(*,I) = net field at location of ith particle
C   POT(I)   = potential at location of ith particle
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION WKSP(1),XA(1),YA(1),POT(1),POD,POP
	DOUBLE COMPLEX QA(1),FIE(1),ZI,ZJ,Z2,Z1,ZAP,ZP,FID,FIP,FI1,FI2
	DOUBLE COMPLEX QTOT,CZERO
	DATA DONE/1.0/,BOXSIZ/128.0/,DZERO/0.0/
	DATA CZERO/(0.0,0.0)/
C
C ----- Set potentials and fields to zero
C       before direct calculation
C
	IF ((IFLAG7.EQ.1).OR.(IFLAG7.EQ.3)) IFLAG=1
	IF (IFLAG7.EQ.2) IFLAG = 2
	IF (IFLAG .EQ. 1) THEN
	   DO 800 I=1,IMAX
	      FIE(I) = CZERO
 800       CONTINUE
	ELSE
	   DO 900 I=1,IMAX
	      POT(I) = DZERO
	      FIE(I) = CZERO
 900       CONTINUE
	END IF
C
C ----- Center and scale the positions to BOXSIZ (128.)
C
	NXA = 1
	NYA = NXA + NAT
	NEWADR = NYA + NAT
	NEWSP = NSP - NEWADR + 1
	RSCAL = 0.0
	RMINX =  1.0E30
	RMAXX = -1.0E30
	RMINY =  1.0E30
	RMAXY = -1.0E30
	DO 100 I=1,NAT
	   IF (XA(I) .GT. RMAXX) RMAXX = XA(I)
	   IF (YA(I) .GT. RMAXY) RMAXY = YA(I)
	   IF (XA(I) .LT. RMINX) RMINX = XA(I)
	   IF (YA(I) .LT. RMINY) RMINY = YA(I)
 100    CONTINUE
	RSX = RMAXX-RMINX
	RSY = RMAXY-RMINY
	XONEW = (RMAXX+RMINX)/2.0
	YONEW = (RMAXY+RMINY)/2.0
	RSCAL = DMAX1(RSX,RSY)
	RSCAL = BOXSIZ/RSCAL*0.99
	EPS = EPS7 * RSCAL
CCC        EPS = EPS7
	QTOT = 0.0
	DO 110 I=1,NAT
	   I1 = I-1
	   WKSP(NXA+I1) = (XA(I)-XONEW) *  RSCAL
	   WKSP(NYA+I1) = (YA(I)-YONEW) *  RSCAL
	   QTOT  = QTOT + QA(I)
 110    CONTINUE
	CALL DADIRF (IFLAG,NAT,IMAX,QTOT,RSCAL,XA,YA,
     *        WKSP(NXA),WKSP(NYA),QA,POT,FIE,EPS,CLOSE)
	DO 3000 I=1,IMAX
	   FIE(I) = DCONJG(FIE(I))
 3000   CONTINUE
C
C ----- Rescale potentials and fields
C
	COEF = RSCAL
	COEP = DLOG(RSCAL)
	IF (IFLAG7 .EQ. 1)  THEN
	   DO 3200 I=1,IMAX
	      FIE(I) = FIE(I)*COEF
 3200      CONTINUE
	END IF
	IF (IFLAG7 .EQ. 2)  THEN
	   DO 3210 I=1,IMAX
	      FIE(I) = FIE(I)*COEF
	      POT(I) = POT(I)-(QTOT-QA(I))*COEP
 3210       CONTINUE
	END IF
	IF (IFLAG7 .EQ. 3)  THEN
	   DO 3220 I=1,IMAX
	      FIE(I) = DCONJG(FIE(I))*COEF
 3220      CONTINUE
	END IF
	RETURN
C
C**********************************************************************
C
	ENTRY DADIRP (IFLAG7,NAT,XAP,YAP,XA,YA,QA,POD,FID,EPS7,CLOSEP)
C
C   *** DESCRIPTION:
C
C       This subroutine computes the potentials and electric
C   fields by direct calculation in two dimensions for an
C   arbitrary point in the computationnal cell.
C
C
C   *** OTHER ENTRY THAT CAN BE CALLED:
C
C       DADIRP computes the fields and/or potentials at
C       an arbitrary point in the computationnal cell.
C       IMPORTANT: DADIRE must be called before DADIRP.
C
C
C   *** INPUT PARAMETERS:
C
C   IFLAG7     = 1 for computing fields only
C                2 for computing potential and electric fields
C                 3 for Hilbert Matrix
C   NAT       = number of particles
C   XAP       = first coordinate of point
C   YAP       = second coordinate of point
C   QA(i)     = charge of a particle
C   EPS       = near approach distance.
C               this number is provided by the user to check for very
C               near approaches. If the point and a particle are within
C               EPS of each other, the program executes a near approach
C               subroutine CLOSEP which the user must declare as
C               EXTERNAL in the calling program.
C   CLOSEP    = is the near approach subroutine to be provided
C               by the user. The calling sequence should be
C                 CALL CLOSEP(EPS,I,QI,XAP,YAP,XI,YI,FI,PO)
C               In many applications, when particles are very close,
C               the nature of the force law changes. For example,
C               the charges may be viewed as continuous distributions
C               about their centers instead of as point charges.
C
C   *** OUTPUT PARAMETERS:
C
C   FID      = net field at location of point
C   POD      = potential at location of point
C
C**********************************************************************
C
	IF ((IFLAG7.EQ.3).OR.(IFLAG7.EQ.1)) IFLAG = 1
	IF (IFLAG7.EQ.2) IFLAG = 2
C
C ----- Scale particles positions
C
	XAP1 = (XAP-XONEW) * RSCAL
	YAP1 = (YAP-YONEW) * RSCAL
	EPS = EPS7 * RSCAL
	CALL DADIRQ (IFLAG,NAT,WKSP(NXA),WKSP(NYA),XAP,YAP,XAP1,YAP1,
     *               XA,YA,QA,POD,FID,EPS,RSCAL,QTOT,CLOSEP)
	FID = DCONJG(FID)
C
C ----- Rescale potentials and fields
C
	IF (IFLAG7 .EQ. 1)  THEN
	   FID = FID * COEF
	END IF
	IF (IFLAG7 .EQ. 2)  THEN
	   FID = FID * COEF
	   POD = POD - QTOT * COEP
	END IF
	IF (IFLAG7 .EQ. 3)  THEN
	   FID = DCONJG(FID) * COEF
	END IF
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DADIRF (IFLAG,NAT,IMAX,QTOT,RSCAL,XAU,YAU,XA,
     *                      YA,QA,POT,FIE,EPS,CLOSE)
C
C       See information in DADIRE
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION XAU(1),YAU(1),XA(1),YA(1),POT(1),POD,POP
	COMPLEX *16 QTOT,QA(1),FIE(1),ZI,ZJ,Z2,Z1,ZAP,ZP,FID,FIP,FI1,FI2
	DATA DONE/1.0/,BOXSIZ/128.0/
	COEF = RSCAL
	COEP = DLOG(RSCAL)
	EPS1 = EPS/RSCAL
	EPS2 = EPS**2
C
C ----- Compute the potential and or field directly
C
	IF (IFLAG .EQ. 1)  THEN
	   DO 2001 I=1,IMAX
	      DO 1001 J=I+1,NAT
		 ZI = DCMPLX(XA(I),YA(I))
		 ZJ = DCMPLX(XA(J),YA(J))
		 R = (XA(I)-XA(J))**2+(YA(I)-YA(J))**2
		 IF (R .LE. EPS2) THEN
		    CALL CLOSE(EPS1,I,J,QA(I),QA(J),XAU(I),YAU(I),
     *                         XAU(J),YAU(J),FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
		    FI1 = FI1/COEF
		    FI2 = FI2/COEF
		    FIE(I) = FIE(I)+FI2
		    IF(J.LE.IMAX) FIE(J) = FIE(J)-FI1
		 ELSE
		    Z2 = DONE/(ZI-ZJ)
		    FIE(I)=FIE(I)+QA(J)*Z2
		    IF(J.LE.IMAX) FIE(J)=FIE(J)-QA(I)*Z2
		 END IF
 1001         CONTINUE
 2001      CONTINUE
	ELSE
	   DO 2002 I=1,IMAX
	      DO 1002 J=I+1,NAT
		 ZI = DCMPLX(XA(I),YA(I))
		 ZJ = DCMPLX(XA(J),YA(J))
		 R = (XA(I)-XA(J))**2+(YA(I)-YA(J))**2
		 IF (R .LE. EPS2) THEN
		    CALL CLOSE(EPS1,I,J,QA(I),QA(J),XAU(I),YAU(I),
     *                         XAU(J),YAU(J),FI1,FI2,POT1,POT2)
C
C ----- Scale result from close subroutine
C
		      FI1 = FI1/COEF
		      FI2 = FI2/COEF
		      FIE(I) = FIE(I)+FI2
		      IF(J.LE.IMAX) FIE(J) = FIE(J)-FI1
		      POT1 = POT1 + (QTOT-QA(J))*COEP
		      POT2 = POT2 + (QTOT-QA(I))*COEP
		      POT(I)=POT(I)+POT2
		      IF(J.LE.IMAX) POT(J)=POT(J)+POT1
		 ELSE
		    Z2 = DONE/(ZI-ZJ)
		    Z1 = CDLOG(ZI-ZJ)
		    FIE(I)=FIE(I)+QA(J)*Z2
		    IF(J.LE.IMAX) FIE(J)=FIE(J)-QA(I)*Z2
		    POT(I)=POT(I)+QA(J)*Z1
		    IF(J.LE.IMAX) POT(J)=POT(J)+QA(I)*Z1
		 END IF
 1002         CONTINUE
 2002      CONTINUE
	END IF
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DADIRQ (IFLAG,NAT,XA,YA,XAPU,YAPU,XAP,YAP,XAU,YAU,
     *               QA,POD,FID,EPS,RSCAL,QTOT,CLOSEP)
C
C       See entry DADIRP for information
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION XA(1),YA(1),XAU(1),YAU(1)
	DOUBLE COMPLEX QTOT,QA(1),Z2,Z1,ZAP,ZP,FID,FIP
	DATA DONE/1.0/
	COEF = RSCAL
	COEP = DLOG(RSCAL)
	EPS1 = EPS/RSCAL
	EPS2 = EPS**2
C
C ----- Compute the potential and or field directly
C
	ZAP = DCMPLX(XAP,YAP)
	IF (IFLAG .EQ. 1)  THEN
	   DO 4001 I=1,NAT
	      ZP = DCMPLX(XA(I),YA(I))
	      Z2 = ZAP-ZP
	      RR = CDABS(Z2)**2
	      IF (RR.LE.EPS2) THEN
	      CALL CLOSEP(EPS1,I,QA(I),XAPU,YAPU,XAU(I),YAU(I),FIP,POP)
C
C ----- Scale result from close subroutine
C
		 FIP = FIP/COEF
		 FID = FID + FIP
	      ELSE
		 Z2 = DONE/Z2
		 FID = FID + QA(I)*Z2
	      END IF
 4001      CONTINUE
	ELSE
	   DO 5002 I=1,NAT
	      ZP = DCMPLX(XA(I),YA(I))
	      Z2 = ZAP-ZP
	      RR = CDABS(Z2)**2
	      IF (RR.LT.EPS2) THEN
		 CALL CLOSEP(EPS1,I,QA(I),XAP,YAP,XA(I),YA(I),FIP,POP)
C
C ----- Scale result from close subroutine
C
		 FIP = FIP/COEF
		 POP = POP + QTOT*COEP
		 FID = FID + FIP
		 POD = POD + POP
	      ELSE
		 Z2 = DONE/Z2
		 Z1 = CDLOG(ZAP-ZP)
		 FID = FID + QA(I)*Z2
		 POD = POD + QA(I)*Z1
	      END IF
 5002      CONTINUE
	END IF
	RETURN
	END
C
C####################  END CODE : DIRECT METHOD  ######################
C################   BEGIN CODE : FAST MULTIPOLE METHOD ################
C
C**********************************************************************
C
	SUBROUTINE DAPIF2 (IOUT,IFLAG7,NAT,NAPB,NINIRE,MEX,IERR,INFORM,
     *                    TOL,EPS7,XA,YA,QA,POTEN,FIELD,WKSP,NSP,CLOSE)
C
C**********************************************************************
C
C   *** INFORMATION:
C
C       Fast multipole method.
C
C
C       A detailed descritpion is given in:
C       `A Fast Adaptive Algorithm for Particle Simulation'.
C       J. Carrier, L. Greengard, V. Rokhlin.
C       Yale U. Comp. Sci. Dept. RR # 496. Sep.86, Revised Jan.87
C   _____________________________________________________________
C   |-----------------------------------------------------------|
C   |    WARNING:                                               |
C   |                                                           |
C   |    The numbers of the lists given in the present          |
C   |    description (comment of the program)  differ from      |
C   |    those they have in the reference above. Following      |
C   |    is a correspondance table.                             |
C   |                                                           |
C   |             -------------------------------               |
C   |             | Program  |  YALEU CS RR#496 |               |
C   |             |----------|------------------|               |
C   |             | LIST 1   |  LIST 2          |               |
C   |             | LIST 2   |  Colleagues      |               |
C   |             | LIST 3   |  LIST 1          |               |
C   |             | LIST 4   |  LIST 3          |               |
C   |             -------------------------------               |
C   |                                                           |
C   |    List 4 and 5 of the RR # 496 were introduced to        |
C   |    help the description of the algorithm, but are         |
C   |    not used in the actual implementation of the           |
C   |    method.                                                |
C   |-----------------------------------------------------------|
C   -------------------------------------------------------------
C
C       See also the abreviated manual for more information
C       a detailed description and the precautions to be taken
C       while using this subroutine.
C
C   *** DESCRIPTION:
C
C
C       This subroutine computes the potentials and electric fields
C   due to all pairwise interactions in a system of NAT charged
C   particles in two dimensions. The vector of fields is also the
C   conjugate of the result of the product of the Hilbert matrix
C   of positions by the vector of charges.
C
C   *** OTHER ENTRY THAT CAN BE CALLED:
C
C       DAPIP2  computes the fields and/or potentials at
C       an arbitrary point in the computationnal cell.
C       IMPORTANT: DAPIF2  must be called before DAPIP2 .
C
C
C   *** INPUT PARAMETERS:
C
C   IOUT       = the numbers of the output files
C   IFLAG7     = 1 for computing fields only
C                2 for computing potential and electric fields
C                3 for Hilbert matrix
C   NAT        = number of particles
C   NAPB       = maximum number of particles in a box
C   NINIRE     = number of integer in a real number
C                1 for single precision
C                2 for double precision
C   MEX        = maximum decimal exponent allowed by the
C                computer system. Exemples:
C                 38 on a Vax in double precision
C                 76 on IBM in double precision
C                307 on a Vax in g-floating
C
C   TOL        = desired precision of calculation
C                The precision which can be achieved by a
C                double-precision version of this program can
C                not be set higher than about 1.0D-14. If too
C                high a precision is requested, the subroutine
C                returns with an error code described below.
C   EPS7       = near approach distance.
C                this number is provided by the user to check for very
C                near approaches. If two particles are within EPS of
C                each other, the program executes a near approach
C                subroutine CLOSE which the user must declare as
C                EXTERNAL in the calling program.
C   XA(i)      = first coordinate of i-th particle
C   YA(i)      = second coordinate of i-th particle
C   QA(i)      = charge of i-th particle
C   _____________________________________________________________
C   |-----------------------------------------------------------|
C   |    WARNING:                                               |
C   |    For fields and/or potentials the charges qa(i) have    |
C   |    to be real numbers, for the hilbert matrix vector      |
C   |    product qa(i) can be complex.                          |
C   |    This subroutine returns potentials  and  the product   |
C   |    of the matrix a = [aij] = [1/(zi-zj)] with the vector  |
C   |    qa. This product is also the conjugate of the fields.  |
C   |___________________________________________________________|
C   -------------------------------------------------------------
C   WKSP       = workspace array
C               the amount of workspace required is approximately:
C
C                   16*NAT+1000
C
C   NSP       = length of work array
C   CLOSE     = is the near approach subroutine to be provided by the
C               user. The calling sequence should be
C               CLOSE(EPS,I1,I2,Q1,Q2,X1,Y1,X2,Y2,FI1,FI2,PO1,PO2)
C               which is called when two particles are within EPS of
C               each other.
C               In many applications, when particles are very close,
C               the nature of the force law changes. For example,
C               the charges may be viewed as continuous distributions
C               about their centers instead of as point charges.
C
C   *** OUTPUT PARAMETERS:
C
C   INFORM   = information returned to user .
C
C            INFORM(1) is the number of terms used in the
C               multipole expansions to achieve desired accuracy.
C            INFORM(2) is the total number of boxes used in the
C               calculations.
C            INFORM(3) is the number of grid levels that actually
C               refine the computational cell.
C            INFORM(4) is the amount of workspace actually used
C               in the calculations.
C   FIELD(*,I) = net field at location of ith particle
C                or ith component of Hilbert matrix vector product.
C   POTEN(I)   = potential at location of ith particle
C   IERR       = error code array
C
C              IERR(1) = 0  no errors encountered
C              IERR(1) = 1  Size of the smallest boxes smaller than
C                           near approach distance EPS.
C              IERR(1) = 4  TOL set too high. Too many terms
C                           in multipole expansions needed.
C                           IERR(2) returns the number of terms
C                           in expansion required to satisfy TOL,
C                           which should not exceed 60.
C              IERR(1) = 8  Insufficient workspace allotted.
C                          IERR(2) returns amount of workspace needed
C              IERR(1) = 16 Total number of particles less than maximum
C                          number of particles allowed in a box.
C              IERR(1) = 32 the total number of grid levels is too high.
C                        This number L depends on the required accuracy TOL
C                        and the maximum capacity of the machine M by
C                        the following formula:
C                        N(L-6.5)+1 = log_2(M)
C                        Where  N= int(log_2|TOL|+3)
C
C
C   *** LOCAL VARIABLES:
C
C   N     = number of terms used in expansions
C
C   *** SUBROUTINES CALLED:
C
C   DAEOUT
C   DAASSI
C   DACOMP
C   DAIJAT
C   DAINDB
C   DALI12
C   DALI34
C   DAUPWD
C   DAFLEX
C   DABARR
C   DADNWD
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER IOUT(1),IERR(1),INFORM(1),NL(50),XB,YB,SHIFT,LO,TA
	DOUBLE PRECISION WKSP(1),XA(1),YA(1),POTEN(1)
	DOUBLE COMPLEX FIELD(1),QA(1),FI,QTOT,CZERO
	INTEGER *4 ZEROES(10)
	REAL *8 INZERO
	EQUIVALENCE (ZEROES(1),INZERO)
	DATA BOXSIZ/128.0/,DZERO/0.0/
	DATA CZERO/(0.0,0.0)/
C
C ----- Print out parameters
C
	IOUT1 = IOUT(1)
	IOUT2 = IOUT(2)
	CALL PRINI(IOUT1,IOUT2)
ccc	CALL PRINF('IN THIS CALCULATION IFLAG =*',IFLAG7,1)
ccc	IF (IFLAG7.EQ.1)
ccc     *     CALL PRINF(' FIELD CALCULATION *',I,0)
ccc	IF (IFLAG7.EQ.2)
ccc     *     CALL PRINF(' POTENTIAL AND FIELD CALCULATION*',I,0)
ccc	IF (IFLAG7.EQ.3)
ccc     *     CALL PRINF(' HILBERT MATRIX - VECTOR PRODUCT *',I,0)
ccc	CALL PRINF('   *',I,0)
ccc	CALL PRINF('MAXIMUM NUMBER OF PARTICLES PER BOX =*',NAPB,1)
ccc	CALL PRINF('   *',I,0)
ccc	CALL PRIN2('REQUIRED PRECISION =*',TOL,1)
ccc	call PRIN2('NEAR APPROACH DISTANCE =*',EPS7,1)
ccc	CALL PRINF('   *',I,0)
ccc	CALL PRINF('FIRST  OUTPUT FILE NUMBER =*',IOUT1,1)
ccc	CALL PRINF('SECOND OUTPUT FILE NUMBER =*',IOUT2,1)
ccc	CALL PRINF('   *',I,0)
ccc	CALL PRINF('NUMBER OF PARTICLES =*',NAT,1)
	IF ((IFLAG7.EQ.3).OR.(IFLAG7.EQ.1)) IFLAG = 1
	IF (IFLAG7.EQ.2) IFLAG = 2
C
C       INITIALIZE ANALYTICAL ROUTINES
C
	 CALL DSTTIN(JKER)
	 CALL DSTABI(JKER)
C
C       SET THE USER-PROVIDED WORK ARRAY TO INTEGER ZERO
C
	DO 1080 I=1,10
	ZEROES(I)=0
 1080 CONTINUE
C
	DO 1090 I=1,NSP
	WKSP(I)=INZERO
 1090 CONTINUE
C
C ----- Determine number of multipole expansion terms needed for
C       specified precision
C
	IERR(1) = 0
	IERR(3) = NSP
	N      = -NINT( DLOG(TOL)/DLOG(2.0D0) ) + 3
	INFORM(1) = N
	IF (N .GT. 60) THEN
	   IERR(1) = 4
	   IERR(2) = N
	   CALL DAEOUT(IERR)
	   RETURN
	END IF
C
C ----- Set maximum number of levels allowed
C

	C102 = DLOG(10.0D0)/DLOG(2.0D0)
	REX2 = MEX*C102-1
CCCC      CALL PRIN2(' REX2=*',REX2,1)
	LMX = REX2/N+6
CCC     LMX = ( REX2-1)/N
CCCC      CALL PRINF(' LMX =*',LMX,1)
C
C ----- Check if number of particles NAT is less than  maximum
C       number of particles allowed per box.
C
	IF (NAT.LE.NAPB) THEN
	   IERR(1) = 16
	   IERR(2) = NAT
	   IERR(3) = NAPB
	   CALL DAEOUT(IERR)
	   RETURN
	END IF
C
C ----- Set potentials and fields to zero
C
	IF (IFLAG .EQ. 1) THEN
	   DO 500 I=1,NAT
	      FIELD(I) = CZERO
 500       CONTINUE
	ELSE
	   DO 600 I=1,NAT
	      POTEN(I) = DZERO
	      FIELD(I) = CZERO
 600       CONTINUE
	END IF
C
C ----- Center and scale the positions to BOXSIZ (128.)
C
	NXA = 1
	NYA = NXA + NAT
	NEWADR = NYA + NAT
	NEWSP = NSP - NEWADR + 1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
	RSCAL = 0.0
	RMINX =  1.0E30
	RMAXX = -1.0E30
	RMINY =  1.0E30
	RMAXY = -1.0E30
	DO 100 I=1,NAT
	   IF (XA(I) .GT. RMAXX) RMAXX = XA(I)
	   IF (YA(I) .GT. RMAXY) RMAXY = YA(I)
	   IF (XA(I) .LT. RMINX) RMINX = XA(I)
	   IF (YA(I) .LT. RMINY) RMINY = YA(I)
 100    CONTINUE
	RSX = RMAXX-RMINX
	RSY = RMAXY-RMINY
	XONEW = (RMAXX+RMINX)/2.0
	YONEW = (RMAXY+RMINY)/2.0
	RSCAL = DMAX1(RSX,RSY)
	RSCAL = BOXSIZ/RSCAL*0.99
	EPS = EPS7*RSCAL
	QTOT = 0.0
	DO 110 I=1,NAT
	   I1 = I-1
	   WKSP(NXA+I1) = (XA(I)-XONEW) *  RSCAL
	   WKSP(NYA+I1) = (YA(I)-YONEW) *  RSCAL
	   QTOT  = QTOT + QA(I)
 110    CONTINUE
C
C ----- allocate arrays ibox,xb,yb,nbat,npar,nchi,shift from workspace
C
	IBOX = NEWADR
	NEWADR = IBOX + (NAT + NINIRE-1)/NINIRE
	NEWSP = NSP - NEWADR + 1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
	NBOXMX =  NEWSP*NINIRE / (2*NINIRE+4)
	XB      = NEWADR
	YB      = XB + NBOXMX
	NBAT    = YB + NBOXMX
	NPAR    = NBAT + (NBOXMX + NINIRE-1)/NINIRE
	NCHI    = NPAR + (NBOXMX + NINIRE-1)/NINIRE
	SHIFT   = NCHI + (NBOXMX + NINIRE-1)/NINIRE
C
C ----- Form tables ibox, nbat, npar, nchi, nl
C
CCC        ITT = 0
CCC        ITIME = MRUN(0)
	CALL DAASSI (EPS7,NAPB,NAT,NBOX,NBOXMX,IERR,L,LMX,NL,
     *            WKSP(NXA),WKSP(NYA),QA,WKSP(IBOX),WKSP(XB),
     *            WKSP(YB),WKSP(NBAT),WKSP(NPAR), WKSP(NCHI))
        CALL PRINF('WKSP(NCHI) AFTER DAASSI IS*',WKSP(NCHI),NBOX)
        call PRIN2('WKSP(XB) AFTER DAASSI IS*',WKSP(XB),10)
	IF (IERR(1).NE.0) THEN
	   CALL DAEOUT(IERR)
	   IF (IERR(1).GT.1)  RETURN
	END IF
C
C ----- Compress tables by eliminating empty boxes
C
	CALL DACOMP (NAT,NBOX,NBOXOL,L,NL,WKSP(IBOX),WKSP(XB),
     *              WKSP(YB),WKSP(NBAT),WKSP(NPAR),
     *              WKSP(NCHI),WKSP(SHIFT))
        CALL PRINF('WKSP(NCHI) AFTER COMPRESSION*',WKSP(NCHI),NBOX)
C
C ----- Move tables to free workspace ...
C
	INFORM(2) = NBOX
	INFORM(3) = L-1
	NSHFT = NBOXMX - NBOX
C
C ----- Move YB
C
	NYB = YB - NSHFT
	DO 120 I=1,NBOX
	   I1 = I-1
	   WKSP(NYB+I1) = WKSP(YB+I1)
 120    CONTINUE
C
C ----- Move NBAT
C
	NBOXI= (NBOX + NINIRE - 1)/NINIRE
	NNBAT = NYB + NBOX
	DO 130 I=1,NBOXI
	   I1 = I-1
	   WKSP(NNBAT+I1) = WKSP(NBAT+I1)
 130    CONTINUE
C
C ----- Move NPAR
C
	NNPAR = NNBAT + NBOXI+1
	DO 140 I=1,NBOXI
	   I1 = I-1
	   WKSP(NNPAR+I1) = WKSP(NPAR+I1)
 140    CONTINUE
C
C ----- Move NCHI
C
	NNCHI = NNPAR + NBOXI + 1
	DO 150 I=1,NBOXI
	   I1 = I-1
	   WKSP(NNCHI+I1) = WKSP(NCHI+I1)
 150    CONTINUE
CCCC    CALL PRINF('MOVED NCHI IS*',WKSP(NNCHI),NBOX)
C
C ----- Compute new address of free workspace
C
	NEWADR = NNCHI + NBOXI + 1
C
C ----- Allocate memory for INDB.
C
	INDB = NEWADR
	NEWADR = INDB + NBOXI + 1
	NEWSP = NSP - NEWADR + 1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
C
C ----- Form index of childless and parent boxes  INDB.
C
	CALL DAINDB(NAPB,NAT,MAXP,NBOX,L,NL,WKSP(NNBAT),
     *                    WKSP(NNCHI),WKSP(INDB))
C
C ----- Allocate memory for particle lists
C
	IAT = NEWADR
	JAT = IAT + (NBOX+1 + NINIRE-1)/NINIRE
	NEWADR = JAT + (NBOX+1+NAT + NINIRE-1)/NINIRE
	NEWSP = NSP - NEWADR+1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
C
C ----- Form list of particles for childless boxes
C
	CALL DAIJAT (NAPB,NAT,NBOX,L,NL,WKSP(IBOX),WKSP(NNBAT),
     *               WKSP(INDB),WKSP(IAT),WKSP(JAT))
C
C ----- Allocate memory  for list of colleagues list
C       and first interaction list
C
	L1 = NEWADR
	L2 = L1 + NBOXI+1
	JL = L2 + NBOXI + 1
	NEWADR = JL + 36*(NBOXI+1)
	NEWSP = NSP - NEWADR+1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
C
C ----- Form interaction list (#1) and list of colleagues (#2)
C       #1: separated children of parent's colleagues
C       #2: Adjacent boxes of the same size
C
	CALL DALI12 (NBOX,LEN,L,NL,WKSP(XB),WKSP(NYB),
     *               WKSP(NNPAR),WKSP(NNCHI),WKSP(L1),WKSP(L2),
     *               WKSP(JL))
C
C ----- Compute new address of free workspace
C
	NEWADR = JL + (LEN + NINIRE - 1)/NINIRE
C
C ----- Allocate memory for third and fourth lists
C
	L3 = NEWADR
	L4 = L3 + NBOXI+1
	IDESC = L4 + NBOXI+1
CCC        NEWADR = IDESC + 20000
	NEWADR = IDESC + NBOXI + 1
	IHALF = (NSP - NEWADR + 1)/2
	JL3 = NEWADR
	JL4 = JL3 + IHALF
C
C ----- Form lists #3 and #4 for childless boxes
C       #3: childless adjacent descendents of colleagues
C       #4: separated descendants of colleagues
C
       CALL DALI34 (NAPB,NBOX,LEN3,LEN4,L,NL,WKSP(XB),
     *              WKSP(NYB),WKSP(NNBAT),WKSP(NNPAR),WKSP(NNCHI),
     *              WKSP(INDB),
     *              WKSP(L1),WKSP(L2),WKSP(JL),WKSP(L3),WKSP(L4),
     *              WKSP(IDESC),WKSP(JL3),WKSP(JL4))
C
C ----- Compress JL3 and JL4
C
CCC        ISH3     = 20000
	ISH3     = NBOXI + 1
	ISH4     = IHALF + ISH3 - (LEN3 + NINIRE-1)/NINIRE
	NJL3    = IDESC
	NJL4    = NJL3 + (LEN3 + NINIRE-1)/NINIRE
	NEWADR  = NJL4 + (LEN4 + NINIRE-1)/NINIRE
	JL3P    = JL3 + (LEN3-1 + NINIRE-1)/NINIRE
	JL4P    = JL4 + (LEN4-1 + NINIRE-1)/NINIRE
	NEWSP   = NSP - NEWADR+1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
	DO 160 I=JL3,JL3P
	   WKSP(I-ISH3) = WKSP(I)
 160    CONTINUE
	DO 170 I=JL4,JL4P
	   WKSP(I-ISH4) = WKSP(I)
 170    CONTINUE
C
C ----- All the necessary tables are now available
C
ccc        JTIME = MRUN(0)-ITIME
ccc        ITT = ITT + JTIME
ccc        WRITE (IOUT1,*) ' TIME IN BKKPNG =',JTIME
C
C ----- Allocate memory for multipole and local expansions
C
	LO  =  NEWADR
	TA  =  LO + (NBOX+1)*N*2
	IZ = TA + (NBOX+1)*N*2
	ICH = IZ + MAXP*2
	IBAK = ICH + MAXP*2
	N2 =IBAK + (MAXP + NINIRE -1)/NINIRE
	NEWADR = N2 + (MAXP + NINIRE -1)/NINIRE
	INFORM(4) = NEWADR
	NEWSP   = NSP - NEWADR+1
	IF (NEWSP .LT. 0)  THEN
	   IERR(1) = 8
	   IERR(2) = -NEWSP+NSP
	   CALL DAEOUT(IERR)
	   RETURN
	ENDIF
C
C ----- SET WORKSPACE TO ZERO
C
	DO 180 I=LO,NEWADR
	   WKSP(I) = DZERO
 180    CONTINUE
C
C ----- Print informations
C
ccc	CALL PRINF('NTERMS = NB OF TERMS IN EXPANSIONS*',INFORM,0)
ccc	CALL PRINF('NBOXES = NB OF BOXES*',INFORM,0)
ccc	CALL PRINF('NLEVLS = NB OF LEVELS*',INFORM,0)
ccc	CALL PRINF('NSPACE = AMOUNT OF SPACE USED*',INFORM,0)
	CALL PRINF(' NTERMS  NBOXES  NLEVLS  NSPACE *',INFORM,4)
	CALL PRINF('   *',IFLAG,0)
C
C ----- Upward pass: form multipole expansions for childless boxes
C       and merge them at all coarser mesh levels
C
ccc        ITIME = MRUN(0)
	CALL DAUPWD (N,NAPB,IFLAG,NAT,NBOX,L,NL,WKSP(NXA),
     *               WKSP(NYA),QA,WKSP(IBOX),WKSP(XB),
     *               WKSP(NYB),WKSP(NNBAT),WKSP(NNPAR),WKSP(NNCHI),
     *               WKSP(INDB),WKSP(IAT),WKSP(JAT), WKSP(LO),
     *               WKSP(IZ),WKSP(ICH))
ccc        JTIME = MRUN(0)-ITIME
ccc        ITT = ITT + JTIME
ccc        WRITE (IOUT1,*) ' TIME IN DAUPWD =',JTIME
CCCC         CALL PRINF('AFTER DAUPWD, IFLAG IS*',IFLAG,1)
C
C ----- At each level, convert multipole expansion of box to local
C       expansion about the centers of boxes in fist list, and add
C       these to initial local expansions obtained from parent level.
C
ccc        ITIME = MRUN(0)
	CALL DAFLEX (N,NAPB,IFLAG,NAT,NBOX,L,NL,WKSP(NXA),
     *               WKSP(NYA),WKSP(IBOX),WKSP(XB),
     *               WKSP(NYB),WKSP(NNBAT),WKSP(NNPAR),WKSP(NNCHI),
     *               WKSP(IAT),WKSP(JAT),WKSP(L1),WKSP(L2),
     *               WKSP(JL),WKSP(LO),WKSP(TA),POTEN,FIELD)
ccc        JTIME = MRUN(0)-ITIME
ccc        ITT = ITT + JTIME
ccc        WRITE (IOUT1,*) ' TIME IN DAFLEX =',JTIME
CCCC         CALL PRINF('AFTER DAFLEX, IFLAG IS*',IFLAG,1)
C
C ----- For childless boxes compute interactions with second/third
C       list and fourth list
C
CCC        ITIME = MRUN(0)
	CALL DABARR (N,NAPB,IFLAG,NAT,NBOX,L,NL,QTOT,RSCAL,XA,YA,
     *              WKSP(NXA),WKSP(NYA),QA,WKSP(IBOX),WKSP(XB),
     *              WKSP(NYB),WKSP(NNBAT),WKSP(NNPAR),WKSP(NNCHI),
     *              WKSP(INDB),
     *              WKSP(IAT),WKSP(JAT),WKSP(L1),WKSP(L2),WKSP(JL),
     *              WKSP(L3),WKSP(L4),WKSP(NJL3),WKSP(NJL4),WKSP(LO),
     *              WKSP(TA),WKSP(IZ),WKSP(ICH),WKSP(IBAK),WKSP(N2),
     *              POTEN,FIELD,EPS,CLOSE)
CCC        JTIME = MRUN(0)-ITIME
CCC        ITT = ITT + JTIME
CCC        WRITE (IOUT1,*) ' TIME IN DABARR =',JTIME
CCCC         CALL PRINF('AFTER DABARR, IFLAG IS*',IFLAG,1)
C
C ----- Shift the local expansions from parents to children,
C ----- evaluate local expansions at particles' position
C
ccc        ITIME = MRUN(0)
	CALL DADNWD (N,IFLAG,NAT,NBOX,L,NL,WKSP(NXA),WKSP(NYA),QA,
     *               WKSP(IBOX),WKSP(XB),WKSP(NYB),WKSP(NNBAT),
     *               WKSP(NNPAR),WKSP(NNCHI),WKSP(TA),POTEN,FIELD)
ccc        JTIME = MRUN(0)-ITIME
ccc        ITT = ITT + JTIME
ccc        WRITE (IOUT1,*) ' TIME IN DADNWD =',JTIME
CCCC         CALL PRINF('AFTER DADNWD, IFLAG IS*',IFLAG,1)
C
C ----- Rescale potentials and fields
C
	COEF = RSCAL
	COEP = DLOG(RSCAL)
	IF (IFLAG7 .EQ. 1)  THEN
	   DO 200 I=1,NAT
	      FIELD(I) = FIELD(I)*COEF
 200       CONTINUE
	END IF
	IF (IFLAG7 .EQ. 2) THEN
	   DO 210 I=1,NAT
	      FIELD(I) = FIELD(I)*COEF
	      POTEN(I) = POTEN(I)-(QTOT-QA(I))*COEP
 210       CONTINUE
	END IF
	IF (IFLAG7 .EQ. 3)  THEN
	   DO 220 I=1,NAT
	      FIELD(I) = DCONJG(FIELD(I))*COEF
 220       CONTINUE
	END IF
	RETURN
C
C**********************************************************************
C
       ENTRY DAPIP2 (NAPB,IFLAG7,NAT,XA,YA,QA,
     *               WKSP,NSP,XAP,YAP,PO,FI,EPS7,CLOSEP)
C
C   *** INFORMATION :
C
C       Second entry of the subroutine DAPIF2 .  For a detailed
C       description see the abbreviated manual.
C
C   *** DESCRIPTION :
C
C
C       Compute the potential and/or field at an arbitrary point
C       in the computationnal cell.
C
C   *** INPUT PARAMETERS:
C
C   NAPB      = maximum number of particles in a box
C   IFLAG7    = 1 for computing fields only
C               2 for computing potential and electric fields
C               3 for Hilbert Matrix
C   NAT       = number of particles
C   QA(i)     = charge of a particle
C   WKSP      = workspace array
C   NSP       = length of work array
C   XAP       = first coordinate of the point
C   YAP       = second coordinate of the point
C   EPS       = near approach distance.
C               this number is provided by the user to check for very
C               near approaches. If the point and a particle are within
C               EPS ofeach other, the program executes a near approach
C               subroutine CLOSEP which the user must declare as
C               EXTERNAL in the calling program.
C   CLOSEP    = is the near approach subroutine to be provided by the
C               user. The calling sequence should be
C                 CALL CLOSEP(EPS,I,QI,XAP,YAP,XI,YI,FI,PO)
C               In many applications, when particles are very close,
C               the nature of the force law changes. For example,
C               the charges may be viewed as continuous distributions
C               about their centers instead of as point charges.
C
C   *** OUTPUT PARAMETERS:
C
C   FI         = net field at location of the point
C   PO         = potential at location of the point
C
C   *** SUBROUTINE CALLED:
C
C       DAAUPO
C
C
C**********************************************************************
C
	IF ((IFLAG7.EQ.3).OR.(IFLAG7.EQ.1)) IFLAG = 1
	IF (IFLAG7.EQ.2) IFLAG = 2
C
C ----- Scale arbitrary point position
C
	XAP1 = (XAP-XONEW) * RSCAL
	YAP1 = (YAP-YONEW) * RSCAL
	EPS = EPS7 * RSCAL
	CALL DAAUPO (N,NAPB,IFLAG,NAT,NBOX,L,NL,
     *       WKSP(NXA), WKSP(NYA),XA,YA,QA,WKSP(IBOX),WKSP(XB),
     *       WKSP(NYB), WKSP(NNBAT),WKSP(NNPAR),WKSP(NNCHI),
     *       WKSP(INDB),WKSP(IAT),
     *       WKSP(JAT),WKSP(LO),XAP,YAP,XAP1,YAP1,PO,FI,EPS,RSCAL,
     *       QTOT,CLOSEP)
       FI = DCONJG(FI)
C
C ----- Scale potentials and fields back
C
	IF (IFLAG7 .EQ. 1)  THEN
	   FI = FI * COEF
	END IF
	IF (IFLAG7 .EQ. 2)  THEN
	   FI = FI * COEF
	   PO = PO - QTOT * COEP
	END IF
	IF (IFLAG7 .EQ. 3)  THEN
	   FI = DCONJG(FI) * COEF
	END IF
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAASSI (EPS,NAPB,NAT,NBOX,NBOXMX,IERR,L,LMX,NL,
     *                     XA,YA,QA,IBOX,XB,YB,NBAT,NPAR,NCHI)
C
C**********************************************************************
C
C   *** DESCRIPTION:
C
C       Form tables to describe irregular grid.
C
C
C   *** INPUT PARAMETERS:
C
C   NAPB      =  maximum number of particles per box
C   NAT       =  number of particles
C   NBOXMX    =  maximum number of boxes
C   LMX       =  maximum number of levels
C   XA(I)     =  first coordinate of ith particle
C   YA(I)     =  second coordinate of ith particle
C   QA(I)     =  charge of ith particle
C
C   *** OUTPUT PARAMETERS:
C
C   NBOX      =  number of boxes
C   L         =  number of levels plus one
C   NL(I)     =  pointer to first box of ith level
C   IBOX(I)   =  address of the box containing ith particle
C   XB(I)     =  first coordinate of the center of ith box
C   YB(I)     =  second coordinate of the center of ith box
C   NBAT(I)   =  number of particles in ith box
C   NPAR(I)   =  address of ith box's parent (0 for boxes at 1st level)
C   NCHI(I)   =  address of ith box's first child (0 for a childless box)
C
C   *** SUBROUTINES CALLED ***
C
C       DAINIA
C       DAIBNB
C       DAFAMI
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER IBOX(1),IERR(1),NL(1),NBAT(1),NPAR(1),NCHI(1)
	DOUBLE PRECISION XA(1),YA(1),XB(1),YB(1)
	DOUBLE COMPLEX QA(1)
C
C ----- Initializations
C
	CALL DAINIA (NAT,L,MAX,NBOX,NBOXMX,NL,
     *               NBAT,NPAR,NCHI,XB,YB,IBOX)
C
C ----- Divide only boxes  with more than NAPB particles
C
	DO 10 INTY=1,100000000
C
C ----- Check if number of levels does not exceed maximum allowed
C
	IF (L.GT.LMX) THEN
	   IERR(1) = 32
	   IERR(2) = L-1
	   IERR(3) = LMX-1
	   RETURN
	END IF
C
C ----- Compare size of smallest box with EPS.
C       If smallest box smallest box is smaller than EPS,
C       stop the refinement
C
	SSB =  2.0**(-(L-1))
	IF (SSB.LE.EPS) THEN
	   IERR(1) = 1
	   RETURN
	END IF
	   IF (MAX .LE. NAPB) GOTO 11
C
C ----- Form IBOX and NBAT
C
	CALL DAIBNB (NAPB,NAT,NBOX,L,NL,IBOX,NBAT,NPAR,NCHI,XB,YB,
     *               XA,YA,QA)
C
C ----- Form parents and children tables: NPAR and NCHI
C
	CALL DAFAMI (NAPB,NAT,NBOX,NBOXMX,IERR,L,NL,NBAT,
     *               NPAR,NCHI,XB,YB,MAX)
	IF (IERR(1) .NE. 0) RETURN
C
C ----- End of the loop
C
 10     CONTINUE
 11     CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAINIA (NAT,L,MAX,NBOX,NBOXMX,NL,
     *                   NBAT,NPAR,NCHI,XB,YB,IBOX)
C
C**********************************************************************
C
C
C   *** DESCRIPTION :
C
C       Initializations for the grid at the first level
C
C   *** PARAMETERS :
C
C       See subroutine DAASSI
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER NL(1),IBOX(1),NBAT(1),NPAR(1),NCHI(1)
	DOUBLE PRECISION XB(1),YB(1)
	DATA BOXSIZ/128.0/
	DO 1 I=1,4
	   NBAT(I) = 0
	   NPAR(I) = 0
 1      CONTINUE
	DO 3 I=1,NAT
	   IBOX(I) = 0
 3      CONTINUE
	NL(1) = 1
	NL(2) = 5
	NBOX = 4
	D = 0.25*BOXSIZ
	XB(1) = -D
	YB(1) =  D
	XB(2) =  D
	YB(2) =  D
	XB(3) = -D
	YB(3) = -D
	XB(4) =  D
	YB(4) = -D
	L = 1
	MAX = NAT
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAIBNB (NAPB,NAT,NBOX, L,NL,IBOX,NBAT,NPAR,NCHI,
     *                   XB,YB,XA,YA,QA)
C
C**********************************************************************
C
C
C   *** DESCRIPTION :
C
C       Form IBOX(I)  adress of box containing ith particle
C    and NBAT(I) number of particles in ith box
C
C    *** PARAMETERS :
C
C       See subroutine DAASSI
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER  IBOX(1), NBAT(1), NPAR(1), NCHI(1),NL(1)
	DOUBLE PRECISION XB(1),YB(1),XA(1),YA(1)
	DOUBLE COMPLEX QA(1)
C
C ----- Set number of particles to zero for all boxes at level L
C
	DO 1 I=NL(L),NL(L+1)
	   NBAT(I) = 0.0
 1      CONTINUE
C
C ----- Loop on particles, update ibox and nbat
C
	DO 100 I=1,NAT
	   IB = IBOX(I)
	   X  = XA(I)
	   Y  = YA(I)
	   IF (IB .EQ. 0) THEN
	      NB = NAT
	      XP = 0.0
	      YP = 0.0
	      NC = 1
	   ELSE
	      NB = NBAT(IB)
	      IF (NB .LE. NAPB) GOTO 100
	      XP = XB(IB)
	      YP = YB(IB)
	      NC = NCHI(IB)
	   ENDIF
	   IF (X .LT. XP) THEN
	      IF (Y .LT. YP) THEN
		 IBOX(I) = NC+2
		 NBAT(NC+2) = NBAT(NC+2)+1
	      ELSE
		 IBOX(I) = NC
		 NBAT(NC) = NBAT(NC)+1
	      ENDIF
	   ELSE
	      IF (Y .LT. YP) THEN
		 IBOX(I) = NC+3
		 NBAT(NC+3) = NBAT(NC+3)+1
	      ELSE
		 IBOX(I) = NC+1
		 NBAT(NC+1) = NBAT(NC+1)+1
	      END IF
	   ENDIF
 100    CONTINUE
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAFAMI (NAPB,NAT,NBOX,NBOXMX,IERR,L,NL,NBAT,
     *                     NPAR,NCHI,XB,YB,MAX)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Form tables of parents (NPAR) and children (NCHI) and compute
C   the coordinates of the centers of boxes (XB and YB).
C
C   *** ARGUMENTS :
C
C       See subroutine DAASSI
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER IERR(1),NL(1),NBAT(1),NPAR(1),NCHI(1)
	DOUBLE PRECISION XB(1),YB(1)
	DATA BOXSIZ/128.0/
	NBL = NBOX
	NLL = NL(L)
	RL2 = 2.0**L
	D = 0.25/RL2 *BOXSIZ
	MAX = 0
C
C ----- Loop on the boxes of finest level
C
	DO 100 I=NLL,NBL
	   INB = NBAT(I)
	   IF (INB .GT. MAX) MAX = INB
C
C ----- Form table of children
C
	   NCHI(I)=0
	   IF (INB .LE. NAPB) GOTO 100
	   NCHI(I) = NBOX+1
	   DO 10 J=1,4
	      JJ = NBOX+J
C
C ----- Form table of parents
C
	      NPAR(JJ) = I
C
C ----- Compute coordinates of centers
C
	      IF (J .EQ. 1) THEN
		 XB(JJ) = XB(I) - D
		 YB(JJ) = YB(I) + D
	      END IF
	      IF (J .EQ. 2) THEN
		 XB(JJ) = XB(I) + D
		 YB(JJ) = YB(I) + D
	      END IF
	      IF (J .EQ. 3) THEN
		 XB(JJ) = XB(I) - D
		 YB(JJ) = YB(I) - D
	      END IF
	      IF (J .EQ. 4) THEN
		 XB(JJ) = XB(I) + D
		 YB(JJ) = YB(I) - D
	      END IF
 10        CONTINUE
	   NBOX = NBOX +4
	   IF (NBOX .GT. NBOXMX) THEN
	      IERR(1) = 8
	      IERR(2) = (NBOXMX-NBOX)*36+NSP
	      RETURN
	   END IF
 100    CONTINUE
C
C ----- Update level L and table NL
C
	L = L+1
	NL(L+1) = NBOX + 1
CCC        CALL PRINF(' L = *',L,1)
CCC        CALL PRINF(' NL = *',NL,L+1)
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DACOMP (NAT,NBOX,NBOXOL,L,NL,IBOX,XB,
     *                    YB,NBAT,NPAR,NCHI,SHIFT)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Compress tables XB,YB,NBAT,NPAR,NCHI,IBOX,NL
C   by eliminating empty boxes
C
C   *** INPUT ARGUMENTS :
C
C   NAT       =  number of particles
C   NBOXOL   =  number of boxes before elimination of empty boxes
C   L         =  number of level plus one
C
C   *** OUTPUT ARGUMENTS :
C
C   NBOX      =  number of boxes after elimination of empty boxes
C   SHIFT(I)  =  number of empty boxes with address less than i
C
C   *** INPUT/OUTPUT ARGUMENTS :
C
C   NL(I)     = pointer to first box of ith level
C   IBOX(I)   = address of the box containing ith particle
C   XB(I)     = first coordinate of the center of ith box
C   YB(I)     = second coordinate of the center of ith box
C   NBAT(I)   = number of particles in ith box
C   NPAR(I)   = address of ith box's parent (0 for boxes at first level)
C   NCHI(I)   = address of ith box's first child (0 for a childless box)
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER NL(1),IBOX(1),NBAT(1),NPAR(1),NCHI(1),SHIFT(1)
	DOUBLE PRECISION XB(1),YB(1)
C
C ----- Create table shift
C
	NZ = 0
	DO 1 I=1,NBOX
	   NB = NBAT(I)
	   SHIFT(I) = NZ
	   IF (NB.EQ.0) NZ = NZ+1
 1      CONTINUE
	SHIFT(NBOX+1)=NZ
CCC        CALL PRINF(' TABLE SHIFT =*',SHIFT,NBOX+1)
C
C ----- Shift tables: first XB,YB,NBAT,NPAR and NCHI
C
	DO 2 I=1,NBOX
	   NSH = SHIFT(I)
	   NB = NBAT(I)
	   IF (NB .GT. 0) THEN
	      NEWI = I - NSH
	      XB(NEWI)    = XB(I)
	      YB(NEWI)    = YB(I)
	      NBAT(NEWI)  = NBAT(I)
	      NPAR(NEWI)  = NPAR(I) - SHIFT(NPAR(I))
	      NCHI(NEWI)  = NCHI(I) - SHIFT(NCHI(I))
	   ENDIF
 2      CONTINUE
C
C ----- Now shift IBOX and NL
C
	DO 3 I=1,NAT
	   IBOX(I) = IBOX(I) - SHIFT(IBOX(I))
 3      CONTINUE
	DO 4 I=1,L+1
	   NL(I) = NL(I) - SHIFT(NL(I))
 4      CONTINUE
C
C ----- Compute number of (non-empty) boxes
C
	NBOXOL = NBOX
	NBOX = NBOX - NZ
	NPAR(NBOX+1) = 0
CCC        CALL PRINF( 'IBOX =*',IBOX,NAT)
CCC        CALL PRINF( 'NBAT =*',NBAT,NBOX)
CCC        CALL PRINF( 'NPAR =*',NPAR,NBOX)
CCC        CALL PRINF( 'NCHI =*',NCHI,NBOX)
CCC        CALL PRINF( 'NL =*',NL,L)
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAINDB (NAPB,NAT,MAXP,NBOX,L,NL,NBAT,NCHI,INDB)
C
C**********************************************************************
C
C   *** DESCRIPTION:
C
C       This  subroutine forms the array INDB. INDB(i) is equal
C       to 1 if box i is childless and equal to 0 if  box  i is
C       a parent box.
C       This subroutine also returns MAXP, actual maximum number
C       of particles in any box.
C
C
C   *** INPUT PARAMETERS:
C
C   NAPB      = maximum number of particles in a box
C   NAT       = number of particles
C   NBOX      = number of boxes in the subdivision of the
C                computational cell.
C   L         =  number of levels plus one
C   NL(I)     =  pointer to first box of ith level
C   NBAT(I)   =  number of particles in ith box
C   NCHI(I)   =  address of ith box's first child (0 for a childless box)
C
C   *** OUTPUT PARAMETERS:
C
C   MAXP      = actual maximum number of particles in any box
C   INDB(i)   = 1 for a childless box, 0 for a parent box.
C
C
C**********************************************************************
C
C
	INTEGER NL(1),NBAT(1),NCHI(1),INDB(1)
C
C ----- Build a table : 1 if box is childless, 0 if box is parent.
C
	MAXP = 0
	DO 200 LEV =1,L-1
	   NLD = NL(LEV)
	   NLU = NL(LEV+1)-1
C
C ----- All boxes at finest  level are childless
C
	   IF (LEV.EQ.L-1) THEN
	      DO 100 IB=NLD,NLU
		 INDB(IB) = 1
		 NCHI(IB) = 0
		 NB = NBAT(IB)
		 IF (NB.GT.MAXP) MAXP = NB
 100          CONTINUE
C
C ----- At coarser levels a box is childless if it contains less
C       than NAPB particles.
C
	   ELSE
	      DO 150 IB =NLD,NLU
		 NB = NBAT(IB)
		 IF (NB .LE. NAPB) THEN
		    INDB(IB) = 1
		    IF (NB.GT.MAXP) MAXP = NB
		 ELSE
		    INDB(IB) = 0
		 END IF
 150          CONTINUE
	   END IF
 200    CONTINUE
CCC        CALL PRINF(' TABLE NL =*',NL,L)
CCC        CALL PRINF(' TABLE NBAT =*',NBAT,NBOX)
CCC        CALL PRINF(' TABLE NCHI =*',NCHI,NBOX)
CCC        CALL PRINF(' TABLE INDB =*',INDB,NBOX)
	RETURN
	END
C
C**********************************************************************
C
	SUBROUTINE DAIJAT (NAPB,NAT,NBOX,L,NL,IBOX,NBAT,INDB,
     *             IAT,JAT)
C
C**********************************************************************
C
C   *** DESCRIPTION :
C
C       Form the list of particles in each box JAT
C   and the pointer to this list IAT
C
C   *** INPUT PARAMETERS :
C
C   NAPB        = maximum number of particles per box
C   NAT         = number of particles
C   NBOX        = number of boxes
C   L           = number of levels
C   NL(I)       = address of the first box of ith levels
C   IBOX(I)     = address of box containing ith particle
C   NBAT(I)     = number of particles in ith box
C   INDB(I)     = 0 if i is a parent box, 1 if it is a childless box
C
C   *** OUTPUT PARAMETERS :
C
C   IAT(I)      = pointer to the list of particles in ith box
C                 negative for boxes that are not childless
C   JAT(IAT(I)) = is  equal to NBAT(I)
C   JAT(J)      = address of particle contained in ith box
C                 if J in [ IAT(I)+1 , IAT(I) + NBAT(I) ]
C
C**********************************************************************
C
	INTEGER NL(1),IBOX(1),NBAT(1),INDB(1),IAT(1),JAT(1)
	NACC = 1
	NC   = 0
C
C  ----- Loop on the boxes: form pointer IAT
C
	DO 1 I=1,NBOX
	      NB = NBAT(I)
	   IF (INDB(I) .EQ. 0) THEN
CCC           IF (NB .GT. NAPB)  THEN
	      IAT(I) = -NB
	   ELSE
	      IAT(I) = NACC + NC
	      NC = NB + 1
	      NACC = IAT(I)
	      JAT(NACC) = 0
	   ENDIF
 1      CONTINUE
C
C ----- Loop on particles: form list JAT
C
	DO 2 I=1,NAT
	   IB = IBOX(I)
	   IA = IAT(IB)
	   JAT(IA) = JAT(IA) + 1
	   JAT(IA+JAT(IA)) = I
 2      CONTINUE
CCC        CALL PRINF(' IAT=*',IAT,NBOX+1)
CCC        CALL PRINF(' JAT=*',JAT,NAT+NBOX+1)
	RETURN
	END
C
