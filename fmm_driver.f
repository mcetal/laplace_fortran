      program FMM_DRIVER
c
c  This program demonstrates how to use the 2d FMM code for the Coulomb
c  Potential
c
      implicit real*8 (a-h,o-z)
      parameter (nmax = 1000000, nsp = 16*nmax+1000)
c
c Source Points:
      dimension xat(nmax), yat(nmax)
c
c Source Strengths:
      dimension randvec(nmax)
      complex*16 qa(nmax)
c
C Target points:
	integer m
	dimension randvec1(nmax),randvec2(nmax),xp(nmax),yp(nmax)
c
c Exact potenial and field
      dimension poten_ex(nmax)
      complex*16 cfield_ex(nmax), zi, zj, eye, zdis
c
c Exact potential and field at target points
	dimension ptar_ex(nmax)
	complex*16 cftar_ex(nmax)
c
c FMM potential and field at target points
	dimension ptar(nmax)
	complex*16 cftar(nmax)
	
c Arrays returned by FMM
      complex*16 cfield(nmax)
      dimension poten(nmax)
c
c Fast Multipole Work Arrays
      dimension wksp(nsp)
      integer*4 iout(2), inform(10), ierr(10)
c
c For cpu times
      REAL*4 TIMEP(2), ETIME
c
c define complex I and pi
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
c initialize plotting routine
         call PRINI (6,13)
c
c for lack of a better idea, set up source points on a grid
         icount = 0
         do ix = -5, 5
            do iy = -5, 5
               xat(icount+1) = float(ix)
               yat(icount+1) = float(iy)
               icount = icount+1
            end do
         end do
         nat = icount
         call PRINF (' NUMBER OF SOURCES = *', nat, 1)
         call prin2 (' X COORD OF SOURCES = *', xat, nat)
         call prin2 (' Y COORD OF SOURCES = *', yat, nat)
c
c Calculate random source densities on (0,1)
         call RANLUX (randvec,nat)
         do i = 1, nat
            qa(i) = randvec(i)
         end do
         call PRIN2 (' SOURCE STRENGTHS = *', qa, nat)
c
c Assign random target points
	 m = 10
	 call RANLUX(randvec1,m)
	 call RANLUX(randvec2,m)
	 do i = 1,m
		xp(i) = randvec1(i) + i*2.5d0
		yp(i) = randvec2(i) + i*2.5d0
	 end do
	call PRIN2("TARGET POINTS' X COORDINATES = *", XP,M)
	call PRIN2("TARGET POINTS' Y COORDINATES = *",YP,M)
	 
c Calculate exact potential and field
c let z_i = (x_i,y_i) or x_i + I y_i
c     potential at z_i = sum_j qa log|z_i - z_j|
c     field at z_i = sum_j qa (z_i - z_j)/|z_i-z_j|^2
        tbeg = etime(timep)
	  call DIPOTF(XAT,YAT,QA,NAT,XAT,YAT,NAT,POTEN_EX,CFIELD_EX)
	  call DIPOTF(XAT,YAT,QA,NAT,XP,YP,M,PTAR_EX,CFTAR_EX) 
        tend = etime(timep)
        call PRIN2 (' Time for Direct Calculation = *', tend-tbeg, 1)
c
c Call FMM to calculate field and compare to exact answer
c First set up input arguments
c      plotting - setting iout(1) = 6 will write FMM info to screen
         iout(1) = 0
         iout(2) = 13
c
c      iflag7 = 2 for computing potentials and fields, see DAPIF2 for 
c      description of other variables
         iflag7 = 2
         napb = 20
         ninire = 2
         mex = 300
         eps7 = 1.d-14
c
c      tolerance for FMM
         tol = 1.d-14
         tbeg = etime(timep)
         call DAPIF2 (iout, iflag7, nat, napb, ninire, mex, ierr, 
     1                inform, tol, eps7, xat, yat, qa, poten,  
     2                cfield, wksp, nsp, xp, yp, m,ptar,cftar, CLOSE)
 
         call PRINI (6, 13)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
         tend = etime(timep)
         call PRIN2 (' Time for FMM Calculation = *', tend-tbeg, 1)
c
c Check error in calculation
        err_pot = 0.d0
        err_field = 0.d0
	  epot_tar = 0.d0
	  ecfield_tar = 0.d0
        do i = 1, nat
           err_pot = max(err_pot,dabs(poten(i)-poten_ex(i)))
           err_field = max(err_field,cdabs(cfield(i)-cfield_ex(i)))
        end do
	  do i = 1, m
		epot_tar = max(epot_tar,dabs(ptar(i)-ptar_ex(i)))
		ecfield_tar = max(ecfield_tar,cdabs(cftar(i)-cftar_ex(i)))
	  end do
        call prin2 (' error in potential = *', err_pot, 1)
        call prin2 (' error in field = *', err_field, 1)
        call PRIN2 (' Exact Potential = *', poten_ex, nat)
        call PRIN2 (' FMM Potential = *', poten, nat) 
        call PRIN2 (' Exact Field = *', cfield_ex, 2*nat)
        call PRIN2 (' FMM Field = *', cfield, 2*nat)
	  call PRIN2 (' error in target potential = *', epot_tar,1)
	  call PRIN2 (' error in target field = *', ecfield_tar,1)
	  call PRIN2 (' Exact Potential at target points = *', ptar_ex, m)
	  call PRIN2 (' Exact Field at target points = *',cftar_ex, m)
	  call PRIN2 (' FMM potential at target points = *',ptar,m)
	  call PRIN2 (' FMM field at target points = *',cftar,m)
c
      stop
      end
         

