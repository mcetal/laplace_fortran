      PROGRAM LAPLACE_INTEQN
c-----------------------------------------------------------------------
c
c  Solves integral equation for Laplace's Equation in a bounded
c  simply-connected domain
c
c  I've included some mechanisms to allow for multiply-connectedness
c  in future
c     kmax = number of boundary curves
c            1 for simply connected, >1 for multiply connected
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      parameter (kmax = 1, npmax = 2048, nmax = kmax*npmax)
c
c Geometry of Domain
      dimension ak(kmax), bk(kmax), ncyc(kmax)
      complex*16 zk(kmax)
c
c Geometry variables for discretization 
      dimension x(nmax), y(nmax), dsdth(nmax), rkappa(nmax)
      complex*16 z(nmax), dz(nmax)
c
c  Grid variables
      parameter (nx_max = 1000, ny_max = 1000, ng_max = nx_max*ny_max)
      dimension igrid(ng_max), u_gr(ng_max), x_gr(ng_max), y_gr(ng_max)
c
c target points are used to check accuracy of final solution
      parameter (ntar=100)
      dimension xtar(ntar), ytar(ntar), u_tar(ntar)
      complex*16 ztar(ntar)    
c
c system matrix
      dimension amat((nmax+kmax)**2)      
c
c Density 
      dimension density(nmax)
c
c  Matrix equation variables for GMRES
c  MAXL is the maximum nubmer of GMRES iterations performed
c       before restarting.
c  LRWORK is the dimension of a real workspace needed by DGMRES.
c  LIWORK is the dimension of an integer workspace needed by DGMRES.
c  GMWORK and IWORK are work arrays used by DGMRES
c
      parameter (maxl = 50,liwork=30,  
     1           lrwork=10+(nmax+kmax)*(maxl+6)+maxl*(maxl+3))
      dimension gmwork(lrwork), igwork(liwork)
      dimension rhs(nmax+kmax), soln(nmax+kmax)
c
c Fast Multipole Arrays
      parameter (nsp = 20*nmax + 20*ng_max)
      dimension x_zeta(nmax+ng_max), y_zeta(nmax+ng_max)
      complex*16 qa(nmax+ng_max), cfield(nmax+ng_max)
      dimension poten(nmax+ng_max), wksp(nsp)
c
c FFT arrays
      dimension wsave(4*npmax+15)
      complex*16 zf(npmax)
c
c Other arrays
      dimension alpha(nmax), w(nmax)
      complex*16 zeta_k(kmax)
      REAL*4 TIMEP(2), ETIME
c
c common blocks
      common /geometry/ x, y, zk, z, dz, dsdth, rkappa
      common /sys_size/ k, nd, nbk
c
c Read in domain geometry, system size and initialize various
c things
         call INITIALIZE (k, nd, nbk, nx, ny, ak, bk, zk, ncyc, wsave)
c
c Construct contour geometry 
         call MAKE_GEO (k, nd, nbk, ak, bk, zk, ncyc, x, y, z, dz, 
     1                  dsdth, rkappa)
c
c Construct target points to check solution
         call GET_TARGET (ntar, xtar, ytar, ztar)
c
c Get the Dirichlet BCs to form RHS of integral equation
         call GET_BCS (k, nd, nbk, zk, z, rhs)
ccc         call PRIN2 (' rhs = *', rhs, nbk)
c
c Solve!
         call SOLVE (nd, k, nbk, maxl, rhs, soln, density, gmwork,  
     1               lrwork, igwork, liwork)
c
c Calculate solution at target points to check accuracy
         call GET_SOL_TAR (k, nd, nbk, ntar, density, zk, z, dz, ztar, 
     1                     u_tar)
c
      stop
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine INITIALIZE (k, nd, nbk, nx, ny, ak, bk, zk, ncyc,
     1                       wsave)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(*), bk(*), ncyc(*), wsave(*)
      complex*16 zk(*)
c
c initialize PRINT output locations to screen and fort.13
         call PRINI (6,13)
c
c read in initial data
c typically, ak, bk are major/minor axes of ellipse
c (cx,cy) is the centre
c ncyc is added if we want star shapes
         open (unit = 12, file = 'input.data')
         read (12,*) k, nd
         read (12,*) nx, ny
         do kbod = 1, k
            read(12,*) ak(kbod), bk(kbod), cx, cy, ncyc(kbod)
            zk(kbod) = dcmplx(cx,cy)
         end do
         close (12)
         nbk = k*nd
         call PRINF ('Number of countours = *', k, 1)
         call PRINF ('Number of points per contour = *', nd, 1)
         call PRINF ('Total system size = *', nbk, 1)
         call PRIN2 ('   ak = *', ak, k)
         call PRIN2 ('   bk = *', bk, k)
         call PRIN2 ('   zk = *', zk, 2*k)
         call PRINF ('   ncyc = *', ncyc, k)
         call PRINF ('Number of grid points in x = *', nx, 1)
         call PRINF ('Number of grid points in y = *', ny, 1)
c
c initialize fft workspace array
         call DCFFTI (nd, wsave)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CURVE_PARAM (th, ai, bi, thetax, ncyc, xp, yp, xdot,
     *                        ydot, xddot, yddot)
c---------------
c
c  The initial curve parameterization is given in term of theta
c  This routine returns curve points and derivative data given
c  theta
c
      implicit real*8 (a-h, o-z)
c
         if (ncyc.eq.0) then
            cs = dcos(th-thetax)
            sn = dsin(th-thetax)
            a2 = (ai)**2
            b2 = (bi)**2
            rnorm = dsqrt(b2*cs**2 + a2*sn**2)
            radius = ai*bi/rnorm
            rdot = -ai*bi*cs*sn*(-b2+a2)/rnorm**3
            rddot =  -ai*bi*(2.d0*a2*b2*cs**2*sn**2
     *                     +a2**2*cs**4 + a2*b2-a2**2+b2**2*cs**4
     *                     -2.d0*b2**2*cs**2)/rnorm**5
            xp = radius*dcos(th)
            yp = radius*dsin(th)
c
            cs = dcos(th)
            sn = dsin(th)
            xdot = rdot*cs - radius*sn
            ydot = rdot*sn + radius*cs
            xddot = rddot*cs - 2.d0*rdot*sn - radius*cs
            yddot = rddot*sn + 2.d0*rdot*cs - radius*sn
           elseif (ncyc.lt.0) then
            R = ai
            anu = bi
            den = (1.d0-2.d0*anu*dcos(2.d0*th)+anu**2)
     *           *dsqrt(1.d0+anu**2)
            ddot = 4.d0*anu*dsin(2.d0*th)*dsqrt(1.d0+anu**2)
            dddot = 8.d0*anu*dcos(2.d0*th)*dsqrt(1.d0+anu**2)
            xp = (1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dcos(th)/den
            yp = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dsin(th)/den
            xdot = -(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dsin(th)/den
     *             -xp*ddot/den            
            xddot = -xp 
     *              +(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)
     *                            *dsin(th)*ddot/den**2
     *              - xdot*ddot/den - xp*dddot/den 
     *              + xp*ddot**2/den**2
            ydot = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dcos(th)/den
     *             -yp*ddot/den
            yddot = -yp 
     *              -(1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)
     *                            *dcos(th)*ddot/den**2
     *              - ydot*ddot/den - yp*dddot/den 
     *              + yp*ddot**2/den**2 
           else
            eps = bi
            cs = dcos(th)
            sn = dsin(th)
            snn = dsin(NCYC*th)
            csn = dcos(NCYC*th)
            xp = ai*cs + bi*csn
            yp = ai*sn - bi*snn
            xdot = -ai*sn - NCYC*bi*snn
            xddot = -ai*cs - NCYC**2*bi*csn
            ydot = ai*cs - NCYC*bi*csn
            yddot = -ai*sn + NCYC**2*bi*snn
         end if
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine MAKE_GEO (k, nd, nbk, ak, bk, zk, ncyc, x, y, z,  
     1                     dz, dsdth, rkappa)
c---------------
c constructs contour geometry
c
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), ncyc(k)
      dimension x(nbk), y(nbk), dsdth(nbk), rkappa(nbk)
      complex*16 zk(k), z(nbk), dz(nbk)
c
         pi = 4.d0*datan(1.d0)
c
c open matlab files for plotting contour geometry
         open (unit = 21, file = 'geometry.m')
         open (unit = 22, file = 'curvature.m')
c
         dth = 2.d0*pi/nd
         istart = 0
         do kbod = 1, k
            if (ncyc(kbod).eq.0) then
               area = pi*ak(kbod)*bk(kbod)
              else
               area = pi*(ak(kbod)**2 - ncyc(kbod)*bk(kbod)**2)
            end if
            call PRINF ('+++++  NBOD = *',kbod,1)
            call PRIN2 ('   AREA = *',area,1)
            do i = 1, nd
               th = dth*(i-1.d0)
c
c            thetax allows for rotation of ellipses, ignore for now
               thetax = 0.d0
               call CURVE_PARAM (th, ak(kbod), bk(kbod), thetax,
     1                           ncyc(kbod), xp, yp, xdot, ydot,
     2                           xddot, yddot)
c
               z(istart+i) = zk(kbod) + dcmplx(xp,yp)
               x(istart+i) = dreal(z(istart+i))
               y(istart+i) = dimag(z(istart+i))
               dsdth(istart+i) = dsqrt((xdot)**2 + (ydot)**2)
c
c Curve oriention, 
c kbod = 1 is assumed to be the bounding curve, so it should be 
c         traversed counter clockwise
c kbod > 1 is assumed to be an interior curve, so it should be 
c         traversed clockwise
c note that this orientation affects dz and rkappa only. For convenience
c z will always have points listed in a counter clockwise direction
c     
               if (kbod.eq.1) then
                  dz(istart+i) = dcmplx(xdot,ydot)
                  rkappa(istart+i) = (xdot*yddot-ydot*xddot)/
     *                                dsdth(istart+i)**3
                else
                  dz(istart+i) = -dcmplx(xdot,ydot)
                  rkappa(istart+i) = -(xdot*yddot-ydot*xddot)/
     *                                 dsdth(istart+i)**3
               endif 
            end do
            call RSCPLOT (z(istart+1), nd, kbod, 21)
            call RS1PLOT (rkappa(istart+1), nd, kbod, 22)
            istart = istart + nd
         end do
c
         close (21)
         close (22)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GET_TARGET (ntar, xtar, ytar, ztar)
c---------------
c Get target points where we will check exact solution
c
      implicit real*8 (a-h,o-z)
      dimension xtar(ntar), ytar(ntar)
      complex*16 ztar(ntar), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         dth = 2*pi/ntar
c
         do i = 1, ntar
            th = (i-1)*dth
            ztar(i) = 0.25d0*cdexp(eye*th)
            xtar(i) = dreal(ztar(i))
            ytar(i) = dimag(ztar(i))
         end do
c
         open (unit = 21, file = 'targets.m')
         call RSCPLOT (ztar, ntar, 1, 21)
         close(1)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine POLAR_COORD (z, theta, r)      
c---------------
c
c  Given a point z = x + i y, calculate polar coords r and theta
c
      implicit real*8 (a-h, o-z)
      complex*16 z
c
        pi = 4.d0*datan(1.d0)
c
        r = cdabs(z)
        th = dacos(dreal(z)/r)
        if (dimag(z).lt.0.d0) th = 2.d0*pi - th
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      double precision function U_EXACT (k, zk, z)     
c---------------
c  This is an exact test function. It should be harmonic and well 
c  defined in the domain
c  It's used to generate boundary conditions and to test accuracy of 
c  solutions
      implicit real*8 (a-h,o-z)
      complex*16 zk(k), z, zsrc, zdis
c
c  itest is a flag that can be used to switch test functions
         itest = 1
c
         if (itest.eq.1) then 
c
c     u(x,y) = log|z-zsrc|
c     Define a source point - this should be out of the domain
            zsrc = dcmplx(2.d0,2.d0)
            zdis = z-zsrc
            u = cdlog(zdis)
           elseif (itest.eq.2) then
c
c     u(x,y) = r^n cos(n theta), n should be positive integer for
c     bounded domains
            zdis = z - zk(1)
            n = 2
            call POLAR_COORD (zdis, theta, r)
            u = (r**n)*dcos(theta*n)
         end if
         U_EXACT = u
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GET_BCS (k, nd, nbk, zk, z, rhs)
c---------------
c Get Dirichlet boundary conditions according to U_EXACT
c
      implicit real*8 (a-h,o-z)
      dimension rhs(*)
      complex*16 zk(k), z(nbk)
c
         pi = 4.d0*datan(1.d0)
c         
         do i = 1, nbk
	    rhs(i) = U_EXACT (k, zk, z(i))
         end do
c 
      return
      end
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine SOLVE (nd, k, nbk, maxl, rhs, soln, u, gmwork, lrwork, 
     1                  igwork, liwork)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC, MSOLVE
c
c  System
c
      dimension soln(*), rhs(*), u(*)
c
c  DGMRES work arrays
c
      dimension gmwork(lrwork), igwork(liwork)
c
c  Timings
c
      real*4 timep(2), etime
c
         pi = 4.d0*datan(1.d0)
c
c  solve linear system using GMRES.
c
c
c     parameters for DGMRES - see dgmres.f to see what they mean
         itol = 0
         tol = 1.0d-12
         isym = 0
         do i=2,liwork
            igwork(i) = 0
         enddo
         norder = nbk
c
c  Preconditioner flag 
c
         igwork(4) = 0
c
c  Restart flag, will restart after iwork(5) iterations. This value 
c  should be less than or equal to maxl
         igwork(5) = maxl      
c
c  provide initial guess soln
         do i=1,norder
            soln(i) = rhs(i)
         enddo
c
         t0 = etime(timep)
         call DGMRES (norder, rhs, soln, nelt, ia, ja, a, isym,
     1                MATVEC, MSOLVE, itol, tol, itmax, iter, err,  
     1                ierr, 6, sb, sx, gmwork, lrwork, igwork, 
     1                liwork, rw, iw)
         call PRINI (6,13)
         call PRINF ('  # GMRES ITERATIONS = *',iter,1)
         if (ierr.gt.2) then
            call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
            call PRINF ('  iwork = *',iwork,10)
            stop
           elseif (ierr.ge.0) then
            t1 = etime(timep)
            tsec = t1 - t0
            call PRIN2 (' TIME TAKEN IN GMRES = *', tsec, 1)
         end if
c
c  unpack RHS into U
c
            do i = 1,nbk
               u(i) = soln(i)
            end do 
ccc            call PRIN2 (' u = *', u, nbk) 
c
      return
      end
c
c
c---------------
      subroutine FASMVP (k, nd, nbk, zk, x, y, z, dz, rkappa, dsdth, 
     1                   u, w, nsp, charge, poten, cfield, wksp)
c---------------
c  Calculates the matrix vector product
c     u is the current guess for the density
c     w is (0.5I + K) u
c
      implicit real*8 (a-h,o-z)
      dimension x(nbk), y(nbk), rkappa(nbk), dsdth(nbk), u(*), w(*)
      complex*16 zk(k), z(nbk), dz(nbk), eye, delz, zn, zcauchy
c FMM

	integer nsource,ntarget,ifpot,ifgrad,ifcharge,ifdipole,ifhess, 
     1	     ifpottarg,ifgradtarg,ifhesstarg,iprec,ier,i,j,k,l
	real (kind=8) source(2,1000),dipvec(2,1000), 
     1			     targ(2,1000)		
	complex*16  charge(nbk),dipstr(1000),pot(1000), 		
     1		   grad(2,1000),hess(3,1000),pottarg(1000), 
     1		   gradtarg(2,1000),hesstarg(3,1000) 
	
c local variables
	complex*16 ztar(1000)
	dimension xtar(1000),ytar(1000)

c
c
       pi = 4.d0*datan(1.d0)
	 eye = dcmplx(0.d0,1.d0)
	 z2pii = 1.d0/(2.d0*pi*eye)
	 h = 2*pi/nd
c 
c set density for fmm call
       istart = 0
       do kbod = 1, k
		do i = 1, nd
            charge(istart+i) = u(istart+i)*dz(istart+i)*h
            end do
		istart = istart + nd
       end do

	nsource = nsp
	
c Getting target points
	call GET_TARGET(ntarget,xtar,ytar,ztar) 
	do i = 1,ntarget
		targ(1,i) = xtar(i)
		targ(2,i) = ytar(i)
	end do

C     set parameters for FMM routine DAPIF2

	
	 iprec = 3
	 ifcharge =1
	 ifdipole = 0
	 ifpot = 1
	 ifgrad = 1
	 ifhess = 0
	 ifpottarg = 1
	 ifgradtarg = 1
	 ifhesstarg = 0
	 j = 1
	 l = 1 
		

		
C      set source points
	 nsource = nbk	
	 do j = 1,nsource
			 source(1,j) = x(j)
			 source(2,j) = y(j)
	 end do
		


		
	 call lfmm2dparttarg(ier,iprec,nsource,source,ifcharge,charge, 
     1                      ifdipole,dipstr,dipvec,ifpot,pot,ifgrad, 
     1                      grad,ifhess,hess,ntarget,targ, 
     1                      ifpottarg,pottarg,ifgradtarg,gradtarg,
     1			    ifhesstarg,hesstarg)

	 print *,"Number of targets:",ntarget
	 print *,"Potential at target 10", pottarg(9)	
	 call prin2("Potential at Target points:",pottarg,ntarget)
	
C
         
C	 call DAPIF2 (IOUT,IFLAG7,NBK,NAPB,NINIRE,MEX,IERR,INFORM,
C    *                TOL,EPS7,X,Y,QA,POTEN,CFIELD,WKSP,NSP,CLOSE)
C         call PRINI(6,13)
         if (ier.ne.0) then
            write (6,*) '  IERR FROM FMM = ',IER
            stop
         end if
c
c  discrete integral operator
c         istart = 0
c         do kbod = 1, k
c	    do i = 1, nd
c		 self = 0.25d0*h*rkappa(istart+i)*dsdth(istart+i)/pi
c               zcauchy = self*u(istart+i) -
c     1                          dreal(z2pii*cfield(istart+i))
c               w(istart+i) = 0.5d0*u(istart+i) + dreal(zcauchy)
c            end do
c            istart = istart + nd
c        end do
c
      return
      end
c
c
c---------------
      subroutine MATVEC (N, XX, YY, NELT, IA, JA, A, ISYM)
c---------------
c
c  Required by DGMRES with this precise calling sequence.
c  We ignore most of the parameters except N, XX and YY
c  and call the fast multipole routine FASMVP to do the actual
c  work.
c
      implicit double precision (a-h,o-z)
      dimension xx(n), yy(n)
c
c The following parameters need to be the same size as in the main 
c routine to insure common block arrays are dimensioned correctly
      parameter (kmax = 1, npmax = 2048, nmax = kmax*npmax)
c
c common blocks
      common /geometry/ x, y, zk, z, dz, dsdth, rkappa
      common /sys_size/ k, nd, nbk
c
      dimension dsdth(nmax), rkappa(nmax), x(nmax), y(nmax)
      complex*16 zk(kmax), dz(nmax), z(nmax)
c
c  FMM workspace arrays 
      parameter (nsp = 16*nmax+1000)
      COMPLEX*16 QA(NMAX), CFIELD(NMAX), WKSP(NSP)
      DIMENSION POTEN(NMAX)
c
c  local work arrays
      real*4 timep(2), etime
c
         t0 = etime(timep)
c         
         call FASMVP (k, nd, nbk, zk, x, y, z, dz, rkappa, dsdth, 
     1                xx, yy, nsp, qa, poten, cfield, wksp)
         t1 = etime(timep)
         tsec = t1 - t0
ccc         WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
ccc         WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
c
      RETURN
      END
c
c
c---------------
      subroutine MSOLVE(N, R, U, NELT, IA, JA, A, ISYM, RWORK, IWORK)
c---------------
c
c  Another routine required by DGMRES. It allows the use of a
c  preconditioner.
c
      implicit double precision (A-H,O-Z)
c
          write (6,*) 'in msolve'
c
      RETURN
      END
c
c
c---------------
      subroutine GET_SOL_TAR (k, nd, nbk, ntar, density, zk, z, dz,  
     1                        ztar, u_tar)
c---------------
c
c  Calculate solution at target points and check accuracy
c
      implicit double precision (a-h, o-z)
      dimension u_tar(ntar), density(nbk)
      complex*16 z(nbk), dz(nbk), ztar(ntar), eye, zcauchy, z2pii, 
     1           zk(k)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         z2pii = 1.d0/(2.d0*pi*eye)
         dth = 2.d0*pi/nd
c
c
         err = 0.d0
         do itar = 1, ntar
            u_tar(itar) = 0.d0
            istart = 0
            do i = 1, nbk
               zcauchy = density(i)*dz(i)/(z(i) - ztar(itar))
               zcauchy = dth*zcauchy*z2pii
               u_tar(itar) = u_tar(itar) + dreal(zcauchy)
            end do
            u_ex = U_EXACT(k, zk, ztar(itar))
            err = max(err,dabs(u_ex-u_tar(itar)))
         end do
         call PRIN2 (' MAX ERROR AT TARGET POINTS = *', err, 1)
c
      return
      end 
      
                  
