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
      complex*16 qa(nmax)
c
c Exact potenial and field
      dimension poten_ex(nmax)
      complex*16 cfield_ex(nmax), zi, zj, eye, zdis
c
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
         do ix = -50, 50
            do iy = -50, 50
               xat(icount+1) = float(ix)
               yat(icount+1) = float(iy)
               icount = icount+1
            end do
         end do
         nat = icount
         call PRINF (' NUMBER OF SOURCES = *', nat, 1)
ccc         call prin2 (' X COORD OF SOURCES = *', xat, nat)
ccc         call prin2 (' Y COORD OF SOURCES = *', yat, nat)
c
c Again, for lack of a better idea, define the source densities to be
c the x-coordinate value
         do i = 1, nat
            qa(i) = xat(i)
         end do
ccc         call PRIN2 (' SOURCE STRENGTHS = *', qa, nat)
c
c Calculate exact potential and field
c let z_i = (x_i,y_i) or x_i + I y_i
c     potential at z_i = sum_j qa log|z_i - z_j|
c     field at z_i = sum_j qa (z_i - z_j)/|z_i-z_j|^2
         tbeg = etime(timep)
         do i = 1, nat
            poten_ex(i) = 0.d0
            cfield_ex(i) = 0.d0
            zi = dcmplx(xat(i),yat(i))
            do j = 1, nat
               zj = dcmplx(xat(j),yat(j))
               zdis = zi - zj
               dis = cdabs(zdis)
               if (j.ne.i) then
                  poten_ex(i) = poten_ex(i) + qa(j)*dlog(dis)
                  cfield_ex(i) = cfield_ex(i) + qa(j)*zdis/dis**2
               end if
            end do
         end do
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
     2                cfield, wksp, nsp, CLOSE)
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
        do i = 1, nat
           err_pot = max(err_pot,dabs(poten(i)-poten_ex(i)))
           err_field = max(err_field,cdabs(cfield(i)-cfield_ex(i)))
        end do
        call prin2 (' error in potential = *', err_pot, 1)
        call prin2 (' error in field = *', err_field, 1)
ccc        call PRIN2 (' Exact Potential = *', poten_ex, nat)
ccc        call PRIN2 (' FMM Potential = *', poten, nat) 
ccc        call PRIN2 (' Exact Field = *', cfield_ex, 2*nat)
ccc        call PRIN2 (' FMM Field = *', cfield, 2*nat)
c
      stop
      end
         

