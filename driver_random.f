C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PROGRAM LUXTST
C         Exercise for the RANLUX Pseudorandom number generator.
C
      DIMENSION RVEC(1000)
      DIMENSION ISDEXT(25)
C
C         check that we get the right numbers (machine-indep.)
      WRITE (6,'(/A)')  '  CALL RANLUX(RVEC,100)'
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers 101-105:',
     +    (RVEC(L),L=1,5)
C
      WRITE (6,'(/A)')  ' CALL RLUXGO(0,0,0,0)'
      CALL RLUXGO(0,0,0,0)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0,   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0, 101-105:',
     +    (RVEC(L),L=1,5)
C
      WRITE (6,'(/A)')  '   CALL RLUXGO(389,1,0,0)'
      CALL RLUXGO(389,1,0,0)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389,   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389, 101-105:',
     +    (RVEC(L),L=1,5)
C
      WRITE (6,'(/A)')  '  CALL RLUXGO(75,0,0,0)'
      CALL RLUXGO(75,0,0,0)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75,   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75, 101-105:',
     +    (RVEC(L),L=1,5)
C
      WRITE (6,'(/A)')  '  test restarting from the full vector'
      CALL RLUXUT(ISDEXT)
      WRITE (6,'(/A/(1X,5I14))') '  current RANLUX status saved:',ISDEXT
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',
     +    (RVEC(L),L=1,5)
C
      WRITE (6,'(/A)')   '   previous RANLUX status will be restored'
      CALL RLUXIN(ISDEXT)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',
     +    (RVEC(L),L=1,5)
C
      WRITE (6,'(/A)')  '     test the restarting by skipping'
      CALL RLUXGO(4,7674985,0,0)
      CALL RLUXAT(I1,I2,I3,I4)
      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4
      DO 150 LI= 1, 10
  150 CALL RANLUX(RVEC,1000)
      CALL RLUXAT(I1,I2,I3,I4)
      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4
      CALL RANLUX(RVEC,200)
      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',
     +                             RVEC(1), RVEC(200)
      CALL RLUXGO(I1,I2,I3,I4)
      CALL RANLUX(RVEC,200)
      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',
     +                             RVEC(1), RVEC(200)
C
      WRITE (6,'(/A)') ' The following should provoke an error message'
      CALL RLUXGO(4,11111,31,0)
      STOP
C
C   OUTPUT FROM THE ABOVE TEST PROGRAM SHOULD BE:
C   --------------------------------------------
C  CALL RANLUX(RVEC,100)
C RANLUX DEFAULT INITIALIZATION:    314159265
C RANLUX DEFAULT LUXURY LEVEL =   3      p = 223
C RANLUX default numbers   1-  5:
C           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
C RANLUX default numbers 101-105:
C           0.43156743  0.03774416  0.24897110  0.00147784  0.90274453
C
C  CALL RLUXGO(0,0,0,0)
C RANLUX LUXURY LEVEL SET BY RLUXGO : 0     P=  24
C RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
C RANLUX luxury level 0,   1-  5:
C           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
C RANLUX luxury level 0, 101-105:
C           0.41538775  0.05330932  0.58195311  0.91397446  0.67034441
C
C   CALL RLUXGO(389,1,0,0)
C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
C RANLUX INITIALIZED BY RLUXGO FROM SEEDS           1           0           0
C RANLUX luxury p=389,   1-  5:
C           0.94589490  0.47347850  0.95152789  0.42971975  0.09127384
C RANLUX luxury p=389, 101-105:
C           0.02618265  0.03775346  0.97274780  0.13302165  0.43126065
C
C  CALL RLUXGO(75,0,0,0)
C RANLUX P-VALUE SET BY RLUXGO TO:   75
C RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
C RANLUX luxury p= 75,   1-  5:
C           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
C RANLUX luxury p= 75, 101-105:
C           0.25600731  0.23443210  0.59164381  0.59035838  0.07011414
C
C  test restarting from the full vector
C
C  current RANLUX status saved:
C       16156027      16534309      15243811       2751687       6002207
C        7979506       1301976       4567313       4305996       5872599
C       12003090       2146823      12606367       4111505       5979640
C       12739666      10489318      14036909      11729352       8061448
C        7832659       6069758       3197719       1832730      75080216
C RANLUX numbers 1- 5:
C           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
C RANLUX numbers 101-105:
C           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
C
C   previous RANLUX status will be restored
C FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:
C         16156027    16534309    15243811     2751687     6002207
C          7979506     1301976     4567313     4305996     5872599
C         12003090     2146823    12606367     4111505     5979640
C         12739666    10489318    14036909    11729352     8061448
C          7832659     6069758     3197719     1832730    75080216
C RANLUX P-VALUE SET BY RLUXIN TO:   75
C RANLUX numbers 1- 5:
C           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
C RANLUX numbers 101-105:
C           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
C
C     test the restarting by skipping
C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
C RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985           0           0
C  RLUXAT values =         4   7674985         0         0
C  RLUXAT values =         4   7674985    161840         0
C  Next and 200th numbers are:  0.019648  0.590586
C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
C RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985      161840           0
C  Next and 200th numbers are:  0.019648  0.590586
C
C The following should provoke an error message
C RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
C RANLUX INITIALIZED BY RLUXGO FROM SEEDS       11111          31           0
C  Error in RESTARTING with RLUXGO:
C  The values      11111         31          0 cannot occur at luxury level    4
      END

************************************************************************


