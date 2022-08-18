      module newuoa_mod
      use model_mod, only: CALFUN
      use eco_params, only: nparams_opt
      use io, only: results_un    
      
      contains  

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newuoa.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(nparams_opt),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+11)*(NPT+N)+N*(3*N+11)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IGQ=IXP+N*NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
      CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),
     2  NDIM,W(ID),W(IVL),W(IW))
   20 RETURN
      END SUBROUTINE

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newuob.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
     1  XOPT,XNEW,XPT,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(nparams_opt) 
      DIMENSION XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),GQ(*),
     1  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),VLAG(*),W(*)
      integer :: niter=1 !Luo
C
C     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
C       to the corresponding arguments in SUBROUTINE NEWUOA.
C     XBASE will hold a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XOPT will be set to the displacement from XBASE of the vector of
C       variables that provides the least calculated F so far.
C     XNEW will be set to the displacement from XBASE of the vector of
C       variables for the current calculation of F.
C     XPT will contain the interpolation point coordinates relative to XBASE.
C     GQ will hold the gradient of the quadratic model at XBASE.
C     HQ will hold the explicit second derivatives of the quadratic model.
C     PQ will contain the parameters of the implicit second derivatives of
C       the quadratic model.
C     BMAT will hold the last N columns of H.
C     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
C       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
C       the elements of DZ are plus or minus one, as specified by IDZ.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     D is reserved for trial steps from XOPT.
C     VLAG will contain the values of the Lagrange functions at a new point X.
C       They are part of a product that requires VLAG to be of length NDIM.
C     The array W will be used for working space. Its length must be at least
C       10*NDIM = 10*(NPT+N).
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NFTEST=MAX0(MAXFUN,1)
C
C     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
C
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 J=1,(N*NP)/2
   30 HQ(J)=ZERO
      DO 40 K=1,NPT
      PQ(K)=ZERO
      DO 40 J=1,NPTM
   40 ZMAT(K,J)=ZERO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF,.).
C
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      NF=0
   50 NFM=NF
      NFMM=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              XPT(NF,NFM)=RHOBEG
          ELSE IF (NFM .GT. N) THEN
              XPT(NF,NFMM)=-RHOBEG
          END IF
      ELSE
          ITEMP=(NFMM-1)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT .GT. N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XIPT=RHOBEG
          IF (W(IPT+NP) .LT. W(IPT+1)) XIPT=-XIPT
          XJPT=RHOBEG
          IF (W(JPT+NP) .LT. W(JPT+1)) XJPT=-XJPT
          XPT(NF,IPT)=XIPT
          XPT(NF,JPT)=XJPT
      END IF
C
C     Calculate the next value of F, label 70 being reached immediately
C     after this calculation. The least function value so far and its index
C     are required.
C
      DO 60 J=1,N
   60 X(J)=XPT(NF,J)+XBASE(J)
      GOTO 310
   70 W(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          FOPT=F
          KOPT=1
      ELSE IF (F .LT. FOPT) THEN
          FOPT=F
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in
C     the cases when NF is at most 2*N+1.
C
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              GQ(NFM)=(F-FBEG)/RHOBEG
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/RHOBEG
                  BMAT(NF,NFM)=ONE/RHOBEG
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NFM .GT. N) THEN
              BMAT(NF-N,NFMM)=HALF/RHOBEG
              BMAT(NF,NFMM)=-HALF/RHOBEG
              ZMAT(1,NFMM)=-RECIQ-RECIQ
              ZMAT(NF-N,NFMM)=RECIQ
              ZMAT(NF,NFMM)=RECIQ
              IH=(NFMM*(NFMM+1))/2
              TEMP=(FBEG-F)/RHOBEG
              HQ(IH)=(GQ(NFMM)-TEMP)/RHOBEG
              GQ(NFMM)=HALF*(GQ(NFMM)+TEMP)
          END IF
C
C     Set the off-diagonal second derivatives of the Lagrange functions and
C     the initial quadratic model.
C
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          IF (XIPT .LT. ZERO) IPT=IPT+N
          IF (XJPT .LT. ZERO) JPT=JPT+N
          ZMAT(1,NFMM)=RECIP
          ZMAT(NF,NFMM)=RECIP
          ZMAT(IPT+1,NFMM)=-RECIP
          ZMAT(JPT+1,NFMM)=-RECIP
          HQ(IH)=(FBEG-W(IPT+1)-W(JPT+1)+F)/(XIPT*XJPT)
      END IF
      IF (NF .LT. NPT) GOTO 50
C
C     Begin the iterative procedure, because the initial model is complete.
C
      RHO=RHOBEG
      DELTA=RHO
      IDZ=1
      DIFFA=ZERO
      DIFFB=ZERO
      XOPTSQ=ZERO
      DO 80 I=1,N
      XOPT(I)=XPT(KOPT,I)
   80 XOPTSQ=XOPTSQ+XOPT(I)**2
   90 NFSAV=NF
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve the model.
C
  100 KNEW=0
      CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN)
      DSQ=ZERO
      DO 110 I=1,N
  110 DSQ=DSQ+D(I)**2
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          KNEW=-1
          DELTA=TENTH*DELTA
          RATIO=-1.0D0
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          IF (NF .LE. NFSAV+2) GOTO 460
          TEMP=0.125D0*CRVMIN*RHO*RHO
          IF (TEMP .LE. DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
          GOTO 490
      END IF
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C     to BMAT that do not depend on ZMAT.
C
  120 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          TEMPQ=0.25D0*XOPTSQ
          DO 140 K=1,NPT
          SUM=ZERO
          DO 130 I=1,N
  130     SUM=SUM+XPT(K,I)*XOPT(I)
          TEMP=PQ(K)*SUM
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          DO 140 I=1,N
          GQ(I)=GQ(I)+TEMP*XPT(K,I)
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          VLAG(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
          IP=NPT+I
          DO 140 J=1,I
  140     BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 180 K=1,NPTM
          SUMZ=ZERO
          DO 150 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
  150     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 170 J=1,N
          SUM=TEMPQ*SUMZ*XOPT(J)
          DO 160 I=1,NPT
  160     SUM=SUM+W(I)*XPT(I,J)
          VLAG(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 170 I=1,NPT
  170     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 180 I=1,N
          IP=I+NPT
          TEMP=VLAG(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 180 J=1,I
  180     BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
C
C     The following instructions complete the shift of XBASE, including
C     the changes to the parameters of the quadratic model.
C
          IH=0
          DO 200 J=1,N
          W(J)=ZERO
          DO 190 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  190     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 200 I=1,J
          IH=IH+1
          IF (I .LT. J) GQ(J)=GQ(J)+HQ(IH)*XOPT(I)
          GQ(I)=GQ(I)+HQ(IH)*XOPT(J)
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  200     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 210 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
  210     XOPT(J)=ZERO
          XOPTSQ=ZERO
      END IF
C
C     Pick the model step if KNEW is positive. A different choice of D
C     may be made later, if the choice of D by BIGLAG causes substantial
C     cancellation in DENOM.
C
      IF (KNEW .GT. 0) THEN
          CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,KNEW,
     1      DSTEP,D,ALPHA,W,W(NP),W(NDIM+1))
      END IF
C
C     Calculate VLAG and BETA for the current choice of D. The first NPT
C     components of W_check will be held in W.
C
      DO 230 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 220 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  220 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
  230 VLAG(K)=SUM
      BETA=ZERO
      DO 250 K=1,NPTM
      SUM=ZERO
      DO 240 I=1,NPT
  240 SUM=SUM+ZMAT(I,K)*W(I)
      IF (K .LT. IDZ) THEN
          BETA=BETA+SUM*SUM
          SUM=-SUM
      ELSE
          BETA=BETA-SUM*SUM
      END IF
      DO 250 I=1,NPT
  250 VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      DO 280 J=1,N
      SUM=ZERO
      DO 260 I=1,NPT
  260 SUM=SUM+W(I)*BMAT(I,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 270 K=1,N
  270 SUM=SUM+BMAT(JP,K)*D(K)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  280 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is positive and if the cancellation in DENOM is unacceptable,
C     then BIGDEN calculates an alternative model step, XNEW being used for
C     working space.
C
      IF (KNEW .GT. 0) THEN
          TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
          IF (DABS(TEMP) .LE. 0.8D0) THEN
              CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1          KNEW,DSTEP,D,VLAG,BETA,XNEW,W,W(5*NDIM+1),W)
          END IF
      END IF
C
C     Calculate the next value of the objective function.
C
  290 DO 300 I=1,N
      XNEW(I)=XOPT(I)+D(I)
  300 X(I)=XBASE(I)+XNEW(I)
      NF=NF+1
  310 IF (NF .GT. NFTEST) THEN
          NF=NF-1
          IF (IPRINT .GT. 0) PRINT 320
  320     FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 530
      END IF
      CALL CALFUN(N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 330, NF,F,(X(I),I=1,N)
  330      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .LE. NPT) GOTO 70
      IF (KNEW .EQ. -1) GOTO 530
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and set DIFF to the error of this prediction.
C
      VQUAD=ZERO
      IH=0
      DO 340 J=1,N
      VQUAD=VQUAD+D(J)*GQ(J)
      DO 340 I=1,J
      IH=IH+1
      TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
      IF (I .EQ. J) TEMP=HALF*TEMP
  340 VQUAD=VQUAD+TEMP*HQ(IH)
      DO 350 K=1,NPT
  350 VQUAD=VQUAD+PQ(K)*W(K)
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. The branch when KNEW is positive occurs if D is not
C     a trust region step.
C
      FSAVE=FOPT
      IF (F .LT. FOPT) THEN
          FOPT=F
          XOPTSQ=ZERO
          DO 360 I=1,N
          XOPT(I)=XNEW(I)
  360     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      KSAVE=KNEW
      IF (KNEW .GT. 0) GOTO 410
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD .GE. ZERO) THEN
          IF (IPRINT .GT. 0) PRINT 370
  370     FORMAT (/4X,'Return from NEWUOA because a trust',
     1      ' region step has failed to reduce Q.')
          GOTO 530
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO .LE. TENTH) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
      END IF
      IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
      KTEMP=0
      DETRAT=ZERO
      IF (F .GE. FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO 400 K=1,NPT
      HDIAG=ZERO
      DO 380 J=1,NPTM
      TEMP=ONE
      IF (J .LT. IDZ) TEMP=-ONE
  380 HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
      TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
      DISTSQ=ZERO
      DO 390 J=1,N
  390 DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
      IF (DISTSQ .GT. RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)**3
      IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
          DETRAT=TEMP
          KNEW=K
      END IF
  400 CONTINUE
      IF (KNEW .EQ. 0) GOTO 460
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C     can be moved. Then make this move, and update the quadratic model,
C
  410 CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      IH=0
      DO 420 I=1,N
      TEMP=PQ(KNEW)*XPT(KNEW,I)
      DO 420 J=1,I
      IH=IH+1
  420 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      PQ(KNEW)=ZERO
      DO 440 K=1,NPTM
      IF (ZMAT(KNEW,K) .NE. ZERO) THEN
          TEMP=DIFF*ZMAT(KNEW,K)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 430 J=1,NPT
  430     PQ(J)=PQ(J)+TEMP*ZMAT(J,K)
      END IF
  440 CONTINUE
      DO 450 I=1,N
      XPT(KNEW,I)=XNEW(I)
  450 GQ(I)=GQ(I)+DIFF*BMAT(KNEW,I)
      IF (F .LT. FSAVE) KOPT=KNEW
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case KSAVE>0 occurs
C     when the new function value was calculated by a model step.
C
      IF (F .LE. FSAVE+TENTH*VQUAD) GOTO 100
      IF (KSAVE .GT. 0) GOTO 100
C
C     Alternatively, find out if the interpolation points are close enough
C     to the best point so far.
C
      KNEW=0
  460 DISTSQ=4.0D0*DELTA*DELTA
      DO 480 K=1,NPT
      SUM=ZERO
      DO 470 J=1,N
  470 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  480 CONTINUE
C
C     If KNEW is positive, then set DSTEP, and branch back for the next
C     iteration, which will generate a "model step".
C
      IF (KNEW .GT. 0) THEN
          DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
          DSQ=DSTEP*DSTEP
          GOTO 120
      END IF
      IF (RATIO .GT. ZERO) GOTO 100
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 100
C
C     The calculations with the current value of RHO are complete. Pick the
C     next values of RHO and DELTA.
C
      write(results_un,'(102(1x,1PG19.12,:))') niter,F,X !Luo
      niter = niter + 1 !Luo
      
  490 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 500
  500         FORMAT (5X)
              PRINT 510, RHO,NF
  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 520, FOPT,(XBASE(I)+XOPT(I),I=1,N)
  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 90
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C     it is too short to have been tried before.
C
      IF (KNEW .EQ. -1) GOTO 290
  530 IF (FOPT .LE. F) THEN
          DO 540 I=1,N
  540     X(I)=XBASE(I)+XOPT(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
C         WRITE (7,550) NF
          PRINT 550, NF
  550     FORMAT (/4X,'At the return from NEWUOA',5X,
     1      'Number of function values =',I6)
C         WRITE (7,520) F,(X(I),I=1,N)
          PRINT 520, F,(X(I),I=1,N)
      END IF
      RETURN
      END SUBROUTINE

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bigden.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,DELTA,D,VLAG,BETA,GW,W,WVEC,PROD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  VLAG(*),GW(*),W(*),WVEC(NDIM,*),PROD(NDIM,*)
      DIMENSION DEN(9),DENEX(9),WW(9)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KOPT is the index of the optimal interpolation point.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DELTA is the current trust region bound.
C     D will be set to the step from XOPT to the new point.
C     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
C     BETA will be set to the value that will occur in the updating formula
C       when the KNEW-th interpolation point is moved to its new position.
C     GW, W, WVEC, PROD and the private arrays DEN, DENEX and WW will be
C       used for working space, but on return W will be set to W_check.
C       The space for W may be part of the space for PROD.
C
C     D is calculated in a way that should provide a denominator with a large
C     modulus in the updating formula when the KNEW-th interpolation point is
C     shifted to the new position XOPT+D.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      QUART=0.25D0
      TWO=2.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      NPTM=NPT-N-1
      DSQ=DELTA*DELTA
C
C     Set the initial D and G, where G is the gradient of the KNEW-th
C     Lagrange function at the trust region centre. The gradient of this
C     function at XPT(KNEW,.) is set in W. The array VLAG will hold the
C     second derivative coefficients of the Lagrange function.
C
      DO 10 I=1,N
      D(I)=XPT(KNEW,I)-XOPT(I)
      GW(I)=BMAT(KNEW,I)
   10 W(I)=BMAT(KNEW,I)
      DO 20 K=1,NPT
   20 VLAG(K)=ZERO
      DO 30 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 30 K=1,NPT
   30 VLAG(K)=VLAG(K)+TEMP*ZMAT(K,J)
      DO 50 K=1,NPT
      SUM=ZERO
      TEMPB=ZERO
      DO 40 J=1,N
      SUM=SUM+XPT(K,J)*XOPT(J)
   40 TEMPB=TEMPB+XPT(K,J)*XPT(KNEW,J)
      IF (K .EQ. KOPT) XOPTSQ=SUM
      TEMPA=VLAG(K)*SUM
      TEMPB=VLAG(K)*TEMPB
      DO 50 I=1,N
      GW(I)=GW(I)+TEMPA*XPT(K,I)
   50 W(I)=W(I)+TEMPB*XPT(K,I)
      ALPHA=VLAG(KNEW)
C
C     Revise G if its modulus seems to be unusually small.
C
      TEMP=ZERO
      TEMPA=ZERO
      TEMPB=ZERO
      DO 60 I=1,N
      TEMP=TEMP+D(I)**2
      TEMPA=TEMPA+GW(I)**2
   60 TEMPB=TEMPB+W(I)**2
      IF (TEMPB*DSQ .GT. 1.0D4*TEMP*TEMPA) THEN
          DO 70 I=1,N
   70     GW(I)=W(I)
      END IF
C
C     Begin the iteration by making D and G orthogonal and of length DELTA.
C
      ITERC=0
   80 ITERC=ITERC+1
      DD=ZERO
      DG=ZERO
      DO 90 I=1,N
      DD=DD+D(I)**2
   90 DG=DG+D(I)*GW(I)
      GG=ZERO
      DO 100 I=1,N
      GW(I)=DD*GW(I)-DG*D(I)
  100 GG=GG+GW(I)**2
      TEMPD=DELTA/DSQRT(DD)
      TEMPG=DELTA/DSQRT(GG)
      XOPTD=ZERO
      XOPTG=ZERO
      DO 110 I=1,N
      D(I)=TEMPD*D(I)
      GW(I)=TEMPG*GW(I)
      XOPTD=XOPTD+XOPT(I)*D(I)
  110 XOPTG=XOPTG+XOPT(I)*GW(I)
C
C     Set the coefficients of the first two terms of BETA.
C
      TEMPA=HALF*XOPTD*XOPTD
      TEMPB=HALF*XOPTG*XOPTG
      DEN(1)=DSQ*(XOPTSQ+HALF*DSQ)+TEMPA+TEMPB
      DEN(2)=TWO*XOPTD*DSQ
      DEN(3)=TWO*XOPTG*DSQ
      DEN(4)=TEMPA-TEMPB
      DEN(5)=XOPTD*XOPTG
      DO 120 I=6,9
  120 DEN(I)=ZERO
C
C     Put the coefficients of Wcheck in WVEC.
C
      DO 140 K=1,NPT
      TEMPA=ZERO
      TEMPB=ZERO
      TEMPC=ZERO
      DO 130 I=1,N
      TEMPA=TEMPA+XPT(K,I)*D(I)
      TEMPB=TEMPB+XPT(K,I)*GW(I)
  130 TEMPC=TEMPC+XPT(K,I)*XOPT(I)
      WVEC(K,1)=QUART*(TEMPA*TEMPA+TEMPB*TEMPB)
      WVEC(K,2)=TEMPA*TEMPC
      WVEC(K,3)=TEMPB*TEMPC
      WVEC(K,4)=QUART*(TEMPA*TEMPA-TEMPB*TEMPB)
  140 WVEC(K,5)=HALF*TEMPA*TEMPB
      DO 150 I=1,N
      IP=I+NPT
      WVEC(IP,1)=ZERO
      WVEC(IP,2)=D(I)
      WVEC(IP,3)=GW(I)
      WVEC(IP,4)=ZERO
  150 WVEC(IP,5)=ZERO
C
C     Put the coefficents of THETA*Wcheck in PROD.
C
      DO 220 JC=1,5
      KU=NPT
      IF (JC .EQ. 2 .OR. JC .EQ. 3) KU=NDIM
      DO 160 I=1,NPT
  160 PROD(I,JC)=ZERO
      DO 180 J=1,NPTM
      SUM=ZERO
      DO 170 I=1,NPT
  170 SUM=SUM+ZMAT(I,J)*WVEC(I,JC)
      IF (J .LT. IDZ) SUM=-SUM
      DO 180 I=1,NPT
  180 PROD(I,JC)=PROD(I,JC)+SUM*ZMAT(I,J)
      IF (KU .EQ. NDIM) THEN
          DO 200 I=1,NPT
          SUM=ZERO
          DO 190 J=1,N
  190     SUM=SUM+BMAT(I,J)*WVEC(NPT+J,JC)
  200     PROD(I,JC)=PROD(I,JC)+SUM
      END IF
      DO 220 I=1,N
      SUM=ZERO
      DO 210 K=1,KU
  210 SUM=SUM+BMAT(K,I)*WVEC(K,JC)
  220 PROD(NPT+I,JC)=SUM
C
C     Include in DEN the part of BETA that depends on THETA.
C
      DO 240 K=1,NDIM
      SUM=ZERO
      DO 230 I=1,5
      WW(I)=HALF*PROD(K,I)*WVEC(K,I)
  230 SUM=SUM+WW(I)
      DEN(1)=DEN(1)-WW(1)-SUM
      TEMPA=PROD(K,1)*WVEC(K,2)+PROD(K,2)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,4)+PROD(K,4)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,5)+PROD(K,5)*WVEC(K,3)
      DEN(2)=DEN(2)-TEMPA-HALF*(TEMPB+TEMPC)
      DEN(6)=DEN(6)-HALF*(TEMPB-TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,3)+PROD(K,3)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,5)+PROD(K,5)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,4)+PROD(K,4)*WVEC(K,3)
      DEN(3)=DEN(3)-TEMPA-HALF*(TEMPB-TEMPC)
      DEN(7)=DEN(7)-HALF*(TEMPB+TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,4)+PROD(K,4)*WVEC(K,1)
      DEN(4)=DEN(4)-TEMPA-WW(2)+WW(3)
      TEMPA=PROD(K,1)*WVEC(K,5)+PROD(K,5)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,3)+PROD(K,3)*WVEC(K,2)
      DEN(5)=DEN(5)-TEMPA-HALF*TEMPB
      DEN(8)=DEN(8)-WW(4)+WW(5)
      TEMPA=PROD(K,4)*WVEC(K,5)+PROD(K,5)*WVEC(K,4)
  240 DEN(9)=DEN(9)-HALF*TEMPA
C
C     Extend DEN so that it holds all the coefficients of DENOM.
C
      SUM=ZERO
      DO 250 I=1,5
      WW(I)=HALF*PROD(KNEW,I)**2
  250 SUM=SUM+WW(I)
      DENEX(1)=ALPHA*DEN(1)+WW(1)+SUM
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,2)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,4)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,5)
      DENEX(2)=ALPHA*DEN(2)+TEMPA+TEMPB+TEMPC
      DENEX(6)=ALPHA*DEN(6)+TEMPB-TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,3)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,5)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,4)
      DENEX(3)=ALPHA*DEN(3)+TEMPA+TEMPB-TEMPC
      DENEX(7)=ALPHA*DEN(7)+TEMPB+TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,4)
      DENEX(4)=ALPHA*DEN(4)+TEMPA+WW(2)-WW(3)
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,5)
      DENEX(5)=ALPHA*DEN(5)+TEMPA+PROD(KNEW,2)*PROD(KNEW,3)
      DENEX(8)=ALPHA*DEN(8)+WW(4)-WW(5)
      DENEX(9)=ALPHA*DEN(9)+PROD(KNEW,4)*PROD(KNEW,5)
C
C     Seek the value of the angle that maximizes the modulus of DENOM.
C
      SUM=DENEX(1)+DENEX(2)+DENEX(4)+DENEX(6)+DENEX(8)
      DENOLD=SUM
      DENMAX=SUM
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      WW(1)=ONE
      DO 280 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      WW(2)=DCOS(ANGLE)
      WW(3)=DSIN(ANGLE)
      DO 260 J=4,8,2
      WW(J)=WW(2)*WW(J-2)-WW(3)*WW(J-1)
  260 WW(J+1)=WW(2)*WW(J-1)+WW(3)*WW(J-2)
      SUMOLD=SUM
      SUM=ZERO
      DO 270 J=1,9
  270 SUM=SUM+DENEX(J)*WW(J)
      IF (DABS(SUM) .GT. DABS(DENMAX)) THEN
          DENMAX=SUM
          ISAVE=I
          TEMPA=SUMOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=SUM
      END IF
  280 CONTINUE
      IF (ISAVE .EQ. 0) TEMPA=SUM
      IF (ISAVE .EQ. IU) TEMPB=DENOLD
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-DENMAX
          TEMPB=TEMPB-DENMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new D and test for convergence.
C
      WW(2)=DCOS(ANGLE)
      WW(3)=DSIN(ANGLE)
      DO 290 I=1,N
  290 D(I)=WW(2)*D(I)+WW(3)*GW(I)
      DO 300 J=4,8,2
      WW(J)=WW(2)*WW(J-2)-WW(3)*WW(J-1)
  300 WW(J+1)=WW(2)*WW(J-1)+WW(3)*WW(J-2)
      BETA=ZERO
      DENMAX=ZERO
      TAU=ZERO
      DO 310 J=1,9
      BETA=BETA+DEN(J)*WW(J)
      DENMAX=DENMAX+DENEX(J)*WW(J)
  310 IF (J .LE. 5) TAU=TAU+PROD(KNEW,J)*WW(J)
      IF (ITERC .GE. N) GOTO 390
      IF (ITERC .EQ. 1) DENOLD=ZERO
      IF (DABS(DENMAX) .LE. 1.2D0*DABS(DENOLD)) GOTO 390
C
C     Set G to half the gradient of DENOM with respect to D. Then branch
C     for the next iteration.
C
      DO 320 I=1,N
  320 GW(I)=TAU*BMAT(KNEW,I)
      DO 370 K=1,NDIM
      PRVAL=ZERO
      DO 330 J=1,5
  330 PRVAL=PRVAL+PROD(K,J)*WW(J)
      IF (K .LE. NPT) THEN
          SUM=ZERO
          DO 340 I=1,N
  340     SUM=SUM+XPT(K,I)*(XOPT(I)+D(I))
          IF (K .EQ. KOPT) THEN
              TEMPA=ALPHA*(SUM-XOPTSQ+DSQ)
              TEMPB=TEMPA+ALPHA*SUM
              DO 350 I=1,N
  350         GW(I)=GW(I)+TEMPA*XOPT(I)+TEMPB*D(I)
          END IF
          TEMP=(TAU*VLAG(K)-ALPHA*PRVAL)*SUM
          DO 360 I=1,N
  360     GW(I)=GW(I)+TEMP*XPT(K,I)
      ELSE
          GW(K-NPT)=GW(K-NPT)-ALPHA*PRVAL
      END IF
  370 CONTINUE
      GG=ZERO
      DG=ZERO
      DO 380 I=1,N
      GG=GG+GW(I)**2
  380 DG=DG+D(I)*GW(I)
      TEMP=DG*DG/(DSQ*GG)
      IF (TEMP .LE. ONE-1.0D-8) GOTO 80
C
C     Set the vector VLAG before the RETURN from the subroutine.
C
  390 DO 410 K=1,NDIM
      VLAG(K)=ZERO
      SUM=ZERO
      DO 400 J=1,5
      VLAG(K)=VLAG(K)+PROD(K,J)*WW(J)
  400 SUM=SUM+WVEC(K,J)*WW(J)
  410 W(K)=SUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
      RETURN
      END SUBROUTINE

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% biglag.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,DELTA,D,ALPHA,GW,HCOL,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  GW(*),HCOL(*),W(*)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KOPT is the index of the optimal interpolation point.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DELTA is the current trust region bound.
C     D will be set to the step from XOPT to the new point.
C     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
C     GW, HCOL and W will be used for working space.
C
C     The step D is calculated in a way that attempts to maximize the modulus
C     of TAU=VLAG(KNEW), subject to the bound ||D|| .LE. DELTA.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      DSQ=DELTA*DELTA
      NPTM=NPT-N-1
C
C     Set the first D and GW, where GW is the gradient of the KNEW-th Lagrange
C     function at the trust region centre. Moreover, W(NPT+.) will hold the
C     gradient of this function at XPT(KNEW,.). The first NPT components of
C     HCOL and W will be set to the leading elements of the KNEW-th column of
C     H and the scalar products (XPT(K,.),XOPT), K=1,2,...,NPT, respectively.
C
      DO 10 I=1,N
      D(I)=XPT(KNEW,I)-XOPT(I)
      GW(I)=BMAT(KNEW,I)
   10 W(NPT+I)=BMAT(KNEW,I)
      DO 20 K=1,NPT
   20 HCOL(K)=ZERO
      DO 30 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 30 K=1,NPT
   30 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      DO 50 K=1,NPT
      W(K)=ZERO
      TEMP=ZERO
      DO 40 J=1,N
      W(K)=W(K)+XPT(K,J)*XOPT(J)
   40 TEMP=TEMP+XPT(K,J)*XPT(KNEW,J)
      TEMPA=HCOL(K)*W(K)
      TEMPB=HCOL(K)*TEMP
      DO 50 I=1,N
      GW(I)=GW(I)+TEMPA*XPT(K,I)
   50 W(NPT+I)=W(NPT+I)+TEMPB*XPT(K,I)
      ALPHA=HCOL(KNEW)
C
C     Revise GW if its modulus seems to be unusually small.
C
      TEMP=ZERO
      TEMPA=ZERO
      TEMPB=ZERO
      DO 60 I=1,N
      TEMP=TEMP+D(I)**2
      TEMPA=TEMPA+GW(I)**2
   60 TEMPB=TEMPB+W(NPT+I)**2
      IF (TEMPB*DSQ .GT. 1.0D4*TEMP*TEMPA) THEN
          DO 70 I=1,N
   70     GW(I)=W(NPT+I)
      END IF
C
C     Begin the iteration by making D and GW orthogonal and of length DELTA.
C
      ITERC=0
   80 ITERC=ITERC+1
      DD=ZERO
      DG=ZERO
      DO 90 I=1,N
      DD=DD+D(I)**2
   90 DG=DG+D(I)*GW(I)
      GG=ZERO
      DO 100 I=1,N
      GW(I)=DD*GW(I)-DG*D(I)
  100 GG=GG+GW(I)**2
      TEMPD=DELTA/DSQRT(DD)
      TEMPG=DELTA/DSQRT(GG)
      DO 110 I=1,N
      D(I)=TEMPD*D(I)
  110 GW(I)=TEMPG*GW(I)
C
C     Find the elements of W_check, and accumulate their contributions to
C     the coefficients of TAU.
C
      CF1=ZERO
      CF2=ZERO
      CF3=ZERO
      CF4=ZERO
      CF5=ZERO
      DO 130 K=1,NPT
      TEMPA=ZERO
      TEMPB=ZERO
      DO 120 I=1,N
      TEMPA=TEMPA+XPT(K,I)*D(I)
  120 TEMPB=TEMPB+XPT(K,I)*GW(I)
      TMPA=TEMPA*HCOL(K)
      TMPB=TEMPB*HCOL(K)
      CF1=CF1+HALF*TMPB*TEMPB
      CF2=CF2+TMPA*W(K)
      CF3=CF3+TMPB*W(K)
      CF4=CF4+HALF*(TMPA*TEMPA-TMPB*TEMPB)
  130 CF5=CF5+TMPA*TEMPB
      DO 140 I=1,N
      TEMP=BMAT(KNEW,I)
      CF2=CF2+TEMP*D(I)
  140 CF3=CF3+TEMP*GW(I)
C
C     Seek the value of the angle that maximizes the modulus of TAU.
C
      TAUBEG=CF1+CF2+CF4
      TAUMAX=TAUBEG
      TAUOLD=TAUBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 150 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      IF (DABS(TAU) .GT. DABS(TAUMAX)) THEN
          TAUMAX=TAU
          ISAVE=I
          TEMPA=TAUOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=TAU
      END IF
  150 TAUOLD=TAU
      IF (ISAVE .EQ. 0) TEMPA=TAU
      IF (ISAVE .EQ. IU) TEMPB=TAUBEG
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-TAUMAX
          TEMPB=TEMPB-TAUMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new D and test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      DO 160 I=1,N
  160 D(I)=CTH*D(I)+STH*GW(I)
      IF (ITERC .GE. N) GOTO 210
      IF (ITERC .EQ. 1) TAUBEG=ZERO
      IF (DABS(TAU) .LE. 1.1D0*DABS(TAUBEG)) GOTO 210
C
C     Set GW to the gradient of TAU with respect to D. Then branch for the
C     next iteration unless GW and D are nearly parallel.
C
      DO 170 I=1,N
  170 GW(I)=BMAT(KNEW,I)
      DO 190 K=1,NPT
      SUM=W(K)
      DO 180 I=1,N
  180 SUM=SUM+D(I)*XPT(K,I)
      TEMP=HCOL(K)*SUM
      DO 190 I=1,N
  190 GW(I)=GW(I)+TEMP*XPT(K,I)
      GG=ZERO
      DG=ZERO
      DO 200 I=1,N
      GG=GG+GW(I)**2
  200 DG=DG+D(I)*GW(I)
      TEMP=DG*DG/(DSQ*GG)
      IF (TEMP .LE. ONE-1.0D-8) GOTO 80
  210 RETURN
      END SUBROUTINE

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trsapp.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
     1  D,G,HD,HS,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*),
     1  D(*),G(*),HD(*),HS(*)
C
C     N is the number of variables of a quadratic objective function, Q say.
C     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
C       in order to define the current quadratic model Q.
C     DELTA is the trust region radius, and has to be positive.
C     STEP will be set to the calculated trial step.
C     The arrays D, G, HD and HS will be used for working space.
C     CRVMIN will be set to the least curvature of H along the conjugate
C       directions that occur, except that it is set to zero if STEP goes
C       all the way to the trust region boundary.
C
C     The calculation of STEP begins with the truncated conjugate gradient
C     method. If the boundary of the trust region is reached, then further
C     changes to STEP may be made, each one being in the 2D space spanned
C     by the current STEP and the corresponding gradient of Q. Thus STEP
C     should provide a substantial reduction to Q within the trust region.
C
C     Initialization, which includes setting HD to H times XOPT.
C
      HALF=0.5D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(1.0D0)
      DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      DO 10 I=1,N
   10 D(I)=XOPT(I)
      GOTO 170
C
C     Prepare for the first line search.
C
   20 QRED=ZERO
      DD=ZERO
      DO 30 I=1,N
      STEP(I)=ZERO
      HS(I)=ZERO
      G(I)=GQ(I)+HD(I)
      D(I)=-G(I)
   30 DD=DD+D(I)**2
      CRVMIN=ZERO
      IF (DD .EQ. ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
C
C     Calculate the step to the trust region boundary and the product HD.
C
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
      GOTO 170
   50 DHD=ZERO
      DO 60 J=1,N
   60 DHD=DHD+D(J)*HD(J)
C
C     Update CRVMIN and set the step-length ALPHA.
C
      ALPHA=BSTEP
      IF (DHD .GT. ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC .EQ. 1) CRVMIN=TEMP
          CRVMIN=DMIN1(CRVMIN,TEMP)
          ALPHA=DMIN1(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
C
C     Update STEP and HS.
C
      GGSAV=GG
      GG=ZERO
      DO 70 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      HS(I)=HS(I)+ALPHA*HD(I)
   70 GG=GG+(G(I)+HS(I))**2
      IF (ITERC .EQ. ITERMAX) GOTO 160
C
C     Begin another conjugate direction iteration if required.
C
      IF (ALPHA .LT. BSTEP) THEN
          IF (QADD .LE. 0.01D0*QRED) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO 80 I=1,N
          D(I)=TEMP*D(I)-G(I)-HS(I)
          DD=DD+D(I)**2
          DS=DS+D(I)*STEP(I)
   80     SS=SS+STEP(I)**2
          IF (SS .LT. DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
C
C     Test whether an alternative iteration is required.
C
   90 IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
      SG=ZERO
      SHS=ZERO
      DO 100 I=1,N
      SG=SG+STEP(I)*G(I)
  100 SHS=SHS+STEP(I)*HS(I)
      SGK=SG+SHS
      ANGTEST=SGK/DSQRT(GG*DELSQ)
      IF (ANGTEST .LE. -0.99D0) GOTO 160
C
C     Begin the alternative iteration by calculating D and HD and some
C     scalar products.
C
      ITERC=ITERC+1
      TEMP=DSQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO 110 I=1,N
  110 D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO 130 I=1,N
      DG=DG+D(I)*G(I)
      DHD=DHD+HD(I)*D(I)
  130 DHS=DHS+HD(I)*STEP(I)
C
C     Seek the value of the angle that minimizes Q.
C
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
      IF (QNEW .LT. QMIN) THEN
          QMIN=QNEW
          ISAVE=I
          TEMPA=QSAV
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=QNEW
      END IF
  140 QSAV=QNEW
      IF (ISAVE .EQ. ZERO) TEMPA=QNEW
      IF (ISAVE .EQ. IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+ANGLE)
C
C     Calculate the new STEP and HS. Then test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO 150 I=1,N
      STEP(I)=CTH*STEP(I)+STH*D(I)
      HS(I)=CTH*HS(I)+STH*HD(I)
  150 GG=GG+(G(I)+HS(I))**2
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01D0) GOTO 90
  160 RETURN
C
C     The following instructions act as a subroutine for setting the vector
C     HD to the vector D multiplied by the second derivative matrix of Q.
C     They are called from three different places, which are distinguished
C     by the value of ITERC.
C
  170 DO 180 I=1,N
  180 HD(I)=ZERO
      DO 200 K=1,NPT
      TEMP=ZERO
      DO 190 J=1,N
  190 TEMP=TEMP+XPT(K,J)*D(J)
      TEMP=TEMP*PQ(K)
      DO 200 I=1,N
  200 HD(I)=HD(I)+TEMP*XPT(K,I)
      IH=0
      DO 210 J=1,N
      DO 210 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
  210 HD(I)=HD(I)+HQ(IH)*D(J)
      IF (ITERC .EQ. 0) GOTO 20
      IF (ITERC .LE. ITERSW) GOTO 50
      GOTO 120
      END SUBROUTINE

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
C
C     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
C     interpolation point that has index KNEW. On entry, VLAG contains the
C     components of the vector Theta*Wcheck+e_b of the updating formula
C     (6.11), and BETA holds the value of the parameter that has this name.
C     The vector W is used for working space.
C
C     Set some constants.
C
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      DO 20 J=2,NPTM
      IF (J .EQ. IDZ) THEN
          JL=IDZ
      ELSE IF (ZMAT(KNEW,J) .NE. ZERO) THEN
          TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,JL)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 10 I=1,NPT
          TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL)
   10     ZMAT(I,JL)=TEMP
          ZMAT(KNEW,J)=ZERO
      END IF
   20 CONTINUE
C
C     Put the first NPT components of the KNEW-th column of HLAG into W,
C     and calculate the parameters of the updating formula.
C
      TEMPA=ZMAT(KNEW,1)
      IF (IDZ .GE. 2) TEMPA=-TEMPA
      IF (JL .GT. 1) TEMPB=ZMAT(KNEW,JL)
      DO 30 I=1,NPT
      W(I)=TEMPA*ZMAT(I,1)
      IF (JL .GT. 1) W(I)=W(I)+TEMPB*ZMAT(I,JL)
   30 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      TAUSQ=TAU*TAU
      DENOM=ALPHA*BETA+TAUSQ
      VLAG(KNEW)=VLAG(KNEW)-ONE
C
C     Complete the updating of ZMAT when there is only one nonzero element
C     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
C     then the first column of ZMAT will be exchanged with another one later.
C
      IFLAG=0
      IF (JL .EQ. 1) THEN
          TEMP=ALPHA/DENOM
          TEMPA=DSQRT(DABS(TEMP))
          TEMPB=TEMPA*TAU/ALPHA
          DO 40 I=1,NPT
   40     ZMAT(I,1)=TEMPA*VLAG(I)-TEMPB*W(I)
          IF (IDZ .EQ. 1 .AND. TEMP .LT. ZERO) IDZ=2
          IF (IDZ .GE. 2 .AND. TEMP .GE. ZERO) IFLAG=1
      ELSE
C
C     Complete the updating of ZMAT in the alternative case.
C
          JA=1
          IF (BETA .GE. ZERO) JA=JL
          JB=JL+1-JA
          TEMP=ZMAT(KNEW,JB)/DENOM
          TEMPA=TEMP*BETA
          TEMPB=TEMP*TAU
          TEMP=ZMAT(KNEW,JA)
          SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ)
          SCALB=SCALA*DSQRT(DABS(DENOM))
          DO 50 I=1,NPT
          ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
   50     ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
          IF (DENOM .LE. ZERO) THEN
              IF (BETA .LT. ZERO) IDZ=IDZ+1
              IF (BETA .GE. ZERO) IFLAG=1
          END IF
      END IF
C
C     IDZ is reduced in the following case, and usually the first column
C     of ZMAT is exchanged with a later one.
C
      IF (IFLAG .EQ. 1) THEN
          IDZ=IDZ-1
          DO 60 I=1,NPT
          TEMP=ZMAT(I,1)
          ZMAT(I,1)=ZMAT(I,IDZ)
   60     ZMAT(I,IDZ)=TEMP
      END IF
C
C     Finally, update the matrix BMAT.
C
      DO 70 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 70 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   70 CONTINUE
      RETURN
      END SUBROUTINE

      end module newuoa_mod