C.HR MARKOV.FOR
C--------------------------------------------------------------------
C MAIN PROGRAM - MARKOV.FOR
C
C REFFERENCE:
C
C    TAUCHEN, GEORGE, 1987, "QUADRATURE-BASED METHODS FOR OBTAINING
C    APPROXIMATE SOLUTIONS TO NONLINEAR ASSET PRICING MODELS,"
C    DUKE UNIVERSITY.
C
C THIS PROGRAM CALCULATES A DISCRETE APPROXIMATION TO A CONTINUOUS 
C GAUSSIAN VECTOR AUTREGRESSION
C
C Y(T) = B + A(1)*Y(T-1) + A(2)*Y(T-2) + ... + A(NLAG)*Y(T-NLAG) + E(T)
C
C WHERE
C
C    Y(T) IS AN (NVAR X 1) VECTOR
C    B IS AN (NVAR X 1) VECTOR
C    A(I) ARE (NVAR X NVAR) MATRICES
C    E(T) IS AN (NVAR X 1) VECTOR OF I.I.D. MULTIVARIATE NORMAL 
C       RANDOM VARIABLES WITH MEAN VECTOR ZERO AND VARIANCE SIGMAM 
C 
C THREE FILES USED FOR I/O ARE DEFINED WITH OPEN STATEMENTS 
C 
C     OPEN(3,FILE='CON') 
C     OPEN(10,FILE='PARMS') 
C     OPEN(11,FILE='GHQUAD.DAT') 
C 
C UNIT 11 CONTAINS DISCRETE VALUES AND WEIGHTS FOR N-POINT 
C    GAUSSIAN QUADRATURE RULES UP TO 20 POINTS.  
C
C UNIT 10 CONTAINS THE FOLLOWING INPUT PARAMETER VALUES 
C 
C    NVAR     - NUMBER OF VARIABLES IN THE VAR 
C               FORMAT (I5) 
C    NLAG     - NUMBER OF LAGS IN THE VAR 
C               FORMAT (I5) 
C    NVAL     - NUMBER OF DISCRETE POINTS FOR EACH VARIABLE 
C               NVAR ROWS OF FORMAT (I5) 
C    B        - VECTOR OF INTERCEPTS IN THE VAR 
C               NVAR ROWS OF FORMAT (F10.0) 
C    AUTO     - HORIZONTALLY CONCATENATED MATRICES A(I) 
C               [ A(1) | A(2) | ... | A(NLAG) ]
C               THIS MATRIX IS READ COLUMNWISE
C               NVAR*NVAR*NLAG ROWS OF FORMAT (F10.0)
C    SIGMAM    - VARIANCE MATRIX OF E(T) READ COLUMNWISE
C               NVAR*NVAR ROWS OF FORMAT (F10.0)
C
C UNIT 3 CONTAINS ALL OF THE OUTPUT.  THE VARIABLES NS AND 
C    NSTATE ARE DEFINED AS FOLLOWS
C
C    NS       - NUMBER OF DISCRETE VALUES THAT Y(T) CAN TAKE ON;          
C               EQUAL TO THE PRODUCT OF THE ELEMENTS OF THE MATRIX 
C               NVAL.  
C 
C    NSTATE   - NUMBER OF STATES IN THE SYSTEM; EQUAL TO NS**NLAG.
C               THE STATE IS DEFINED BY THE DISCRETE VALUES OF
C               Y(T-1), Y(T-2), ...,Y(T-NLAG); SINCE EACH LAG CAN
C               TAKE ON NS DIFFERENT VALUES, THERE ARE NS**NLAG
C               STATES IN THE SYSTEM.
C
C  EXAMPLE:  IF THERE ARE TWO VARIABLES (Y1 AND Y2) AND TWO
C            LAGS IN THE VAR AND TWO DISCRETE POINTS ARE USED FOR 
C            EACH VARIABLE, THE STATES ARE ARRANGED AS FOLLOWS:
C
C                                   DISCRETE VALUES
C                      
C            STATE NO.        TIME T-1            TIME T-2
C            ---------     ---------------     ---------------
C                           
C                1         Y1(1)     Y2(1)     Y1(1)     Y2(1)
C                2         Y1(2)     Y2(1)     Y1(1)     Y2(1)
C                3         Y1(1)     Y2(2)     Y1(1)     Y2(1)
C                4         Y1(2)     Y2(2)     Y1(1)     Y2(1)
C                5         Y1(1)     Y2(1)     Y1(2)     Y2(1)
C                6         Y1(2)     Y2(1)     Y1(2)     Y2(1)
C                7         Y1(1)     Y2(2)     Y1(2)     Y2(1)
C                8         Y1(2)     Y2(2)     Y1(2)     Y2(1)
C                9         Y1(1)     Y2(1)     Y1(1)     Y2(2)
C               10         Y1(2)     Y2(1)     Y1(1)     Y2(2)
C               11         Y1(1)     Y2(2)     Y1(1)     Y2(2)
C               12         Y1(2)     Y2(2)     Y1(1)     Y2(2)
C               13         Y1(1)     Y2(1)     Y1(2)     Y2(2)
C               14         Y1(2)     Y2(1)     Y1(2)     Y2(2)
C               15         Y1(1)     Y2(2)     Y1(2)     Y2(2)
C               16         Y1(2)     Y2(2)     Y1(2)     Y2(2)
C
C THREE OUTPUT MATRICES ARE RELEVANT
C
C    YMAT     - (NS X NVAR) MATRIX OF DISCRETE VARIABLE VALUES
C
C    PISTA    - (NSTATE X 1) VECTOR THAT IS THE STATIONARY
C               PROBABILITY DISTRIBUTION OF DISCRETE STATES
C
C    PIMAT    - (NSTATE X NS) MATRIX OF TRANSITION PROBABILITIES.  
C               IF NLAG>1, THIS MATRIX IS NOT SQUARE BECAUSE FROM
C               ANY GIVEN STATE, ONLY NS OTHER STATES ARE REACHABLE
C               IN THE NEXT PERIOD.
C
C ALSO AVAILABLE IN THE SUBROUTINE MCHAIN IS THE (NSTATE X NS)
C MATRIX IREACH WHICH CONTAINS THE REACHABLE STATE NUMBERS.  THUS
C ROW I OF IREACH CONTAINS THE NUMBERS OF THE NS STATES THAT ARE 
C REACHABLE IN THE NEXT PERIOD WHEN THE STATE NUMBER IN THE CURRENT
C PERIOD IS I.
C--------------------------------------------------------------------
	subroutine markov
	use params


      IMPLICIT REAL(prec) (A-H,O-Z)
      INTEGER ZZNVAR,ZZNLAG,ZZNS,ZZTHET,ZZPMAT,ZZIW,ZZIW2,ZZNST,
     &   ZZIW3,ZZW
      PARAMETER(NRQUAD=210)
      PARAMETER(ZZNVAR=2)
      PARAMETER(ZZNLAG=1)
      PARAMETER(ZZNS=50)
      PARAMETER(ZZTHET=ZZNVAR+(ZZNVAR**2)*(ZZNLAG+2)+1)
      PARAMETER(ZZPMAT=(ZZNS**ZZNLAG)*ZZNS)
      PARAMETER(ZZIW=2*ZZNS+2*ZZNLAG)
      PARAMETER(ZZNST=ZZNS**ZZNLAG)
      PARAMETER(ZZW=4*ZZNS+2*ZZNVAR+4*ZZNS*ZZNVAR+2*ZZNVAR*ZZNLAG)
      PARAMETER(ZZIW2=ZZNST*ZZNLAG)
      PARAMETER(ZZIW3=ZZNST*ZZNS)
      INTEGER*4 IWORK(ZZIW),IWORK1(ZZNST),IWORK2(ZZIW2)
      INTEGER*4 IWORK3(ZZIW3)
      REAL(prec) PIMAT(ZZPMAT),GHMAT(NRQUAD,3),WORK(ZZW)
      REAL(prec) WORK1(ZZNST),WORK2(ZZNST)




C
C OPEN FILES
C
      OPEN(3,FILE='mres.txt')
      OPEN(99,FILE='GHQUAD.DAT')
C
C POINTERS TO THETA
C
      IB     = 1
      IAUTO  = IB               + NVAR
      ISIG   = IAUTO            + (NVAR**2)*NLAG
      ISIGNV = ISIG+NVAR**2
      IDETSG = ISIGNV+NVAR**2

      CALL HEADER('/B//_')
      CALL DGMPNT(THETA(IB),NVAR,1)
      CALL HEADER('/AUTO//_')
	
      CALL DGMPNT(THETA(IAUTO),NVAR,NVAR*NLAG)
      CALL HEADER('/SIG//_')
      CALL DGMPNT(THETA(ISIG),NVAR,NVAR)

	
C
C CALCULATE THE NUMBER OF STATES
C this is taken from the module parameters
C
C ERROR TRAPS
C
      IF (NVAR .GT. ZZNVAR) THEN
         WRITE(3,1081) NVAR
         STOP
      ENDIF
      IF (NLAG .GT. ZZNLAG) THEN
         WRITE(3,1082) NLAG
         STOP
      ENDIF
      IF (nsm .GT. ZZNS) THEN
         WRITE(3,1083) nsm
         STOP
      ENDIF
      IF ( (4*(NVAR**2)+4*NVAR) .GT. ZZW) THEN
         IFIX=(4*(NVAR**2)+4*NVAR)-ZZW
         WRITE(3,1084) IFIX
         STOP
      ENDIF
      IF (ZZW .GT. 8000) THEN
         WRITE(3,1085)
         STOP
      ENDIF
C
C READ QUADRATURE DATA POINTS
C	
   
      DO 20 I=1,NRQUAD
   20 READ(99,1101) (GHMAT(I,J), J=1,3)

C
C COMPUTE THE INVERSE OF SIGMAM
C
      DO 30 I=1,NVAR**2
   30 THETA(ISIGNV+I-1)=THETA(ISIG+I-1)

      CALL DSWEEP(THETA(ISIGNV),NVAR,1.D-13,IER)
      CALL HEADER('/SIGINV//_')
      CALL DGMPNT(THETA(ISIGNV),NVAR,NVAR)
      IF (IER .GT. 0) THEN
         WRITE(3,1090) IER
         STOP
      ENDIF
C
C POINTERS TO WORK FOR CALL TO MCHAIN
C
      MR1=1
      MR2=MR1+nsm*NVAR
      MR3=MR2+nsm
      MR4=MR3+NVAR
      MR5=MR4+nsm*NVAR
C
C POINTERS TO IWORK FOR CALL TO MCHAIN
C
      MI1=1
      MI2=MI1+NLAG
      MI3=MI2+nsm
      CALL MCHAIN(THETA(IAUTO),THETA(IB),THETA(ISIG),NLAG,NVAL,
     &  nsm,GHMAT,NRQUAD,NVAR,WORK(MR1),WORK(MR2),IWORK3,
     &  IWORK(MI1),WORK(MR3),IWORK(MI2),WORK(MR4),IWORK1,
     &  PIMAT,WORK1,THETA(ISIGNV),THETA(IDETSG),WORK2,WORK(MR5),
     &  IWORK2,IWORK(MI3))


      CLOSE(99)
      CLOSE(3)

     
 1001 FORMAT(I5)
 1002 FORMAT(F10.0)
 1101 FORMAT(F10.0,2F20.0)
 1081 FORMAT(1X,'ERROR:  INCREASE ZZNVAR TO ',I5)
 1082 FORMAT(1X,'ERROR:  INCREASE ZZNLAG TO ',I5)
 1083 FORMAT(1X,'ERROR:  INCREASE ZZNS TO ',I5)
 1084 FORMAT(1X,'ERROR:  INCREASE ZZW BY ',I5)
 1085 FORMAT(1X,'ERROR:  DECREASE ZZNS OR COMPILE UNDER HUGE MODULE')
 1090 FORMAT(1X,'ERROR:  SIGMAM IS SINGULAR:  IER = ',I5)
      END
C.HR SUBROUTINE MCHAIN 
C@
C--------------------------------------------------------------------
C SUBROUTINE MCHAIN                                                  
C--------------------------------------------------------------------
      SUBROUTINE MCHAIN(AUTO,B,SIGMAM,NLAGd,NVALd,NNS,GHMAT,NRGH,NVARd,
     &   ZMAT,WVEC,IREACH,IS,EXVAL,IR,YMAT,ISEQA,PIMAT,PISTA,SIGINV,
     &   DETSIG,WORK1,WORK,IWORK1,IWORK)
	use params

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IS(NLAG),IR(nsm),IREACH(nsm**NLAG,nsm)
      INTEGER ISEQA(nsm**NLAG),IWORK1(1),IWORK(1)
      REAL(prec) AUTO(NVAR,NVAR*NLAG),B(NVAR),SIGMAM(NVAR,NVAR)
      REAL(prec) GHMAT(NRGH,3),ZMAT(nsm,NVAR),WVEC(nsm),EXVAL(NVAR)
      REAL(prec) YMAT(nsm,NVAR),PIMAT(nsm**NLAG,nsm),pista(nsm*nlag)
      REAL(prec) SIGINV(NVAR,NVAR),WORK1(1),WORK(1)
C--------------------------------------------------------------------
C SUBROUTINE TO CALCULATE DISCRETE MARKOV CHAIN THAT APPROXIMATES VAR
C--------------------------------------------------------------------
C SET QUADRATURE ARRAYS                                               
C--------------------------------------------------------------------
      CALL QUAD(NVAL,NVAR,nsm,ZMAT,WVEC,GHMAT,NRGH)
	


      CALL HEADER('/ZMAT//_')
      CALL DGMPNT(ZMAT,nsm,NVAR)
      CALL HEADER('/WVEC//_')
      CALL DGMPNT(WVEC,nsm,1)
C--------------------------------------------------------------------
C SET ADDITIONAL CONSTANTS AND ARRAYS OF PARAMETERS                  
C--------------------------------------------------------------------
      NSTATE = nsm**NLAG
      CALL SPACER(2)
      WRITE(3,100) NVAR,NLAG,nsm,NSTATE
C--------------------------------------------------------------------
C FORM THE REACHABLE STATES AND PUT THEM IN IREACH                    
C--------------------------------------------------------------------
      CALL SPACER(1)
      WRITE(3,110)
      DO 20 J=1,NSTATE
         CALL DECODE(J,IS,1,nsm,NLAG)
         CALL REACHR(IS,IR,IWORK1,nsm,NLAG)
         DO 10 K=1,nsm
   10    IREACH(J,K)= IR(K)
         ISEQA(J)=J
   20 CONTINUE
      CALL SPACER(1)
      WRITE(3,120)
C--------------------------------------------------------------------
C CALCULATE THE MEAN OF THE PROCESS                                  
C--------------------------------------------------------------------
      IR1=1
      IR2=IR1+NVAR**2
      IR3=IR2+NVAR**2
      CALL ARMEAN(AUTO,B,EXVAL,WORK(IR1),WORK(IR2),WORK(IR3),NVAR,NLAG)
      CALL HEADER('/EXVAL//_')
      CALL DGMPNT(EXVAL,1,NVAR)          
C--------------------------------------------------------------------
C TRANSFORM TO THE Y'S AND ADD IN THE MEANS                          
C--------------------------------------------------------------------
      CALL HEADER('/SIGMA//_')
      CALL DGMPNT(SIGMAM,NVAR,NVAR)
      CALL ZTOY(ZMAT,YMAT,WORK(IR1),EXVAL,SIGMAM,NVAR,DETSIG,nsm,
     &          WORK(IR2))
      CALL HEADER('/YMAT//_')
      CALL DGMPNT(YMAT,nsm,NVAR)         
      do 11 j=1,NVAR
        do 11 i=1,nsm
          write(3,*) YMAT(i,j)
 11   continue
C--------------------------------------------------------------------
C GET TRANSITION PROBABILITIES                                       
C--------------------------------------------------------------------
      CALL SPACER(2)
      WRITE(3,130)
      IR2=IR1+NVAR*NLAG
      IR3=IR2+nsm
      CALL PISUB(ISEQA,NSTATE,YMAT,PIMAT,WORK(IR1),IWORK1,
     &           WORK(IR2),nsm,NVAR,NLAG,WVEC,EXVAL,SIGINV,DETSIG,AUTO,
     &           B,WORK(IR3))




C--------------------------------------------------------------------
C Storing ymat and pimat                                      
C--------------------------------------------------------------------

	mstates=ymat(1:nsm,1)
	mprobs=pimat





      CALL SPACER(1)
      WRITE(3,140)
      CALL HEADER('/PIMAT//_')
      CALL DGMPNT(PIMAT,NSTATE,nsm)         
      do 21 j=1,nsm
        do 21 i=1,NSTATE
          write(3,*) PIMAT(i,j)
 21   continue
C--------------------------------------------------------------------
C COMPUTE THE CHAIN'S STATIONARY DISTRIBUTION                        
C--------------------------------------------------------------------
      CALL SPACER(1)
      WRITE(3,150)
      CALL SPACER(1)
      II1=1
      II2=II1+NLAG
      CALL MARKST(PIMAT,PISTA,WORK1,IWORK(II1),IWORK1,IWORK(II2),
     &            NSTATE,nsm,NLAG)

	mstae=pista ! storing the stationary distriibution
      CALL SPACER(1)
      WRITE(3,160)
      CALL HEADER('/PISTA//_')
      CALL DGMPNT(PISTA,NSTATE,1)         
      RETURN
  100 FORMAT(1X,'NVAR =',I5,'   NLAG =',I5,'   nsm =',I5,
     &          '   NSTATE =',I5)      
  110 FORMAT(1X,'FORMING THE REACHABLE STATES')
  120 FORMAT(1X,'REACHABLE STATES FORMED')
  130 FORMAT(1X,'STARTING TO COMPUTE TRANSITION PROBABILITIES') 
  140 FORMAT(1X,'TRANSITION PROBABILITIES COMPUTED')
  150 FORMAT(1X,'STARTING TO COMPUTE STATIONARY DISTRIBUTION')
  160 FORMAT(1X,'STATIONARY DISTRIBUTION COMPUTED')
      END  
C.HR SUBROUTINE PHIMAT 
C.@
      SUBROUTINE PHIMAT(XMAT,MU,SIGINV,FZMAT,XCTR,XSIG,DETSIG,N,NVAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMAT(N,NVAR),MU(NVAR),SIGINV(NVAR,NVAR)
      REAL*8 XCTR(N,NVAR),XSIG(N,NVAR),FZMAT(N)
      DATA PI /3.141592653589798/
C--------------------------------------------------------------------
C SUBROUTINE FOR ROW-WISE MULTIVARIATE NORMAL DISTRIBUTION                  
C
C    NOTES: THE ELEMENTS OF RETURNED VECTOR ARE THE MULTIVARIATE     
C           NORMAL DENSITY N(MU,SIGMAM) EVALUATED AT EACH ROW OF XMAT 
C--------------------------------------------------------------------
      DO 10 I=1,N
      DO 10 J=1,NVAR
   10 XCTR(I,J)=XMAT(I,J)-MU(J)
      CALL DGMPRD(XCTR,SIGINV,XSIG,N,NVAR,NVAR)
      DO 30 I=1,N
         FZMAT(I)=0
         DO 20 J=1,NVAR
   20    FZMAT(I)=FZMAT(I)+XSIG(I,J)*XCTR(I,J)
         FZMAT(I)=EXP(-0.5*FZMAT(I))/((2*PI*DETSIG)**.5)
   30 CONTINUE
      RETURN
      END
C.HR SUBROUTINE CODE 
C@
C--------------------------------------------------------------------
C SUBROUTINE CODE                                                    
C--------------------------------------------------------------------
      SUBROUTINE CODE(IMAT,ICODE,NR,NC,NBASE)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IMAT(NR,NC),ICODE(NR)
C--------------------------------------------------------------------
C SUBROUTINE TO CODE EACH ROW OF IMAT AS                             
C                                                                    
C        X = N1 + (N2-1)*NBASE + ... + (NL-1)*NBASE(L-1)               
C                                                                    
C WHERE (N1,N2,...,NL) IS A ROW OF IMAT AND L = COLS(IMAT).          
C--------------------------------------------------------------------
C CHECK FOR ERRORS                                                   
C--------------------------------------------------------------------
      IF (NBASE .LT. 1) THEN
        WRITE(3,100) NBASE
        STOP
      ENDIF
      DO 20 I=1,NR
       DO 10 J=1,NC
        IF (IMAT(I,J) .GT. NBASE .OR. IMAT(I,J) .LE. 0) THEN
         WRITE(3,110) 
         STOP
        ENDIF
   10  CONTINUE
   20 CONTINUE
C--------------------------------------------------------------------
C START MAIN PART OF SUBROUTINE                                      
C--------------------------------------------------------------------
      DO 40 I=1,NR
         ICODE(I)=IMAT(I,1)
         DO 30 J=2,NC
   30    ICODE(I)=ICODE(I)+(IMAT(I,J)-1)*(NBASE**(J-1))
   40 CONTINUE
      RETURN
  100 FORMAT(1X,'*** CODE ERROR: BAD NBASE',I5)
  110 FORMAT(1X,'*** CODE ERROR: AN ELEMENT OF THE INPUT ARRAY EITHER',
     &       1X,'  EXCEEDS THE NBASE OR IS LESS THAN OR EQUAL TO ZERO')
      END
C.HR SUBROUTINE DECODE 
C@
      SUBROUTINE DECODE(IX,ISTATE,NR,NBASE,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IX(NR),ISTATE(NR,L)
C--------------------------------------------------------------------
C SUBROUTINE TO DECODE THE ROWS OF IX WHERE A TYPICAL ROW IS OF THE 
C FORM
C                                                                    
C       IX(K,.) =  N1 + (N2-1)*NBASE**1 + ... (NL-1)*NBASE**(L-1)        
C--------------------------------------------------------------------
C CHECK FOR ERRORS                                                   
C--------------------------------------------------------------------
      IF (NBASE .LE. 0 .OR. L .LT. 1) THEN
         WRITE(3,100) NBASE,L
         STOP
      ENDIF
      DO 10 I=1,NR
         IF (IX(I) .LE. 0) THEN
            WRITE(3,110)
         ENDIF  
   10 CONTINUE
C--------------------------------------------------------------------
C START SUBROUTINE                                                   
C--------------------------------------------------------------------
      DO 20 I=1,NR
      DO 20 J=1,L
   20 ISTATE(I,J)=INT(MOD(AINT(REAL(IX(I)-1)/REAL(NBASE**(J-1))),
     &                          REAL(NBASE)))+1
      RETURN
  100 FORMAT(1X,'*** DECODE ERROR BAD NBASE = ',I5,
     &          '  OR BAD LENGTH = ',I5)
  110 FORMAT(1X,'*** DECODE WARNING: AT LEAST ONE ELEMENT OF INPUT',
     &       1X,'    VECTOR IS LESS THAN OR EQUAL TO ZERO')
      END
C.HR SUBROUTINE REACHR 
C@
C--------------------------------------------------------------------
C              SUBROUTINE REACHR                                     
C--------------------------------------------------------------------
      SUBROUTINE REACHR(IVEC,IY,IX,NBASE,NLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IVEC(NLAG),IX(NBASE,NLAG),IY(NBASE)
      DO 20 I=1,NBASE
         IX(I,1)=I
         DO 10 J=2,NLAG
   10    IX(I,J)=IVEC(J-1)
   20 CONTINUE
      CALL CODE(IX,IY,NBASE,NLAG,NBASE)
      RETURN
      END
C.HR SUBROUTINE PISUB
C.@
C--------------------------------------------------------------------
C SUBROUTINE FOR SELECTED ROWS OF PROBABILITY TRANSITION MATRIX      
C--------------------------------------------------------------------
      SUBROUTINE PISUB(JVEC,NRJ,YMAT,PIOUT,X,IST,PIOUTT,nsm,NVAR,NLAG,
     &           WVEC,EXVAL,SIGINV,DETSIG,AUTO,B,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 JVEC(NRJ),IST(NRJ,NLAG)
      REAL*8 YMAT(nsm,NVAR),PIOUT(NRJ,nsm),X(NVAR*NLAG)
      REAL*8 PIOUTT(1,nsm),WVEC(NVAR**NLAG),EXVAL(NVAR)
      REAL*8 SIGINV(NVAR,NVAR),AUTO(NVAR,NVAR*NLAG),B(NVAR),WORK(1)
      I1=1
      I2=I1+nsm
      I3=I2+nsm
      I4=I3+NVAR*NLAG
      I5=I4+NVAR
      CALL DECODE(JVEC,IST,NRJ,nsm,NLAG)
      DO 40 I=1,NRJ
         DO 10 L=1,NLAG
         DO 10 K=1,NVAR
   10    X((L-1)*NVAR+K)=YMAT(IST(I,L),K)
         CALL PIFN(X,1,PIOUTT,YMAT,nsm,NVAR,WVEC,EXVAL,SIGINV,
     &             DETSIG,AUTO,B,WORK(I1),WORK(I2),WORK(I3),
     &             WORK(I4),NLAG,WORK(I5))
         VSUM=0
         DO 20 J=1,nsm
   20    VSUM=VSUM+PIOUTT(1,J)
         DO 30 J=1,nsm
   30    PIOUT(I,J)=PIOUTT(1,J)/VSUM
   40 CONTINUE 
      RETURN
      END 



	  




C.HR SUBROUTINE ARMEAN 
C@
      SUBROUTINE ARMEAN(AUTO,B,MEAN,SUMAUT,SUMINV,WORKSP,NVAR,NLAG)
      REAL*8 AUTO(NVAR,NLAG*NVAR),B(NVAR),MEAN(NVAR)
      REAL*8 SUMAUT(NVAR,NVAR),SUMINV(NVAR,NVAR)
      REAL*8 WORKSP(2*NVAR**2+4*NVAR)
C--------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE MEAN OF AN AUTOREGRESSIVE PROCESS      
C--------------------------------------------------------------------
C START MAIN PART OF SUBROUTINE                                      
C--------------------------------------------------------------------
      DO 30 L=1,NLAG
         DO 20 I=1,NVAR
            DO 10 J=1,NVAR
               IF (L .EQ. 1) SUMAUT(I,J)=0
               SUMAUT(I,J)= SUMAUT(I,J)+AUTO(I,(J+(L-1)*NVAR))
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      DO 50 I=1,NVAR
         DO 40 J=1,NVAR
            SUMAUT(I,J)=-SUMAUT(I,J)
            IF (I .EQ. J) SUMAUT(I,J)=1+SUMAUT(I,J)
   40    CONTINUE
   50 CONTINUE
      CALL DMPINV(SUMAUT,NVAR,NVAR,SUMINV,IRANK,WORKSP)
      CALL DGMPRD(SUMINV,B,MEAN,NVAR,NVAR,1)
      RETURN
      END
C.HR SUBROUTINE QUAD 
C@
C--------------------------------------------------------------------
C SUBROUTINE QUAD                                                    
C--------------------------------------------------------------------
      SUBROUTINE QUAD(NVAL,NVAR,nsm,ZMAT,WVEC,QUADAT,NRQUAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 NVAL(NVAR)
      REAL*8 QUADAT(NRQUAD,3),ZMAT(nsm,NVAR),WVEC(nsm)
C--------------------------------------------------------------------
C SUBROUTINE TO PUT ALL POSSIBLE COMBINATIONS OF QUADRATURE ABSCISSA 
C    AND CORRESPONDING PRODUCTS OF QUADRATURE WEIGHTS IN             
C    RETURNED MATRIX                                                 
C--------------------------------------------------------------------
C START MAIN SUBROUTINE                                              
C--------------------------------------------------------------------

	  
	  
	  NPCUM=1
      DO 10 I=1,nsm
   10 WVEC(I)=1
      DO 50 I=1,NVAR
       NP=NVAL(I)
       DO 40 J=1,NP
        Z=QUADAT(((NP-1)*NP)/2+J,2)
        W=QUADAT(((NP-1)*NP)/2+J,3)
        DO 30 K=1,(nsm/NP)/NPCUM
         DO 20 L=1,NPCUM
          ZMAT((K-1)*NPCUM*NP+(J-1)*NPCUM+L,I)=Z
          WVEC((K-1)*NPCUM*NP+(J-1)*NPCUM+L)=WVEC((K-1)*NPCUM*NP+
     &       (J-1)*NPCUM+L)*W
   20    CONTINUE
   30   CONTINUE
   40  CONTINUE
       NPCUM=NPCUM*NP
   50 CONTINUE
      RETURN
      END 


C.HR SUBROUTINE CHOL
C@
      SUBROUTINE CHOL(S,M,XX,CSTAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 S(M,M),XX((M*(M+1))/2),CSTAT(M,M)
C
C     SUBROUTINE TO FIND THE CHOLESKY FACTORIZATION OF A MATRIX S
C
C     SUBROUTINE CALLED:  DMFSD
C
      DO 10 J=1,M
      DO 10 I=1,J
   10 XX(I+(J*(J-1))/2)=S(I,J)
      CALL DMFSD(XX,M,1.D-13,IER)
      IF (IER.LT.0) THEN
      WRITE(NOUT,100) IER
      STOP
      ENDIF
      DO 30 J=1,M
      DO 20 I=1,M
      CSTAT(I,J)=0
      IF (I .LE. J) CSTAT(I,J)=XX(I+(J*(J-1))/2)
   20 CONTINUE
   30 CONTINUE
      RETURN
  100 FORMAT(1X,'ERROR, CALL TO DMFSD, IER =',I5)
      END
C.HR SUBROUTINE ZTOY 
C@
      SUBROUTINE ZTOY(ZMAT,YMAT,S,EXVAL,SIGMAM,NVAR,DETSIG,nsm,WORKSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ZMAT(nsm,NVAR),YMAT(nsm,NVAR),S(NVAR,NVAR),EXVAL(NVAR)
      REAL*8 SIGMAM(NVAR,NVAR),WORKSP(1)
C--------------------------------------------------------------------
C SUBROUTINE TO TRANSFORM THE Z'S TO THE Y'S                         
C--------------------------------------------------------------------
C START MAIN PART OF SUBROUTINE                                      
C--------------------------------------------------------------------
      CALL CHOL(SIGMAM,NVAR,WORKSP,S)
      CALL DGMPRD(ZMAT,S,YMAT,nsm,NVAR,NVAR)
      DO 20 I=1,nsm
         DO 10 J=1,NVAR
            YMAT(I,J)=YMAT(I,J)+EXVAL(J)
   10    CONTINUE
   20 CONTINUE
      DETSIG=1
      DO 30 J=1,NVAR
         DETSIG=DETSIG*S(J,J)
   30 CONTINUE
      RETURN    
      END
C.HR SUBROUTINE PIFN 
C.@
C--------------------------------------------------------------------
C SUBROUTINE FOR TRANSITION PROBABILITES AS A FUNCTION OF X                
C--------------------------------------------------------------------
      SUBROUTINE PIFN(XMAT,NRXMAT,PIVAL,YMAT,nsm,NVAR,WVEC,EXVAL,
     &           SIGINV,DETSIG,AUTO,B,TEMP,FZMAT,X,CONEXP,
     &           NLAG,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMAT(NRXMAT,NVAR*NLAG),PIVAL(NRXMAT,nsm),YMAT(nsm,NVAR)
      REAL*8 WVEC(nsm),EXVAL(NVAR),SIGINV(NVAR,NVAR)
      REAL*8 AUTO(NVAR,NVAR*NLAG),B(NVAR),TEMP(nsm),FZMAT(nsm)
      REAL*8 X(NVAR*NLAG),CONEXP(NVAR),WORK(1)
      I1=1
      I2=I1+nsm*NVAR
      CALL PHIMAT(YMAT,EXVAL,SIGINV,FZMAT,WORK(I1),WORK(I2),DETSIG,
     &            nsm,NVAR)
      DO 10 I=1,nsm
   10 TEMP(I)= WVEC(I)/FZMAT(I)
      DO 60 J=1,NRXMAT
         DO 20 K=1,NVAR*NLAG
   20    X(K)= XMAT(J,K)
         DO 40 L=1,NVAR
            CONEXP(L)=0
            DO 30 K=1,NVAR*NLAG
   30       CONEXP(L)=CONEXP(L) + X(K)*AUTO(L,K)
            CONEXP(L)=CONEXP(L)+B(L)
   40    CONTINUE
         CALL PHIMAT(YMAT,CONEXP,SIGINV,FZMAT,WORK(I1),WORK(I2),
     &               DETSIG,nsm,NVAR)
         DO 50 I=1,nsm
   50    PIVAL(J,I)=FZMAT(I)*TEMP(I)
   60 CONTINUE
      RETURN
      END
C.HR SUBROUTINE MARKST 
C@
      SUBROUTINE MARKST(PIMAT,PISTA,TEMP,IST,ISTOLD,JOLD,NSTATE,
     &                  nsm,NLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IST(NLAG),ISTOLD(nsm,NLAG),JOLD(nsm)
      REAL*8 PIMAT(NSTATE,nsm),PISTA(NSTATE),TEMP(NSTATE)
C--------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE STATIONARY DISTRIBUTION OF THE MARKOV 
C CHAIN  
C--------------------------------------------------------------------
C INTITIALIZE                                                        
C--------------------------------------------------------------------
      DATA EPSTA /.0001/
      DATA MAXIT /50/
      CRIT=EPSTA+1     
      NBASE = nsm
      DO 10 J=1,NSTATE
   10 PISTA(J)=1/REAL(NSTATE)
C--------------------------------------------------------------------
C MAIN LOOP TO COMPUTE THE CHAIN'S STATIONARY DISTRIBUTION           
C--------------------------------------------------------------------
      DO 90 NUMIT=1,MAXIT
         IF (CRIT .GT. EPSTA) THEN
            IF (NLAG .EQ. 1) THEN
               DO 30 J=1,NSTATE
                  TEMP(J)=0
                  DO 20 K=1,NSTATE
   20             TEMP(J)=TEMP(J) + PISTA(K)*PIMAT(K,J)
   30          CONTINUE
            ELSE
               DO 70 J=1,NSTATE
                  CALL DECODE(J,IST,1,NBASE,NLAG)
                  DO 50 K=1,nsm
                     DO 40 L=1,NLAG-1
   40                ISTOLD(K,L)=IST(L+1)
                     ISTOLD(K,NLAG)=K
   50             CONTINUE
                  CALL CODE(ISTOLD,JOLD,nsm,NLAG,NBASE)
                  SUM=0
                  DO 60 K=1,nsm
   60             SUM=SUM+PIMAT(JOLD(K),IST(1))*PISTA(JOLD(K))
                  TEMP(J)=SUM
   70          CONTINUE
            ENDIF
            CRIT=0
            DO 80 J=1,NSTATE
               DELTA=ABS(TEMP(J)-PISTA(J))
               IF (DELTA .GT. CRIT) CRIT=DELTA
               PISTA(J)=TEMP(J)
   80       CONTINUE
            WRITE(3,100) NUMIT,CRIT
         ENDIF
         IF ((NUMIT .EQ. MAXIT) .AND. (CRIT .GT. EPSTA)) THEN
            WRITE(3,110) NUMIT,DELTA
            CALL DGMPNT(PISTA,NSTATE,1)
         ENDIF
   90 CONTINUE
      RETURN
  100 FORMAT(1X,'NUMIT = ',I5,'  CRIT = ',F7.6)
  110 FORMAT(1X,'*** UNSUCCESSFUL ITERATIONS FOR STATIONARY ',
     &          'DISTRIBUTION',1X,'NUMIT = ',I5,'  DELTA = ',F10.6)
      END 


