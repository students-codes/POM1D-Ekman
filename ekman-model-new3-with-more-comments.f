      PROGRAM MAIN
C**********************************************************************
C
C      THIS IS A 1-D VERSION OF THE MELLOR PRINCETON OCEAN MODEL
C      WITH CONSTANT TURBULENT EDDY VISCOSITY
C
C***********************************************************************
      INCLUDE 'ekman-model.h'
C
C  working arrays
C
      DIMENSION ZZnew(110000),uanalytical(110000),vanalytical(110000)
      DATA PI/3.14159/
C
      DATA ISTART/1/,IPRINT/96/ 
C
C --- NAMELIST CONTAINS THE USER INPUT PARAMETERS
C
      NAMELIST /LIST1/LEVELS,DT1,days,H,sal,temp,WINDX,WINDY,
     1      VERTMIX,ALAT
      READ(*,LIST1)
      WRITE(6,*) 'Parametri di input'
      WRITE(6,*) 'vertical LEVELS      : ',LEVELS   
      WRITE(6,*) '  DT1   (seconds)    : ',DT1
      WRITE(6,*) ' integration days    : ',days
      WRITE(6,*) '  H      (meters)    : ',H
      WRITE(6,*) '  sal     (PSU)      : ',sal
      WRITE(6,*) '  temp    (deg C)    : ',temp
      WRITE(6,*) '  WINDX  (m/s)       : ',WINDX
      WRITE(6,*) '  WINDY  (m/s)       : ',WINDY
      WRITE(6,*) '  VERTMIX(m**2/s)    : ',VERTMIX
      WRITE(6,*) 'latitude for Coriolis: ',ALAT
C
C --- TRANSFORM IN MODEL VARIABLES THE INPUT VALUES
C
      KB=LEVELS
      KBM1=KB-1
      rhoair=1.25
      CD=1.0e-03
      rho0=1025.
C
C -- POM requires that the surface forcing is given with the opposite sign and divided by rho0
c
      WUSURF=-rhoair*CD*SQRT(WINDX**2+WINDY**2)*WINDX/rho0
      WVSURF=-rhoair*CD*SQRT(WINDX**2+WINDY**2)*WINDY/rho0
C
C --- true value to be output
C
      wwx=-wusurf*rho0
      wwy=-wvsurf*rho0
      WRITE(6,504) Wwx,Wwy
 504  FORMAT(/,'X component WIND STRESS ', 1PE10.2,
     1     ' Y component WIND STRESS ', 1PE10.2,/)

C*********************************************************************
C        
C        SUBROUTINE DEPTH ESTABLISHES THE VERTICAL NUMERICAL GRID
C         Z,ZZ, DZ AND DZZ ARE DIMENSIONLESS
C
C *********************************************************************
C
      CALL DEPTH(Z,ZZ,DZ,DZZ,KB,H)
C
C --- MODEL CONSTANT INITIALIZATION
C --- COR   CORIOLIS PARAMETER
C --- GRAV  GRAVITY
C --- SMOTH  NUMERICAL PARAMETER TO MIX THE TIME MARCHING STEPS 
C --- MOLECULAR VISCOSITY
C
      TIME=0.                   
      COR=2.*7.29e-5*sin(alat*.01745329)
      write(6,*)' coriolis parameter ', cor
      GRAV=9.806
      SMOTH=1.
      UMOL=0.0
C
      ISTART=1
      IEND=days*86400./DT1
      IPRINT=IEND/5
C

      Write(6,*)
      WRITE(6,15) ISTART,IEND,DT1
15    FORMAT(//' MODEL START TIME STEP ISTART =',I4,
     1  ' MODEL END TIME STEP IEND =',I6,
     1  '    TIME STEP IN SECONDS =',F5.1,/)
C
C
C****INSERT USER CHOICE OF INITIAL CONDITIONS*********************
C
      DO 5 K=1,KB
      UB(K)=0.
      U(K)=UB(K)
      UF(K)=U(K)
      VB(K)=0.
      V(K)=VB(K)
      VF(K)=V(K)
      TB(k)=temp
      T(K)=TB(K)
      TF(K)=T(K)
      SB(K)=sal
      S(K)=SB(K)
      SF(K)=S(K)
C
C --- VERTICALLY CONSTANT VALUE OF TURBULENT VISCOSITY
C
      KM(K)=VERTMIX
   5  CONTINUE
C
      DT2=2.*DT1
      DT4=4.*DT1
      DAYI=DT1/86400.
C
C
C
C***********************************************************************
C*                                                                     *
C*      BEGIN  time integration                                        *
C*                                                                     *
C***********************************************************************
C

      DO 9000 IINT=ISTART,IEND
C
      TIME=TIME+DT1/86400.
C
      DO 380 K=1,KB-1
      UF(K)=UB(K)+DT2*COR*V(K)
      VF(K)=VB(K)-DT2*COR*U(K)
 380  CONTINUE
      CALL PROFU(DT2,KB)
      CALL PROFV(DT2,KB)
C
C--- time filter 
C
      DO 382 K=1,KB
      U(K)=U(K)+.5*SMOTH*(UF(K)+UB(K)-2.*U(K))
      V(K)=V(K)+.5*SMOTH*(VF(K)+VB(K)-2.*V(K))
      UB(K)=U(K)
      U(K)=UF(K)
      VB(K)=V(K)
 382  V(K)=VF(K)
C
 9000 CONTINUE
C
C
C***********************************************************************
C*                                                                     *
C*      END   time integration                                        *
C*                                                                     *
C***********************************************************************
C
C
c    now start to print out
c
      WRITE(6,'('' END OF SIMULATION, TIME (DAYS) ='',F10.2)') TIME
      write(6,*)'  simulation results      '
      WRITE(6,501)
 501  FORMAT(1X,'           DEPTH     U         V         KM')
      DO 550 K=1,KBm1
      ZZD=ZZ(K)*H
      WRITE(6,502) IEND,ZZD,UB(K),VB(K),KM(K)
      write(97,*) ub(K)
      write(98,*) vb(K)
      write(99,*) ZZD
 550  CONTINUE
 502  FORMAT(1X,'solution at',I10,3X,'time step',3X,F7.1,3(1PE10.2))
C
C ********analytical solution case for meridional constant wind stress
C
C
      write(6,*)' analytical solution values   '
      km0=vertmix
      tau0=-wvsurf*rho0
      DE=PI*sqrt(2.*km0/cor)
      V0=(1.4142*PI*TAU0)/(DE*cor*RHO0)
      write(6,*)' parameters of the analytical solution'
      write(6,*)' meridional wind stress amplitude ', tau0
      write(6,*)' vertical mixing coefficient ',km0
      write(6,*)' density constant value ',rho0
      write(6,*)' e-folding EKMAN depth  ',de/pi
      write(6,*)' velocity amplitude    ', V0
      write(6,*)' depth      u         v      UDIFF    VDIFF'
c
c
      do 9100 k=1,kbm1
      theta=(pi/4)+((pi/DE)*(zz(k)*h))
      uanalytical(k)=V0*cos(theta)*exp(pi*zz(k)*h/de)
      vanalytical(k)=V0*sin(theta)*exp(pi*zz(k)*h/de)
      UDIFF=uanalytical(k)-UB(K)
      VDIFF=Vanalytical(k)-VB(K)
      WRITE(6,503)ZZ(K)*H,Uanalytical(K),Vanalytical(K),UDIFF,VDIFF
 9100 continue
c
c *************** write analytical solution
c

      do 9200 k=1,kbm1
      zzd=zz(k)*h
      write(87,*)uanalytical(k)
      write(88,*)vanalytical(k)
      write(89,*)zzd
 9200 continue
c
 503  FORMAT(1X,F7.1,4(1PE10.2))
c
      STOP
      END
      
C
      SUBROUTINE DEPTH(Z,ZZ,DZ,DZZ,KB,H)
      DIMENSION Z(KB),ZZ(KB),DZ(KB),DZZ(KB)
      DIMENSION ZZnew(KB)
c
C***********************************************************************
C   THIS SUBROUTINE ESTABLISHES THE VERTICAL NUMERICAL GRID 
c   depths and layers
C***********************************************************************
C
      DO 10 K=1,KB
   10 Z(K)=-FLOAT(K-1)/FLOAT(KB-1)
C
      DO 11 K=1,KB-1
      ZZ(K)=.5*(Z(K)+Z(K+1))
      DZ(K)=Z(K)-Z(K+1)
      ZZnew(K)=ZZ(K)*H
   11 CONTINUE
C
      ZZ(KB)=ZZ(KB-1)-dz(KB-1)
      ZZnew(KB)=ZZ(KB)*H
C
      DO 6 K=1,KB-1
      DZZ(K)=ZZ(K)-ZZ(K+1)
6     CONTINUE
C
C
      WRITE(6,*) 'Model levels in meters'
      WRITE(6,*) '    K       Z' 
      WRITE(6,'(I5,F14.7 )') (K,ZZnew(K),K=1,KB)
      RETURN
      END

C
      SUBROUTINE PROFU(DT2,KB)
        INCLUDE 'ekman-model.h'            
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C         DT2*(KM*U')' - U= -UB                                        *
C                                                                      *
C***********************************************************************
 85   DH=H
      DO 100 K=2,KB-1
      A(K-1)=-DT2*(KM(K)+UMOL)/(DZ(K-1)*DZZ(K-1)*DH
     1     *DH)
      C(K)=-DT2*(KM(K)+UMOL)/(DZ(K)*DZZ(K-1)*DH
     1     *DH)
 100  CONTINUE
      EE(1)=A(1)/(A(1)-1.)
      GG(1)=(-DT2*WUSURF/(-DZ(1)*DH)-UF(1))
     1   /(A(1)-1.)
      DO 101 K=2,KB-2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K)
      GG(K)=(C(K)*GG(K-1)-UF(K))*GG(K)
 101  CONTINUE
      CBC=AMAX1(.0025,.16/ALOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
      CBC=CBC*SQRT(UB(KB-1)**2+(.25*(VB(KB-1)
     1     +VB(KB-1)+VB(KB-1)+VB(KB-1)))**2)
C********************************************************************
C        TO RESTORE BOTTOM B.L. DELETE NEXT LINE
C*********************************************************************
      CBC=0.
      UF(KB-1)=(C(KB-1)*GG(KB-2)-UF(KB-1))/(CBC
     1     *DT2/(-DZ(KB-1)*DH)-1.-(EE(KB-2)-1.)*C(KB-1))
      DO 103 K=2,KB-1
      KI=KB-K
      UF(KI)=EE(KI)*UF(KI+1)+GG(KI)
 103  CONTINUE
 92   WUBOT=-CBC*UF(KB-1)
      RETURN
      END
C
      SUBROUTINE PROFV(DT2,KB)
        INCLUDE 'ekman-model.h'           
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C         DT2*(KM*V')' -V= -VB                                         *
C                                                                      *
C***********************************************************************
      DH=H
      DO 100 K=2,KB-1
      A(K-1)=-DT2*(KM(K)+UMOL  )/(DZ(K-1)*DZZ(K-1)*DH
     1     *DH)
      C(K)=-DT2*(KM(K)+UMOL  )/(DZ(K)*DZZ(K-1)*DH
     1     *DH)
 100  CONTINUE
      EE(1)=A(1)/(A(1)-1.)
      GG(1)=(-DT2*WVSURF/(-DZ(1)*DH)-VF(1))
     1   /(A(1)-1.)
 98   CONTINUE
      DO 101 K=2,KB-2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K)
      GG(K)=(C(K)*GG(K-1)-VF(K))*GG(K)
 101  CONTINUE
      CBC=AMAX1(.0025,.16/ALOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
      CBC=CBC*SQRT(UB(KB-1)**2+(.25*(VB(KB-1)
     1     +VB(KB-1)+VB(KB-1)+VB(KB-1)))**2)
C********************************************************************
C        TO RESTORE BOTTOM B.L. DELETE NEXT LINE
C*********************************************************************
      CBC=0.
      VF(KB-1)=(C(KB-1)*GG(KB-2)-VF(KB-1))/(CBC
     1     *DT2/(-DZ(KB-1)*DH)-1.-(EE(KB-2)-1.)*C(KB-1))
      DO 103 K=2,KB-1
      KI=KB-K
      VF(KI)=EE(KI)*VF(KI+1)+GG(KI)
 103  CONTINUE
 92   WVBOT=-CBC*VF(KB-1)
      RETURN
      END
C********************************************************************
C                   END OF FILE                       
C*********************************************************************
	    
