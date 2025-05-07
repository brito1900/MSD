        PROGRAM SPHERE
	implicit integer*8 (I-N)

        COMMON / BLOCK1 / RX, RY, VX, VY
        COMMON / BLOCK2 / COLTIM, PARTNR
        COMMON / MEZCLA / sigmas,mass,type
	COMMON / NUMBER / NN
	Common / parms / ss1,ss2,mm1,mm2,al11,al12,al22,del11,del12,del22
	REAL*8 ss1,ss2,mm1,mm2,al11,al12,al22,del11,del12,del22
        

C    *******************************************************************
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                    NUMBER OF ATOMS                  **
C    ** REAL    RX(N),RY(N)          ATOM POSITIONS                   **
C    ** REAL    VX(N),VY(N)          ATOM VELOCITIES                  **
C    ** REAL    COLTIM(N)            TIME TO NEXT COLLISION           **
C    ** INTEGER PARTNR(N)            COLLISION PARTNER                **
C    ** REAL    SIGMA                ATOM DIAMETER                    **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE )                                  **
C    **    READS IN CONFIGURATION                                     **
C    ** SUBROUTINE CHECK ( SIGMA, OVRLAP, E )                         **
C    **    CHECKS CONFIGURATION AND CALCULATES ENERGY                 **
C    ** SUBROUTINE UPLIST ( SIGMA, I )                                **
C    **    SEEKS COLLISIONS WITH J>I                                  **
C    ** SUBROUTINE DNLIST ( SIGMA, I )                                **
C    **    SEEKS COLLISIONS WITH J<I                                  **
C    ** SUBROUTINE BUMP ( SIGMA, I, J, e, lx, b, alpha, delta )       **
C    **    DOES COLLISION DYNAMICS AND CALCULATES COLLISION VIRIAL    **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT CONFIGURATION                                   **
C    **                                                               **
C    ** OUTPUT FILES                                                  **
C    ** fort.20  time, collisions, energy                             **
C    ** fort.30  configurations, x(i),y(i),vx(i),vy(i)                **
C    ** fort.40  cumulants of the velocity distr. function (empty)    **
C    ** fort.50  histrograms (empty)                                  **
C    ** fort.21  t delta n(k,t) (complex number) for k=2*pi/L         **
C    ** fort.22  the same as fort.21 for k=2pi/L * 2                  **
C    ** fort.3X  t,delta v_perp(k,t) (complex number) for k=2*pi/L*X  **
C    ** fort.4X  t,delta v_par (k,t) for k=2*pi/L*X                   **
C    **                                                               **
C    **                                                               **
C    ** Delta parameters have to be corrected with Lx,                **
C    ** as they have velocity dimensions                              **
C    *******************************************************************

        INTEGER*8     N,NN,irate,ik
        PARAMETER ( N = 50 000 )
	
        REAL*8      TIMBIG, EPS, ALPHA, DELTA
        PARAMETER ( TIMBIG = 1.0E10 )
	

        REAL*8      RX(N), RY(N), VX(N), VY(N)
        REAL*8      sigmas(N),mass(N)
	integer*8   type(N)
        REAL*8      COLTIM(N)
        INTEGER*8     PARTNR(N)
        REAL*8      SIGMA
	
	real*8      b
	integer*8   indexb, nbin,ninel
	real*8      histo(0:100)
	real*8      RXbis, RYbis
	real*8      tau,EFINAL
	real*16     Etracer,Colltracer

        INTEGER*8     I, J, K, COLL, irs
        REAL*8      TIJ, T, LY,LX, E, kt, elost
        CHARACTER   CNFILE*30
        LOGICAL     OVRLAP


        INTEGER*8     inumcol(8)
        INTEGER*8       iturn 

        INTEGER*8     maxik,xturn,yturn
        PARAMETER ( maxik=8 )
	
	COMPLEX*16  deltan,deltav,deltal
        real*8      kx,ky,PI,nk(maxik)
        
	COMPLEX, PARAMETER:: ICOMP=(0,1)
	PI=4.0d0*datan(1.0d0)
	xturnn=0
	yturnn=0
	xturn=0
	yturn=0
	Etracer=0.0d0
	colltracer=0.0d0


	do iiii=1,maxik
	nk(iiii)=3*iiii-2
	enddo

	

        iturn=1
        inumcol(1)=250000
        inumcol(2)=500000
	inumcol(3)=750000
	inumcol(4)=1000000
	inumcol(5)=ninel/5
	inumcol(6)=ninel/6
	inumcol(7)=ninel/7
	inumcol(8)=ninel/8

c    elost is the percentage of energy lost between sucessive 
c          writings of configurations. 

	irate=1 

	call readwrite(LX,LY,ALPHA,DELTA,NELAS,NINEL,KT,IRS)

	IF (NN.GT.N) THEN 
	     print*,'TOO MANY PARTICLES, ABORTING.....'
	     stop
	     endif

C    ** READ IN CONFIGURATION **

15       CALL INIT(kt, irs, lx) 
	
	do iii=1,nn
	sigmas(iii)= sigmas(iii) / LX
	enddo

C    ** CHECK FOR PARTICLE OVERLAPS **
C    ** CALCULATE ENERGY            **

        CALL CHECK ( SIGMA, OVRLAP )

         IF ( OVRLAP ) goto 15

C    ** SET UP INITIAL COLLISION LISTS COLTIM AND PARTNR **

        DO 10 I = 1, NN
           COLTIM(I) = TIMBIG
           PARTNR(I) = NN
10      CONTINUE

        DO 20 I = 1, NN
           CALL UPLIST ( SIGMA, I )
20      CONTINUE

c
c This Tmp includes an Lx factor as it was renormalized, by dimensionality condition to 1/lx
	Temp= pi*al22**2/(4*(1-al22**2)**2)*(1+sqrt(1+4*(1-al22**2)/
     @   (pi*al22**2)))**2*(del22*Lx)**2*mm2

        WRITE(*,'(//'' **** START OF DYNAMICS **** '')')
	write(*,'(''TEMP   '',1p3e12.4)') Temp

C    *******************************************************************
C    ** ELASTIC LOOP BEGINS                                           **
C    *******************************************************************
	nbin = 40
	indexb = 0
        DO 1000 COLL = 1, NELAS
C       ** LOCATE MINIMUM COLLISION TIME **

           TIJ = TIMBIG

           DO 200 K = 1, NN
              IF ( COLTIM(K) .LT. TIJ ) THEN
                 TIJ = COLTIM(K)
                 I   = K
              ENDIF
200        CONTINUE

           J = PARTNR(I)

c	if(TIJ.le.0.0d0) print*, " ERROR, T NEGAtIVO!!",tij, coll,i,j 
c	print*, "    "
c	print*, " T, coll, i,j",tij, coll,i,j 

C       ** MOVE PARTICLES FORWARD BY TIME TIJ AND REDUCE COLLISION TIMES **
C       ** APPLY PERIODIC BOUNDARIES          **

           T = T + TIJ

           DO 300 K = 1, NN
              COLTIM(K) = COLTIM(K) - TIJ
              RX(K) = RX(K) + VX(K) * TIJ
              RY(K) = RY(K) + VY(K) * TIJ
              RX(K) = RX(K) - ANINT ( RX(K) )
              RY(K) = RY(K) - ANINT ( RY(K) )
300        CONTINUE


C       ** COMPUTE COLLISION DYNAMICS **

           CALL BUMP ( SIGMA, I, J, E, lx, b, alpha, delta) 
	tau = 2.*coll/NN


C       ** RESET COLLISION LISTS FOR THOSE PARTICLES WHICH NEED IT **

           DO 400 K = 1, NN
              IF ( ( K .EQ. I ) .OR. ( PARTNR(K) .EQ. I ) .OR.
     :             ( K .EQ. J ) .OR. ( PARTNR(K) .EQ. J )     ) THEN
                 CALL UPLIST ( SIGMA, K )
              ENDIF
400        CONTINUE

           CALL DNLIST ( SIGMA, I )
           CALL DNLIST ( SIGMA, J )

1000    CONTINUE
C    *******************************************************************
C    ** ELASTIC LOOP ENDS.                                            **
C    *******************************************************************
c	CALL maxwell(lx)
        T = 0.0d0
	E=0.0d0
	do 1030 i=1,nn
	E=E+vx(i)*vx(i)+vy(i)*vy(i)
 1030	continue
	E=E/2.0d0

 2020   format(1pe20.10,i15,1pe20.10)
	write(30,2024) T,0,E*lx*lx
       	do ip =1,nn
	WRITE(30,198)RX(ip)*LX+lx/2,RY(ip)*LX+lx/2,VX(ip)*lx,VY(ip)*lx
	enddo    

C    *******************************************************************
C    ** INELASTIC LOOP BEGINS                                         **
C    *******************************************************************
	print*, 'INELASTIC LOOP BEGINS       **'
	EFINAL=0
	do ii=1,nn
c	print*, (0.5+rx(ii))*lx,(0.5+ry(ii))*lx,sigmas(ii)*lx/2
	write(88,*) vx(ii)*lx,vy(ii)*lx
	EFINAL=EFINAL+vx(ii)**2+vy(ii)**2
 	enddo
c	print*, "E inicial:",EFINAL*lx**2/2,"  TInicial:", EFINAL/2/nn*Lx**2
	print*, "E inicial:",EFINAL*lx**2/2,"  TInicial:", EFINAL/2/nn*LX**2


        DO 2000 ICOLL = 1, NINEL-1, IRATE
c	print*,icoll-1,t, e*lx*lx

6666 	format ( 1pe16.7,'  (',1pe16.7,','1pe16.7,')')
	
C	write(20,2020) T,icoll-1,E*lx*lx

	do 2222 COLL = 1, IRATE
C       ** LOCATE MINIMUM COLLISION TIME **
           TIJ = TIMBIG

           DO 2200 K = 1, NN
              IF ( COLTIM(K) .LT. TIJ ) THEN
                 TIJ = COLTIM(K)
                 I   = K
              ENDIF
 2200      CONTINUE

           J = PARTNR(I)
c	print*, "    "
c	print*, " T, coll, i,j",tij, coll,i,j 

C       ** MOVE PARTICLES FORWARD BY TIME TIJ AND REDUCE COLLISION TIMES **
C       ** APPLY PERIODIC BOUNDARIES          **

           T = T + TIJ

           DO 2300 K = 1, NN
              COLTIM(K) = COLTIM(K) - TIJ
              RX(K) = RX(K) + VX(K) * TIJ
              RY(K) = RY(K) + VY(K) * TIJ
		if(k.eq.1) then
			if(rx(1).gt.0.5d0) xturn=xturn+1
			if(rx(1).lt.-0.5d0) xturn=xturn-1
			if(ry(1).gt.0.5d0) yturn=yturn+1
			if(ry(1).lt.-0.5d0) yturn=yturn-1
		endif
		if(k.eq.NN) then
                        if(rx(NN).gt.0.5d0) xturnn=xturnn+1
                        if(rx(NN).lt.-0.5d0) xturnn=xturnn-1
                        if(ry(NN).gt.0.5d0) yturnn=yturnn+1
                        if(ry(NN).lt.-0.5d0) yturnn=yturnn-1
                endif
              RX(K) = RX(K) - ANINT ( RX(K) )
              RY(K) = RY(K) - ANINT ( RY(K) )
 2300      CONTINUE

C       ** COMPUTE COLLISION DYNAMICS **

           CALL BUMP ( SIGMA, I, J, E, lx, b, alpha, delta)

	tau = 2.*(icoll+coll-1)/NN

C       ** RESET COLLISION LISTS FOR THOSE PARTICLES WHICH NEED IT **

           DO 2400 K = 1, NN
              IF ( ( K .EQ. I ) .OR. ( PARTNR(K) .EQ. I ) .OR.
     :             ( K .EQ. J ) .OR. ( PARTNR(K) .EQ. J )     ) THEN
                 CALL UPLIST ( SIGMA, K )
              ENDIF
 2400      CONTINUE

           CALL DNLIST ( SIGMA, I )
           CALL DNLIST ( SIGMA, J )

2222    CONTINUE

C ** CONDITION TO WRITE THE CONFIGURATION. 

        if ((icoll+irate-1).eq.inumcol(iturn)) then
		iturn=iturn+1
c	        call maxwell(lx)
		write(30,2024) T,icoll+irate-1,E*lx*lx
		write(*,2024) T,icoll+irate-1,(icoll+irate-1.0)/nn,E*lx*lx
        	do ip =1,nn
		WRITE(30,198)RX(ip)*LX+lx/2,RY(ip)*LX+lx/2,VX(ip)*lx,VY(ip)*lx
		enddo    
	endif 

 198 	format(1p4e13.5)
 199	format(1p4e18.10)
 2024   format('# time=',1pe15.8,'   Coll= ',i8 ,
     @       ' Coll/NN=',1pe12.5,' Energy=',1pe18.10)

	write(40,*) T,lx*(rx(1)+xturn),lx*(ry(1)+yturn),icoll
	write(50,*) T,lx*(rx(nn)+xturnn),lx*(ry(nn)+yturnn),icoll
		Etracer=Etracer + vx(1)**2+vy(1)**2
		colltracer=colltracer+1

2000    CONTINUE

 	
	EFINAL=0
	do ii=1,nn
c	print*, (0.5+rx(ii))*lx,(0.5+ry(ii))*lx,sigmas(ii)*lx/2
	write(88,*) vx(ii),vy(ii)
	EFINAL=EFINAL+vx(ii)**2+vy(ii)**2
 	enddo
	print*, "E final:",EFINAL*Lx**2/2,"  TFinal:", EFINAL/2/nn*Lx**2,
     @   " TfinalTracer=", ETracer/2*Lx**2/Colltracer

C    *******************************************************************
C    ** INELASTIC LOOP ENDS.                                          **
C    *******************************************************************

	write(40,333) E*lx*lx/NN,NINEL, ETracer/2*Lx**2/Colltracer
	write(50,333) E*lx*lx/NN,NINEL, ETracer/2*Lx**2/Colltracer
 333    format('# Energy, colls, Ttracer=',1pe14.6,i20,1pe14.6)
        WRITE(*,'(//'' **** END OF DYNAMICS **** '')')

C    ** CHECK FOR PARTICLE OVERLAPS **

c        CALL CHECK ( SIGMA, OVRLAP )

c        IF ( OVRLAP ) THEN
c           WRITE(*,'('' PARTICLE OVERLAP IN FINAL CONFIGURATION '')')
c        ENDIF

        STOP
        END

	SUBROUTINE READWRITE(LX,LY,ALPHA,DELTA,NELAS,NINEL,KT,IRS)
	implicit integer*8 (I-N)
	COMMON /NUMBER / NN
        COMMON / MEZCLA / sigmas,mass,type
	Common/parms/ ss1,ss2,mm1,mm2,al11,al12,al22,del11,del12,del22
	REAL*8 ss1,ss2,mm1,mm2,al11,al12,al22,del11,del12,del22
        PARAMETER ( N = 50 000 )

	integer*8 NN,NELAS,NINEL,IRS
	REAL*8 LX,LY,ALPHA,DELTA,KT
        REAL*8      sigmas(N),mass(N)
	integer*8   type(N)
	integer   IRSS    


C    *******************************************************************

        WRITE(*,'(//  '' MOLECULAR DYNAMICS OF HARD SPHERES     '')')

C    ** READ IN BASIC SIMULATION PARAMETERS **

        WRITE(*,'('' ENTER NUMBER OF PARTICLES of each type     '')')
        READ (*,*) NN1,nn2
        WRITE(*,'('' ENTER SIGMAS                               '')')
        READ (*,*) SS1,SS2
        WRITE(*,'('' ENTER MASSES                               '')')
        READ (*,*) MM1,MM2
        WRITE(*,'('' ENTER ALPHA11,ALPHA12,ALPHA22              '')')
        READ (*,*) al11,al12,al22
        WRITE(*,'('' ENTER DELTA11,DELTA12,DELTA22              '')')
        READ (*,*) del11,del12,del22
        WRITE(*,'('' ENTER BOX LENGTH                           '')')
        READ (*,*) LX
        WRITE(*,'('' ENTER NUMBER OF ELASTIC COLLISIONS         '')')
        READ (*,*) NELAS
        WRITE(*,'('' ENTER NUMBER OF INELASTIC COLLISIONS       '')')
        READ (*,*) NINEL
        WRITE(*,'('' ENTER TEMPERATURE kT                       '')')
        READ (*,*) kt   
        WRITE(*,'('' ENTER SEED OF RANDOM NUMBER GENERATOR      '')')
        READ (*,*) IRS
	IRSS=IRS
C  Srand iniciaiza el gen de aleatorios
	call srand(IRSS)
	xxx=rand(IRSS)

	do i=20,90,10
        WRITE(i,'(''#NUMBER OF PARTICLES    '',2I7)'    ) NN1,nn2
        WRITE(i,'(''#SIGMAS                 '',2f15.5)'    ) SS1,ss2
        WRITE(i,'(''#MASSES                 '',2f15.5)'    ) mm1,mm2
        WRITE(i,'(''#ALPHAS                 '',3f15.5)'    ) al11,al12,al22
        WRITE(i,'(''#DELTAS                 '',3f15.5)'    ) del11,del12,del22
        WRITE(i,'(''#BOX LENGHT             '',F15.5)'  ) LX  
        WRITE(i,'(''#ELASTIC COLLISIONS       '',I15)'  ) NELAS
        WRITE(i,'(''#INELASTIC COLLISIONS       '',I15)'  ) NINEL 
        WRITE(i,'(''#TEMPERATURE kT         '',f15.5)'  ) kt    
        WRITE(i,'(''#SEED OF RANDOM NUMBER  '',I15)'  ) IRS   
	enddo

c	do i=20,49 
c        WRITE(i ,'(''#MOLECULAR DYNAMICS-DISSIPATIVE HARD SPHERES'')')
c        WRITE(i ,'(''#'')')
c        WRITE(i ,'(''#NUMBER OF PARTICLES....'',I15)'    ) NN
c        WRITE(i ,'(''#BOX LENGHT.............'',F15.5)'  ) LX   
c        WRITE(i ,'(''#ELASTIC COLLISIONS.....'',I15)'  ) NELAS
c        WRITE(i ,'(''#INELASTIC COLLISIONS.....'',I15)'  ) NINEL 
c        WRITE(i ,'(''#ALPHA AND DELTA PARAMETERS   '',2f13.5)') 
c     @              ALPHA, DELTA
c        WRITE(i ,'(''#TEMPERATURE kT.........'',f15.5)'  ) kt    
c        WRITE(i ,'(''#SEED OF RANDOM NUMBER..'',I15 )'  ) IRS   
c        WRITE(i ,'(''#'')')
c	enddo
C HEADERS FOR INDIVIDUAL FILES
	write(20,2022)
 2022   format('#     Time                 Collision        Energy    ')
	do i=1,9
	write(20+i,'(''#    time     delta n (k) with k= 2pi/L *'',i2)')i
	write(30+i,'(''#    time     delta vT(k) with k= 2pi/L *'',i2)')i
	write(40+i,'(''#    time     delta vL(k) with k= 2pi/L *'',i2)')i
	enddo

c	define matrices
	nn=nn1+nn2
	do i=1,nn1
	sigmas(i)=ss1
	mass(i)=mm1
	type(i)=1
	enddo
	do i=1,nn2
	sigmas(nn1+i)=ss2
	mass(nn1+i)=mm2
	type(nn1+i)=2
	enddo
	del11=del11/Lx
	del12=del12/Lx
	del22=del22/Lx
	
	
	RETURN 
	END

        SUBROUTINE CHECK ( SIGMA, OVRLAP)
	implicit integer*8 (I-N)

        COMMON / BLOCK1 / RX, RY, VX, VY
        COMMON / BLOCK2 / COLTIM, PARTNR
        COMMON / MEZCLA / sigmas,mass,type
	COMMON / NUMBER / NN

C    *******************************************************************
C    ** TESTS FOR PAIR OVERLAPS                                       **
C    *******************************************************************

        INTEGER*8     N,NN
        PARAMETER ( N = 50 000 )

        REAL*8      RX(N), RY(N), VX(N), VY(N)
        REAL*8      COLTIM(N)
        REAL*8      sigmas(N),mass(N)
	integer*8   type(N)
        INTEGER*8     PARTNR(N)

        REAL*8      SIGMA
        LOGICAL     OVRLAP

        INTEGER*8     I, J
        REAL*8      RXI, RYI, RXIJ, RYIJ, RIJSQ, SIGSQ, RIJ
        REAL*8      TOL
        PARAMETER ( TOL = 1.0E-8 )

C    *******************************************************************

        OVRLAP = .FALSE.

        DO 100 I = 1, NN - 1

           RXI = RX(I)
           RYI = RY(I)

           DO 99 J = I + 1, NN
c		print*,i,j
        SIGSQ  = ((sigmas(i)+sigmas(j))/2)** 2

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RXIJ = RXIJ - ANINT ( RXIJ )
              RYIJ = RYIJ - ANINT ( RYIJ )
              RIJSQ = RXIJ ** 2 + RYIJ ** 2

              IF ( RIJSQ .LT. SIGSQ ) THEN

                 RIJ = SQRT ( RIJSQ / SIGSQ )
                 WRITE(*,'('' I,J,RIJ/SIGMA = '',2I5,F15.8)')
     :              I, J, RIJ
                 WRITE(*,'(''POSITIONS = '',4F15.8)')
     :              rx(i),rx(j),ry(i),ry(j) 

                 IF ( ( 1.0 - RIJ ) .GT. TOL ) OVRLAP = .TRUE.

              ENDIF
99         CONTINUE

100     CONTINUE

        RETURN
        END

        SUBROUTINE UPLIST ( SIGMA, I )
	implicit integer*8 (I-N)

        COMMON / BLOCK1 / RX, RY, VX, VY
        COMMON / BLOCK2 / COLTIM, PARTNR
        COMMON / MEZCLA / sigmas,mass,type
	COMMON / NUMBER / NN

C    *******************************************************************
C    ** LOOKS FOR COLLISIONS WITH ATOMS J > I                         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 50 000 )

        REAL*8      TIMBIG
        PARAMETER ( TIMBIG = 1.0E10 )

        INTEGER*8     I
        REAL*8      SIGMA
        REAL*8      RX(N), RY(N), VX(N), VY(N)
        REAL*8      sigmas(N),mass(N)
	integer*8   type(N)
        REAL*8      COLTIM(N)
        INTEGER*8     PARTNR(N)

        INTEGER     J
        REAL*8      RXI, RYI, RXIJ, RYIJ
        REAL*8      VXI, VYI, VXIJ, VYIJ
        REAL*8      RIJSQ, VIJSQ, BIJ, TIJ, DISCR, SIGSQ

C    *******************************************************************

        IF ( I .EQ. NN ) RETURN

C        SIGSQ = SIGMA ** 2
        COLTIM(I) = TIMBIG
        RXI = RX(I)
        RYI = RY(I)
        VXI = VX(I)
        VYI = VY(I)

        DO 100 J = I + 1, NN
        SIGSQ = ((sigmas(i)+sigmas(j))/2) ** 2

           RXIJ = RXI - RX(J)
           RYIJ = RYI - RY(J)
c		print*, "RXIJ, RYIJ", RXIJ,RYIJ
           RXIJ = RXIJ - ANINT ( RXIJ )
           RYIJ = RYIJ - ANINT ( RYIJ )
           VXIJ = VXI - VX(J)
           VYIJ = VYI - VY(J)
c		print*, "VXIJ, VYIJ", VXIJ,VYIJ
           BIJ  = RXIJ * VXIJ + RYIJ * VYIJ
c		print*, "BIJ",BIJ

           IF ( BIJ .LT. 0.0 ) THEN

              RIJSQ = RXIJ ** 2 + RYIJ ** 2
              VIJSQ = VXIJ ** 2 + VYIJ ** 2
              DISCR = BIJ ** 2 - VIJSQ * ( RIJSQ - SIGSQ )

              IF ( DISCR .GT. 0.0 ) THEN
                 TIJ = ( -BIJ - SQRT ( DISCR ) ) / VIJSQ
                 IF ( TIJ .LT. COLTIM(I) ) THEN
                    COLTIM(I) = TIJ
                    PARTNR(I) = J
                 ENDIF
              ENDIF

           ENDIF
C	print*, TIJ
100     CONTINUE
        RETURN
        END



        SUBROUTINE DNLIST ( SIGMA, J )
	implicit integer*8 (I-N)

        COMMON / BLOCK1 / RX, RY, VX, VY
        COMMON / BLOCK2 / COLTIM, PARTNR
        COMMON / MEZCLA / sigmas,mass,type
	COMMON / NUMBER / NN

C    *******************************************************************
C    ** LOOKS FOR COLLISIONS WITH ATOMS I < J                         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 50 000 )

        REAL*8      TIMBIG
        PARAMETER ( TIMBIG = 1.E10 )

        INTEGER*8     J
        REAL*8      SIGMA
        REAL*8      RX(N), RY(N), VX(N), VY(N)
        REAL*8      COLTIM(N)
        INTEGER*8     PARTNR(N)
        REAL*8      sigmas(N),mass(N)
	integer*8   type(N)

        INTEGER     I
        REAL*8      RXJ, RYJ, RXIJ, RYIJ
        REAL*8      VXJ, VYJ, VXIJ, VYIJ
        REAL*8      RIJSQ, VIJSQ, BIJ, TIJ, DISCR, SIGSQ

C    *******************************************************************

        IF ( J .EQ. 1 ) RETURN

c        SIGSQ = SIGMA ** 2
        RXJ = RX(J)
        RYJ = RY(J)
        VXJ = VX(J)
        VYJ = VY(J)

        DO 100 I = 1, J - 1
        SIGSQ = ((sigmas(i)+sigmas(j))/2) ** 2

           RXIJ = RX(I) - RXJ
           RYIJ = RY(I) - RYJ
           RXIJ = RXIJ - ANINT ( RXIJ )
           RYIJ = RYIJ - ANINT ( RYIJ )
           VXIJ = VX(I) - VXJ
           VYIJ = VY(I) - VYJ
           BIJ  = RXIJ * VXIJ + RYIJ * VYIJ 

           IF ( BIJ .LT. 0.0 ) THEN

              RIJSQ = RXIJ ** 2 + RYIJ ** 2
              VIJSQ = VXIJ ** 2 + VYIJ ** 2
              DISCR = BIJ ** 2 - VIJSQ * ( RIJSQ - SIGSQ )

              IF ( DISCR .GT. 0.0 ) THEN
                 TIJ = ( - BIJ - SQRT ( DISCR ) ) / VIJSQ
                 IF ( TIJ .LT. COLTIM(I) ) THEN
                    COLTIM(I) = TIJ
                    PARTNR(I) = J
                 ENDIF
              ENDIF
           ENDIF
100     CONTINUE

        RETURN
        END

        SUBROUTINE BUMP ( SIGMA, I, J, E, lx, b, alpha, delta)
	implicit integer*8 (I-N)

        COMMON / BLOCK1 / RX, RY, VX, VY
        COMMON / MEZCLA / sigmas,mass,type
	COMMON / NUMBER / NN
	Common/parms/ ss1,ss2,mm1,mm2,al11,al12,al22,del11,del12,del22
	REAL*8 ss1,ss2,mm1,mm2,al11,al12,al22,del11,del12,del22
	REAL*8  m12,m21,a12,d12

C    *******************************************************************
C    ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
C    ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **  
C    ** Also computes collision parameter, b                          **
C    *******************************************************************

        INTEGER*8     N,NN
        PARAMETER ( N = 50 000 )

        INTEGER*8     I, J
        REAL*8      SIGMA, EPS, E, b, lx
        REAL*8      RX(N), RY(N), VX(N), VY(N)

        REAL*8      RXIJ, RYIJ, VXIJ, VYIJ, FACTOR
        REAL*8      DELVX, DELVY, SIGSQ
        REAL*8      DELVX1, DELVY1, DELVX2, DELVY2
        REAL*8      alpha, delta, delta1,delta2, xx
        REAL*8      sigmas(N),mass(N)
        REAL*8      sigmag,sisj
	integer*8   type(N)

C define the collision alpha, delta ans reduced masses
	a12=al22
	d12=del22
	if (type(i).ne.type(j)) then
		a12=al12 
		d12=del12
		else 
		if (type(i).eq.1) a12=al11
		if (type(i).eq.1) d12=del11
		endif

	m12=mass(i)/(mass(i)+mass(j))
	m21=mass(j)/(mass(i)+mass(j))
C define sigma, relative velocities, etc...
c	print*, "in BUMP, i,j, a12,d12,m12,m21", i,j,a12,d12,m12,m21
c	print*, "in BUMP, i,j, a12,d12,m12,m21", i,j,a12,d12,m12,m21

        SIGSQ = SIGMA * SIGMA
        RXIJ = RX(I) - RX(J)
        RYIJ = RY(I) - RY(J)
        RXIJ = RXIJ - DNINT ( RXIJ )
        RYIJ = RYIJ - DNINT ( RYIJ )
	VXIJ = VX(I) - VX(J)
	VYIJ = VY(I) - VY(J)
c	print*, "Posc Part i,j,",i,j ,rx(i),ry(i),rx(j),ry(j)
c	print*, "Vels Part ",i,j ,vx(i),vy(i),vx(j),vy(j)

	xx=vx(i)**2+vx(j)**2+vy(i)**2+vy(j)**2

C COMPUTES IMPACT PARAMETER
c	b=lx*(RYIJ*VXIJ - RXIJ*VYIJ)/SQRT(VXIJ**2+VYIJ**2)
c	if ( (b .lt. -1.0) .or. (b .gt. 1.0) ) then
c	WRITE (50,1234) b
c	write (50,1235) i,j,rxij, ryij, vxij, vyij
c	endif
1234	FORMAT (F9.6)
1235	format (i5, i5, f12.8, f12.8, f9.6, f9.6)
C END OF IMPACT PARAMETER

C  sigma media: sisj
	sisj =(sigmas(i)+sigmas(j))/2.0d0

         sigmag =  ( RXIJ * VXIJ + RYIJ * VYIJ ) / sisj
        FACTOR = - ( RXIJ * VXIJ + RYIJ * VYIJ ) / SIGSQ *
     @              (1.0d0+alpha)*0.5d0
c	print*, "In Bump, sigmag=", sigmag
	sigmag=-sigmag
        DELVX1 =  RXIJ/sisj * ( -m21*(1+a12)*sigmag-2*m21*d12)
        DELVY1 =  RYIJ/sisj * ( -m21*(1+a12)*sigmag-2*m21*d12)
        DELVX2 =  RXIJ/sisj * ( +m12*(1+a12)*sigmag+2*m12*d12)
        DELVY2 =  RYIJ/sisj * ( +m12*(1+a12)*sigmag+2*m12*d12)
C	print*,"BumpDistance:",rxij**2+ryij**2,((sigmas(i)+sigmas(j))/2)**2
c	print*, "In Bump", delvx1,delvy1,delvx2,delvy2

        VX(I) = VX(I) - DELVX1
        VX(J) = VX(J) - DELVX2
        VY(I) = VY(I) - DELVY1
        VY(J) = VY(J) - DELVY2
c	print*, "Vels PostColl Part ",i,j ,vx(i),vy(i),vx(j),vy(j)

	e=e-(xx-(vx(i)**2+vx(j)**2+vy(i)**2+vy(j)**2))/2

        RETURN
        END

	subroutine INIT(kt,irs,lx)
	common /block1/ rx, ry, vx, vy
        COMMON / MEZCLA / sigmas,mass,type
	common /number/ nn

        INTEGER*8     N,NN, irs
        PARAMETER ( N = 50 000 )

        REAL*8      kt, RX(N), RY(N), VX(N), VY(N), lx, vxt, vyt
        REAL*8      sigmas(N),mass(N)
	real*8 distance , ET,Tt
	integer*8   type(N)

c	print*, sigmas(1),mass(1),type(1)

	call ranint(irs)

	skt=kt**0.5
	
	vxt=0.0d0
	vyt=0.0d0
	do 1 i=1,nn
	aux=(-2*log(1-rand(0)))**0.5
	raux=2*3.1415927*ranf(0)
	vx(i)=aux*cos(raux)*skt 
	vy(i)=aux*sin(raux)*skt 
	vxt=vxt+vx(i)
	vyt=vyt+vy(i)
 1	continue 

	do 8 i=1,nn
	vx(i)=vx(i)-vxt/nn
	vy(i)=vy(i)-vyt/nn
 8	continue 


	ET=0
        do ii=1,nn
        ET=ET+vx(ii)**2+vy(ii)**2
        enddo
        Tt=ET/2/nn
	do ii=1,nn
	vx(ii)=vx(ii)*Sqrt(kt/Tt)
	vy(ii)=vy(ii)*Sqrt(kt/Tt)
	enddo


	icount=1
	rx(1)=lx*rand(0)
	ry(1)=lx*rand(0)
	do 2 i=2,nn
 3  	icount=icount+1
	if (icount.gt.1000*nn) then 
		print*, ' Problems with the initial condition. Aborting...'
		print*, ' Check the density or change the seed.'
		stop
		endif
 	rx(i)=lx*rand(0)
 	ry(i)=lx*rand(0)
 	rx(i)=lx*rand(0)
 	ry(i)=lx*rand(0)
	iover=0
		do j=1,i-1
		distance=((sigmas(i)+sigmas(j))/2)**2
		d=(rx(i)-rx(j))**2+(ry(i)-ry(j))**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j)+lx)**2+(ry(i)-ry(j))**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j)-lx)**2+(ry(i)-ry(j))**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j))**2+(ry(i)-ry(j)+lx)**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j))**2+(ry(i)-ry(j)-lx)**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j)+lx)**2+(ry(i)-ry(j)+lx)**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j)+lx)**2+(ry(i)-ry(j)-lx)**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j)-lx)**2+(ry(i)-ry(j)+lx)**2
		if (d.lt.distance) iover=1 
		d=(rx(i)-rx(j)-lx)**2+(ry(i)-ry(j)-lx)**2
		if (d.lt.distance) iover=1 
		enddo
	if(iover.eq.1) goto 3
 2	continue

c	print*, "*** HERE *** " 

c	do ip=1,nn
c 	print*, lx,ip,rand(0),    rx(ip),ry(ip),sigmas(ip)/2
c	write(90,*) mass(ip)
c	enddo


	do 4 ip=1, nn
	rx(ip)=rx(ip)/lx-0.50d0
	ry(ip)=ry(ip)/lx-0.50d0
	vx(ip)=vx(ip)/lx
	vy(ip)=vy(ip)/lx
 4  	continue

	return
	end 
	
	subroutine maxwell( lx )
	common /block1/ rx, ry, vx, vy
	common /number/ nn

        INTEGER*8     N,NN, nmax
        PARAMETER ( N = 50 000 , nmax=20 )

        REAL*8      RX(N), RY(N), VX(N), VY(N), lx, v2, v4, sum,v6
	real*8	    v2x, v2y, v4x, v4y, sumx, sumy
	integer  maxw(0:nmax)

	do 1 i=0,nmax
	maxw(i)=0
 1	continue

	vmax=0
	v2=0.0d0
	v2x=0.0d0
	v2y=0.0d0
	v4=0.0d0
	v4x=0.0d0
	v4y=0.0d0
	v6=0.0d0
	do 2 i=1,nn
	xx=vx(i)**2+vy(i)**2
	if (xx.gt.vmax) vmax=xx
	sum=vx(i)*vx(i)+vy(i)*vy(i)
	sumx=vx(i)*vx(i)
	sumy=vy(i)*vy(i)
	v2=v2+sum
	v2x=v2x+sumx
	v2y=v2y+sumy
	v4=v4+sum*sum
	v4x=v4x+sumx*sumx
	v4y=v4y+sumy*sumy
	v6=v6+sum*sum*sum
 2	continue
	vmax=vmax**0.5
c	write(40,50) lx*lx*v2/nn,lx*lx*lx*lx*v4/nn,nn*v4/v2/v2,
c     @    nn*v4x/v2x/v2x, nn*v4y/v2y/v2y, 1.0d0*nn*nn*v6/v2/v2/v2 
 50	format(e14.5,e14.5,e14.5,e14.5,e14.5,e14.5)

c	do 3 i=1,nn
c	v=(vx(i)**2+vy(i)**2)**0.5 /vmax * nmax
c	maxw(v)=maxw(v)+1
c 3	continue

c	do 4 i=0,nmax
c	write(40,5)(i+0.5)*vmax/nmax*lx,maxw(i)/float(nn)/(vmax/nmax*lx)
c 4	continue
c 5	format(1p2e16.8)

	return
	end

c------------------------------------------------------------
c------------------------------------------------------------
c   From here on is the random number generator.        
c------------------------------------------------------------
c------------------------------------------------------------

      function ranf(iran)

c     berkeley random number generator
c     range changed to 0 < 1

      implicit real*8 (a-h,o-z)
      real*4 ranf
      parameter (m3=647,m2=1442,m1=3707,m0= 373)
      parameter (t4=2.0**48,t3=2.0**36,t2=2.0**24,t1=2.0**12)
      parameter (mm=4096)

c     mm = 2 ** 12
      common /rjran/ i3,i2,i1,i0

c     the random number is:
c     (i3*t3+i2*t2+i1*t1+i0)/2.0**48
c     the multiplier is:
c     (m3*t3+m2*t2+m1*t1+m0)

      ranf = float(i3)/t1 + float(i2)/t2 + float(i1)/t3 + float(i0)/t4
      if(ranf.ge.0.9999999)  ranf = 0.0

c     multiply i's and m's

      j0 = m0 * i0
      j1 = m0 * i1 + m1 * i0
      j2 = m0 * i2 + m1 * i1 + m2 * i0
      j3 = m0 * i3 + m1 * i2 + m2 * i1 + m3 * i0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      i0 = mod(k0,mm)
      i1 = mod(k1,mm)
      i2 = mod(k2,mm)
      i3 = mod(k3,mm)
      return
      end


      subroutine ranint(istart)

c     initialize random number generator

      implicit real*8 (a-h,o-z)
      dimension iran(4)
	integer*8 istart

      iran(1)=12345 + istart*1000
      iran(2)=12345 + istart*1000
      iran(3)=12345 + istart*1000
      iran(4)=12345 + istart*1000
      call ranset(iran)
      call ranget(iran)

      return
      end


      subroutine ranset(iran)

      implicit real*8 (a-h,o-z)
      parameter (mm = 4096, nn = 100000)
      dimension iran(4)
      common /rjran/ ii3,ii2,ii1,ii0
      i3 = iran(1)
      i2 = iran(2)
      i1 = iran(3)
      i0 = iran(4)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii0)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii1)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii2)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii3)
      return
      end


      subroutine ranget(iran)

      implicit real*8 (a-h,o-z)
      parameter (mm = 100000)
      parameter (m10 = 4096)
      parameter (m21 =  167, m20 = 77216)
      parameter (m32 =    6, m31 = 87194, m30 = 76736)
      dimension iran(4)
      common /rjran/ ii3,ii2,ii1,ii0
      j0 = ii0 + m10 * ii1 + m20 * ii2 + m30 * ii3
      j1 =                   m21 * ii2 + m31 * ii3
      j2 =                               m32 * ii3
      j3 =                                       0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      iran(4) = mod(k0,mm)
      iran(3) = mod(k1,mm)
      iran(2) = mod(k2,mm)
      iran(1) = mod(k3,mm)
      return
      end


      subroutine divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,id,ir)

c     given the integer i = i0 + nn * (i1 + nn * (i2 + nn * i3))
c     this routine calculates j = i / id and ir = mod(i, id)
c     j is expressed as i, ir is just an integer

      implicit real*8 (a-h,o-z)

      j3 = i3 / id
      k3 = mod(i3, id)
      k2 = k3 * nn + i2
      j2 = k2 / id
      k2 = mod(k2, id)
      k1 = k2 * nn + i1
      j1 = k1 / id
      k1 = mod(k1, id)
      k0 = k1 * nn + i0
      j0 = k0 / id
      ir = mod(k0, id)
      return
      end


 

