
	integer*8 it,NINEL,i,j
	parameter (it=50000000)
	real*8 rx(it),ry(it),t(-1:it)
	real*8 vprom(0:it),tprom(-1:it)
	real*8 mm1,mm2,lx,kt,dens,nu,DE,temp,collfreq
	integer num(0:it)
	character aa*20


        WRITE(*,'(''# ENTER NUMBER OF PARTICLES of each type     '')')
        READ (*,*) aa,aa,aa,NN1,nn2
        WRITE(*,'(''# PARTICLES of each type '' ,2I7)'    )nn1,nn2
        WRITE(*,'(''# ENTER SIGMAS                               '')')
        READ (*,*) aa,SS1,SS2
        WRITE(*,'(''#  SIGMAS  '',1p2e12.4)') ss1,ss2
        WRITE(*,'(''# ENTER MASSES                               '')')
        READ (*,*) aa,MM1,MM2
        WRITE(*,'(''# MASSES  '',1p2e12.4)')mm1,mm2
        WRITE(*,'(''# ENTER ALPHA11,ALPHA12,ALPHA22              '')')
        READ (*,*) aa,al11,al12,al22
        WRITE(*,'(''# ALPHA11,ALPHA12,ALPHA22  '',1p3e12.4)') al11,al12,al22
        WRITE(*,'(''# ENTER DELTA11,DELTA12,DELTA22              '')')
        READ (*,*) aa,del11,del12,del22
        WRITE(*,'(''# DELTA11,DELTA12,DELTA22  '',1p3e12.4)') DEL11,DEL12,DEL22
        WRITE(*,'(''# ENTER BOX LENGTH                           '')')
        READ (*,*) aa,aa,LX
        WRITE(*,'(''# BOX LENGTH '',1pe12.4)') LX
        WRITE(*,'(''# ENTER NUMBER OF ELASTIC COLLISIONS         '')')
        READ (*,*) aa,aa,NELAS
        WRITE(*,'(''# ENTER NUMBER OF INELASTIC COLLISIONS       '')')
        READ (*,*) aa,aa,NINEL
        WRITE(*,'(''# ENTER TEMPERATURE kT                       '')')
        READ (*,*) aa,aa,kt
        WRITE(*,'(''# ENTER SEED OF RANDOM NUMBER GENERATOR      '')')
        READ (*,*) aa
C
C  DE is the difussion coefficient - García Rojo, Luding, Brey, PRE74, 061305 (2006)
C

	t(0)=0
	t(-1)=0

	do i=1, NINEL-1
	read(*,*) t(i),rx(i),ry(i)
c 	print*, i,t(i)
	enddo

	do i=0,NINEL-1
	vprom(i)=0
	tprom(i)=0
	num(i)=0
	enddo
	read (*,*) aa,aa,aa,aa,Temp
	WRITE(*,'(''# Temperature '',1pe14.6)') Temp
	
C	print*, "finished reading"
	dens=(NN1+NN2)/Lx/LX
	nu=dens*3.1415926536/4
	DE=(Temp/3.1415927*mm2)**0.5/(2*dens*ss2*(1-7*nu/16.))*(1-nu)**2
	collfreq=(2*3.1415927)**0.5*dens*ss2*(1-7*nu/16.)/(1-nu)**2*(2*Temp)**0.5
	WRITE(*,'(''# Elastic Coll freq:  '',1p3e12.4)') collfreq
	WRITE(*,'(''# Elastic Diff Coeff: '',1p3e12.4)') DE

	collfreq=2.0d0*NINEL/t(NINEL-1)/(nn1+nn2)
	WRITE(*,'(''# Calculated Coll freq:  '',1p3e12.4)') collfreq

	do i=1, NINEL-1
	if (t(i).gt.50.0d0/collfreq) then
	istep=i
	goto 3333
	endif
	enddo
 3333	continue


C I use the real coll freuqnecy and the istep calculated above..
c	collfreq=(2*3.1415927)**0.5*dens*ss2*(1-7*nu/16.)/(1-nu)**2*(2*Temp)**0.5
c	istep= 50.0 *(nn1+nn2)/collfreq

	write(*,666) collfreq, istep
 666    format('#  collfreq=', 1pe12.6,'    istep=', i8)

	write(*,'(''#        TIME         Mean SQ Displacemt     Stats'')')

 	lag=100
	do i=1,NINEL-2*istep ,lag
c	print*, i
	do j=0,istep
	despl=(rx(i)-rx(i+j))**2 + (ry(i)-ry(i+j))**2
	vprom(j)=vprom(j)+despl
	tprom(j)=tprom(j)+t(i+j)-t(i)
	num(j)=num(j)+1
	enddo
	enddo

	do i=0,istep
	tprom(i)=tprom(i)/num(i)
	vprom(i)=vprom(i)/num(i)
	write(*,44) tprom(i),vprom(i),num(i)
	enddo
c	do i=1,NINEL/2,istep
c	print*, tprom(i),vprom(i),num(i)
c	enddo

  44 	format(1p2e19.8,i10)

	tmin=10.0d0/collfreq
	tmax=20.0d0/collfreq
	xy=0
	x2=0
	y2=0
	xm=0
	ym=0
	nn=0
	do i=0, istep
	if ( (tprom(i).gt.tmin).and.(tprom(i).lt.tmax)) then 
		nn=nn+1
		xy=xy+tprom(i)*vprom(i)
		x2=x2+tprom(i)**2
		y2=y2+vprom(i)**2
		xm=xm+tprom(i)	
		ym=ym+vprom(i)	
		endif
	enddo

	D=x2-xm*xm/nn
	E=xy-xm*ym/nn
	ord=ym/nn-e/d*xm/nn
	
	D0S=e/d/4.0d0*dens*(2.0d0/Temp)**0.5

	write(*,54) tmin, tmax, nn
  54	format('#  tmin:', 1pe14.5,'   tmax:',1pe14.5,'  number of points:',i8)
	write(*,55) 
  55 	Format('#   density      Temperature       slope       intercept      D0S   ')
	write(*,56)  dens, temp, e/d, ord, D0S
  56	format('#', 1p7e14.5)


	

	stop
	end
