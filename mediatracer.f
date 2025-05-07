	character a 
	real*8 x,xm,xm2


	xm=0.d0
	xm2=0.d0
	do i=1,10
	read(*,*) a,x,x,x,x,x
	xm=xm+x
	xm2=xm2+x*x
	enddo
	x1=xm/10.d0
	y1=((xm2-xm*xm/10)/9)**0.5

        xm=0.d0
        xm2=0.d0
        do i=1,10
        read(*,*) a,x,x,x,x,x
        xm=xm+x
        xm2=xm2+x*x
        enddo
        x2=xm/10.d0
        y2=((xm2-xm*xm/10)/9)**0.5


	print*, x1,y1,x2,y2

	stop
	end
