program nbody
IMPLICIT NONE
!three body system of SUN-JUPITER-EARTH
!Writing a simple second order code
double precision :: r(1:3,1:8), v(1:3,1:8), a(1:3,1:8), a0(1:3,1:8), m(1:8), dt, tcurr, tout, G, dr(1:3), dsun(1:8)
!r,v,a are the positions, velocities and accelerations of the bodies.
!dt and tcurr are the timestep and current time in the simulations.
!G is the gravitational constant.
double precision :: rad2, AU, yr, com(1:3), cov(1:3), mtot, KE, PE, Ei, Ecurr, delE, rad, tcount, techeck, tdcount, td
!AU is the astronomical unit and rad2 is the sum of the squares of distances in x,y,z directions
integer :: i,j,k,n
!i,j and k are the indexes for iteration.

!-----------------------------------------------------------------------------
!number of objects
n=8
!Mass of SUN-JUPITER-EARTH-VENUS-MARS-SATURN-URANUS-NEPTUNE(in Kg)
m(1:n) = (/1.989e30, 1898.e24, 5.97e24, 4.87e24, 0.642e24, 568.e24, 86.8e24, 102.e24/) 
!Astronomical unit
AU = 1.5e11 !(in m)
!Orbital distances of SUN-JUPITER-EARTH-VENUS-MARS-SATURN-URANUS-NEPTUNE (in m)
r(1,1:n) = (/0., 778.6e9, 149.6e9, 108.2e9, 227.9e9, 1433.5e9, 2872.5e9, 4495.1e9/)
r(2:3,1:n) = 0. !forcing the objects to remain in the XY plane
!Orbital velocities of SUN-JUPITER-EARTH-VENUS-MARS-SATURN-URANUS-NEPTUNE in the solar system (in m/s)
v(2,1:n) = (/0., 13.1e3, 29.8e3, 35.0e3, 24.1e3, 9.7e3, 6.8e3, 5.4e3/)
v(1,1:n) = 0. 
v(3,1:n) = 0. !allowing the object to move in circular paths on the XY plane.
!dsun
dsun(1:n) = 0.
!Gravitational constant (in N m^2 Kg^-2)
G = 6.67e-11 
!timestep, dt (in sec)
dt = 10.
!Storing a current time
tcurr = 0.
!Time for the simulations, tout = 1000.0 yrs
tout = 1000.*365.*24.*3600.
!conversion from yr to seconds
yr = 365.*24.*3600.

!time after which the dsun is stored; storing dsun after every 1000th iteration (1000*1000s)
!implies taking measurement after 12 days
td = 1.e6 
!Time after which to do energy check (approx. 6 months)
techeck = 6.*30.*24.*3600.
!--------------------------------------------------------------------------------
!setting the initial conditions
dr(1:3) = 0. !separation between the objects in the three directions
com(1:3) = 0.
cov(1:3) = 0.
mtot = 0.
tcount = 0. !counter for energy check
tdcount = 0. !for the dsun calculation
Ei = 0.
Ecurr = 0.
delE = 0.

!Creating a file for storing the x and y position values
!OPEN(3,file='pcorr.txt',status='new')

!Creating a file for storing the distance from the sun and time values
!OPEN(4,file='Dsunnew.txt',status='new')

!------------------------------------------------------
!Doing the centre of mass/velocity correction
do i=1,n
	com(1:3) = com(1:3) + m(i)*r(1:3,i)
	cov(1:3) = cov(1:3) + m(i)*v(1:3,i)
	mtot = mtot + m(i)
end do
com = com/mtot
cov = cov/mtot
do i=1,n
	r(1:3,i) = r(1:3,i) - com(1:3) !shifts the positions and velocities of the three bodies basesd on the centre of mass and velocity 
	v(1:3,i) = v(1:3,i) - cov(1:3)
!We are in the COM frame
end do
!Write(3,*) tcurr*(1/yr), r(1:2,1:3)/AU !writing the initial position
!------------------------------------------------------
!Finding Initial Kinetic and Potential energy
KE = 0.
do i=1,n
	KE = KE + 0.5*m(i)*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2)
end do

PE = 0.
rad = 0.
do i=1,n-1
	do j=i+1,n
		dr(1:3) = r(1:3,i) - r(1:3,j)
		rad = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
		PE = PE - G*m(i)*m(j)/rad
	end do
end do

Ei = PE + KE !initial energy
!-------------------------------------------------------
!Finding initial acceleration,a0 of the SUN-JUPITER-EARTH-VENUS-MARS-SATURN
dr(1:3) = 0.
do i=1,n
	a0(1:3,i) = 0.
	do j=1,n
		if(i==j) cycle
		do k=1,3
		dr(k) = r(k,i) - r(k,j)
		end do
		
		rad2 = dr(1)**2 + dr(2)**2 + dr(3)**2
        
		do k=1,3
			a0(k,i) = a0(k,i) - G*m(j)*(dr(k))/(sqrt(rad2)*rad2) !dr(k) represents the direction of the position vector
		end do
		
	end do
end do 

!Creating a file for storing current time and fractional energies
!OPEN(3,file='Frnew.txt',status='new')


!Introducing a time loop to iterate for each time step
do 
	!get the new position, r1
	do i=1,n
		r(1:3,i) = r(1:3,i) + v(1:3,i)*dt + 0.5*a0(1:3,i)*dt**2 
	end do
	!if(tcurr>0.50*tout) then
	!write(6,*) 'rsun=',sqrt(r(1,1)**2 + r(2,1)**2)/AU
	!end if
	!Finding new acceleration,a of the SUN-JUPITER-EARTH
	do i=1,n
		a(1:3,i) = 0.
		do j=1,n
			if(i==j) cycle
			do k=1,3
				dr(k) = r(k,i) - r(k,j) !this uses new position r1.
			end do
			
			rad2 = dr(1)**2 + dr(2)**2 + dr(3)**2
			
			do k=1,3
				!get the new acceleration, a1 
				a(k,i) = a(k,i) - G*m(j)*(dr(k))/(sqrt(rad2)*rad2) !dr(k) represents the direction of the position vector
			end do
			
		end do
	end do
	
	!get the new velocity, v1 using the second order expansion
	do i=1,n
		v(1:3,i) = v(1:3,i) + 0.5*(a0(1:3,i) + a(1:3,i))*dt
	end do
	a0 = a !store a into a0 for the next timestep
	
	
	
	!Distance from the sun against the current time
	do i=1,n
		dsun(i) = sqrt(r(1,i)**2 + r(2,i)**2)
	end do
	
	
	
	!Energy check
	!getting the current total energy
	tcount = tcount + dt
	if(tcount>techeck) Then
		KE = 0.
		do i=1,n
			KE = KE + 0.5*m(i)*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2)
		end do

		PE = 0.
		rad = 0.
		do i=1,n-1
			do j=i+1,n
				dr(1:3) = r(1:3,i) - r(1:3,j)
				rad = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
				PE = PE - G*m(i)*m(j)/rad
			end do
		end do
		
		Ecurr = KE + PE
		!if(tcurrdt) then
		!write(6,*) Ei,Ecurr
		!stop 'stopped'
		!end if
	
		!write(6,*) '---------'
		!write(6,*) 'A,KE=',KE
		!write(6,*) 'B,PE=',PE
		!if(tcurr>(0.5*tout)) then
		!	STOP  'stopped'
		!end if
		
		!fractional energy loss,delE
		delE = ((Ei - Ecurr))/Ei
		!write(6,*) '-------'
		!write(6,*) 'delE=', delE
		!if(tcurr>(0.5*tout)) then
		!	stop 'stopped'
		!end if
		
		tcount = 0. !reinitialising to zero for the next calculation of delE
		!STOP 'I AM HERE delE'
		!Write the current time(in yr) and the fractional energy
		!write(3,*) (tcurr*(1/yr)),delE
		
		!if(Ei==0) then
		!	Stop 'Stopped'
		!else
		!Write(6,*) 'FREAKY'
		!end if
		
	end if
	tcurr = tcurr + dt
	
	!for storing the dsun value for every 1000th iteration
	tdcount = tdcount + dt
	!if(tdcount==td) then
	!	!write(6,*) tcurr
	!	!STOP 'stopped'
	!	write(4,*) dsun*(1/AU), tcurr*(1/yr)
	!	tdcount = 0.
	!end if
	if(tcurr>(1e7*dt)) then 
	Write(6,*) 'A:',delE,Ei
	end if
	!write(3,*) tcurr*(1/yr),r(1:2,1:3)*(1/AU) !writing the new position
	if(tcurr>tout) exit
end do
		
		
	
!CLOSE(3)
!CLOSE(4)
!write(6,*) 'rsun_jup=',(sqrt(r(1,2)**2 + r(2,2)**2)-sqrt(r(1,1)**2 + r(2,1)**2))/AU
!write(6,*) 'rsun_earth=',(sqrt(r(1,3)**2 + r(2,3)**2)-sqrt(r(1,1)**2 + r(2,1)**2))/AU
!write(6,*) 'v_earth=',sqrt(v(1,3)**2 + v(2,3)**2)
!write(6,*) 'ajup_sun=',sqrt(a(1,2)**2 + a(2,2)**2)
!CLOSE(3)
!WRITE(6,*) 'shape of r array',SHAPE(r)
!WRITE (6,*) ' final positions of S-J-E:\n', r(1,1:3)
!WRITE (6,*) ' Orbital velocity of earth: ', v(2,3)
!WRITE (6,*) 'ACCELERATION:', a(3,1:3)
end program nbody

!DO THE FRACTIONAL ENERGY CHECKS
!THEN SAVE THE CODE AND DO THE BOOTSTRAPPING PART.
            
            
            