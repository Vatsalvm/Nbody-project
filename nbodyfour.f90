program nbodyfour
implicit none
!fourth order predictor-corrector method used here

!declaration
double precision r(1:3,1:3,-8:8),v(1:3,1:3,-8:8),a(1:3,1:3,-8:8),a_temp(1:3,1:3,0:8),m(1:3),dt,tcurr,tout,G
double precision r_pf(1:3,1:3,-8:8), v_pf(1:3,1:3,-8:8),a_pf(1:3,1:3,-8:8) !r_pf, v_pf, a_pf store the past information and the future positions
double precision dr(1:3), rad2, AU, yr, com(1:3), cov(1:3), mtot, KE, PE, Ei, Ecurr, delE, rad, tcount, dsun(1:3,-8:1)
double precision techeck, tdsun, day
integer i,j,k,n,p

!----------------------------------------
!number of objects
n = 3
!Year (in s) and a astronomical unit (in m)
yr = 365.*24.*3600.
AU = 1.5e11
!Orbital distances for the three-body system of the Sun-Jupiter-Earth (in m)
!r_past, v_past and a_past arrays store the past information
r(1,1:n,0) = (/0., 778.6e9, 149.6e9/)
r(2:3,1:n,0) = 0.
!Mass of the bodies (in Kg)
m(1:n) = (/1.989e30, 1.898e27, 5.97e24/)
!Total mass of the system (in Kg)
mtot = 0.
!Orbital velocities of the bodies in the system (in m/s)
v(2,1:n,0) = (/0., 13.1e3, 29.8e3/)
v(1,1:n,0) = 0. !forcing the objects to orbit in the XY plane
v(3,1:n,0) = 0.
!Gravitational constant (in N m^2 Kg^-2)
G = 6.67e-11
!Timestep (in s)
dt = 1000.
!Current time (in s)
tcurr = 0.
!Time for all simulations to run (in s)- setting 10. s so can get enough future steps which then will be shifted to past steps
tout = 10.*yr
!Current distance from the sun (in m)
dsun(1:n,-8:1) = 0.
!Time after which fractional energy is determined (in s). Fractional energy determined after every 6 months
techeck = 0.5 * yr
!Time after which separation from sun measured (in s)
tdsun = 30. * 24. * 3600.
!Counter for the simulations (in s)
tcount = 0.

!--------------------------------------
!Initial conditions
dr(1:3) = 0.
com(1:3) = 0.
cov(1:3) = 0.
!Initial energy (in J)
Ei = 0.
!Kinetic and Potential energy (in J)
KE = 0.
PE = 0.
!Distance between each object (in m)
rad = 0.
!Distance squared between each object (in m2)
rad2 = 0.
!Current energy (in J)
Ecurr = 0.
!Fractional energy (in J)
delE = 0.

!Position,velocity and acceleration arrays which store past information and the future values
r_pf(1:3,1:n,-8:8) = 0.
v_pf(1:3,1:n,-8:8) = 0.
a_pf(1:3,1:n,-8:8) = 0.

!-----------------------------------
!Nbody units


!-----------------------------------
!Centre of mass correction
do i=1,n
    com(1:3) = com(1:3) + m(i)*r(1:3,i,0)
    cov(1:3) = cov(1:3) + m(i)*v(1:3,i,0)
    mtot = mtot + m(i)
end do

com = com/mtot
cov = cov/mtot

do i=1,n
    r(1:3,i,0) = r(1:3,i,0) - com(1:3)  !considering t=0 time for position and velocity of the objects which implies to initial values.
    v(1:3,i,0) = v(1:3,i,0) - cov(1:3)  !shifting the position and velocity of the objects based on the centre-of-mass position and velocity
!Now the objects are in the CoM frame
end do

!----------------------------------
!Initial energy determination
do i=1,n
    KE = KE + 0.5*m(i)*(v(1,i,0)*v(1,i,0) + v(2,i,0)*v(2,i,0) + v(3,i,0)*v(3,i,0))
end do

do i=1,n-1
    do j=i+1,n
        dr(1:3) = r(1:3,i,0) - r(1:3,j,0)  !separation between each body in x,y,z directions
        rad = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)) 
        PE = PE - (G*m(i)*m(j)/rad)
    end do
end do

Ei = PE + KE                          !Initial Energy
!write(6,*) 'Initial energy=',Ei
!----------------------------------    
!Initial acceleration calculation at time t=0 s
dr(1:3) = 0.
do i=1,n
    a(1:3,i,0) = 0.                  !Initial acceleration at t=0
    do j=1,n
        if(i==j) cycle                    !to avoid the case of finding separation between the object and itself
        do k=1,3
            dr(k) = r(k,i,0) - r(k,j,0)
        end do
        rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)                  !r^2 = x^2 + y^2 + z^2
    
        do k=1,3                                                        !runs over the all directions
            a(k,i,0) = a(k,i,0) - (G*m(j)*dr(k)/(sqrt(rad2)*rad2))    !dr(k) is the unit vector indicating the direction of the force from
                                                                        !objects 'j' to 'i'
        end do
    end do
end do
!write(6,*) 'Initial acceleration=',a(1,1:n,0)
!a(1:3,1:n,0) = 0.
!a(1:3,1:n,0) = a_temp(1:3,1:n,0)      !stores initial acceleration in main acceleration array 'a'. This will save a0 from getting overwritten.

!!!!!!!!!!!!!!shifting the t=0 step to t=-8 step
r(1:3,1:n,-8) = r(1:3,1:n,0)
v(1:3,1:n,-8) = v(1:3,1:n,0)
a(1:3,1:n,-8) = a(1:3,1:n,0)
!write(6,*) '%%%%%%%%%%%%%%%%%%%%%'
!write(6,*) 'Initial acceleration of jupiter:A=',a(1:3,2,-8)
!write(6,*) '%%%%%%%%%%%%%%%%%%%%%%'
!----------------------------------
!Bootstrapping the second order code- to get the past steps required for predictor-corrector method
!Time loop created which iterates over 'p'
do p=-7,0            !keeping track of the steps in time
    !Future position obtained from the initial position,velocity and acceleration
    do i=1,n
        r(1:3,i,p) = r(1:3,i,p-1) + v(1:3,i,p-1)*dt + 0.5*a(1:3,i,p-1)*dt*dt      !p-1 implies t=0 step.
    end do
    !write(6,*) '---------------------------'
    !write(6,*) 'A:New x position of Jupiter=',r(1,2,p)/AU
    !write(6,*) '---------------------------'
    
    !Future acceleration loop
    do i=1,n
        a(1:3,i,p) = 0.
        do j=1,n
            if(i==j) cycle
            do k=1,3
                dr(k) = r(k,i,p) - r(k,j,p)
            end do
            rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
            !write(6,*) 'rad2=',rad2 
            !STOP
            do k=1,3
                a(k,i,p) = a(k,i,p) - (G*m(j)*dr(k)/(sqrt(rad2)*rad2))        !get the future acceleration
            end do
        end do
    end do
    !write(6,*) '$$$$$$$$$$$$$$$$$$$$$$$'
    !write(6,*) 'B:New y acceleration of Jupiter=',a(1:3,2,p)
    !write(6,*) '$$$$$$$$$$$$$$$$$$$$$$$'
    
    !Future velocity
    do i=1,n
        v(1:3,i,p) = v(1:3,i,p-1) + 0.5*(a(1:3,i,p-1) + a(1:3,i,p))*dt
    end do
    
    !a0(1:3,1:3,p) = a(1:3,1:3,p)    !stores the future acceleration 'a' into 'a0'. Here a0 is the temporary array used in this loop.
    !write(6,*) '##########################'
    !write(6,*) 'C:New y velocity of Jupiter=',v(2,2,p)
    !write(6,*) '##########################'
    !Now to shift the current and future steps to the past positions
    !like shifting t=0 to t=-8 step
    !r_pf(1:3,1:n,-p) = r(1:3,1:n,p) !for p=1 iteration, (-1) step taken the value of (+1) step
    !v_pf(1:3,1:n,-p) = v(1:3,1:n,p) !?????????
    !a_pf(1:3,1:n,-p) = a(1:3,1:n,p)

    !-------------------
    !Distance calculated from the Sun
    do i=1,n
        dsun(i,p) = sqrt(r(1,i,p)*r(1,i,p) + r(2,i,p)*r(2,i,p) + r(3,i,p)*r(3,i,p))
    end do
    !---------------------
    !Energy conservation check
    !determining the fractional energy
    tcount = tcount + dt
    if(tcount==techeck) then
        KE = 0.
        do i=1,n
            KE = KE + 0.5*m(i)*(v(1,i,p)*v(1,i,p) + v(2,i,p)*v(2,i,p) + v(3,i,p)*v(3,i,p))
        end do

        PE = 0.
        do i=1,n-1
            do j=i+1,n
                dr(1:3) = r(1:3,i,p) - r(1:3,j,p)
                rad = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                PE = PE - (G*m(i)*m(j)/rad)
            end do
        end do
        
        Ecurr = PE + KE  !Total energy
        delE = (Ecurr-Ei)/Ei
        tcount = 0.
        
        
    end if
    tcurr = tcurr + dt
    !if(tcurr>tout) exit
end do

!Predictor-corrector method

end program nbodyfour

        
        





