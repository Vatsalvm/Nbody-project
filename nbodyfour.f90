program nbodyfour
implicit none
!fourth order predictor-corrector method used here

!declaration
double precision r(1:3,1:8,-8:1),v(1:3,1:8,-8:1),a(1:3,1:8,-8:1),m(1:8),dt,tcurr,textracount,tout,G
double precision r_newpast(1:3,1:8,-3:-1), v_newpast(1:3,1:8,-3:-1), a_newpast(1:3,1:8,-3:-1)  !these arrays store the new position, velocity and acceleration produced when the timestep is halved. They can be considered as t=-1/2 step
!double precision r_threehalf(1:3,1:3,-1), v_threehalf(1:3,1:3,-1), a_threehalf(1:3,1:3,-1) !Similar to the above arrays, but these are taken as t=-3/2 step.
double precision r_c(1:3,1:8,1), v_c(1:3,1:8,1),a_c(1:3,1:8,1) !r_c, v_c, a_c store the future corrected values
double precision dr(1:3), rad2, AU, yr, com(1:3), cov(1:3), mtot, KE, PE, Ei, Ecurr, delE, rad, tcount, dsun(1:8)
double precision techeck, tdsun, day, relerr, small, err_fact_r(1:3,1:8), err_fact_v(1:3,1:8), err_max_r(1:8), err_max_v(1:8)
double precision err_r, err_v
integer i,j,k,n,p,q

!----------------------------------------
!counter for xy position
q = 1
!number of objects
n = 8
!a day
day = 24.*3600.
!Year (in s) and a astronomical unit (in m)
yr = 365.*24.*3600.
AU = 1.5e11
!Orbital distances for Sun-Jupiter-Earth-Venus-Mars-Saturn-Uranus-Neptune (in m)
r(1,1:n,0) = (/0., 778.6e9, 149.6e9, 108.2e9, 227.9e9, 1433.5e9, 2872.5e9, 4495.1e9/)
r(2:3,1:n,0) = 0.
!Mass of the bodies (in Kg)
m(1:n) = (/1.989e30, 1.898e27, 5.97e24, 4.87e24, 0.642e24, 568.e24, 86.8e24, 102.e24/)
!Total mass of the system (in Kg)
mtot = 0.
!Orbital velocities of the bodies in the system (in m/s)
v(2,1:n,0) = (/0., 13.1e3, 29.8e3, 35.0e3, 24.1e3, 9.7e3, 6.8e3, 5.4e3/)
v(1,1:n,0) = 0. !forcing the objects to orbit in the XY plane
v(3,1:n,0) = 0.
!Gravitational constant (in N m^2 Kg^-2)
G = 6.67e-11
!Timestep (in s)
dt = 500.
!Current time (in s)
tcurr = 0.
textracount = 0.
!Time for all simulations to run (in s)- setting 10. s so can get enough future steps which then will be shifted to past steps
tout = 1000.*yr
!Current distance from the sun (in m)
dsun(1:n) = 0.
!Time after which fractional energy is determined (in s). Fractional energy determined after every 6 months
techeck = 0.5*yr
!Time after which the separation from sun is written to the file (in s)
tdsun = 30. * 24. * 3600.
!Counter for the simulations (in s)
tcount = 0.
!relative error criterion for predictor-corrector (from the Adams-Bashforth-Moulton paper). This is determines whether to half or double the step-size, dt
relerr = 5.e-10
!'Small', an offset used in the error calculation from the predicted and corrected values
small = 1.0e-5
!err_fact- will be used in the error calculation using the predicted and corrected values
err_fact_r(1:3,1:n) = 0.     !position values are used here
err_fact_v(1:3,1:n) = 0.     !velocity values are used here
err_max_r(1:n) = 0.         !This stores the maximum error for the three objects among the errors for each x,y,z directions in position
err_max_v(1:n) = 0.         !This stores the maximum error for the three objects among the errors for each x,y,z directions in velocity
err_r = 0.
err_v = 0.
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

!Predicted position,velocity and acceleration stored in the arrays below
r(1:3,1:n,1) = 0.
v(1:3,1:n,1) = 0.
a(1:3,1:n,1) = 0.
!Position,velocity and acceleration arrays which store the corrected position,velocity and acceleration
r_c(1:3,1:n,1) = 0.
v_c(1:3,1:n,1) = 0.
a_c(1:3,1:n,1) = 0.
!Below store the new positions,velocities,accelerations produced when the timestep is halved.
r_newpast(1:3,1:n,-3:-1) = 0.
v_newpast(1:3,1:n,-3:-1) = 0.
a_newpast(1:3,1:n,-3:-1) = 0.
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
!write(6,*) 'A: INITIAL KE=',KE
!WRITE(6,*) 'B:INITIAL PE=',PE
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

!!!!!!!!!!!!!!shifting the t=0 step to t=-8 step
r(1:3,1:n,-8) = r(1:3,1:n,0)
v(1:3,1:n,-8) = v(1:3,1:n,0)
a(1:3,1:n,-8) = a(1:3,1:n,0)
!----------------------------------
!Bootstrapping the second order code- to get the past steps required for predictor-corrector method
!loop created which iterates over 'p'. This produces the past steps.
do p=-7,0
    !Future position obtained from the initial position,velocity and acceleration
    do i=1,n
        r(1:3,i,p) = r(1:3,i,p-1) + v(1:3,i,p-1)*dt + 0.5*a(1:3,i,p-1)*dt*dt
    end do
    !write(6,*) 'xy position of Earth:',r(1:2,3,p)
    !acceleration loop
    do i=1,n
        a(1:3,i,p) = 0.
        do j=1,n
            if(i==j) cycle
            do k=1,3
                dr(k) = r(k,i,p) - r(k,j,p)
            end do
            rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

            do k=1,3
                a(k,i,p) = a(k,i,p) - (G*m(j)*dr(k)/(sqrt(rad2)*rad2))
            end do
        end do
    end do

    !Future velocity
    do i=1,n
        v(1:3,i,p) = v(1:3,i,p-1) + 0.5*(a(1:3,i,p-1) + a(1:3,i,p))*dt
    end do
end do
!write(6,*) '----------------------------'
!write(6,*) 'past xy position of earth:',r(1:2,3,-8:0)/AU
!----------------------------------

!--------Fourth order Predictor-Corrector method----------
!Uses four past steps, for example, for t+1 position, t=0,-1,-2,-3 past positions are required
!---------Adams-Bashforth-Moulton Method---------
!OPEN(3,file='xypos_2.txt',status='new')
!open(5,file='encheck_sixobj_dt500.txt',status='new')
!open(3,file='dsun_1.txt',status='new')
do
    !x--------PREDICTOR---------x
    !Predicted position, r(t+1)
    do i=1,n
        r(1:3,i,1) = r(1:3,i,0) + dt*(-9.*v(1:3,i,-3) + 37.*v(1:3,i,-2) - 59.*v(1:3,i,-1) + 55.*v(1:3,i,0))/24. !coefficients are taken from the PREDICTOR-CORRECTOR METHOD PAPER
    !Predicted velocity, v(t+1)
        v(1:3,i,1) = v(1:3,i,0) + dt*(-9.*a(1:3,i,-3) + 37.*a(1:3,i,-2) - 59.*a(1:3,i,-1) + 55.*a(1:3,i,0))/24.
    end do
!    if(tcurr<0.1*tout) Then
!        write(6,*) '-----------------------'
!        write(6,*) 'predicted y velocity of Earth:',v(2,3,1)/1000.
!    end if


    !For the predicted acceleration
    dr(1:3) = 0.
    do i=1,n
        a(1:3,i,1) = 0.                  !Initial acceleration at t=0
        do j=1,n
            if(i==j) cycle                    !to avoid the case of finding separation between the object and itself
            do k=1,3
                dr(k) = r(k,i,1) - r(k,j,1)
            end do
            rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)                  !r^2 = x^2 + y^2 + z^2

            do k=1,3                                                        !runs over the all directions
                a(k,i,1) = a(k,i,1) - (G*m(j)*dr(k)/(sqrt(rad2)*rad2))    !dr(k) is the unit vector indicating the direction of the force from
                                                                            !objects 'j' to 'i'
            end do
        end do
    end do

    !x--------CORRECTOR---------x
    !Now, the predicted value at t+1 is used alongwith the t=0,-1,-2 positions for determining the corrected value at t+1.
    do i=1,n
        !Corrected position, r_c(t+1)
        r_c(1:3,i,1) = r(1:3,i,0) + dt*(v(1:3,i,-2) - 5.*v(1:3,i,-1) + 19.*v(1:3,i,0) + 9.*v(1:3,i,1))/24.
        !Corrected velocity, v_c(t+1)
        v_c(1:3,i,1) = v(1:3,i,0) + dt*(a(1:3,i,-2) - 5.*a(1:3,i,-1) + 19.*a(1:3,i,0) + 9.*a(1:3,i,1))/24.

    end do

!    textracount = textracount + dt
!    if(textracount==0.5)
!    write(6,*) '-------------------------'
!    write(6,*) 'y corrected velocity of earth:',v_c(2,3,1)/1000.
!
    !Corrected acceleration, a_c(t+1)
    dr(1:3) = 0.
    rad2 = 0.
    do i=1,n
        a_c(1:3,i,1) = 0.
        do j=1,n
            if(i==j) cycle                    !to avoid the case of finding separation between the object and itself
            do k=1,3
                dr(k) = r_c(k,i,1) - r_c(k,j,1)
            end do
            rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)                  !r^2 = x^2 + y^2 + z^2

            do k=1,3                                                        !runs over the all directions
                a_c(k,i,1) = a_c(k,i,1) - (G*m(j)*dr(k)/(sqrt(rad2)*rad2))    !dr(k) is the unit vector indicating the direction of the force from
                                                                            !objects 'j' to 'i'
            end do
        end do
    end do
    !write(6,*) 'Corrected acceleration of earth:',a_c(1:2,3,1)/AU
    tcurr = tcurr + dt
    textracount = textracount + dt !counter for noting the xy positions of the objects
    !write(3,*) tcurr/yr,r_c(1:2,1:n,1)/AU  !writing the xy position (corrected) for the objects to a file
    !write(6,*) 'textracount=',textracount/yr

!    if(tcurr<(500.*yr)) then !a reading every 6 months
!        write(6,*) 'xy pos',dt
!        !write(3,*) tcurr/yr,r_c(1:2,1:n,1)/AU  !writing the xy position (corrected) for the objects to a file
!        q = q+1
!        !textracount = 0.
!    end if
!
    !-------------------------------------------------------
    !DOWNSHIFTING arrays of position,velocity and acceleration
    !The final past step (t=-8) is shifted out of the array, and the corrected value at t=+1 is shifted to t=0.
    !Done for position, velocity and acceleration
    do p=-8,0
        do i=1,n
            if(p==0) Then
                r(1:3,i,p) = r_c(1:3,i,p+1)  !This IF condition is introduced so that the corrected value at t=+1 is shifted to t=0 in the main position array.
                v(1:3,i,p) = v_c(1:3,i,p+1)
                a(1:3,i,p) = a_c(1:3,i,p+1)
            end if
            r(1:3,i,p) = r(1:3,i,p+1)   !the past steps are downshifted so that the corrected value at t=+1 can occupy t=0 step
            v(1:3,i,p) = v(1:3,i,p+1)
            a(1:3,i,p) = a(1:3,i,p+1)
        end do
    end do

    !------------Fractional energy calculation---------------
    !Energy check
    !getting the current total energy
    tcount = tcount + dt
    if(tcount==(0.5*yr)) Then
    !write(6,*) 'true'
        KE = 0.
        do i=1,n
            KE = KE + 0.5*m(i)*(v_c(1,i,1)**2 + v_c(2,i,1)**2 + v_c(3,i,1)**2)
        end do

        PE = 0.
        rad = 0.
        dr(1:3) = 0.
        do i=1,n-1
            do j=i+1,n
                dr(1:3) = r_c(1:3,i,1) - r_c(1:3,j,1)
                rad = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                PE = PE - G*m(i)*m(j)/rad
            end do
        end do

        Ecurr = KE + PE
!        WRITE(6,*) '########'
!        WRITE(6,*) 'a: KE=',KE
!        WRITE(6,*) 'b: PE=',PE
        !fractional energy loss,delE
        delE = ((Ei - Ecurr))/Ei
!        write(6,*) '----------------------------'
!        write(6,*) 'initial energy=',Ei
!        write(6,*) 'Current energy=',Ecurr
!        write(6,*) 'delE=',delE
!        write(6,*) '--------------'
        !write(5,*) tcurr/yr,delE
        tcount = 0. !reinitialising to zero for the next calculation of delE
    end if

    !----------Measuring the distance from the sun-------------
    do i=1,n
        dsun(i) = sqrt(r_c(1,i,1)*r_c(1,i,1) + r_c(2,i,1)*r_c(2,i,1))
    end do
    !write(3,*) tcurr/yr,dsun(3)/AU
    if(tcurr<300.*yr) Then
        !write(6,*) 'true'
        write(6,*) 'orbital distance of earth=',dsun(3)/AU
        !tcount=0.
    end if

    !--------Adaptive timestep check------------
    !Changing the timestep based on the error found between the predicted and corrected values
    do i=1,n
        err_fact_r(1:3,i) = (19./270.)*(ABS(r_c(1:3,i,1) - r(1:3,i,1))/(ABS(r_c(1:3,i,1)) + small))
        err_fact_v(1:3,i) = (19./270.)*(ABS(v_c(1:3,i,1) - v(1:3,i,1))/(ABS(v_c(1:3,i,1)) + small))

        err_max_r(i) = maxval(err_fact_r(1:3,i)) !this finds the maximum error among the x,y,z components of the position and velocity for each object.
        err_max_v(i) = maxval(err_fact_v(1:3,i))
    end do

    !The maximum error among the objects is compared against the relative error.
    !Maximum error found for position and velocity values.
    err_r = maxval(err_max_r(1:n))
    err_v = maxval(err_max_v(1:n))  !Stores the maximum error among the all the objects
!    write(6,*) '------------------------------'
!    write(6,*) 'max(err_r,err_v)=',max(err_r,err_v)
    if(max(err_r,err_v)>relerr) Then  !finds the maximum error contribution from the position or velocity
        dt = dt*0.5
        !New past steps are introduced.

        do i=1,n
            !This is the t=-1/2 step.
            !r_newpast(1:3,1:n,-1) = (-5.*r(1:3,1:n,-4) + 28.*r(1:3,1:n,-3) - 70.*r(1:3,1:n,-2) + 140.*r(1:3,1:n,-1) + 35.*r(1:3,1:n,0))/128.
            v_newpast(1:3,i,-1) = (-5.*v(1:3,i,-4) + 28.*v(1:3,i,-3) - 70.*v(1:3,i,-2) + 140.*v(1:3,i,-1) + 35.*v(1:3,i,0))/128.
            a_newpast(1:3,i,-1) = (-5.*a(1:3,i,-4) + 28.*a(1:3,i,-3) - 70.*a(1:3,i,-2) + 140.*a(1:3,i,-1) + 35.*a(1:3,i,0))/128.

            !This is the t=-3/2 step
            !r_newpast(1:3,1:n,-3) = (3.*r(1:3,1:n,-4) - 20.*r(1:3,1:n,-3) + 90.*r(1:3,1:n,-2) + 60.*r(1:3,1:n,-1) - 5.*r(1:3,1:n,0))/128.
            v_newpast(1:3,i,-3) = (3.*v(1:3,i,-4) - 20.*v(1:3,i,-3) + 90.*v(1:3,i,-2) + 60.*v(1:3,i,-1) - 5.*v(1:3,i,0))/128.
            a_newpast(1:3,i,-3) = (3.*a(1:3,i,-4) - 20.*a(1:3,i,-3) + 90.*a(1:3,i,-2) + 60.*a(1:3,i,-1) - 5.*a(1:3,i,0))/128.
        end do
        !So the first four past steps taken in by the predictor will be: ???????????????DOUBT HERE
        r(1:3,1:n,-2) = r(1:3,1:n,-1)     !the -1 step becomes the new -2 paststep
        v(1:3,1:n,-2) = v(1:3,1:n,-1)
        a(1:3,1:n,-2) = a(1:3,1:n,-1)

        r(1:3,1:n,-1) = r_newpast(1:3,1:n,-1)  !the -1/2 step going into -1 step.
        v(1:3,1:n,-1) = v_newpast(1:3,1:n,-1)
        a(1:3,1:n,-1) = a_newpast(1:3,1:n,-1)

        r(1:3,1:n,-3) = r_newpast(1:3,1:n,-3)   !the -3/2 step going into the -3 paststep
        v(1:3,1:n,-3) = v_newpast(1:3,1:n,-3)
        a(1:3,1:n,-3) = a_newpast(1:3,1:n,-3)

        err_max_r(1:n) = 0.
        err_max_v(1:n) = 0.

!        do p=-3,-1
!            if(p==-2) Then
!                r(1:3,1:n,p) = r(1:3,1:n,p+1)
!            end if
!            r(1:3,1:n,p) = r_newpast(1:3,1:n,p)
!            v(1:3,1:n,p) = v_newpast(1:3,1:n,p)
!            a(1:3,1:n,p) = a_newpast(1:3,1:n,p)
!        end do

    else if(max(err_r,err_v)<(0.01*relerr)) then
        dt = dt*2.

        r(1:3,1:n,-1) = r(1:3,1:n,-2)   !the -2 step becomes the new -1 paststep
        v(1:3,1:n,-1) = v(1:3,1:n,-2)
        a(1:3,1:n,-1) = a(1:3,1:n,-2)

        r(1:3,1:n,-2) = r(1:3,1:n,-4)   !the -4 step becomes the new -2 paststep
        v(1:3,1:n,-2) = v(1:3,1:n,-4)
        a(1:3,1:n,-2) = a(1:3,1:n,-4)

        r(1:3,1:n,-3) = r(1:3,1:n,-6)   !the -6 step becomes the new -3 paststep
        v(1:3,1:n,-3) = v(1:3,1:n,-6)
        a(1:3,1:n,-3) = a(1:3,1:n,-6)

        err_max_r(1:n) = 0.
        err_max_v(1:n) = 0.
    end if

    !write(6,*) 'dt=',dt
    if(tcurr>tout) exit
end do
!CLOSE(3)
end program nbodyfour

        
        





