program nbody
IMPLICIT NONE
!three body system of SUN-EARTH-JUPITER"
!Simulating the motions of the objects"
double precision :: r(1:3,1:3), v(1:3,1:3), a(1:3,1:3), a0, ai(1:3,1:3), da(1:3,1:3), m(1:3), dt, tcurr, G, Dse, Dsj, Dej, dr(1:3)
!r,v,a are the positions, velocities and accelerations of the bodies.
!dt and tcurr are the timestep and current time in the simulations."
!G is the gravitational constant."
double precision :: rad2, AU
integer :: i,j,k,n
!i,j and k are the indexes for iteration."

!Mass of SUN-JUPITER-EARTH (in Kg)
m(1:3) = (/1.989e30, 1.898e27, 5.972e24/)
!mj = 1.898*10**27
!ms = 1.989*10**30
!me = 5.972*10**24 !approximately, Ms = 1000Mj = 1000000Me 
!Astronomical unit
AU = 1.5e11 !(in Km)
!Orbital distances of SUN-JUPITER-EARTH (in m)
r(1,1:3) = (/0., 5.2*1.5e11, 1.0*1.5e11/)
!Orbital velocities of SUN-JUPITER-EARTH in the solar system (in m/s)
v(2,1:3) = (/0., 13.7*1000 ,30.*1000/)
!Gravitational constant (in N m^2 Kg^-2)
G = 6.67e-11 !try to set G=1 to avoid large exponentials
!timestep, dt (in sec)
dt = 1000.0
!number of objects
n=3

!setting the initial conditions
!initial velocity and position
v(1,1:n) = 0. 
v(3,1:n) = 0. !allowing the object to move in circular paths on the XY plane.
r(2:3,1:n) = 0. !forcing the objects to remain in the XY plane
dr(1:3) = 0. !sepaaration between the objects in the three directions


!Finding acceleration of the SUN-JUPITER-EARTH
!Need to store the initial acceleration for calculating the new velocity
!ai = 0.
a = 0.
do i=1,n
    a(1:3,i) = 0.
    do j=1,n
        if(i==j) cycle
        do k=1,3
            dr(k) = r(k,i) - r(k,j)
        end do
        rad2 = dr(1)**2 + dr(2)**2 + dr(3)**2
        
        do k=1,3
            a(k,i) = a(k,i) - G*m(j)*(dr(k)*AU)/(sqrt(rad2)*rad2*(AU**3)) !dr(k) represents the direction of the position vector
            !get the new acceleration 
        end do
    end do
end do
!WRITE (6,*) 'sqrt(rad2) vector', sqrt(rad2)
da = a - ai !'a' being the final acceleration
ai = a
tcurr = 0. !current time
!determining the position and velocity using the second order expansion
!1st iteration
do i=1,n
    r(1:3,i) = r(1:3,i) + v(1:3,i)*dt + 0.5*a(1:3,i)*dt**2
    v(1:3,i) = v(1:3,i) + a(1:3,i)*dt + 0.5*da(1:3,i)*dt
end do

!WRITE (6,*) ' final positions of S-J-E:\n', r(1,1:3)
!WRITE (6,*) ' Orbital velocity of earth: ', v(2,3)
!WRITE (6,*) 'ACCELERATION:', a(3,1:3)
end program nbody
 

            
            
            