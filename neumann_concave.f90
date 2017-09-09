!
!  This file contains a primitive fortran code for solving the Neumann
!  boundary value problem
!
!        \Delta u        =  0      in \Omega^c
!             D_\nu u    =  f      on \partial\Omega,
!
!  where \Omega is a curve with a single corner point with a user-specified
!  interior angle > pi and D_\nu denotes differentiation with respect to the 
!  outward-pointing unit norm.
!
!   It operates by discretizing the integral equation
!
!    (1/2 I + K^* ) \sigma =  f,
!
!  where K is the adjoint of the double layer operator, and inverting the 
!  resulting linear system using an LAPACK routine.
!
!  When the ifscale flag is set, the method used is equivalent to a 
!  Galerkin method, the convergence of which is established (for instance) in
!
!      M. Costabel, "Boundary Integral Operators On Lipschitz Domains"
!      SIAM Journal on Mathematical Analysis, 19 (1988), pg. 613-626
!
!  When it is not set, the method used is a standard Nystrom method.  This results
!  in large condition numbers.
!

module neumann_subroutines

implicit double precision (a-h,o-z)

double precision :: theta,beta

data pi / 3.14159265358979323846264338327950288d0  /
contains


subroutine curve(t,x,y,dx,dy)
implicit double precision (a-h,o-z)

!
! This subroutine supplies a counter-clockwise parameterization
! of a curve with a single corner of angle theta located at the 
! origin.
!

if (t .lt. 0) then

x    = 2*sin(t*3/2)
y    = -beta*sin(t)

dx   = 2*cos(t*3/2)*3/2
dy   = -beta*cos(t)

else

x    = -2*sin(t*3/2)
y    = -beta*sin(t)

dx   = -2*cos(t*3/2)*3/2
dy   = -beta*cos(t)

endif

end subroutine



subroutine kernel(s,t,val)
implicit double precision (a-h,o-z)

!
!  Evaluate the double layer kernel; care is taken to ensure it is
!  accurately done when s is close to t and when s and t are both close
!  to the corner point
!


val = 0
if (s .gt. 0 .AND. t .lt. 0) then

if (s .lt. 1.0d-2 .AND. t .gt. -1.0d-2) then
val =  (beta*(0.32d4*s**3-0.104d4*s**5+0.6d2*s**7+s**4*t*(1940-0.39d3*t**2+  &
0.27d2*t**4)-0.6d1*s**2*t*(2080-480*t**2+39*t**4)+0.4d1*t*(3840-  &
1040*t**2+97*t**4)))/(0.512d4*Pi*((0.3d1*s-(9*s**3)/0.8d1+  &
(81*s**5)/0.64d3+0.3d1*t-(9*t**3)/0.8d1+(81*t**5)/0.64d3)**2+  &
(beta**2*(0.12d3*s-0.2d2*s**3+s**5-0.1d1*t*(120-20*t**2+  &
t**4))**2)/0.144d5))
return
endif

val = (beta*((6+0.4d1*Cos(s))*Sin(s/0.2d1)**3+0.3d1*Cos((3*s)/0.2d1)*Sin(t)+&
0.2d1*Cos(s)*Sin((3*t)/0.2d1)))/(0.2d1*Pi*(4-0.2d1*Cos(3*s)-  &
0.2d1*Cos(3*t)+beta**2*(Sin(s)-0.1d1*Sin(t))**2+  &
0.8d1*Sin((3*s)/0.2d1)*Sin((3*t)/0.2d1)))

return
endif

if (s .lt. 0 .AND. t .gt. 0) then


if (s .gt. -1.0d-2 .AND. t .lt. 1.0d-2) then
val =  -(beta*(0.32d4*s**3-0.104d4*s**5+0.6d2*s**7+s**4*t*(1940-0.39d3*t**2+&
0.27d2*t**4)-0.6d1*s**2*t*(2080-480*t**2+39*t**4)+0.4d1*t*(3840-  &
1040*t**2+97*t**4)))/(0.512d4*Pi*((0.3d1*s-(9*s**3)/0.8d1+  &
(81*s**5)/0.64d3+0.3d1*t-(9*t**3)/0.8d1+(81*t**5)/0.64d3)**2+  &
(beta**2*(0.12d3*s-0.2d2*s**3+s**5-0.1d1*t*(120-20*t**2+  &
t**4))**2)/0.144d5))
return
endif


val = (beta*(0.3d1*Cos((3*s)/0.2d1)*(Sin(s)-Sin(t))-  &
0.2d1*Cos(s)*(Sin((3*s)/0.2d1)+  &
Sin((3*t)/0.2d1))))/(0.2d1*Pi*(beta**2*Sin(s)**2+  &
0.4d1*Sin((3*s)/0.2d1)**2-0.2d1*beta**2*Sin(s)*Sin(t)+  &
beta**2*Sin(t)**2+0.8d1*Sin((3*s)/0.2d1)*Sin((3*t)/0.2d1)+  &
0.4d1*Sin((3*t)/0.2d1)**2))
return
endif

if (s .lt. 0 .AND. t .lt. 0) then

!
!  Both points close to the corner
!
if (s .gt. -1.0d-3 .AND. t .gt. -1.0d-3) then
delta = t-s
val = delta**3*((-423*beta)/(0.256d3*(9+beta**2)**2*Pi)-  &
(41*beta**3)/(0.768d3*(9+beta**2)**2*Pi))-(5*beta*delta)/(0.16d2*(9+  &
beta**2)*Pi)+(delta**2*((-2655*beta)/(0.256d3*(9+beta**2)**2*Pi)-  &
(95*beta**3)/(0.256d3*(9+beta**2)**2*Pi))-(15*beta)/(0.16d2*(9+  &
beta**2)*Pi))*s+((-5*beta*(531+0.19d2*beta**2)*delta)/(0.128d3*(9+  &
beta**2)**2*Pi)-(5*beta*(269487+0.13878d5*beta**2+  &
0.415d3*beta**4)*delta**3)/(0.6144d4*(9+beta**2)**3*Pi))*s**2+  &
((-5*beta*(369+0.11d2*beta**2))/(0.128d3*(9+beta**2)**2*Pi)-  &
(5*beta*(532737+0.29106d5*beta**2+  &
0.857d3*beta**4)*delta**2)/(0.6144d4*(9+beta**2)**3*Pi))*s**3
return
endif

if (abs(s-t) .lt. 1.0d-3) then
delta = t-s
val   = (-3*beta*(-0.2d1*Cos((3*s)/0.2d1)*Sin(s)+  &
0.3d1*Cos(s)*Sin((3*s)/0.2d1)))/(0.8d1*Pi*(beta**2*Cos(s)**2+  &
0.9d1*Cos((3*s)/0.2d1)**2))-  &
(beta*delta**3*(0.41d2*beta**6*Cos(s)**7*Cos((3*s)/0.2d1)+  &
0.2007d4*beta**4*Cos(s)**5*Cos((3*s)/0.2d1)**3+  &
0.26163d5*beta**2*Cos(s)**3*Cos((3*s)/0.2d1)**5+  &
0.102789d6*Cos(s)*Cos((3*s)/0.2d1)**7-  &
0.12d2*beta**6*Cos(s)**5*Cos((3*s)/0.2d1)*Sin(s)**2-  &
0.6696d4*beta**4*Cos(s)**3*Cos((3*s)/0.2d1)**3*Sin(s)**2-  &
0.59292d5*beta**2*Cos(s)*Cos((3*s)/0.2d1)**5*Sin(s)**2-  &
0.288d3*beta**6*Cos(s)**3*Cos((3*s)/0.2d1)*Sin(s)**4+  &
0.2592d4*beta**4*Cos(s)*Cos((3*s)/0.2d1)**3*Sin(s)**4+  &
0.198d3*beta**6*Cos(s)**6*Sin(s)*Sin((3*s)/0.2d1)+  &
0.16362d5*beta**4*Cos(s)**4*Cos((3*s)/0.2d1)**2*Sin(s)*Sin((3*s)/0.2d1)  +  &
0.71442d5*beta**2*Cos(s)**2*Cos((3*s)/0.2d1)**4*Sin(s)*Sin((3*s)/0.2d1)  &
-0.538002d6*Cos((3*s)/0.2d1)**6*Sin(s)*Sin((3*s)/0.2d1)+  &
0.432d3*beta**6*Cos(s)**4*Sin(s)**3*Sin((3*s)/0.2d1)-  &
0.23328d5*beta**4*Cos(s)**2*Cos((3*s)/0.2d1)**2*Sin(s)**3*Sin((3*s)/0.2d1)  &
+0.34992d5*beta**2*Cos((3*s)/0.2d1)**4*Sin(s)**3*Sin((3*s)/0.2d1)-  &
0.2187d4*beta**4*Cos(s)**5*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)**2+  &
0.91854d5*beta**2*Cos(s)**3*Cos((3*s)/0.2d1)**3*Sin((3*s)/0.2d1)**2+  &
0.1003833d7*Cos(s)*Cos((3*s)/0.2d1)**5*Sin((3*s)/0.2d1)**2+  &
0.34992d5*beta**4*Cos(s)**3*Cos((3*s)/0.2d1)*Sin(s)**2*Sin((3*s)/0.2d1)**2  -  &
0.314928d6*beta**2*Cos(s)*Cos((3*s)/0.2d1)**3*Sin(s)**2*Sin((3*s)/0.2d1)**2  &
-0.8748d4*beta**4*Cos(s)**4*Sin(s)*Sin((3*s)/0.2d1)**3+  &
0.472392d6*beta**2*Cos(s)**2*Cos((3*s)/0.2d1)**2*Sin(s)*Sin((3*s)/0.2d1)**3  &
-0.708588d6*Cos((3*s)/0.2d1)**4*Sin(s)*Sin((3*s)/0.2d1)**3-  &
0.118098d6*beta**2*Cos(s)**3*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)**4+  &
0.1062882d7*Cos(s)*Cos((3*s)/0.2d1)**3*Sin((3*s)/0.2d1)**4))/(0.768d3*Pi*(beta**2*Cos(s)**2 &
 +0.9d1*Cos((3*s)/0.2d1)**2)**4)+  &
(beta*delta*((-5*Cos(s)*Cos((3*s)/0.2d1))/(0.8d1*(beta**2*Cos(s)**2+  &
0.9d1*Cos((3*s)/0.2d1)**2))+(((3*Cos((3*s)/0.2d1)*Sin(s))/0.2d1-  &
(9*Cos(s)*Sin((3*s)/0.2d1))/0.4d1)*(0.2d1*beta**2*Cos(s)*Sin(s)+  &
0.27d2*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)))/(0.2d1*(beta**2*Cos(s)**2+ &
0.9d1*Cos((3*s)/0.2d1)**2)**2)))/(0.2d1*Pi)+  &
(beta*delta**2*((-(Cos((3*s)/0.2d1)*Sin(s))/0.8d1+  &
(27*Cos(s)*Sin((3*s)/0.2d1))/0.64d2)/(beta**2*Cos(s)**2+  &
0.9d1*Cos((3*s)/0.2d1)**2)-  &
(5*Cos(s)*Cos((3*s)/0.2d1)*(0.2d1*beta**2*Cos(s)*Sin(s)+  &
0.27d2*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)))/(0.16d2*(beta**2*Cos(s)**2+&
0.9d1*Cos((3*s)/0.2d1)**2)**2)+(((3*Cos((3*s)/0.2d1)*Sin(s))/0.2d1-  &
(9*Cos(s)*Sin((3*s)/0.2d1))/0.4d1)*(0.16d2*beta**4*Cos(s)**4+  &
0.468d3*beta**2*Cos(s)**2*Cos((3*s)/0.2d1)**2+  &
0.2916d4*Cos((3*s)/0.2d1)**4+0.36d2*beta**4*Cos(s)**2*Sin(s)**2-  &
0.108d3*beta**2*Cos((3*s)/0.2d1)**2*Sin(s)**2+  &
0.1296d4*beta**2*Cos(s)*Cos((3*s)/0.2d1)*Sin(s)*Sin((3*s)/0.2d1)-  &
0.243d3*beta**2*Cos(s)**2*Sin((3*s)/0.2d1)**2+  &
0.6561d4*Cos((3*s)/0.2d1)**2*Sin((3*s)/0.2d1)**2))/(0.48d2*(beta**2*Cos(s)**2  &
+0.9d1*Cos((3*s)/0.2d1)**2)**3)))/(0.2d1*Pi)
return
endif


val = (beta*(-0.2d1*(3+2*Cos(s))*Sin(s/0.2d1)**3-  &
0.3d1*Cos((3*s)/0.2d1)*Sin(t)+  &
0.2d1*Cos(s)*Sin((3*t)/0.2d1)))/(0.2d1*Pi*(4-0.2d1*Cos(3*s)-  &
0.2d1*Cos(3*t)+beta**2*(Sin(s)-0.1d1*Sin(t))**2-  &
0.8d1*Sin((3*s)/0.2d1)*Sin((3*t)/0.2d1)))
return

endif


! s > 0, t > 0

if (s .lt. 1.0d-3 .AND. t .lt. 1.0d-3) then
delta = t-s
val   = (15*beta*s)/(0.16d2*(9+beta**2)*Pi)+((1845*beta)/(0.128d3*(9+  &
beta**2)**2*Pi)+(55*beta**3)/(0.128d3*(9+beta**2)**2*Pi))*s**3+  &
delta*((5*beta)/(0.16d2*(9+beta**2)*Pi)+((2655*beta)/(0.128d3*(9+  &
beta**2)**2*Pi)+(95*beta**3)/(0.128d3*(9+beta**2)**2*Pi))*s**2)+  &
delta**3*((beta*(1269+0.41d2*beta**2))/(0.768d3*(9+beta**2)**2*Pi)+  &
(5*beta*(269487+0.13878d5*beta**2+0.415d3*beta**4)*s**2)/(0.6144d4*(9+&
beta**2)**3*Pi))+delta**2*((5*beta*(531+0.19d2*beta**2)*s)/(0.256d3*(9  &
+beta**2)**2*Pi)+(5*beta*(532737+0.29106d5*beta**2+  &
0.857d3*beta**4)*s**3)/(0.6144d4*(9+beta**2)**3*Pi))
return
endif

if (abs(s-t) .lt. 1.0d-3) then
delta = t-s
val   = (3*beta*(-0.2d1*Cos((3*s)/0.2d1)*Sin(s)+  &
0.3d1*Cos(s)*Sin((3*s)/0.2d1)))/(0.8d1*Pi*(beta**2*Cos(s)**2+  &
0.9d1*Cos((3*s)/0.2d1)**2))+  &
(beta*delta**3*(0.41d2*beta**6*Cos(s)**7*Cos((3*s)/0.2d1)+  &
0.2007d4*beta**4*Cos(s)**5*Cos((3*s)/0.2d1)**3+  &
0.26163d5*beta**2*Cos(s)**3*Cos((3*s)/0.2d1)**5+  &
0.102789d6*Cos(s)*Cos((3*s)/0.2d1)**7-  &
0.12d2*beta**6*Cos(s)**5*Cos((3*s)/0.2d1)*Sin(s)**2-  &
0.6696d4*beta**4*Cos(s)**3*Cos((3*s)/0.2d1)**3*Sin(s)**2-  &
0.59292d5*beta**2*Cos(s)*Cos((3*s)/0.2d1)**5*Sin(s)**2-  &
0.288d3*beta**6*Cos(s)**3*Cos((3*s)/0.2d1)*Sin(s)**4+  &
0.2592d4*beta**4*Cos(s)*Cos((3*s)/0.2d1)**3*Sin(s)**4+  &
0.198d3*beta**6*Cos(s)**6*Sin(s)*Sin((3*s)/0.2d1)+  &
0.16362d5*beta**4*Cos(s)**4*Cos((3*s)/0.2d1)**2*Sin(s)*Sin((3*s)/0.2d1)  +  &
0.71442d5*beta**2*Cos(s)**2*Cos((3*s)/0.2d1)**4*Sin(s)*Sin((3*s)/0.2d1)  &
-0.538002d6*Cos((3*s)/0.2d1)**6*Sin(s)*Sin((3*s)/0.2d1)+  &
0.432d3*beta**6*Cos(s)**4*Sin(s)**3*Sin((3*s)/0.2d1)-  &
0.23328d5*beta**4*Cos(s)**2*Cos((3*s)/0.2d1)**2*Sin(s)**3*Sin((3*s)/0.2d1)  &
+0.34992d5*beta**2*Cos((3*s)/0.2d1)**4*Sin(s)**3*Sin((3*s)/0.2d1)-  &
0.2187d4*beta**4*Cos(s)**5*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)**2+  &
0.91854d5*beta**2*Cos(s)**3*Cos((3*s)/0.2d1)**3*Sin((3*s)/0.2d1)**2+  &
0.1003833d7*Cos(s)*Cos((3*s)/0.2d1)**5*Sin((3*s)/0.2d1)**2+  &
0.34992d5*beta**4*Cos(s)**3*Cos((3*s)/0.2d1)*Sin(s)**2*Sin((3*s)/0.2d1)**2  -  &
0.314928d6*beta**2*Cos(s)*Cos((3*s)/0.2d1)**3*Sin(s)**2*Sin((3*s)/0.2d1)**2  &
-0.8748d4*beta**4*Cos(s)**4*Sin(s)*Sin((3*s)/0.2d1)**3+  &
0.472392d6*beta**2*Cos(s)**2*Cos((3*s)/0.2d1)**2*Sin(s)*Sin((3*s)/0.2d1)**3  &
-0.708588d6*Cos((3*s)/0.2d1)**4*Sin(s)*Sin((3*s)/0.2d1)**3-  &
0.118098d6*beta**2*Cos(s)**3*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)**4+  &
0.1062882d7*Cos(s)*Cos((3*s)/0.2d1)**3*Sin((3*s)/0.2d1)**4))/(0.768d3*Pi*(beta**2*Cos(s)**2  &
+0.9d1*Cos((3*s)/0.2d1)**2)**4)+  &
(beta*delta*((5*Cos(s)*Cos((3*s)/0.2d1))/(0.8d1*(beta**2*Cos(s)**2+  &
0.9d1*Cos((3*s)/0.2d1)**2))+(((-3*Cos((3*s)/0.2d1)*Sin(s))/0.2d1+  &
(9*Cos(s)*Sin((3*s)/0.2d1))/0.4d1)*(0.2d1*beta**2*Cos(s)*Sin(s)+  &
0.27d2*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)))/(0.2d1*(beta**2*Cos(s)**2+ &
0.9d1*Cos((3*s)/0.2d1)**2)**2)))/(0.2d1*Pi)+  &
(beta*delta**2*(((Cos((3*s)/0.2d1)*Sin(s))/0.8d1-  &
(27*Cos(s)*Sin((3*s)/0.2d1))/0.64d2)/(beta**2*Cos(s)**2+  &
0.9d1*Cos((3*s)/0.2d1)**2)+  &
(5*Cos(s)*Cos((3*s)/0.2d1)*(0.2d1*beta**2*Cos(s)*Sin(s)+  &
0.27d2*Cos((3*s)/0.2d1)*Sin((3*s)/0.2d1)))/(0.16d2*(beta**2*Cos(s)**2+&
0.9d1*Cos((3*s)/0.2d1)**2)**2)+(((-3*Cos((3*s)/0.2d1)*Sin(s))/0.2d1+  &
(9*Cos(s)*Sin((3*s)/0.2d1))/0.4d1)*(0.16d2*beta**4*Cos(s)**4+  &
0.468d3*beta**2*Cos(s)**2*Cos((3*s)/0.2d1)**2+  &
0.2916d4*Cos((3*s)/0.2d1)**4+0.36d2*beta**4*Cos(s)**2*Sin(s)**2-  &
0.108d3*beta**2*Cos((3*s)/0.2d1)**2*Sin(s)**2+  &
0.1296d4*beta**2*Cos(s)*Cos((3*s)/0.2d1)*Sin(s)*Sin((3*s)/0.2d1)-  &
0.243d3*beta**2*Cos(s)**2*Sin((3*s)/0.2d1)**2+  &
0.6561d4*Cos((3*s)/0.2d1)**2*Sin((3*s)/0.2d1)**2))/(0.48d2*(beta**2*Cos(s)**2  &
+0.9d1*Cos((3*s)/0.2d1)**2)**3)))/(0.2d1*Pi)
return
endif

val =(beta*(-0.3d1*Cos((3*s)/0.2d1)*(Sin(s)-Sin(t))+  &
0.2d1*Cos(s)*(Sin((3*s)/0.2d1)-  &
Sin((3*t)/0.2d1))))/(0.2d1*Pi*(beta**2*Sin(s)**2+  &
0.4d1*Sin((3*s)/0.2d1)**2-0.2d1*beta**2*Sin(s)*Sin(t)+  &
beta**2*Sin(t)**2-0.8d1*Sin((3*s)/0.2d1)*Sin((3*t)/0.2d1)+  &
0.4d1*Sin((3*t)/0.2d1)**2))
return

end subroutine


subroutine potential(x,y,val,derx,dery)
implicit double precision (a-h,o-z)

!
!  Evaluate a function f(x,y) which satsifies Laplace's equation in
!  the exterior of the domain.
!

x0   =  1.0d0
y0   =  0.0d0

dd   =  (x-x0)**2 + (y-y0)**2 
val  =  0.5d0 * log (dd )
derx = (x-x0)/dd
dery = (y-y0)/dd

end subroutine





subroutine legequad(n,xs,whts)
implicit double precision (a-h,o-z)

double precision, allocatable, intent(out) :: xs(:),whts(:)

!
!  Return the nodes and weights of the n-point Gauss-Legendre quadrature
!  rule on the interval [-1,1].
!

double precision :: pols(0:n)

maxiters = 10
allocate(xs(n),whts(n))

!
!   Use Newton's method and the recurrence relation to find half of the
!   roots --- the others are obtained via symmetry.
!
!   Note that we also store the value of the derivative at each of the obtained
!   roots for use in computing the weights below.
!

pi    = acos(-1.0d0)

ifodd = 0
nn = (n+1)/2
if (nn /= n/2) then
ifodd=1
nn=nn-1
endif

!
!  Use Tricomi's formula to  estimate the roots of P_{n+1}
!

do i =nn+1,n
dk    = i-nn-1
dn    = n
theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi
x0    = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)* &
        (39.0d0-28.0d0/sin(theta)**2)
xs(i) = x0*cos(theta)
enddo


!
!  Perform Newton iterations in order to refine the estimates.
!
do iter = 1,maxiters

!
!  Evaluate the Legendre polynomial of degree n at each point; save
!  the values of the derivative in the whts array.
!
do i=nn+1,n
call legepol(n,xs(i),pols(i),whts(i))
end do

!
!  Perform one Newton iteration
!
pols(nn+1:n) = pols(nn+1:n)/whts(nn+1:n)
xs(nn+1:n)   = xs(nn+1:n) - pols(nn+1:n)

!
!  Check to see if we are done
!

if ( sum(abs(pols(nn+1:n))) < eps0) then
exit
endif

end do

if (iter == maxiters)  then
print *,"legequad failed!"
stop
end if

!
! Compute the weights using the derivatives we stored above.
!
do j=nn+1,n
   x       = xs(j)
   dd      = 2.0d0/(1.0d0-x**2)
   whts(j) = dd/(whts(j)**2)
end do

!
! Reflect the quadrature nodes.
!
do j=1,nn
xs(j)   = -xs(n-j+1)
whts(j) = whts(n-j+1)
end do

!
! Handle the root at 0 if n is odd.
!

if (ifodd .eq. 1) then
x0          = 0
call legepol(n,x0,pol,der)
xs(nn+1)   = x0
whts(nn+1) = 2.0d0/(der**2)
endif

end subroutine



subroutine legepol(n,x,pol,der)
implicit double precision (a-h,o-z)

!
!  Return the value of the Legendre polynomial of degree n and its
!  derivative at a specified point in the interval [-1,1].
!
!  Input parameters:
!    n - the maximum degree of the polynomials to evaluate
!  
!  Output parameters:
!    pol - the value of polynomial
!    der - the value of its derivative

dimension pols0(0:n+1)

!
!  Handle the endpoints of the intervals
!

if ( abs(x-1) .lt. eps0) then
pol      = 1.0d0
der      = n*(n+1.0d0)/2.0d0
return
endif

if ( abs(x+1) .lt. eps0) then
pol      = (-1)**n 
der      = (-1)**(n+1.0d0)*(n+1.0d0)/2.0d0
return
endif

!
!  Compute the values of the polynomials of orders 0 through n+1 via the
!  recurrence relation
!

pols0(0) = 1.0d0
pols0(1) = x

do i = 1,n
pols0(i+1) =  (2*i+1.0d0) / (i+1.0d0) * x * pols0(i) - (i+0.0d0)/(i+1.0d0) * pols0(i-1)
end do

pol = pols0(n)

!
!  Evaluate the derivative
!

if (n .eq. 0) then
der = 0
elseif (n .eq. 1) then
der = 1.0d0
else
der = (n+1.0d0)/(1.0d0-x**2) * ( x*pols0(n) - pols0(n+1)   ) 
endif


end subroutine




subroutine singquad(nquad,xs,whts)
implicit double precision (a-h,o-z)
double precision, allocatable :: xs(:),whts(:)

!
!  Return a quadrature on the interval [-pi,pi] designed to integrate
!  functions with singularities at 0.
!

double precision, allocatable :: xslege(:),whtslege(:),ab(:,:)
double precision, allocatable :: ab0(:,:),xs0(:),whts0(:)

k      = 30
dsub   = 4.0d0
dd     = 2.0d0
rad    = 0.1d0
nints0 = 20
nints1 = 10

call legequad(k,xslege,whtslege)
allocate(ab0(2,nints0+nints1), xs0(k*nints0+k*nints1), whts0(k*nints0+k*nints1))

ab0   = 0
xs0   = 0
whts0 = 0

!
!  First, build a quadrature on [0,rad]
!

rad0 = rad**(1/dsub)
idx  = 0

do int=1,nints0
a = rad0*dd**(-nints0+int-1)
b = rad0*dd**(-nints0+int)

do i=1,k
x0   = xslege(i)   * (b-a)/2  + (b+a)/2
wht0 = whtslege(i) * (b-a)/2
x    = x0**dsub
wht  = x0**(dsub-1)*dsub*wht0

idx        = idx + 1
xs0(idx)   = x
whts0(idx) = wht
end do

end do

do int=1,nints1
a = rad + (1.0d0-rad)/(nints1+0.0d0) * (int-1)
b = rad + (1.0d0-rad)/(nints1+0.0d0) * (int)

do i=1,k
idx        = idx + 1
xs0(idx)   = (b-a)/2 * xslege(i) + (b+a)/2
whts0(idx) = (b-a)/2 * whtslege(i)
end do

end do


nquad = 2*idx
allocate(xs(nquad),whts(nquad))
do i=1,idx

xs(idx-i+1)   = -xs0(i)
whts(idx-i+1) = whts0(i)

xs(idx+i)     = xs0(i)
whts(idx+i)   = whts0(i)
end do

xs   = pi*xs
whts = pi*whts

end subroutine

end module


program neumann
use neumann_subroutines
implicit double precision (a-h,o-z)

double precision, allocatable :: xs(:),whts(:)
double precision, allocatable :: amatr(:,:),rhs(:),sol(:)

double precision, allocatable :: work(:),bmatr(:,:),sing(:)
integer, allocatable          :: ipiv(:)

double precision, allocatable :: xslege(:),whtslege(:)

!
!
!  theta   - angle of the corner
!
!  ifscale - set this to 1 in order to scale everything by square roots,
!            and thereby execute a method equivalent to Galerkin discretization


pi      = acos(-1.0d0)
theta   = 3*pi/2
beta    = tan((2*pi-theta)/2)
ifscale = 1

!
!  Construct a quadrature on the interval [-pi,pi] to use in the discretization
!  of the integral operator
!

call singquad(n,xs,whts)
print *,"# of quadrature nodes                   = ",n

!
!   Evaluate the entries of the matrix discretizing the operator 1/2I + K
! 

allocate(amatr(n,n))

do i=1,n
t    = xs(i)
twht = whts(i)
do j=1,n
s    = xs(j)
swht = whts(j)
call kernel(t,s,val)

if (ifscale .eq. 1) then
amatr(i,j) = val * sqrt(twht)*sqrt(swht)
else
amatr(i,j) = val * swht
endif

end do

amatr(i,i) = amatr(i,i) + 0.5d0
end do

!
!  Compute the condition number of the matrix using an LAPACK routine
!

allocate(bmatr(n,n),sing(n))
bmatr=amatr

lwork = -1
call dgesdd('N',n,n,bmatr,n,sing,u,n,vt,n,ww,lwork,iw,info)
lwork = ww

allocate(work(lwork))
call dgesdd('N',n,n,bmatr,n,sing,u,n,vt,n,work,lwork,iw,info)
if (info .ne. 0) then
print *,"ddegsdd failed with info = ",info
stop
endif

deallocate(work,bmatr)

! uncomment to display the singular values of the matrix
! print *,"singular values = "
! print "(4(D25.16,2X))",sing

print *,"CONDITION NUMBER                        = ",sing(1)/sing(n)

!
!  Form the right hand side f of the integral equation
!
!    1/2 sigma(x) + K[\sigma](x) = f(x)
!

allocate(rhs(n))
do i=1,n
t    = xs(i)
twht = whts(i)

call curve(t,x,y,dx,dy)
call potential(x,y,val,derx,dery)

val = (derx*dy - dery*dx)

if (ifscale .eq. 1) then
rhs(i) = val*sqrt(twht)
else
rhs(i) = val
endif

end do

!
!  Solve the integral equation for \sigma 
!

allocate(ipiv(n),sol(n),bmatr(n,n))
bmatr = amatr
sol   = rhs
nrhs  = 1
call dgesv(n,nrhs,bmatr,n,ipiv,sol,n,info)

if (info .ne. 0) then
print *,"dgesv failed with info = ",info
stop
endif

deallocate(bmatr)

!
!  Evaluate the resulting potential at the point (3,1) and compare
!  the result to the known value of the potential
!
!

x0   = 3.0d0
y0   = 1.0d0
dsum = 0.0d0

do i=1,n

t    = xs(i)
twht = whts(i)

call curve(t,x,y,dx,dy)

dd   = (x-x0)**2 + (y-y0)**2
val  = 0.5d0*log(dd)

if (ifscale .eq. 1) then
dsum = dsum + val*sol(i)*sqrt(twht)
else
dsum = dsum + val*sol(i)*twht
endif

end do

dsum = dsum / (2*pi)
call potential(x0,y0,dsum0,derx,dery)

print *,"true value of potential at (3,1)        = ",dsum
print *,"computed value of potential at (3,1)    = ",dsum0
print *,"relative error                          = ",abs(dsum-dsum0)/abs(dsum0)

end program
