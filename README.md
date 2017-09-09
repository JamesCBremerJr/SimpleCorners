This repository contains primitive codes for solving certain boundary
value problems for Laplace's equation on a domain with a single corner
point.  The file dirichlet.f90 contains a code for solving the interior
Dirichlet problem

	  \Delta u =  0      in \Omega
	         u =  f	     on \partial\Omega,

where \Omega is a domain with a single corner point with an angle
specified by the user.  The file neumann.f90 contains a code
for solving the exterior Neumann problem


          \Delta u = 0       in \Omega^c
	   D_\nu u = f       on \partial\Omega,

where \Omega is a domain with a single corner point with angle
specified by the user and D_\nu denotes the derivative with
respect to the outward-pointing unit normal.

The code for solving the Dirichlet problem can be compiled via the
following command:
    
    gfortran -O3 dirichlet.f90 -llapack

The code for the Neumann problem can be compiled via the following
command:

    gfortran -O3 neumann.f90 -llapack

The Dirichlet problem is equivalent to the integral equation

    (1/2 I + K ) \sigma = f,                                   (1)

where K is the double layer operator


     K[\sigma](x) = 1/(2\pi) \int   (x-y) \cdot \nu(y) / |x-y|^2 \sigma(y) dS(y),
                               \partial\Omega

and the Neumann problem corresponds to the integral equation

	(1/2 I + K^*) \sigma = f.                              (2)

A weighted Nystrom scheme which is equivalent to a standard
Galerkin discretization method is used to discretize (1) and (2).
It differs from standard Nystrom schemes only in that the 
matrices discretizing (1) and (2) are conjugated with a diagonal
matrix whose entries consist of square roots of quadrature weights.
The resulting linear system is sufficiently small that it can be 
inverted with an LAPACK routine.

It is well-known that the operators 1/2I + K and 1/2I + K* are
isomorphisms, each of  which is a compact perturbation of a coercive
operator.  See, for instance,

      M. Costabel, "Boundary Integral Operators On Lipschitz Domains"
      SIAM Journal on Mathematical Analysis, 19 (1988), pg. 613-626.

That Galerkin methods for the discretization of (1) and
(2) are effective follow from this observation.   Indeed, such discretizations 
are "quasi-optimal," meaning that the error in the obtained solution
is proportional to the error with which the discretization scheme 
approximates the true solution of the equation.

There seems to be a great deal of confusion regarding these results,
particularly the fact that the convergence of such methods was
established decades ago.  The only issue is the efficiency
with which the solutions are represented, not the convergence
of discretizations or the conditioning of the resulting
linear system.

Arbitrarily high order algebraic convergence can be easily obtained
using techniques only slightly more complicated than those used in
the primitive codes of this pacakage.  See, for instance,

  J. Bremer, V. Rokhlin and I. Sammis , "Universal quadratures for boundary 
  integral equations on two-dimensional domains with corners."
  Journal of Computational Physics, 229 (2010), pp. 8259-8280

in which a procedure for the construction of efficient  methods for 
the representation of  solutions of boundary value problems for
elliptic partial differential equation on domains with corner
points is described.