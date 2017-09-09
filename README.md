This repository contains primitive codes for solving certain boundary
value problems for Laplace's equation on a domain with a single corner
point.  The file dirichlet.f90 contains a code for solving the interior
Dirichlet problem

	  \Delta u =  0      in \Omega
	         u =  f	     on \partial\Omega,

where \Omega is a convex domain with a single corner point (so the
interior angle  is between 0 and pi).   The file neumann.f90 contains
a code for solving  the exterior Neumann problem


          \Delta u = 0       in \Omega^c
	   D_\nu u = f       on \partial\Omega,

where \Omega is a convex domain with a single corner point 
(so the interior angle is between is and pi) and D_\nu denotes the derivative with
respect to the outward-pointing unit normal.

The file neumann_concave.f90 contains a code for solving the 
Neumann boundary value problem on a concave domain with a single corner
point (so the interior angle is between pi and 2*pi).

The code for solving the Dirichlet problem on a can be compiled via the
following command:
    
    gfortran -o dirichlet.out -O3 dirichlet.f90 -llapack; ./dirichlet.out

The code for the Neumann problem on a convex domain can be compiled
via the following command:

    gfortran -o neumann.out -O3 neumann.f90 -llapack; ./neumann.out

The code for the Neumann problem on a concave domain can be compiled
and executed via the following command:

    gfortran -o neumann_concave.out -O3 neumann_concave.f90 -llapack; ./neumann_concave.out

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
The resulting linear systems are sufficiently small that they can be 
inverted with LAPACK routines for dense matrices (i.e., not
fast solvers are required to solve these problems).

It is well-known that the operators 1/2I + K and 1/2I + K^* are
isomorphisms which are compact perturbations coercive operators.
See, for instance,

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
linear system.  This is true for convex domains as well as concave
domains, provided proper discretization methods are used.

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