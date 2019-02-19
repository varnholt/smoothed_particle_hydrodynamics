# Smoothed Particle Hydrodynamics (SPH)
by Matthias Varnholt a.k.a mueslee^hjb in 2015/2016

This code and my (most likely incorrect) notes below are based on a
webinar held by Alan Heirich


## Why SPH?

We have a bunch of particles we'd like to animate just like real water
would move like. The more physically correct the simulation is the more
convincing the effect looks like. Navier and Stokes developed equations
that serve this very purpose. SPH gives approximations for every term
in the Navier-Stokes equations, i.e. SPH helps us doing the actual math.


## Terms

Examples for fluids: Gasses (e.g. air), Liquids (e.g. water), aerosols  
Examples for liquids: Water, everything that's (mostly) incompressible  
Examples for aerosols:  Smoke  


## Incompressible Navier-Stokes

Mr. Navier and Mr. Stokes developed equations for the motion of fluids.
The 'Incompressible Navier-Stokes' variant can be used to compute the
motion of fluids that are not compressible - like liquids, or water in
particular.

The Incompressible Navier-Stokes equation consist of two elements:
- The equation of motion
- The equation of mass continuity

### The 'equation of motion'
The equation of motion says that changes in velocity are described as
a function of

- gravity g
- (a gradient of) pressure ∇p:  
   fluid flows from high pressure to low pressure regions
- velocity u∇^2v
- viscosity u:
   the fluid stickiness, low: air water, high: honey, mud

I.e. the fluid moves
- under the force of gravity
- under the force of pressure
- under the force of the 'viscous term'

### The 'mass continuity equation'

```      
         rho     ( ∇ dot v ) = 0
         ^         ^
         density   divergence
                   of velocity
```
The equation of mass continuity expresses that mass is neither created
nor destroyed.

For the incompressible equation the density is assumed to be constant.
Therefore for our purpose it does not need to be considered. All the
mass we will start with, we will still have at the end.

Since the mass continuity equations just says that we must not create
nor destroy particles while we're computing (i.e. we are assuming a
constant number of particles), this term does not need to be solved.

So back to the equation of motion.

It looks like this:
```      
         rho     * [dv/dt       + v dot ∇v]      = rho * g -  ∇p +    u∇^2v
```
|rho:|density of the fluid, a scalar variable  
|p:|pressure of the fluid, also a scalar variable  
|g:|gravity, a 3D vector  
|v:|velocity, also a 3D vector  

In words:
```

                    the           the convective
         density * [derivative  + acceleration  ] = gravity - pressure + viscosity
                    of velocity   term                        gradient   term
```
```
      The 'convective acceleration term':

          v dot ∇v = [v  * dv  / dx; v  * dv   / dy; v  * dv  / dz]
                       x     x        y     y         z     z

                           ^              ^               ^
                           partial derivatives of v with respect
                           to x, y and z
```
If fluids have to move through an area of limited space, their
velocity increases. That's because the same volume of fluid is moved
from a larger region to an area with less space. A good example is
water flowing through a hose.

```
      The 'pressure gradient':

          - ∇p = [dp / dx; dp / dy; dp / dz]

                  ^
                  partial derivatives of p with respect
                  to x, y and z
```
A gradient can be understood as a slope in a higher dimension (like
a hillside). The sign of the pressure gradient function is negative.
This is because the system moves from a region of high pressure to a
region of low pressure.

The pressure p is defined as follows:
```
         p = k       ( rho              - rho0 )
             ^         ^                  ^
             some      actual             resting density
             scalar    density of        (the density of the fluid
                       the fluid at       at equilibrium)
                       a given point
```
rho > rho0 => positive pressure
rho < rho0 => negative pressure (suction)

```
      The 'viscosity term':

         u * ∇^2 * v = [ ∇^2 v , ∇^2v , ∇^2 v  ]
                              x      y       z

         ∇^2 is also called the 'Laplacian operator' or the 'diffusion
         operator'. It diffuses whatever quantity it is applied to.
         => It is diffusing velocity.
         => It is diffusing the momentum of the system.
         Diffusing the velocity of the system will result in a system having
         identical velocities at every position sooner or later.
```


### The material derivative

We want to represent the dynamics of the system by a set of
particles and we do that by taking the material derivative. The
material derivative is the derivative along a path with a given
velocity v.
```
            rho * Dv/Dt = rho * g - ∇ p + u ∇^2 v
```

Now we have a simple equation for the motion of a single particle i:
```
            dv
              i        1         u    
            --- = g -  ---- ∇p + ---- ∇^2 v
            dt         rho       rho
                          i         i
```
The derivative of the velocity of particle i with respect to time is
gravity - (the pressure gradient / the density rho) plus (the
viscosity term devided by our density rho).

These are the equations we will actually solve.



## From Navier Stokes to Smoothed Particle Hydrodynamics
We can represent any quantity by the summation of nearby points
multiplied with a weighting function W. W is also called 'smoothing
kernel'. They've been introduced by Monaghan in '92 - *sigh*.

W will give more strength to points that are close to us
- Points that are further away from us have a weaker influence on us.
- For points away more than a certain distance, W will return 0,
   i.e. they don't affect us at all.  
   How far a quantity must be away before it stops interacting with us
   is called the 'interaction radius'

This is evaluating a quantity by sampling a neighborhood of space
and weighting points by how close they are to our sampling points.

Monaghan also introduced approximations to the terms of the
incompressible Navier-Stokes equation.

### Smoothed Particle Hydrodynamics
Here are the approximations for the individual terms of the
incompressible Navier-Stokes equations:

      '~=' : 'approximately equal'

1)
```
         rho  ~= SUM m  W(r - r , h)
            i         j        j
```

I.e. the density is approximately equal to the sum of the masses of
nearby points, weighted by the smoothing kernel W. Density is the
amount of mass at a given point, so this should be fairly reasonable.
So we approximate the density at point i by the summation of
neighboring points j (weighted appropriately).


2)
```
         ∇p               p        p
           i               i        j       
         ----  ~= SUM m ( ------ + ------ )  ∇W (r - r , h)
         rho           j   rho^2    rho^2             j
            i                   i        j
```
The pressure gradient divided by rho i is approximately equal to the
sums of the masses at various points j multiplied by a scalar quantity
of pressure over density. The term is multiplied by the gradient of a
smoothing kernel (∇W) which is a vector expression (since a gradient
is a vector). Therefore the result of this whole expression is a
vector value.

3)
```
                                       v  - v
         u               u              j    i  
         ---- ∇^2 v   ~= ---- SUM m  ( ------- ) ∇^2 W (r - r , h)
         rho       i     rho       j   rho                   j
            i               i             j
```

u is a scalar coefficient that says how viscous is the fluid. If you
choose a small value for u it behaves like water, for larger values
of u its behavior is more like sirup.

This viscosity term is approximately u divided by rho multiplied
with the sum of the masses of several points j nearby multiplied by
the difference in velocity between two points (i and j), divided by
the density of j. Finally the term is multiplied by the Laplacian of
the smoothing kernel ∇^2W. If vj and vi are equal, there's no
viscous interaction between them, if the difference between them is
small, there's a small viscous interaction and for large differences
between vi and vj, there will be a large viscous interaction.

Over time this term has the effect of encouraging particles to travel
together (to move in the same direction).

### Smoothing kernels:

Attributes of smoothing kernels:
- At a distance h, w will drop to 0. I.e. particles that are too far
   away will not interact with the particle currently processed.
- Over a sphere of radius h, w sums to 1.

```
         ||r - r  ||^2: the distance between two points, squared
                b
```
There's not much intuition about these :)

1)
```
                   315
         W =  ------------- ( h^2 - || r - r ||^2 )^3
              64 * PI * h^9                 b
```

2)

```         
                                                    r - r
                -45                                      b
          ∇W =  -------- * (h - ||r - r  ||)^2 * ------------
                PI * h^6               b         || r - r  ||
                                                         b
```
3)

```
                             45
          ∇^2 W(r - r , h) = -------- * (h - ||r - r  ||)
                     b       PI * h^6               b
```

   Algorithm:
```
      dv
        i        1         u    
      --- = g -  ---- ∇p + ---- ∇^2 v
      dt         rho       rho
                    i         i
```
Based on the material derivative, we do the following:

1) Compute an approximation for rho for every particle
   using equation 1.

2) Evaluate equation 2, the pressure gradient. It depends on rho,
   so rho needs to be computed first. It also depends on the pressure
   gradient which is calculated by computing the difference between
   rho and rho0, the resting potential (see above).

3) Evaluate equation 3, the viscous term. It depends on rho and on
   the current velocity of the particle.

4) After having computed all the approximations, we can put them
   together in the material derivative formula. Using a simple numerical
   scheme we can timestep the velocity first and then compute the new
   position of the particle.

In short:

   For each particle
   1) compute the density
   2) compute the pressure from density
   3) compute the pressure gradient
   4) calc the viscosity term
   5) timestep the velocity and calculate the next position


### Optimizations:
Without optimization we'd have to test every particle against every
particle (O(n^2)). The result will still be correct because the smoothing
kernel W will give 0 for as result for all particles that are outside the
interaction radius.

1) Voxels
   Device space into local regions of space (voxels). Those voxels should
   have the size 2h on a side (twice the interactivity radius; quantities
   outside this radius are not interacting with the particle).

   So a particle can only interact with
   - particles in the same voxel and in
   - adjacent voxels

   If a particle was exactly in the center of a voxel, it could only
   interact with particles in the same voxel; if it's not in the perfect
   center, it could interact with particles of 2x2x2 voxels.

2) Subsets
   Choose a random subset of n particles (for example 32) because in
   most cases it's not required to take all the particles into account.  
   Those 32 are good enough.
