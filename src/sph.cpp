// base
#include "sph.h"

// sph
#include "particle.h"

// Qt
#include <QTime>

// cmath
#include <math.h>

// openmp
#include <omp.h>



/*

   Smoothed particle hydrodynamics (SPH),
   by Matthias Varnholt a.k.a mueslee^hjb in 2015/2016

   This code and my (most likely incorrect) notes below are based on a
   webinar held by Alan Heirich


   Why SPH?

      We have a bunch of particles we'd like to animate just like real water
      would move like. The more physically correct the simulation is the more
      convincing the effect looks like. Navier and Stokes developed equations
      that serve this very purpose. SPH gives approximations for every term
      in the Navier-Stokes equations, i.e. SPH helps us doing the actual math.


   Terms

      Examples for fluids: Gasses (e.g. air), Liquids (e.g. water), aerosols
      Examples for liquids: Water, everything that's (mostly) incompressible
      Examples for aerosols: Smoke


   Incompressible Navier-Stokes

      Mr. Navier and Mr. Stokes developed equations for the motion of fluids.
      The 'Incompressible Navier-Stokes' variant can be used to compute the
      motion of fluids that are not compressible - like liquids, or water in
      particular.

      The Incompressible Navier-Stokes equation consist of two elements:
      - The equation of motion
      - The equation of mass continuity


      The 'equation of motion'

         The equation of motion says that changes in velocity are described as
         a function of

            - gravity g
                                       __
            - (a gradient of) pressure \/p:
              fluid flows from high pressure to low pressure regions
                        __
            - velocity u\/^2v

            - viscosity u:
              the fluid stickiness, low: air water, high: honey, mud

         I.e. the fluid moves
         - under the force of gravity
         - under the force of pressure
         - under the force of the 'viscous term'


      The 'mass continuity equation'
                   __
         rho     ( \/ dot v ) = 0
         ^         ^
         density   divergence
                   of velocity

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
                                        __                     __        __
         rho     * [dv/dt       + v dot \/v]      = rho * g -  \/p +    u\/^2v

         rho: density of the fluid, a scalar variable
         p:   pressure of the fluid, also a scalar variable
         g:   gravity, a 3D vector
         v:   velocity, also a 3D vector

      In words:

                    the           the convective
         density * [derivative  + acceleration  ] = gravity - pressure + viscosity
                    of velocity   term                        gradient   term


      The 'convective acceleration term':
               __
         v dot \/v = [v  * dv  / dx; v  * dv   / dy; v  * dv  / dz]
                       x     x        y     y         z     z

                           ^              ^               ^
                           partial derivatives of v with respect
                           to x, y and z

         If fluids have to move through an area of limited space, their
         velocity increases. That's because the same volume of fluid is moved
         from a larger region to an area with less space. A good example is
         water flowing through a hose.


      The 'pressure gradient':
           __
         - \/p = [dp / dx; dp / dy; dp / dz]

                  ^
                  partial derivatives of p with respect
                  to x, y and z

         A gradient can be understood as a slope in a higher dimension (like
         a hillside). The sign of the pressure gradient function is negative.
         This is because the system moves from a region of high pressure to a
         region of low pressure.

         The pressure p is defined as follows:

         p = k       ( rho              - rho0 )
             ^         ^                  ^
             some      actual             resting density
             scalar    density of        (the density of the fluid
                       the fluid at       at equilibrium)
                       a given point

         rho > rho0 => positive pressure
         rho < rho0 => negative pressure (suction)


      The 'viscosity term':
             __           __       __      __
         u * \/^2 * v = [ \/^2 v , \/^2v , \/^2 v  ]
                                x       y        z
         __
         \/^2 is also called the 'Laplacian operator' or the 'diffusion
         operator'. It diffuses whatever quantity it is applied to.
         => It is diffusing velocity.
         => It is diffusing the momentum of the system.
         Diffusing the velocity of the system will result in a system having
         identical velocities at every position sooner or later.



      The material derivative

         We want to represent the dynamics of the system by a set of
         particles and we do that by taking the material derivative. The
         material derivative is the derivative along a path with a given
         velocity v.

                                    __       __
            rho * Dv/Dt = rho * g - \/ p + u \/^2 v


         Now we have a simple equation for the motion of a single particle i:

            dv
              i        1    __    u    __
            --- = g -  ---- \/p + ---- \/^2 v
            dt         rho        rho
                          i          i

         The derivative of the velocity of particle i with respect to time is
         gravity - (the pressure gradient / the density rho) plus (the
         viscosity term devided by our density rho).

         These are the equations we will actually solve.



   From Navier Stokes to Smoothed Particle Hydrodynamics

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



   Smoothed Particle Hydrodynamics

      Here are the approximations for the individual terms of the
      incompressible Navier-Stokes equations:

      '~=' : 'approximately equal'

      1)


         rho  ~= SUM m  W(r - r , h)
            i         j        j


         I.e. the density is approximately equal to the sum of the masses of
         nearby points, weighted by the smoothing kernel W. Density is the
         amount of mass at a given point, so this should be fairly reasonable.
         So we approximate the density at point i by the summation of
         neighboring points j (weighted appropriately).


      2)

         __
         \/p              p        p
            i              i        j        __
         ----  ~= SUM m ( ------ + ------ )  \/ W (r - r , h)
         rho           j   rho^2    rho^2                j
            i                   i        j

         The pressure gradient divided by rho i is approximately equal to the
         sums of the masses at various points j multiplied by a scalar quantity
         of pressure over density. The term is multiplied by the gradient of a
         smoothing kernel (\/W) which is a vector expression (since a gradient
         is a vector). Therefore the result of this whole expression is a
         vector value.


      3)
                                        v  - v
         u    __          u              j    i        __
         ---- \/^2 v   ~= ---- SUM m  ( ------------ ) \/^2 W (r - r , h)
         rho        i     rho       j   rho                         j
            i                i             j

         u is a scalar coefficient that says how viscous is the fluid. If you
         choose a small value for u it behaves like water, for larger values
         of u its behavior is more like sirup.

         This viscosity term is approximately u divided by rho multiplied
         with the sum of the masses of several points j nearby multiplied by
         the difference in velocity between two points (i and j), divided by
         the density of j. Finally the term is multiplied by the Laplacian of
         the smoothing kernel \/^2W. If vj and vi are equal, there's no
         viscous interaction between them, if the difference between them is
         small, there's a small viscous interaction and for large differences
         between vi and vj, there will be a large viscous interaction.

         Over time this term has the effect of encouraging particles to travel
         together (to move in the same direction).



      Smoothing kernels:

         Attributes of smoothing kernels:
         - At a distance h, w will drop to 0. I.e. particles that are too far
           away will not interact with the particle currently processed.
         - Over a sphere of radius h, w sums to 1.


         ||r - r  ||^2: the distance between two points, squared
                b

         There's not much intiution about these :)

         1)

                   315
         W =  ------------- ( h^2 - || r - r ||^2 )^3
              64 * PI * h^9                 b


         2)
                                                    r - r
         __     -45                                      b
         \/W =  -------- * (h - ||r - r  ||)^2 * ------------
                PI * h^6               b         || r - r  ||
                                                         b

         3)

         __                  45
         \/^2 W(r - r , h) = -------- * (h - ||r - r  ||)
                     b       PI * h^6               b




   Algorithm:

      dv
        i        1    __    u    __
      --- = g -  ---- \/p + ---- \/^2 v
      dt         rho        rho
                    i          i

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


   Optimizations:

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



   ----------------------------------------------------------------------------

   temp notes:

      g: gravity
      v: velocity

      m: mass
      V: volume
      density (rho) = m/V

      density of water: 999,97 kg/m³

      number density (n): concentration of particles in our volume
                      n = N (number of particles) / V


   TODO:
   [x] inverse of rhoi and rhoj can be precomputed in one loop
   [x] scale!
   [x] check computePressure, deltaRho seems to be wrong
   [x] change grid size to somewhat like: 32, 16, 16
   [x] mass = 0.001 does not work, 1.0 is better.
       find a suitable value for mass! if rho0 is 1000 is 1000 a good value?
   [ ] missing optimization
       // the distance between two particles ought to be passed to computeDensity
       Particle* SPH::evaluateNeighbor(
   [ ] get rid of sqrt in CFL condition check
   [ ] inverse rho computation in computeDensity does not make sense

*/




SPH::SPH()
 : mParticleCount(0),
   mGridCellCount(0),
   mRho0(0.0f),
   mStopped(false),
   mPaused(false)
{
   // grid
   float h = 3.34f;
   mSimulationScale = 0.004f;
   mSimulationScaleInverse = 1.0f / mSimulationScale;
   mH = h;
   mH2 = pow(h, 2);
   mHTimes2 = h * 2.0f;
   mHTimes2Inv = 1.0f / mHTimes2;
   mHScaled  = h * mSimulationScale;
   mHScaled2 = pow(h * mSimulationScale, 2);
   mHScaled6 = pow(h * mSimulationScale, 6);
   mHScaled9 = pow(h * mSimulationScale, 9);
   mParticleCount = 64 * 1024;
   mGridCellsX = 32;
   mGridCellsY = 32;
   mGridCellsZ = 32;
   mGridCellCount = mGridCellsX * mGridCellsY * mGridCellsZ;
   mCellSize = 2.0f * h;
   mMaxX = mCellSize * mGridCellsX;
   mMaxY = mCellSize * mGridCellsY;
   mMaxZ = mCellSize * mGridCellsZ;

   // physics
   mRho0 = 1000.0f;
   mStiffness = 0.75f;
   mGravity = vec3(0.0f, -9.8f, 0.0f);
   mViscosityScalar = 10.0f;
   mTimeStep = 0.0042f;
   mDamping = 0.75f;

   // float x = (1000.0f / (float)mParticleCount) * 2.0f;
   // printf("%f\n", x);
   float mass = 0.001f;
   mCflLimit = 100.0f;
   mCflLimit2 = mCflLimit * mCflLimit;

   // smoothing kernels
   mKernel1Scaled = 315.0f / (64.0f * M_PI * mHScaled9);
   mKernel2Scaled = -45.0f / (M_PI * mHScaled6);
   mKernel3Scaled = -mKernel2Scaled;

   // we do not examine a particle against more than 32 other particles
   mExamineCount = 8; // 8?

   mSrcParticles = new Particle[mParticleCount];
   mVoxelIds= new int[mParticleCount];
   mVoxelCoords= new vec3i[mParticleCount];

   for (int i = 0; i < mParticleCount; i++)
   {
      mSrcParticles[i].mMass = mass;
   }

   mGrid = new QList<Particle*>[mGridCellCount];

   mNeighbors = new Particle*[mParticleCount*mExamineCount];
   mNeighborDistancesScaled = new float[mParticleCount*mExamineCount];

   // randomize particle start positions
   // initParticlePositionsRandom();
   initParticlePolitionsSphere();
}


SPH::~SPH()
{
   stopSimulation();
   quit();
   wait();
}


bool SPH::isStopped() const
{
   mMutex.lock();
   bool stopped = mStopped;
   mMutex.unlock();

   return stopped;
}


bool SPH::isPaused() const
{
   mMutex.lock();
   bool paused = mPaused;
   mMutex.unlock();

   return paused;
}



void SPH::run()
{
   while(!isStopped())
   {
      if (!isPaused())
      {
         step();
      }
   }
}


void SPH::step()
{
   int timeVoxelize = 0;
   int timeFindNeighbors = 0;
   int timeComputeDensity = 0;
   int timeComputePressure = 0;
   int timeComputeAcceleration = 0;
   int timeIntegrate = 0;
   QTime t;

   // put particles into voxel grid
   t.start();
   voxelizeParticles();
   timeVoxelize = t.elapsed();

   // find neighboring particles
   t.restart();
   #pragma omp parallel for
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];
      const vec3i& voxel= mVoxelCoords[particleIndex];

      // neighbors for this particle
      Particle** neighbors= &mNeighbors[particleIndex*mExamineCount];

      // a) examine a local region of 2x2x2 voxels using the spatial index to
      //    look up the particles that we are going to examine
      //    -> create a neighbor map based on the analaysis
      findNeighbors(particle, particleIndex, neighbors, voxel.x, voxel.y, voxel.z);
   }
   timeFindNeighbors = t.elapsed();

   // compute density
   //    we only compute interactions with 32 particles.
   //    this number is somewhat arbitrary, 8 or 16 would work, too but
   //    it might not look too convincing.
   //    -> compute the interaction and the physics (with these 32 particles)
   t.restart();
   #pragma omp parallel for
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      // neighbors for this particle
      Particle** neighbors= &mNeighbors[particleIndex*mExamineCount];
      float* neighborDistances= &mNeighborDistancesScaled[particleIndex*mExamineCount];

      computeDensity(particle, neighbors, neighborDistances);
   }
   timeComputeDensity = t.elapsed();

   // compute pressure
   t.restart();
   #pragma omp parallel for
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      computePressure(particle);
   }
   timeComputePressure = t.elapsed();

   // compute acceleration
   t.restart();
   #pragma omp parallel for
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      // neighbors for this particle
      Particle** neighbors= &mNeighbors[particleIndex*mExamineCount];
      float* neighborDistances= &mNeighborDistancesScaled[particleIndex*mExamineCount];

      computeAcceleration(particle, neighbors, neighborDistances);
   }
   timeComputeAcceleration = t.elapsed();

   // integrate
   t.restart();
   #pragma omp parallel for
   for (int particleIndex = 0; particleIndex < mParticleCount; particleIndex++)
   {
      Particle* particle = &mSrcParticles[particleIndex];

      integrate(particle);
   }
   timeIntegrate = t.elapsed();

   emit updateElapsed(
      timeVoxelize,
      timeFindNeighbors,
      timeComputeDensity,
      timeComputePressure,
      timeComputeAcceleration,
      timeIntegrate
   );

   emit stepFinished();
}


void SPH::pauseResume()
{
   mMutex.lock();
   mPaused = !mPaused;
   mMutex.unlock();
}


void SPH::stopSimulation()
{
   mMutex.lock();
   mStopped = true;
   mMutex.unlock();
}



void SPH::initParticlePositionsRandom()
{
   qsrand(QTime::currentTime().msec());

   for (int i = 0; i < mParticleCount; i++)
   {
      float x = qrand() / (float)RAND_MAX;
      float y = qrand() / (float)RAND_MAX;
      float z = qrand() / (float)RAND_MAX;

      x *= mGridCellsX * mHTimes2 * 0.1f;
      y *= mGridCellsY * mHTimes2 * 0.75f;
      z *= mGridCellsZ * mHTimes2;

      if (x == (float)mGridCellsX)
         x -= 0.00001f;
      if (y == (float)mGridCellsY)
         y -= 0.00001f;
      if (z == (float)mGridCellsZ)
         z -= 0.00001f;

      mSrcParticles[i].mPosition.set(x, y, z);
   }

   // just set up random directions
   for (int i = 0; i < mParticleCount; i++)
   {
      // have a range from -1 to 1
      float x = ((qrand() / (float)RAND_MAX) * 2.0f) - 1.0f;
      float y = ((qrand() / (float)RAND_MAX) * 2.0f) - 1.0f;
      float z = ((qrand() / (float)RAND_MAX) * 2.0f) - 1.0f;

      mSrcParticles[i].mVelocity.set(x, y, z);
   }
}


void SPH::initParticlePolitionsSphere()
{
   qsrand(QTime::currentTime().msec());

   float dist = 0.0f;

   float x = 0.0f;
   float y = 0.0f;
   float z = 0.0f;

   vec3 sphereCenter;
   sphereCenter.set(
      mMaxX * 0.5f,
      mMaxY * 0.9f,
      mMaxZ * 0.5f
   );

   float radius = 25.0f;

   int a = mParticleCount * 0.95f;

   for (int i = 0; i < a; i++)
   {
      do
      {
         x = qrand() / (float)RAND_MAX;
         y = qrand() / (float)RAND_MAX;
         z = qrand() / (float)RAND_MAX;

         x *= mGridCellsX * mHTimes2;
         y *= mGridCellsY * mHTimes2;
         z *= mGridCellsZ * mHTimes2;

         if (x == (float)mGridCellsX)
            x -= 0.00001f;
         if (y == (float)mGridCellsY)
            y -= 0.00001f;
         if (z == (float)mGridCellsZ)
            z -= 0.00001f;

         dist = (vec3(x,y,z) - sphereCenter).length();
      }
      while (dist > radius);

      mSrcParticles[i].mPosition.set(x, y, z);
      mSrcParticles[i].mVelocity.set(0, 0, 0);
   }

   // ground
   for (int i = a; i < mParticleCount; i++)
   {
      do
      {
         x = qrand() / (float)RAND_MAX;
         y = qrand() / (float)RAND_MAX;
         z = qrand() / (float)RAND_MAX;

         x *= mGridCellsX * mHTimes2;
         y *= mGridCellsY * mHTimes2;
         z *= mGridCellsZ * mHTimes2;

         if (x == (float)mGridCellsX)
            x -= 0.00001f;
         if (y == (float)mGridCellsY)
            y -= 0.00001f;
         if (z == (float)mGridCellsZ)
            z -= 0.00001f;

         dist = (vec3(x,y,z) - sphereCenter).length();
      }
      while (y > 10.0f);

      mSrcParticles[i].mPosition.set(x, y, z);
      mSrcParticles[i].mVelocity.set(0, 0, 0);
   }
}



void SPH::clearGrid()
{
   for (int i = 0; i < mGridCellCount; i++)
   {
      mGrid[i].clear();
   }
}


void SPH::voxelizeParticles()
{
   omp_lock_t writelock;
   omp_init_lock(&writelock);

   clearGrid();

   #pragma omp parallel for
   for (int i = 0; i < mParticleCount; i++)
   {
      Particle* particle = &mSrcParticles[i];

      // compute a scalar voxel id from a position
      vec3 pos = particle->mPosition;


      // the voxel size is 2h * 2h * 2h
      // to get the voxel just divide the coordinate by 2h
      // then convert this x, y, z index to a serial number
      // between 0.. maximum number of voxels
      int voxelX = (int)floor(pos.x * mHTimes2Inv);
      int voxelY = (int)floor(pos.y * mHTimes2Inv);
      int voxelZ = (int)floor(pos.z * mHTimes2Inv);

      // it has been seen the positions can run slightly out of bounds for
      // one solver step. so the positions are temporarily fixed here.
      if (voxelX < 0) voxelX= 0;
      if (voxelY < 0) voxelY= 0;
      if (voxelZ < 0) voxelZ= 0;
      if (voxelX >= mGridCellsX) voxelX= mGridCellsX-1;
      if (voxelY >= mGridCellsY) voxelY= mGridCellsY-1;
      if (voxelZ >= mGridCellsZ) voxelZ= mGridCellsZ-1;

      // don't write into particle but into separate memory
      mVoxelCoords[i].x= voxelX;
      mVoxelCoords[i].y= voxelY;
      mVoxelCoords[i].z= voxelZ;

      int voxelId = computeVoxelId(voxelX, voxelY, voxelZ);

      // the lock performs terribly! (~17ms instead of 3ms)
      //      omp_set_lock(&writelock);
      //      mGrid[voxelId].push_back(particle);
      //      omp_unset_lock(&writelock);

      // instead, remember the voxelId and put the particle in later
      mVoxelIds[i]= voxelId;
   }

   // put each particle into according voxel (sequential)
   for (int i = 0; i < mParticleCount; i++)
   {
       Particle* particle = &mSrcParticles[i];
       mGrid[ mVoxelIds[i] ].push_back(particle);
   }
}


void SPH::findNeighbors(Particle* particle, int particleIndex, Particle** neighbors, int voxelX, int voxelY, int voxelZ)
{
   float xOrientation = 0.0f;
   float yOrientation = 0.0f;
   float zOrientation = 0.0f;

   int x = 0;
   int y = 0;
   int z = 0;

   int particleOffset = 0;
   int particleIterateDirection = 0;
   int neighborIndex = 0;
   bool enoughNeighborsFound = false;

   vec3 pos = particle->mPosition;

   // get voxel orientation within voxel
   //
   //
   // 3d illustration
   //
   //    1 voxel:
   //      ____________
   //     /     /     /|
   //    /     /     / |
   //   +-----+-----+  |  the interaction radius of particle x may
   //   |     |  x  |  |  reach into 7 other voxels
   //   |     |     | /|  => 8 voxels to check including 'self'
   //   +-----+-----+/ |
   //   |     |     |  |
   //   |     |     | /
   //   +-----+-----+/
   //
   //   <---< h >--->
   //
   //
   // 2d illustration of 1 slice
   //
   //    +-----+-----+-----+-----
   //    |     |/////|/////|
   //    |     |/////|/////|
   //    +-----+-----+-----+-----
   //    |     |////x|/////|
   //    |     |/////|/////|
   //    +-----+-----+-----+-----
   //    |     |     |     |
   //    |     |     |     |
   //    +-----+-----+-----+-----
   //    |     |     |     |
   //    |     |     |     |

   // this gives us the relative position; i.e the orientation within a voxel
   xOrientation = pos.x - (voxelX * mHTimes2);
   yOrientation = pos.y - (voxelY * mHTimes2);
   zOrientation = pos.z - (voxelZ * mHTimes2);

   // get neighbour voxels
   x = 0;
   y = 0;
   z = 0;

   // retrieve location within voxel
   // 1 voxel side ^= 2h
   // => check in which half of the voxel a particle is located
   //    btw, we just ignore the fact a particle can be positioned exactly on a
   //    center axis of a voxel
   (xOrientation > mH) ? x++ : x--;
   (yOrientation > mH) ? y++ : y--;
   (zOrientation > mH) ? z++ : z--;

   // neighbour voxels
   int vx[8];
   int vy[8];
   int vz[8];

   // same slice
   vx[0] = voxelX;
   vy[0] = voxelY;
   vz[0] = voxelZ;

   // distance 1
   vx[1] = voxelX + x;
   vy[1] = voxelY;
   vz[1] = voxelZ;

   vx[2] = voxelX;
   vy[2] = voxelY + y;
   vz[2] = voxelZ;

   vx[3] = voxelX;
   vy[3] = voxelY;
   vz[3] = voxelZ + z;

   // distance 2
   vx[3] = voxelX + x;
   vy[3] = voxelY + y;
   vz[3] = voxelZ;

   vx[5] = voxelX + x;
   vy[5] = voxelY;
   vz[5] = voxelZ + z;

   vx[6] = voxelX;
   vy[6] = voxelY + y;
   vz[6] = voxelZ + z;

   // distance 3
   vx[7] = voxelX + x;
   vy[7] = voxelY + y;
   vz[7] = voxelZ + z;

   int vxi;
   int vyi;
   int vzi;

   for (int voxelIndex = 0; voxelIndex < 8; voxelIndex++)
   {
      vxi = vx[voxelIndex];
      vyi = vy[voxelIndex];
      vzi = vz[voxelIndex];

      // check if voxels can be processed
      if (
            vxi > 0 && vxi < mGridCellsX
         && vyi > 0 && vyi < mGridCellsY
         && vzi > 0 && vzi < mGridCellsZ
      )
      {
         // +-----+-----+-----+-----+-----+-----+-----+-----
         // |     |     |     |     |     |     |     |
         // |     |     |    0|   *1|     |     |     |
         // +-----+-----+-----+-----+-----+-----+-----+-----
         // |     |     | *  x|  * *|     |     |     |
         // |     |     |  * 2| * *4|     |     |     |
         // +-----+-----+-----+-----+-----+-----+-----+-----
         // |     |     |     |     |     |     |     |
         // |     |     |     |     |     |     |     |
         // +-----+-----+-----+-----+-----+-----+-----+-----

         const QList<Particle*>& voxel = mGrid[computeVoxelId(vxi, vyi, vzi)];

         if (!voxel.isEmpty())
         {
            // randomize idea:
            //
            // for every voxel n
            //
            //    pick a random start particle within voxel n
            //    randomly go up or down from this start particle
            //    check if that particle is within the interaction radius
            //    increment 'found counter'
            //
            //    if (32 particles found),
            //       break
            //
            //    optimization: store dot product in order to avoid
            //                  computing it twice

            // TODO: that's neither fast no a good idea
            //       if there's only 1 particle nearby, the code below is pretty pointless
            particleOffset = qrand() % voxel.length();
            particleIterateDirection = (particleIndex % 2) ? -1 : 1;

            int i = 0;
            while (true)
            {
               int nextIndex = particleOffset + i * particleIterateDirection;

               // leave if we're out out the voxel's bounds
               if (nextIndex < 0 || nextIndex > voxel.length() - 1)
                  break;

               Particle* neighbor = voxel[nextIndex];
               i++;

               Particle* validNeighbor = evaluateNeighbor(particle, neighbor);

               if (validNeighbor)
               {
                  neighbors[neighborIndex] = validNeighbor;
                  neighborIndex++;
               }

               // leave if we have sufficient neighbor particles
               enoughNeighborsFound = (neighborIndex > mExamineCount - 1);
               if (enoughNeighborsFound)
                  break;
            }
         }
      }

      // no need to process any other voxels
      if (enoughNeighborsFound)
         break;
   }

   particle->mNeighborCount = neighborIndex;
}


Particle* SPH::evaluateNeighbor(
   Particle* current,
   Particle* neighbor
)
{
   Particle* validNeighbor = 0;

   if (current != neighbor)
   {
      // we save the sqrt here and use mInteractionRadius2 (h^2) instead
      // of mInteractionRadius
      //
      // float distance = sqrt(dot);
      // if (distance < mInteractionRadius)

      vec3 dist = current->mPosition - neighbor->mPosition;
      float dot = dist * dist;

      // the dot product is unscaled and so is MH2;
      // so there's no need to add any simulation scale here
      if (dot < mH2)
      {
         validNeighbor = neighbor;
      }
   }

   return validNeighbor;
}



void SPH::computeDensity(Particle* particle, Particle** neighbors, float* neighborDistances)
{
   // rho at a given point i is approximately equal to the sum of the masses
   // at various points j (weighted by the weighting function W)
   //
   //    density = (sum of the masses of nearby points * W(r - r , h) )
   //                                                           j
   //    W: smoothing kernel W
   //
   //    if (distance is greater than h return w = 0)
   //    else
   //
   //              315
   //    w =  ------------- ( h^2 - || r - r ||^2 )^3
   //         64 * PI * h^9                 b
   //
   //    || vector || = length of a vector
   //    || vector ||^2 = length of a vector without sqrt
   //                   = x * x + y * y + z * z
   //                   = dot(v, v)
   //                   = dot(v1, v2) =  v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
   //                   = dot(v, v) = v.x * v.x + v.y * v.y + v.z * v.z
   //                   = length(v) = sqrt( v.x * v.x + v.y * v.y + v.z * v.z )
   //                   => length(v)^2 = dot(v,v)


   float density = 0.0f;
   float mass = 0.0f;
   vec3 pos = particle->mPosition;
   float w = 0.0f;
   float rightPart = 0.0f;

   for (int neighborIndex = 0; neighborIndex < particle->mNeighborCount; neighborIndex++)
   {
      Particle* neighbor = neighbors[neighborIndex];

      if (!neighbor)
         break;

      if (neighbor != particle)
      {
         // add mass of neighbor
         mass = neighbor->mMass;

         // apply smoothing kernel to mass
         vec3 dist = pos - neighbor->mPosition;
         float dot = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
         float distance = sqrt(dot);
         float distanceScaled = distance * mSimulationScale;
         neighborDistances[neighborIndex] = distanceScaled;

         if (distanceScaled > mHScaled)
         {
            w = 0.0f;
         }
         else
         {
            // ( h^2 - ||r - r ||^2 )^3
            //                b
            // the dot product is not used here since it's not properly scaled
            rightPart = (mHScaled2 - (distanceScaled * distanceScaled));
            rightPart = (rightPart * rightPart * rightPart);

            //           315
            // w =  -------------  * rightPart
            //      64 * PI * h^9
            w = mKernel1Scaled * rightPart;

            // apply weighted neighbor mass to our density
            density += (mass * w);
         }
      }
   }

   float inverseDensity = ((density > 0.0f) ? (1.0f / density) : 1.0f);
   particle->mDensity = density;

   // we'll need those for the acceleration evaluation later
   // that's why we precompute them now
   particle->mDensityInverse = inverseDensity;
   particle->mDensityInverse2 = inverseDensity * inverseDensity;
}


void SPH::computePressure(Particle* particle)
{
   // pressure p = k(rho - rho0)
   // pressure is the difference between rho and rho (the resting potential)

   // rho0: resting density
   float deltaRho = particle->mDensity - mRho0;
   float p = mStiffness * deltaRho;
   particle->mPressure = p;
}


void SPH::computeAcceleration(Particle* p, Particle** neighbors, float* neighborDistances)
{
   // pressure gradient
   //
   //    (1)
   //
   //    pressure gradient (divided by rho ) at a given point i is
   //                                     i
   //    approximately equal to
   //                                                     p              p
   //                                                      i              j
   //    the sum of the masses at various points j  *  --------  +    --------
   //                                                   rho  ^2        rho  ^2
   //                                                      i              j
   //
   //    (2)
   //                                                     __
   //    multiplied by the gradient of a smoothing kernel \/W
   //
   //    The gradient is a vector; therefore this whole expression is a
   //    vector value.
   //                                                  r - r
   //       __       -45                                    b
   //       \/W =  -------- * (h - ||r - r  ||)^2 * ------------
   //              PI * h^6               b         || r - r  ||
   //                                                       b


   // viscosity term
   //
   //     u = a scalar coefficient that says how viscous is the fluid
   //         small value => water
   //         large value => sirup
   //
   //     the viscosity term is approximately
   //
   //         u
   //        ----  *  the sum of the masses at various points j
   //        rho      multiplied by the velocity between two points i and j
   //           i     divided by the pressure of j                     __
   //                 multiplied by the laplacian of a smoothin kernel \/^2W
   //
   //
   //     that is:
   //
   //               (3)
   //                | starting here
   //                |
   //                        v    v
   //         u               j -  i   __
   //        ---- * SUM  m * ------- * \/^2W
   //        rho     j    j    rho
   //           i                 j
   //
   //        __         45
   //        \/^2W = -------- * (h - ||r - r  ||)
   //                PI * h^6               b

   Particle* neighbor = 0;
   float distanceToNeighborScaled = 0.0f;

   float pi = p->mPressure;
   float rhoiInv = p->mDensityInverse;
   float rhoiInv2 = p->mDensityInverse2;
   float piDivRhoi2 = pi * rhoiInv2;
   vec3 r = p->mPosition;
   vec3 vi = p->mVelocity;

   float pj = 0.0f;
   float rhoj = 0.0f;
   float rhojInv = 0.0f;
   float rhojInv2 = 0.0f;
   float mj = 0.0f;
   vec3 rj;
   vec3 vj;
   vec3 rMinusRj;
   vec3 rMinusRjScaled;

   // pressure gradient..
   vec3 pressureGradient(0.0f, 0.0f, 0.0f);
   vec3 pressureGradientContribution;

   // ..and viscous term
   vec3 viscousTerm(0.0f, 0.0f, 0.0f);
   vec3 viscousTermContribution;

   // are added to the final acceleration
   vec3 acceleration(0.0f, 0.0f, 0.0f);

   for (int neighborIndex = 0; neighborIndex < p->mNeighborCount; neighborIndex++)
   {
      neighbor = neighbors[neighborIndex];

      pj = neighbor->mPressure;
      rhoj = neighbor->mDensity;
      rhojInv = neighbor->mDensityInverse;
      rhojInv2 = neighbor->mDensityInverse2;
      rj = neighbor->mPosition;
      vj = neighbor->mVelocity;
      mj = neighbor->mMass;

      // pressure gradient
      //
      // (2)
      //
      //    r - r
      //         j
      // ------------
      // || r - r  ||
      //         j
      rMinusRj = (r - rj);
      rMinusRjScaled = rMinusRj * mSimulationScale;
      distanceToNeighborScaled = neighborDistances[neighborIndex];
      pressureGradientContribution = rMinusRjScaled / distanceToNeighborScaled;

      //   -45
      // --------
      // PI * h^6
      pressureGradientContribution *= mKernel2Scaled;

      // (h - ||r - r  ||)^2
      //             j
      float centerPart = (mHScaled - distanceToNeighborScaled);
      centerPart = centerPart * centerPart;
      pressureGradientContribution *= centerPart;

      // (1)
      //      p          p
      //       i          j
      // m  * -------- + --------
      //  j   rho  ^2    rho  ^2
      //         i          j
      float factor = mj * piDivRhoi2 * (pj * rhojInv2);
      pressureGradientContribution *= factor;

      // add pressure gradient contribution to pressure gradient
      pressureGradient += pressureGradientContribution;


      // viscosity
      //
      // (3)
      //
      //     v    v
      //      j -  i
      // m * -------
      //  j    rho
      //          j
      viscousTermContribution = vj - vi;
      viscousTermContribution *= rhojInv;
      viscousTermContribution *= mj;

      // (4)
      //      45
      //   -------- * (h - ||r - r  ||)
      //   PI * h^6               b
      viscousTermContribution *= mKernel3Scaled;
      viscousTermContribution *= (mHScaled - distanceToNeighborScaled);

      // add contribution to viscous term
      viscousTerm += viscousTermContribution;
   }

   viscousTerm *= (mViscosityScalar * rhoiInv);

   // potential optimization:
   // multiplication with rhoiInv could be done here in just one go
   acceleration = viscousTerm - pressureGradient;

   // check CFL condition
   float dot =
        acceleration.x * acceleration.x
      + acceleration.y * acceleration.y
      + acceleration.z * acceleration.z;

   bool limitExceeded = (dot > mCflLimit2);
   if (limitExceeded)
   {
      float length = sqrt(dot);
      float cflScale = mCflLimit / length;
      acceleration *= cflScale;
   }

   // yay. done.
   p->mAcceleration = acceleration;
}


void SPH::integrate(Particle* p)
{
   vec3 acceleration = p->mAcceleration;
   vec3 position = p->mPosition;
   vec3 velocity = p->mVelocity;

   // apply external forces
   acceleration += mGravity;

   // semi-implicit euler integration
   float posTimeStep = mTimeStep * mSimulationScaleInverse;
   vec3 newVelocity = velocity + (acceleration * mTimeStep);
   vec3 newPosition = position + (newVelocity * posTimeStep);

   // let particles bounce back if they collide with the grid boundaries
   handleBoundaryConditions(
      position,
      &newVelocity,
      posTimeStep,
      &newPosition
   );

   p->mVelocity = newVelocity;
   p->mPosition = newPosition;
}


void SPH::handleBoundaryConditions(
   vec3 position,
   vec3* newVelocity,
   float timeStep,
   vec3* newPosition
)
{
   // x coord
   if (newPosition->x < 0.0f)
   {
      vec3 normal(1, 0, 0);
      float intersectionDistance = -position.x / newVelocity->x;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
   else if (newPosition->x > mMaxX)
   {
      vec3 normal(-1, 0, 0);
      float intersectionDistance = (mMaxX - position.x) / newVelocity->x;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }

   // y coord
   if (newPosition->y < 0.0f)
   {
      vec3 normal(0, 1, 0);
      float intersectionDistance = -position.y / newVelocity->y;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
   else if (newPosition->y > mMaxY)
   {
      vec3 normal(0, -1, 0);
      float intersectionDistance = (mMaxY - position.y) / newVelocity->y;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }

   // z coord
   if (newPosition->z < 0.0f)
   {
      vec3 normal(0, 0, 1);
      float intersectionDistance = -position.z / newVelocity->z;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
   else if (newPosition->z > mMaxZ)
   {
      vec3 normal(0, 0, -1);
      float intersectionDistance = (mMaxZ - position.z) / newVelocity->z;

      applyBoundary(
         position,
         timeStep,
         newPosition,
         intersectionDistance,
         normal,
         newVelocity
      );
   }
}


void SPH::applyBoundary(
      vec3 position,
      float timeStep,
      vec3* newPosition,
      float intersectionDistance,
   vec3 normal,
   vec3* newVelocity
)
{
   vec3 intersection = position + (*newVelocity * intersectionDistance);

   float dotProduct =
        newVelocity->x * normal.x
      + newVelocity->y * normal.y
      + newVelocity->z * normal.z;

   vec3 reflection = *newVelocity - (normal * dotProduct * 2.0f);

   float remaining = timeStep - intersectionDistance;

   // apply boundaries
   *newVelocity = reflection;
   *newPosition = intersection + reflection * (remaining * mDamping);
}


int SPH::computeVoxelId(int voxelX, int voxelY, int voxelZ)
{
   return (voxelZ * mGridCellsY + voxelY) * mGridCellsX + voxelX;
}


void SPH::clearNeighbors()
{
   memClear32(mNeighbors, mParticleCount * mExamineCount * sizeof(Particle*));
}


void SPH::memClear32(void* dst, int len)
{
   unsigned int* dst32= (unsigned int*)dst;
   len>>=2;
   while (len--)
      *dst32++= 0;
}


float SPH::getCellSize() const
{
   return mCellSize;
}


Particle* SPH::getParticles()
{
   return mSrcParticles;
}


int SPH::getParticleCount() const
{
   return mParticleCount;
}


void SPH::getGridCellCounts(int &x, int &y, int &z)
{
   x = mGridCellsX;
   y = mGridCellsY;
   z = mGridCellsZ;
}


void SPH::getParticleBounds(float &x, float &y, float &z)
{
   x = mMaxX;
   y = mMaxY;
   z = mMaxZ;
}


float SPH::getInteractionRadius2() const
{
   return mHScaled2;
}


QList<Particle *>* SPH::getGrid()
{
   return mGrid;
}



vec3 SPH::getGravity() const
{
   return mGravity;
}


void SPH::setGravity(const vec3 &gravity)
{
   mGravity = gravity;
}


float SPH::getCflLimit() const
{
   return mCflLimit;
}


void SPH::setCflLimit(float cflLimit)
{
   mCflLimit = cflLimit;
   mCflLimit2 = mCflLimit * mCflLimit;
}


float SPH::getDamping() const
{
   return mDamping;
}


void SPH::setDamping(float damping)
{
   mDamping = damping;
}


float SPH::getTimeStep() const
{
   return mTimeStep;
}


void SPH::setTimeStep(float timeStep)
{
   mTimeStep = timeStep;
}


float SPH::getViscosityScalar() const
{
   return mViscosityScalar;
}


void SPH::setViscosityScalar(float viscosityScalar)
{
   mViscosityScalar = viscosityScalar;
}


float SPH::getStiffness() const
{
   return mStiffness;
}


void SPH::setStiffness(float stiffness)
{
   mStiffness = stiffness;
}

