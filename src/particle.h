#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec3.h"

class Particle
{
public:

   Particle();

   vec3 mPosition;
   vec3 mVelocity;
   vec3 mAcceleration;

   float mMass;
   float mDensity;
   float mDensityInverse;
   float mDensityInverse2;
   float mPressure;

   int mNeighborCount;
};

#endif // PARTICLE_H
