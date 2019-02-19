#ifndef VEC3_H
#define VEC3_H

#include <math.h>

class vec3
{

   public:

      vec3();
      vec3(float x, float y, float z);

      float x;
      float y;
      float z;


      inline void set(float sx, float sy, float sz)
      {
         x = sx;
         y = sy;
         z = sz;
      }

      // dot product
      inline float operator * (const vec3& v) const
      {
         return x * v.x + y * v.y + z * v.z;
      }


      inline vec3 operator + (const vec3 &rhs) const
      {
         return vec3(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z);
      }


      inline vec3 operator - (const vec3 &rhs) const
      {
         return vec3(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z);
      }


      inline vec3 operator * (float rhs) const
      {
         return vec3(this->x * rhs, this->y * rhs, this->z * rhs);
      }


      inline vec3 operator / (float rhs) const
      {
         rhs= 1.0f / rhs;
         return vec3(x*rhs, y*rhs, z*rhs);
      }


      inline void operator += (const vec3 &rhs)
      {
         this->x += rhs.x;
         this->y += rhs.y;
         this->z += rhs.z;
      }


      inline void operator -= (const vec3& rhs)
      {
         this->x -= rhs.x;
         this->y -= rhs.y;
         this->z -= rhs.z;
      }


      inline void operator /= (float rhs)
      {
         rhs= 1.0f / rhs;
         this->x *= rhs;
         this->y *= rhs;
         this->z *= rhs;
      }


      inline void operator *= (float rhs)
      {
         this->x *= rhs;
         this->y *= rhs;
         this->z *= rhs;
      }

      inline bool operator == (const vec3 &rhs) const
      {
         return (this->x == rhs.x && this->y == rhs.y && this->z == rhs.z);
      }


      inline bool operator != (const vec3 &rhs) const
      {
         return this->x != rhs.x || this->y != rhs.y || this->z != rhs.z;
      }

      inline float length() const
      {
         return static_cast<float>(sqrt(x * x + y * y + z * z));
      }

      // returns length()^2
      inline float length2() const
      {
         return x * x + y * y + z * z;
      }
};

#endif // VEC3_H
