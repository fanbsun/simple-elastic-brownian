// This is the particle class
// It lists features of the particle as members and also lists how its position, velocity, etc are updated

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <cmath>
#include <vector>
#include "vector3d.h"
#include "edge.h"

class EDGE;


class PARTICLE
{
  public:

  // members
  int ty;           // type of the particle
  double diameter;	// diameter of the particle
  double m; 		// mass of the particle
  VECTOR3D bx;          //vector of box lengths (lj reduced unit)
  VECTOR3D position;	// position vector of the particle
  VECTOR3D velocity;	// velocity vector of the particle
  VECTOR3D force;	// force vector on the particle
  VECTOR3D sforce;
  double pe;		// potential energy
  double ke;	    // kinetic energy
  double se;
  VECTOR3D fdrag;
  VECTOR3D fran;
  std::vector<EDGE*> itsE;              //its Edges

  
  // member functions
  
  // make a particle
  PARTICLE(int initial_type = 1, double initial_diameter = 0, double initial_mass = 0, VECTOR3D initial_position = VECTOR3D(0,0,0), VECTOR3D initial_velocity = VECTOR3D(0,0,0), VECTOR3D boxlength_vec = VECTOR3D(4,4,4))
  {
    ty = initial_type;
    diameter = initial_diameter;
    m = initial_mass;
    position = initial_position;
    velocity = initial_velocity;
    bx = boxlength_vec;
  }
  
  // the next two functions are central to the velocity-Verlet algorithm:

  // update position of the particle
  void update_position(double dt, double L)		// dt is the time-step, L is the box length
  {
    position = ( position + (velocity ^ dt) );	// position updated to a full time-step
      if (position.x > L/2.0)
          position.x = position.x - L;
      if (position.x < -L/2.0)
          position.x = position.x + L;
      if (position.y > L/2)
          position.y = position.y - L;
      if (position.y < -L/2.0)
          position.y = position.y + L;
      if (position.z > L/2.0)
          position.z = position.z - L;
      if (position.z < -L/2.0)
          position.z = position.z + L;
  }



  // update velocity of the particle
  void update_velocity(double dt)
  {
    velocity = ( velocity + ( (force) ^ ( dt / (2.0*m) ) ));	// notice the half time-step
  }

  void update_velocity_brownian_firsthalf(double half_dt, double reduced_energy, double fric_zeta, double guassian)
  {
      velocity = ( velocity + ( (force) ^ ( half_dt / m ) ) + sqrt(reduced_energy * fric_zeta * half_dt) * guassian);	// notice the half time-step
  }

    // calculate kinetic energy of a particle
  void kinetic_energy()				
  {
    ke = 0.5 * m * velocity.Magnitude() * velocity.Magnitude();
  }

   void update_stretching_energy(double ks);

   void update_stretching_force(double ks);

   VECTOR3D dist(PARTICLE *A, PARTICLE *B) {
        VECTOR3D r_vec; //= (A->pos - B->pos);
        r_vec.x = A->position.x - B->position.x;
        r_vec.y = A->position.y - B->position.y;
        r_vec.z = A->position.z - B->position.z;
        VECTOR3D box = A->bx;
        VECTOR3D hbox = A->bx ^ 0.5;
        if (r_vec.x > hbox.x) r_vec.x -= box.x;
        else if (r_vec.x < -hbox.x) r_vec.x += box.x;
        if (r_vec.y > hbox.y) r_vec.y -= box.y;
        else if (r_vec.y < -hbox.y) r_vec.y += box.y;
        if (r_vec.z > hbox.z) r_vec.z -= box.z;
        else if (r_vec.z < -hbox.z) r_vec.z += box.z;
        return r_vec;
    }
};

#endif