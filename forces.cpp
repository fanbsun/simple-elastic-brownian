// This file contains the routine that computes the force on the particle exerted by all other particles

using namespace std;

#include <vector>
#include "particle.h"
#include "edge.h"
#include<iostream>

void compute_forces(vector<PARTICLE>& particle, vector<EDGE>& edge, double lj_epsilon, double r_c, double L, double ks)
{
  double r;
  VECTOR3D r_vec, flj;

   for (unsigned int i = 0; i < edge.size(); i++){
        edge[i].update_length();
    }
   for (unsigned int i = 0; i < particle.size(); i++){
        particle[i].update_stretching_force(ks);
   }

  
  for (unsigned int i = 0; i < particle.size(); i++)
  {
    flj = VECTOR3D(0,0,0);
    for (unsigned int j = 0; j < particle.size(); j++)
    {
      if (j == i) continue;
      r_vec = particle[i].position - particle[j].position;

      if (r_vec.x > L/2)
          r_vec.x = r_vec.x - L;
      if (r_vec.x < -L/2)
          r_vec.x = r_vec.x + L;
      if (r_vec.y > L/2)
          r_vec.y = r_vec.y - L;
      if (r_vec.y < -L/2)
          r_vec.y = r_vec.y + L;
      if (r_vec.z > L/2)
          r_vec.z = r_vec.z - L;
      if (r_vec.z < -L/2)
          r_vec.z = r_vec.z + L;

      r = r_vec.Magnitude();
      if (r < r_c)
          //flj = flj + (r_vec / r) * (48 * (1 / r) * ((1 / pow(r,12)) - 0.5 * (1 / pow(r,6))));
          flj = flj + (r_vec / r) ^ (18 * (1 / r) * ((1 / pow(r,9)) - (1 / pow(r,6))));
      else
          flj = flj + VECTOR3D(0,0,0);
    }
    particle[i].force =  flj+ particle[i].se;
  }
}