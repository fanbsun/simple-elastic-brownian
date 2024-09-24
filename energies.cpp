// This function has the potential energy evaluations

using namespace std;
#include <vector>
#include "particle.h"
#include<iostream>

double compute_potential_energies(vector<PARTICLE>& particle, double lj_epsilon, double r_c, double L)
{
  // what you need to compute energy
  VECTOR3D r_vec;
  double r;
  
  // energy is computed as pair-wise sums
  for (unsigned int i = 0; i < particle.size(); i++)
  {
    double uljpair = 0.0;
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
      if (r < r_c) {
          double energy_shift = 1 * (2 / pow(r_c, 9) - 3 / pow(r_c, 6));
          uljpair = uljpair + 1 * (2 / pow(r, 9) - 3 / pow(r, 6)) - energy_shift;
      }
      else
          uljpair = uljpair + 0;
    }
    particle[i].pe = uljpair;
  }
  
  double total_pe = 0;
  for (unsigned int i = 0; i < particle.size(); i++) {
      total_pe += particle[i].pe;
  }
  total_pe = 0.5 * total_pe;    // factor of 0.5 to account for double counting

  return total_pe;
}