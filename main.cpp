// This is main.
// This is MD simulation of many particles interacting with a pair potential

using namespace std;

#include "functions.h"
#include "particle.h"
#include "edge.h"
#include <iostream>
#include <chrono>
#include <random>

void compute_forces(vector<PARTICLE>&, vector<EDGE>&, double, double, double, double);
double compute_potential_energies(vector<PARTICLE>&, double, double, double);

int main(int argc, char* argv[]) 
{
  // we begin with defining the Boltzmann constant
  const double kB = 1.38e-23;	// Joules per Kelvin

  cout << "\n---Simulating fluid (Argon)---\n";

  // key descriptors of the many particle system
  double diameter = 3.405e-10;  // what is the size of each particle? (in meters)
  double mass = 6.634e-26;  // what is the mass of each particle? (in kg)
  double ljenergy = 120 * kB; // what is the strength of typical particle-particle collisions? (in Joules)

  // reduced units
  double unitlength = diameter; // unit of length is the simulation (diameter of the particle)
  double unitmass = mass; // unit of mass in the simulation (mass of the particle)
  double unitenergy = ljenergy; // unit of energy in the simulation (characteristic pair potential strength)
  double unittemperature = ljenergy/kB;

  double unittime = sqrt(unitmass * unitlength * unitlength / unitenergy); // Unit of time

  cout << "unit of length is " << unitlength << " meters" << endl;
  cout << "unit of mass is " << unitmass << " kilograms" << endl;
  cout << "unit of energy is " << unitenergy << " Joules" << endl;
  cout << "unit of time is " << unittime << " seconds" << endl;
  cout << "unit of temperature is " << unittemperature << " Kelvin" << endl;
  cout << "\n";

  double reduced_diameter = diameter / unitlength;
  double reduced_mass = mass / unitmass;
  double reduced_ljenergy = ljenergy / unitenergy;
  double reduced_temperature = 1.0; // 120K
  double damp = 1000; //unit of time
    double ks = 5000.0;

  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = 69501;
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution (0.0,1.0);

  std::uniform_real_distribution<> distr(-0.5,0.5);

  vector<PARTICLE> particle;
  vector<EDGE> particle_edge;

  double boxL = 4; // box length
  double distance_cutoff = 2.5;
  //create particle
    //uncomment here for 2 particles with 1 spring
  PARTICLE fresh_particle1 = PARTICLE(1,reduced_diameter,reduced_mass,VECTOR3D(0,0,0),VECTOR3D(1,0,0), VECTOR3D(boxL,boxL,boxL));
  PARTICLE fresh_particle2 = PARTICLE(1,reduced_diameter,reduced_mass,VECTOR3D(1,0,0),VECTOR3D(-1,0,0), VECTOR3D(boxL,boxL,boxL));
  particle.push_back(fresh_particle1);
  particle.push_back(fresh_particle2);


    //uncomment here for 3 particles connected with a ring
    /*

    PARTICLE fresh_particle1 = PARTICLE(1,reduced_diameter,reduced_mass,VECTOR3D(0,0,0),VECTOR3D(1,0,0), VECTOR3D(boxL,boxL,boxL));
    PARTICLE fresh_particle2 = PARTICLE(1,reduced_diameter,reduced_mass,VECTOR3D(1,0,0),VECTOR3D(-1,0,0), VECTOR3D(boxL,boxL,boxL));
    PARTICLE fresh_particle3 = PARTICLE(1,reduced_diameter,reduced_mass,VECTOR3D(0.5,0.866,0),VECTOR3D(0,0,0), VECTOR3D(boxL,boxL,boxL));
    particle.push_back(fresh_particle1);
    particle.push_back(fresh_particle2);
    particle.push_back(fresh_particle3);
    */

  //create edge
    //uncomment here for 2 particles with 1 spring
    EDGE edge1 = EDGE(VECTOR3D(-1,0,0), 0);

    particle_edge.push_back(edge1);
    particle_edge[0].len0 = 1;
    particle_edge[0].length = 1;
    particle_edge[0].itsP.push_back(&particle[0]);
    particle_edge[0].itsP.push_back(&particle[1]);
    particle[0].itsE.push_back(&particle_edge[0]);
    particle[1].itsE.push_back(&particle_edge[0]);

    //uncomment here for 3 particles connected with a ring
    /*
  EDGE edge1 = EDGE(VECTOR3D(-1,0,0), 0);
    EDGE edge2 = EDGE(VECTOR3D(0.5,-0.866,0), 1);
    EDGE edge3 = EDGE(VECTOR3D(0.5,0.866,0), 2);

  particle_edge.push_back(edge1);
    particle_edge.push_back(edge2);
    particle_edge.push_back(edge3);
  particle_edge[0].len0 = 1;
  particle_edge[0].length = 1;
    particle_edge[0].itsP.push_back(&particle[0]);
    particle_edge[0].itsP.push_back(&particle[1]);
    particle[0].itsE.push_back(&particle_edge[0]);
    particle[1].itsE.push_back(&particle_edge[0]);

    particle_edge[1].len0 = 1;
    particle_edge[1].length = 1;
    particle_edge[1].itsP.push_back(&particle[1]);
    particle_edge[1].itsP.push_back(&particle[2]);
    particle[1].itsE.push_back(&particle_edge[1]);
    particle[2].itsE.push_back(&particle_edge[1]);

    particle_edge[2].len0 = 1;
    particle_edge[2].length = 1;
    particle_edge[2].itsP.push_back(&particle[2]);
    particle_edge[2].itsP.push_back(&particle[0]);
    particle[2].itsE.push_back(&particle_edge[2]);
    particle[0].itsE.push_back(&particle_edge[2]);

    */



  // output to screen
  cout << "particle diameter in meters: " << diameter << " and in reduced units: " << reduced_diameter << endl;
  cout << "particle mass in kg: " << mass << " and in reduced units: " << reduced_mass << endl;
  cout << "pair potential energy strength in Joules: " << ljenergy << " and in reduced units: " << reduced_ljenergy << endl;
  cout << "Number of particles inside the box: " << particle.size() << endl;
  cout << "box size is " << boxL << endl;
  cout << "\n";

  // initial energies and forces computation
  double totalke = 0.0;
    double se = 0.0;
    double length = 0.0;

  for (unsigned int i = 0; i < particle.size(); i++)
  {
      particle[i].kinetic_energy();
      totalke += particle[i].ke;
  }

    for (unsigned int i = 0; i < particle_edge.size(); i++){
        particle_edge[i].update_length();
        length = particle_edge[0].length;
    }
    for (unsigned int i = 0; i < particle.size(); i++){
        particle[i].update_stretching_energy(ks);
    }
    for (unsigned int i = 0; i < particle.size(); i++){
        se += particle[i].se;
    }



  double totalpe = compute_potential_energies(particle, reduced_ljenergy, distance_cutoff, boxL);
  compute_forces(particle, particle_edge, reduced_ljenergy, distance_cutoff, boxL, ks);

  // create files for storing movie and energies

  char file_movie[200], file_energy[200];

  // file movie
  sprintf(file_movie, "movie.out");
  ofstream list_propagation(file_movie, ios::out); // create a file to store and visualize 3D data
  make_movie(0,particle,list_propagation, boxL);

  // file ke, pe, and total energies of all particles
  sprintf(file_energy, "energy.out");
  ofstream output_energy(file_energy, ios::out);
    VECTOR3D average_velocity_vector_init = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < particle.size(); i++) {
        //particle[i].update_velocity_brownian_secondhalf(delta_t, reduced_ljenergy, damp, guassian_random);
        average_velocity_vector_init = average_velocity_vector_init + particle[i].velocity;
    }


  output_energy << 0 << "  " << totalke  << "  "  << se << "  "   << length <<  "  " << totalpe << "  " << totalke+totalpe+se<< "  " << average_velocity_vector_init.Magnitude() << endl;

  // print energies
  cout << "initial total kinetic energy of all particles (in reduced units): " << totalke << endl;
    cout << "initial total stretching energy (in reduced units): " << se << endl;
  cout << "initial total potential energy of all particles (in reduced units): " << totalpe << endl;
  cout << "initial total system energy (ke + pe + se) of all particles (in reduced units): " << totalke+totalpe+se << endl;
  cout << "\n";

  double totaltime = 100;
  int steps = 50000*2;		// number of time discretizations (slices)
  double delta_t = 0.002;	// choose steps carefully
  int movie_step = 100;  // movie will be made every movie_step
  int energycalc_step = 100; // energies will be computed every energycalc_step

  cout << "Fluid simulation time (in seconds): " << totaltime*unittime << " and in reduced units: " << totaltime << endl;

  // Molecular Dynamics

//  Brownian Dynamics Equations
//  r1 = -m / damp;
//  r2 = sqrt((24*m*kB)/(damp*dt));
//  frandom = r2 * uniform(-0.5, 0.5);
//  fdrag = r1 * v


  cout << "progress..." << endl;
  for (int num = 1; num <= steps; num++)
  {
      double guassian_random = distribution(generator);

      //compute fdrag
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].fdrag = particle[i].velocity^(-1.0 * particle[i].m / damp);
      }

      // velocity-Verlet
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_velocity(delta_t);
      }

      for (unsigned int i = 0; i < particle.size(); i++) {
          //particle[i].update_position_brownian(delta_t, boxL, reduced_ljenergy, damp, guassian_random);  // update position full timestep
          particle[i].update_position(delta_t, boxL);
      }

      compute_forces(particle, particle_edge,reduced_ljenergy, distance_cutoff, boxL,ks);  // expensive step

      //compute fran
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].fran.x = sqrt((particle[i].m*24.0*reduced_ljenergy)/(damp * delta_t)) * distr(generator);
          particle[i].fran.y = sqrt((particle[i].m*24.0*reduced_ljenergy)/(damp * delta_t)) * distr(generator);
          particle[i].fran.z = sqrt((particle[i].m*24.0*reduced_ljenergy)/(damp * delta_t)) * distr(generator);
      }

      //add fran and fdrag to flj
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].force = particle[i].force + particle[i].fdrag + particle[i].fran;
      }

      for (unsigned int i = 0; i < particle.size(); i++) {
          //particle[i].update_velocity_brownian_secondhalf(delta_t, reduced_ljenergy, damp, guassian_random);
          particle[i].update_velocity(delta_t);
      }
      VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
      for (unsigned int i = 0; i < particle.size(); i++) {
          //particle[i].update_velocity_brownian_secondhalf(delta_t, reduced_ljenergy, damp, guassian_random);
          average_velocity_vector = average_velocity_vector + particle[i].velocity;
      }
       average_velocity_vector = average_velocity_vector*(1.0/particle.size());
       for (unsigned int i = 0; i < particle.size(); i++){
          particle[i].velocity = particle[i].velocity - average_velocity_vector;
       }

      // calling a movie function to get a movie of the simulation every movie_step
      if (num%movie_step == 0)
          make_movie(num,particle,list_propagation,boxL);

      // calculating energies every energycalc_step
      if (num%energycalc_step == 0) {
          totalke = 0.0;
          se = 0.0;
          for (unsigned int i = 0; i < particle.size(); i++) {
              particle[i].kinetic_energy();
              totalke += particle[i].ke;
          }
          for (unsigned int i = 0; i < particle_edge.size(); i++){
              particle_edge[i].update_length();
              length = particle_edge[0].length;
          }
          for (unsigned int i = 0; i < particle.size(); i++){
              particle[i].update_stretching_energy(ks);
          }
          for (unsigned int i = 0; i < particle.size(); i++){
              se += particle[i].se;
          }
          totalpe = compute_potential_energies(particle, reduced_ljenergy, distance_cutoff, boxL);
          // outputting the energy to make sure simulation can be trusted
          output_energy << 0 << "  " << totalke  << "  "  << se << "  "   << length <<  "  " << totalpe << "  " << totalke+totalpe+se<< "  " << average_velocity_vector.Magnitude() << endl;
      }

      // monitoring progress of the simulation
      double progress = ((num)/(double)steps);
      ProgressBar(progress);
  }

  return 0;
} 
// End of main