#ifndef _EDGE_H
#define _EDGE_H

#include <vector>
#include "vector3d.h"
#include "particle.h"

class PARTICLE;

class EDGE {
public:
//member variables
   int id;                                      //identification
   double length;                               //actual length
   double len0;                                 //ideal length
   VECTOR3D lengthvec;                          //actual length vector
   std::vector<PARTICLE*> itsP;                    //vector of particles in the bond

   EDGE(VECTOR3D initial = VECTOR3D(0, 0, 0), int id_i = 0){                
      id = id_i;
      lengthvec = initial;
   }

   PARTICLE *opposite(PARTICLE *theP);
   void update_length();
};

#endif //_EDGE_H
