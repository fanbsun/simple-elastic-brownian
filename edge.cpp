#include <assert.h>
#include <iostream>
#include "edge.h"
#include "particle.h"
#include "functions.h"

using namespace std;

//Finds particle opposite to the one referred in a edge
PARTICLE* EDGE::opposite(PARTICLE *theP){
   assert(itsP.size()==2);
   if (itsP[0]==theP){
      return itsP[1];
   }else if (itsP[1]==theP){
      return itsP[0];
   } return NULL;
}

//update's the length of a edge (PBC)
void EDGE::update_length(){
   lengthvec = itsP[0]->position - itsP[1]->position;
   VECTOR3D box = (itsP[0]->bx);
   VECTOR3D halfbox = (itsP[0]->bx) ^ 0.5;
   if (lengthvec.x > halfbox.x) lengthvec.x -= box.x;
   else if (lengthvec.x < -halfbox.x) lengthvec.x += box.x;
   if (lengthvec.y > halfbox.y) lengthvec.y -= box.y;
   else if (lengthvec.y < -halfbox.y) lengthvec.y += box.y;
   if (lengthvec.z > halfbox.z) lengthvec.z -= box.z;
   else if (lengthvec.z < -halfbox.z) lengthvec.z += box.z;
      
   length = lengthvec.Magnitude(); //edgelength
}