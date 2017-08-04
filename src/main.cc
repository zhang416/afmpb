#include "afmpb.h"

int main(int argc, char **argv) {
  auto parameter = afmpb::init(argc, argv); 

  if (parameter) {
    afmpb::AFMPB system(std::move(parameter)); 

    auto status = system.computePotential(); 

   
    
    //if (system.solve())  
    if (status)
      system.collect(); 


    system.finalize(); 
  } 

  auto err = afmpb::finalize(); 
  assert(err == 0); 

  return 0;
}

