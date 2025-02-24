#include "OpenMM.h"

#ifndef _INTERFACE_H_
#define _INTERFACE_H_

// This is our opaque "handle" class containing all the OpenMM objects that
// must persist from call to call during a simulation. The main program gets 
// a pointer to one of these but sees it as essentially a void* since it 
// doesn't know the definition of this class.
struct MyOpenMMData {
    MyOpenMMData() : system(0), context(0), integrator(0), nonbond(0), custombond(0),cmmr(0), barostat(0) {}
    ~MyOpenMMData() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;
    OpenMM::NonbondedForce* nonbond;
    OpenMM::CustomBondForce* custombond;
    OpenMM::CMMotionRemover* cmmr;
    OpenMM::MonteCarloBarostat* barostat;
};

// This is a struct to wrap all the GCMC paramters and pass them to create a GCMC
// object
struct GcmcParameters {
    //MODELING AND SIMULATION PARAMETERS
    double temperature;
    double beta;
    double mu;
    double boxEdgeLength;
    double cutoffCluster;
    double eps;
    double cationWavelengthCube;
    double anionWavelengthCube;
    
    //FORCE FIELD DATA
    double cationCharge;
    double cationSigma;
    double cationEpsilon;
    double anionCharge;
    double anionSigma;
    double anionEpsilon;

    //GCMC PARAMETERS
    int    M;                     
    int    upperLimit;               
    int    lowerLimit;                
    int    mdSteps;    

    //GCMC Initial State Parameters
    int numReal;
    int numGhost;
    std::vector<int> cationRealList;
    std::vector<int> cationGhostList;
    std::vector<int> anionRealList;
    std::vector<int> anionGhostList;

    //Name Information
    
    //Simulation Object 
    MyOpenMMData *omm;
                  
};


#endif
