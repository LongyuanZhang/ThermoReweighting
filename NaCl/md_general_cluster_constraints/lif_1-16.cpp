/* -----------------------------------------------------------------------------
 *               OpenMM(tm) HelloEthane example in C++ (June 2009)
 * -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <math.h>
#include <typeinfo>
#include <time.h>
#include <chrono>
#include <fstream>
#include "Interface.h"
#include "GcmcIonPair.h"
#include "OpenMM.h"

#include <boost/utility.hpp> // for boost::tie
#include <boost/graph/graphviz.hpp>

#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <set>

//#include <cuda_profiler_api.h>

// -----------------------------------------------------------------------------
//                                 MOCK MD CODE
// -----------------------------------------------------------------------------

//                   MODELING AND SIMULATION PARAMETERS
const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 2000;    // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 40;      // total simulation time (ps)
static const bool   WantEnergy   = true;

const double Temperature         = 298.15;  // Kelvins
const double FrictionInPerPs     = 1.0;     // collisions per picosecond
const double CutoffDistanceInAng = 10.0;    // Angstroms
//const double BoxEdgeLengthInNm   = 2.2;   //nms
//const double BoxEdgeLengthInNm   = 2.2;   //nms
const double BoxEdgeLengthInNm   = 4.1031;   //nms
//const double VolumeInNm3          =pow(BoxEdgeLengthInNm,3); //nm3
const double CutoffClusterInNm   = 0.38;//cutoff for cluster criteria
const double Eps                 = 0.02;

const double R     = 8.3144598;             //ideal gas constant in J/(K*Mol)
const double Beta  = 1.0/(Temperature*R);   //reciprocal temperature 
const double Mu    = -700;             //chemical potential in kilojoule/mole
const double CationWavelength = 0.0210870;   //thermal wavelength of cation in nms
const double AnionWavelength = 0.0169807;    //thermal wavelength of anion in nms
//const double WaveLengthCube = pow(WaveLength,3.0);


//                            FORCE FIELD DATA
const double Coulomb14Scale      = 0.5;
const double LennardJones14Scale = 0.5;

//const double CationCharge        = 0.5000;     //charge for cation 
const double CationEpsilon       = 0.4184;//epsilon for cation in kJ/mol
const double CationSigma         = 0.2584;  //sigma for cation in nanometers
const double CationCharge        = 1;     //charge for cation 
//const double CationEpsilon       = 0.4184;//epsilon for cation in kJ/mol


//const double AnionCharge         = -0.5000;    //charge for anion 
const double AnionEpsilon        = 0.4184; //epsilon for anion in kJ/mol
const double AnionSigma          = 0.4036;  //sigma for anion in nanometers
const double AnionCharge         = -1;    //charge for anion 
//const double AnionEpsilon        = 0.4184; //epsilon for anion in kJ/mol


const double OHDistance          = 0.1; 
const double HHDistance          = 0.16329809;

//                             GCMC PARAMETERS
const int    M     = 10;                    //number of intermediate states + 1
const int    UpperLimit = 8;               //maximum num of atoms
const int    LowerLimit = 1;                //minimum num of atoms
const int    NW = 6909;                         // number of water atoms
//const int    NW = 930;                         // number of water atoms
const int    MdSteps = 1000;                  //md steps
const int    NumGcmc = 40000000;

//                            other related information
const std::string   CationRealName = "Na";
const std::string   CationGhostName = "Nag";
const std::string   CationIntermediateName = "Nai";

const std::string   AnionRealName = "Cl";
const std::string   AnionGhostName = "Clg";
const std::string   AnionIntermediateName = "Cli";


struct AtomType {
    double mass, charge, sigmaInNms, epsilonInkJPerMol;
} atomType[] = {/*0 Na*/ 22.989769282, 1.00, 0.2584, 0.4184,
                /*1 Cl*/ 35.4532, -1.00, 0.4036, 0.4184,
                /*2 O*/ 15.99943, -0.8476, 0.31657195050398818,0.6497752,
                /*3 H*/ 1.007947, 0.4238, 0.065, 0.16628};
//struct AtomType {
//    double mass, charge, sigmaInNms, epsilonInkJPerMol;
//} atomType[] = {/*0 Li*/ 22.989769282, 1.00, 0.2584, 0.4184,
//                /*1 F*/ 35.4532,-1.00, 0.4036, 0.4184,
//                /*2 O*/ 15.99943, -0.8476, 0.31657195050398818,0.6497752,
//                /*3 H*/ 1.007947, 0.4238, 0.065, 0.16628};
const int Na = 0;
const int Cl = 1;
const int OW = 2;
const int HW1 = 3;
const int HW2 = 3;

//                                MOLECULE DATA
const int EndOfList=-1;

struct MyAtomInfo
{   int type; char* pdb; double initPosInAng[3]; double posInAng[3];} 
atoms[UpperLimit*2+NW+1]; 

// -----------------------------------------------------------------------------
//                           INTERFACE TO OpenMM
// -----------------------------------------------------------------------------

struct MyOpenMMData;
static MyOpenMMData* 
myInitializeOpenMM( const MyAtomInfo    atoms[],
                    double              temperature,
                    double              frictionInPerPs,
                    double              stepSizeInFs,
                    double              boxEdgeLengthInNm, 
                    std::string&        platformName); 
static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
static void          myGetOpenMMState(MyOpenMMData*, bool wantEnergy,
                                      double& energy, 
                                      MyAtomInfo atoms[]);
static void          myTerminateOpenMM(MyOpenMMData*);
static void
myWritePDBFrame(GcmcIonPair& gcmc, int numReal, FILE* pdbOutput,
                 const MyAtomInfo atoms[]);


//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
static void
myWritePDBFrame(GcmcIonPair& gcmc, int numReal, FILE* pdbOutput,
                const MyAtomInfo atoms[]) 
{
    const std::vector<OpenMM::Vec3> posInNm = gcmc.simulation->context->getState(OpenMM::State::Positions).getPositions();
    double a = pow(gcmc.simulation->context->getState(0).getPeriodicBoxVolume(),1.0/3.0)*10;
    // Write out in PDB format -- printf is so much more compact than formatted cout.
    fprintf(pdbOutput,"MODEL    %d\n",numReal);
    fprintf(pdbOutput,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",a,a,a,90.0,90.0,90.0);
    for (int i=0; i < numReal ; ++i){
        int n = gcmc.cationRealList[i];
        fprintf(pdbOutput,"%-6s%5d  %-3s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            "ATOM",2*i+1, atoms[n].pdb, "Li+","A",2*i+1, posInNm[n][0]*10, posInNm[n][1]*10, posInNm[n][2]*10,1.00,0.00);
        n = gcmc.anionRealList[i];
        fprintf(pdbOutput,"%-6s%5d  %-3s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            "ATOM",2*i+2, atoms[n].pdb, " F-","A",2*i+2, posInNm[n][0]*10, posInNm[n][1]*10, posInNm[n][2]*10,1.00,0.00);
    }
    for (int n=UpperLimit*2;atoms[n].type != EndOfList ; ++n)
        fprintf(pdbOutput,"%-6s%5d  %-3s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            "ATOM",numReal*2 + (n-UpperLimit*2) + 1, atoms[n].pdb, "HOH","A",numReal*2 + (n-UpperLimit*2)/3 + 1, posInNm[n][0]*10, posInNm[n][1]*10, posInNm[n][2]*10,1.00,0.00);
    fprintf(pdbOutput,"ENDMDL\n");
    fprintf(pdbOutput,"   \n");
}

static void
myWriteGROFrame(GcmcIonPair& gcmc, int numReal, FILE* groOutput,
                const MyAtomInfo atoms[])
{
    const std::vector<OpenMM::Vec3> posInNm = gcmc.simulation->context->getState(OpenMM::State::Positions).getPositions();
    double a = pow(gcmc.simulation->context->getState(0).getPeriodicBoxVolume(),1.0/3.0);
    double timeInPs = gcmc.simulation->context->getState(OpenMM::State::Energy).getTime();
    // no velocity
    fprintf(groOutput, "Time = %.3f ps\n", timeInPs);
    fprintf(groOutput, "%5d\n", UpperLimit*2);

    const char* format = "%8.3f%8.3f%8.3f\n";
    for (int i=0; i < numReal; ++i){
    // keeping the index of atoms unchanged, unlike how GCMC codes output PDB
        fprintf(groOutput, "%5d%-5.5s%5.5s%5d", (2*i+1) % 100000, "Na+", atoms[2*i].pdb, (2*i+1) % 100000);
        fprintf(groOutput, format, posInNm[2*i][0], posInNm[2*i][1], posInNm[2*i][2]);

        fprintf(groOutput, "%5d%-5.5s%5.5s%5d", (2*i+2) % 100000, "Cl-", atoms[2*i+1].pdb, (2*i+2) % 100000);
        fprintf(groOutput, format, posInNm[2*i+1][0], posInNm[2*i+1][1], posInNm[2*i+1][2]);

    }
    // only write down coordinates of ions
    //for (int n=UpperLimit*2;atoms[n].type != EndOfList ; ++n){
    //    fprintf(groOutput, "%5d%-5.5s%5.5s%5d", numReal*2 + (n-UpperLimit*2)/3 + 1, "HOH", atoms[n].pdb, (n+1) % 100000);
    //    fprintf(groOutput, format, posInNm[n][0], posInNm[n][1], posInNm[n][2]);
    //}
    fprintf(groOutput, "%10.5f %9.5f %9.5f\n", a, a, a);
}

// -----------------------------------------------------------------------------
//                           OpenMM-USING CODE
// -----------------------------------------------------------------------------
// The OpenMM API is visible only at this point and below. Normally this would
// be in a separate compilation module; we're including it here for simplicity.
// -----------------------------------------------------------------------------

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
#endif

using OpenMM::Vec3; // so we can just say "Vec3" below



// -----------------------------------------------------------------------------
//                           LiF MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        srand (time(NULL)); 
   
        std::string   platformName;
        //reading pdb file
        std::ifstream pdb("MD150.pdb");
        std::string line;
        for (int i=0; i<4; i++)
            std::getline(pdb, line);
        for (int i=0; i<UpperLimit*2+NW; i++){
            std::getline(pdb,line);
            std::istringstream iss(line);
            std::vector<std::string> atomInfo(std::istream_iterator<std::string>{iss}, 
                                              std::istream_iterator<std::string>());
            if (atomInfo[2] == "OW")
                atoms[i].type = OW;
            else if (atomInfo[2] == "HW1" || atomInfo[2] == "HW2")
                atoms[i].type = HW1;
            else if (atomInfo[2] == "Na")
                atoms[i].type = Na;
            else if (atomInfo[2] == "Cl")
                atoms[i].type = Cl;
            else
                std::cout << line << std::endl;
            atoms[i].pdb = new char[3];
            strcpy(atoms[i].pdb, atomInfo[2].c_str());
            atoms[i].initPosInAng[0] = std::stod(atomInfo[atomInfo.size()-5]);
            atoms[i].initPosInAng[1] = std::stod(atomInfo[atomInfo.size()-4]);
            atoms[i].initPosInAng[2] = std::stod(atomInfo[atomInfo.size()-3]);
            //std::cout << atoms[i].type<<" "<<atoms[i].initPosInAng[0] << " " 
            //<< atoms[i].initPosInAng[1] << " "
            //<< atoms[i].initPosInAng[2] << std::endl;
            atoms[i].posInAng[0] = 0;
            atoms[i].posInAng[1] = 0;
            atoms[i].posInAng[2] = 0;
        }
        pdb.close();
        atoms[UpperLimit*2+NW].type = EndOfList;

        // Set up OpenMM data structures; returns OpenMM Platform name.
        MyOpenMMData* omm = myInitializeOpenMM(atoms, Temperature, 
                            FrictionInPerPs, StepSizeInFs, 
                            BoxEdgeLengthInNm, platformName);

        printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
        // Set up GCMC initial condition
        // (1)set up real atoms and ghost atoms
        // (2)set up nonbonded force for real atoms and ghost atoms
        // (3)set up custombond force for all ion pairs 
        int numReal = UpperLimit;

        std::vector<int> cationRealList, cationGhostList, 
                         anionRealList, anionGhostList;
        // Set one ion pair as interacted ion pairs
        for (int i = 0; i < numReal; i++){
            cationRealList.push_back(2*i);
            anionRealList.push_back(2*i+1);
        }
        for (int i = numReal; i < UpperLimit; i++){
            cationGhostList.push_back(2*i);
            anionGhostList.push_back(2*i+1);
            omm->nonbond->setParticleParameters(2*i,0.0,CationSigma,0.0);
            omm->nonbond->setParticleParameters(2*i+1,0.0,AnionSigma,0.0);
        }
        omm->nonbond->updateParametersInContext(*omm->context);
        std::cout.precision(20);
 
        std::vector<int> cationList, anionList;
        for (int i = 0; atoms[i].type != EndOfList; i++){
            if (atoms[i].type == Na)
                cationList.push_back(i);
            else if (atoms[i].type == Cl)
                anionList.push_back(i);
            else
                ;
        } 

        int numGcmc = 0;
        GcmcParameters parameters;
        parameters.temperature     = Temperature;
        parameters.beta            = Beta;
        parameters.mu              = Mu;
        parameters.boxEdgeLength   = BoxEdgeLengthInNm;
        parameters.cutoffCluster   = CutoffClusterInNm;
        parameters.eps             = Eps;
        parameters.cationWavelengthCube = pow(CationWavelength,3);
        parameters.anionWavelengthCube = pow(AnionWavelength,3);
        parameters.cationCharge    = CationCharge;
        parameters.cationSigma     = CationSigma;
        parameters.cationEpsilon   = CationEpsilon;
        parameters.anionCharge     = AnionCharge;
        parameters.anionSigma      = AnionSigma;
        parameters.anionEpsilon    = AnionEpsilon;
        parameters.M               = M;                     
        parameters.upperLimit      = UpperLimit;               
        parameters.lowerLimit      = LowerLimit;                
        parameters.mdSteps         = MdSteps;    
        parameters.numReal         = numReal;
        parameters.numGhost        = UpperLimit - numReal;
        parameters.cationRealList  = cationRealList;
        parameters.cationGhostList = cationGhostList;
        parameters.anionRealList   = anionRealList;
        parameters.anionGhostList  = anionGhostList;
        parameters.omm             = omm;

        GcmcIonPair gcmc(&parameters);

        const std::vector<OpenMM::Vec3> velocities = gcmc.simulation->context->getState(OpenMM::State::Velocities).getVelocities();
        printf("Initial velocities:\n");
        for (int i=0; i<UpperLimit*2; ++i) {
            printf("%8.3f%8.3f%8.3f\n", velocities[i][0], velocities[i][1], velocities[i][2]);
        }
        
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        //cudaProfilerStart();
        
        //production
        FILE *groOutput;
        groOutput = fopen("out.gro", "w");
        gcmc.f = 0.0;
        int count = 0;
        double potEnergyInKJ;
        std::ofstream write;
        int WriteNucl8 = 0;
        int NumEquRun = 100000; // equ = 400 ps, each step = 2*2 = 4ps
        int NumMDRun = 15000000; // MD = 60 ns
        myWriteGROFrame(gcmc, gcmc.numReal, groOutput, atoms);
        for (int i = 0; i < NumEquRun; i++){
            gcmc.step();
        }
        for (int i = 0; i < NumMDRun; i++){
            gcmc.step();
            if (i % 250 == 0){ // every 1 ps
                double cbond_kJmol = gcmc.simulation->context->getState(OpenMM::State::Energy, false, 4).getPotentialEnergy();
                if (cbond_kJmol == 0) myWriteGROFrame(gcmc, gcmc.numReal, groOutput, atoms);
            }
        }
 
        //cudaProfilerStop();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << elapsed << std::endl;
        // Clean up OpenMM data structures.

        fclose(groOutput); 
	// delete &gcmc;
        //myTerminateOpenMM(omm);

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}

// -----------------------------------------------------------------------------
//                      INITIALIZE OpenMM DATA STRUCTURES
// -----------------------------------------------------------------------------
// We take these actions here:
// (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
// (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
//     in a manner which is opaque to the caller.
// (3) Fill the OpenMM::System with the force field parameters we want to
//     use and the particular set of atoms to be simulated.
// (4) Create an Integrator and a Context associating the Integrator with
//     the System.
// (5) Select the OpenMM platform to be used.
// (6) Return the MyOpenMMData struct and the name of the Platform in use.
//
// Note that this function must understand the calling MD code's molecule and
// force field data structures so will need to be customized for each MD code.
static MyOpenMMData* 
myInitializeOpenMM( const MyAtomInfo    atoms[],
                    double              temperature,
                    double              frictionInPerPs,
                    double              stepSizeInFs,
                    double              boxEdgeLengthInNm, 
                    std::string&        platformName) 
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
       ("/home/lzhang657/Schmidt_work/Nucleation_Ajay/GCMC/openmm/lib/plugins");
    OpenMM::Platform& platform = OpenMM::Platform::findPlatform({});
  
    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
    OpenMM::NonbondedForce&         nonbond     = *(omm->nonbond = new OpenMM::NonbondedForce());
    OpenMM::CustomBondForce&        custombond  = *(omm->custombond = 
                                    new OpenMM::CustomBondForce("lamda*k*max(0,r-rmax)^2"));
    OpenMM::CMMotionRemover&        cmmr = *(omm->cmmr = new OpenMM::CMMotionRemover());
    OpenMM::MonteCarloBarostat&     barostat = *(omm->barostat = new OpenMM::MonteCarloBarostat(1.01325,298.15,25));
    system.addForce(&nonbond);
    nonbond.setForceGroup(0);
    system.addForce(&cmmr);
    cmmr.setForceGroup(1);
    system.addForce(&custombond);
    custombond.setForceGroup(2);
    system.addForce(&barostat);
    barostat.setForceGroup(3);
    //OpenMM::AndersenThermostat&     thermostat  = *new OpenMM::AndersenThermostat(
    //        temperature,      // kelvins
    //        frictionInPerPs); // collision frequency in 1/picoseconds
    //system.addForce(&thermostat);
    //thermostat.setForceGroup(3);
    
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    //  (4) Add custombond interaction as restraint
    std::vector<Vec3> initialPosInNm;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atype.mass);
        nonbond.addParticle(atype.charge,
                            atype.sigmaInNms, 
                            atype.epsilonInkJPerMol); 
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }
    
    std::vector< std::pair<int,int> >   bondPairs;

    //use rigid water model
    for (int i = 0; atoms[i].type != EndOfList; ++i) {
        if (atoms[i].type == OW){
            system.addConstraint(i,i+1,OHDistance);
            system.addConstraint(i,i+2,OHDistance);    
            system.addConstraint(i+1,i+2,HHDistance);
//            std::cout << i << " " << i+1 << " " << i+2 << std::endl;
            bondPairs.push_back(std::make_pair(i, i+1));
            bondPairs.push_back(std::make_pair(i, i+2));
        }
   }
    // Populate nonbonded exclusions
    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    //set custom bond force for constraints
    custombond.addGlobalParameter("k",100000);
    custombond.addGlobalParameter("rmax",CutoffClusterInNm);
    custombond.addPerBondParameter("lamda"); 
    std::vector<int> cationList, anionList;
    for (int i = 0; atoms[i].type != EndOfList; i++){
        if (atoms[i].type == Na)
            cationList.push_back(i);
        else if (atoms[i].type == Cl)
            anionList.push_back(i);
        else
            ;
    } 
    for (int i = 0; i < cationList.size(); i++)
        for (int j = 0; j < anionList.size(); j++)
            custombond.addBond(cationList[i],anionList[j],{0.0});
    custombond.setUsesPeriodicBoundaryConditions(true); 
    //create periodic box
    nonbond.setNonbondedMethod(OpenMM::NonbondedForce::PME);
    nonbond.setEwaldErrorTolerance(0.00013); // correspond to n_mesh = about 32
    //nonbond.setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
    nonbond.setCutoffDistance(CutoffDistanceInAng * OpenMM::NmPerAngstrom);
    nonbond.setUseDispersionCorrection(false);
    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLengthInNm,0,0),
                                        Vec3(0,boxEdgeLengthInNm,0), 
                                        Vec3(0,0,boxEdgeLengthInNm));        

    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    omm->integrator = new OpenMM::LangevinIntegrator(temperature, frictionInPerPs, 
                                                     stepSizeInFs * OpenMM::PsPerFs);
    omm->integrator->setConstraintTolerance(0.00001);
    //omm->integrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform, {{"Precision","mixed"}});
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocitiesToTemperature(Temperature,rand());
    platformName = omm->context->getPlatform().getName();
    //std::cout << omm->context->getPlatform().getPropertyValue(*omm->context, "Precision") << std::endl;
    return omm;
}


// -----------------------------------------------------------------------------
//                     COPY STATE BACK TO CPU FROM OPENMM
// -----------------------------------------------------------------------------
static void
myGetOpenMMState(MyOpenMMData* omm, bool wantEnergy, 
                 double& energyInKJ,
                 MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
        infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    // Forces are also available (and cheap).

    const OpenMM::State state = omm->context->getState(infoMask);
    //timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    //const std::vector<Vec3>& positionsInNm = state.getPositions();
    //for (int i=0; i < (int)positionsInNm.size(); ++i)
    //    for (int j=0; j < 3; ++j)
    //        atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;

    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKJ = 0;
    //std::cout << state.getPotentialEnergy() << "  " << state.getKineticEnergy() << std::endl;
    if (wantEnergy) 
        energyInKJ = state.getPotentialEnergy();
}


// -----------------------------------------------------------------------------
//                     TAKE MULTIPLE STEPS USING OpenMM 
// -----------------------------------------------------------------------------
static void 
myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
    omm->integrator->step(numSteps);
}

// -----------------------------------------------------------------------------
//                     DEALLOCATE OpenMM OBJECTS
// -----------------------------------------------------------------------------
static void 
myTerminateOpenMM(MyOpenMMData* omm) {
    delete omm;
}





