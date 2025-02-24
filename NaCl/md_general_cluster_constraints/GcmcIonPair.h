#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#ifndef _GCMCIONPAIR_H_
#define _GCMCIONPAIR_H_

typedef boost::adjacency_list<boost::vecS, boost::vecS, 
        boost::undirectedS,
        boost::property<boost::vertex_name_t,int>, 
        boost::property<boost::edge_weight_t, double>> Graph;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex; 
typedef boost::property_map<Graph, boost::vertex_name_t>::type NameMap;

class GcmcIonPair
{
public:
    //members
    //MODELING AND SIMULATION PARAMETERS
    double temp;                 
    double beta;                 
    double mu;                   
    double boxEdgeLength;        
    double volume;               
    double cutoffCluster;
    double eps; 
    double cationWavelengthCube; 
    double anionWavelengthCube;      
    double wavelengthCube;

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
    int    numReal;              
    int    numGhost;
    int    cationNumReal;
    int    anionNumReal;
    int    cationNumGhost;
    int    anionNumGhost;             
    std::vector<int> cationRealList;  
    std::vector<int> cationGhostList; 
    std::vector<int> anionRealList;   
    std::vector<int> anionGhostList; 

    MyOpenMMData *simulation;    

    //Other Parameters
    int    cationRandIndex; 
    int    anionRandIndex;  
    int    cationRefIndex;  
    int    anionRefIndex;   
    int    numInsr;    
    int    numDel;     
    int    numMD;
    int    sucInsr;    
    int    sucDel;
    int    sucMD;     
    int    m;               
    int    type; 
    double f;               
    int    cationNIn;       
    int    anionNIn;        
    double vIn;             
    std::vector<OpenMM::Vec3> positions; 
    std::vector<OpenMM::Vec3> oldPositions; 
    std::vector<OpenMM::Vec3> velocities; 
    std::vector<OpenMM::Vec3> oldVelocities; 
    std::vector<OpenMM::Vec3> newVelocities;

    std::normal_distribution<double> gauss;
    std::default_random_engine generator;

    Graph  g;
    Graph g_plus_constraints;
    std::vector<Edge> spanning_tree;
    std::set<std::pair<int, int>> constraints, constraintsToBeRemoved;
    std::vector<std::pair<int, int>> edgesToBeChecked;
    bool gConnectivity;
    NameMap nameMap;
  
    double *eta;
    double *lamda;
    int    *histogram;
    double *energy;
    double *freeEnergy;
    Vertex *index;

    std::random_device rd;
    std::mt19937 gen;

public:
    //constructors   
    GcmcIonPair(GcmcParameters *parameters); 
    ~GcmcIonPair();

public: 
    //functions
    double getDistance(std::vector<OpenMM::Vec3>& positions, int atomA, int atomB);

    void generateGraph(std::vector<OpenMM::Vec3>& positions, Graph& g);

    void generateGraph2(std::vector<OpenMM::Vec3>& positions, Graph& g, std::vector<std::pair<int, int>>& edgesToBeChecked);

    bool checkClusterCriteria(Graph& g);

    void generateMST(Graph&g, std::vector<Edge> &spanning_tree);

    void updateForce(int type, double scaling);

    void addRandRefRestraint(int type);
  
    void removeRandRefRestraint(int type);

    void addMSTRestraints(Graph &g,std::vector<Edge> &spanning_tree);

    void removeMSTRestraints(Graph &g,std::vector<Edge> &spanning_tree);

    // set random positions for insertion
    void setRandomPositions(int type);

    // pick randome atom for deletion
    bool getRandomAtoms(int type);

    // set random velocity for new atom according to Boltzmann distribution
    void setRandomVelocities(int type);

    void insertion();

    void deletion();

    void MDstep();

    void printEnergyComponents();

    void printGraphEdges(Graph& g);

    void printEdgeSet(std::set<std::pair<int, int>>& edgeSet);

    void printEdgeVector(std::vector<std::pair<int, int>>& edgeVec);

    void step();

    double getInsrProbability(); 

    double getDelProbability();
    
    double getMDProbability();

    void rezeroStat();  

    void rezeroHistogram();
   
    void rezeroEnergy();

    void updateWangLaudauFactor();

    // This funtion gets average energy for each state, 
    // it's useless unless for debugging 
    double* getEnergy();

    double* getFreeEnergy();  
};

#endif


