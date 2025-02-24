#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/copy.hpp>
#include "Interface.h"
#include "GcmcIonPair.h"
#include "OpenMM.h"

typedef boost::adjacency_list<boost::vecS, boost::vecS, 
        boost::undirectedS,
        boost::property<boost::vertex_name_t,int>, 
        boost::property<boost::edge_weight_t, double>> Graph;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex; 
typedef boost::property_map<Graph, boost::vertex_name_t>::type NameMap;

GcmcIonPair::GcmcIonPair(GcmcParameters *parameters){

    //MODELING AND SIMULATION PARAMETERS
    temp                 = parameters->temperature;
    beta                 = parameters->beta;             
    mu                   = parameters->mu;              // chemical potential
    boxEdgeLength        = parameters->boxEdgeLength;
    volume               = pow(boxEdgeLength,3);
    cutoffCluster        = parameters->cutoffCluster;
    eps                  = parameters->eps;
    cationWavelengthCube = parameters->cationWavelengthCube;// de brogile wavelength for cation**3
    anionWavelengthCube  = parameters->anionWavelengthCube; // de brogile wavelength for antion**3
    wavelengthCube       = sqrt(cationWavelengthCube*anionWavelengthCube);

    //FORCE FIELD DATA
    cationCharge         = parameters->cationCharge; 
    cationSigma          = parameters->cationSigma; 
    cationEpsilon        = parameters->cationEpsilon;
    anionCharge          = parameters->anionCharge;  
    anionSigma           = parameters->anionSigma;  
    anionEpsilon         = parameters->anionEpsilon;

    //GCMC PARAMETERS
    M                    = parameters->M;                  // number of intermediates 
    upperLimit           = parameters->upperLimit;    
    lowerLimit           = parameters->lowerLimit;
    mdSteps              = parameters->mdSteps;            // MD steps after each GCMC.
                
    //GCMC Initial State Parameters
    numReal              = parameters->numReal;            
    numGhost             = parameters->numGhost;           
    cationNumReal        = numReal;                         
    anionNumReal         = numReal;
    cationNumGhost       = numGhost;
    anionNumGhost        = numGhost;
    cationRealList  = parameters->cationRealList;
    cationGhostList = parameters->cationGhostList;
    anionRealList   = parameters->anionRealList;
    anionGhostList = parameters->anionGhostList;

    simulation    = parameters->omm;

    //Other Parameters
    cationRandIndex      = 1;
    anionRandIndex       = 0;
    cationRefIndex       = 1;
    anionRefIndex        = 0;
    numInsr              = 0;
    numDel               = 0;
    numMD                = 0;
    sucMD                = 0;
    sucInsr              = 0;
    sucDel               = 0;
    m                    = 0;  //index for intermediate state
    type                 = 0;  //denote the insertion/deletion related to cation(0) or anion(1) 
    f                    = 1.0/8192.0;//wang-laudau initial modification factor
    cationNIn            = 0.0;//coordination number of cation
    anionNIn             = 0.0;//coordination number of anion 
    vIn                  = 4.0/3.0*M_PI*pow(cutoffCluster,3);//excluded volume
 
    positions = 
        simulation->context->getState(OpenMM::State::Positions).getPositions();
    oldPositions = 
        simulation->context->getState(OpenMM::State::Positions).getPositions();
    velocities = 
        simulation->context->getState(OpenMM::State::Velocities).getVelocities();
    oldVelocities = 
        simulation->context->getState(OpenMM::State::Velocities).getVelocities();
    newVelocities = 
        simulation->context->getState(OpenMM::State::Velocities).getVelocities();
    gauss = std::normal_distribution<double>(0.0,1.0);     
    generator = std::default_random_engine(time(NULL));  // normal distribution random generator

    g = Graph();
    g_plus_constraints = Graph();                        // trail graph 
    spanning_tree = std::vector<Edge>();
    constraints = std::set<std::pair<int, int>>();
    edgesToBeChecked = std::vector<std::pair<int, int>>();
    constraintsToBeRemoved = std::set<std::pair<int, int>>();
    gConnectivity = false;
    NameMap nameMap;                        // name map for finding the index of the atom in the graph

    eta = new double[upperLimit*M+1];       // wang-landau bias 
    lamda = new double[M/2+1];              // scaling parameter for intermediate states
    histogram = new int[upperLimit*M + 1]; 
    energy = new double[upperLimit*M + 1];
    freeEnergy = new double[upperLimit*M + 1];
    index = new Vertex[upperLimit*2];       // array to store the vertex in the graph for each atom
                                            // opposite of NameMap

    for (int i = 0; i < M/2; i++)
        lamda[i] = i*i*1.0/(M*M/4.0);
    lamda[M/2] = 1.0;
 
    for (int i = 0; i < upperLimit*M + 1; i++){
        histogram[i] = 0;
        eta[i] = 0.0;
        freeEnergy[i] = 0.0;
        energy[i] = 0.0;
    }
   
    for (int i = 0; i < upperLimit*2 ; i++)
        index[i] = 0;

    gen.seed(rd());
}

GcmcIonPair::~GcmcIonPair(){
    delete eta;
    delete lamda;
    delete histogram;
    delete energy;
    delete freeEnergy;
    delete index;
}

double GcmcIonPair::getDistance(std::vector<OpenMM::Vec3>& positions, int atomA, int atomB){
    double dSq = 0.0;
    double d = 0.0;
    OpenMM::Vec3 dVec = positions[atomA] - positions[atomB];
    dSq += pow((dVec[0]-boxEdgeLength*round(dVec[0]/boxEdgeLength)),2);
    dSq += pow((dVec[1]-boxEdgeLength*round(dVec[1]/boxEdgeLength)),2);
    dSq += pow((dVec[2]-boxEdgeLength*round(dVec[2]/boxEdgeLength)),2);
    d = sqrt(dSq);
    return d;
}

void GcmcIonPair::generateGraph(std::vector<OpenMM::Vec3>& positions, Graph& g){
    int count = 0;
    int cationSize = cationRealList.size();
    int anionSize = anionRealList.size();
    bool cluster = false;
    g.clear();
    for (int i = 0; i < cationSize; i++){
        Vertex v = boost::add_vertex(cationRealList[i],g);
        index[cationRealList[i]] = v;
    }
    for (int i = 0; i < anionSize; i++){
        Vertex v = boost::add_vertex(anionRealList[i],g);
        index[anionRealList[i]] = v;
    }
    for (int i = 0; i < cationSize; i++)
        for (int j = 0; j < anionSize; j++){
            double d = getDistance(positions,cationRealList[i],anionRealList[j]);
            if (d < cutoffCluster)
                boost::add_edge(index[cationRealList[i]],index[anionRealList[j]],d,g);
        }
    nameMap = (boost::get(boost::vertex_name, g));
}

void GcmcIonPair::generateGraph2(std::vector<OpenMM::Vec3>& positions, Graph& g, std::vector<std::pair<int, int>>& edgesToBeChecked){
    int count = 0;
    int cationSize = cationRealList.size();
    int anionSize = anionRealList.size();
    bool cluster = false;
    g.clear();
    edgesToBeChecked.clear();
    for (int i = 0; i < cationSize; i++){
        Vertex v = boost::add_vertex(cationRealList[i],g);
        index[cationRealList[i]] = v;
    }
    for (int i = 0; i < anionSize; i++){
        Vertex v = boost::add_vertex(anionRealList[i],g);
        index[anionRealList[i]] = v;
    }
    for (int i = 0; i < cationSize; i++)
        for (int j = 0; j < anionSize; j++){
            double d = getDistance(positions,cationRealList[i],anionRealList[j]);
            if (d < cutoffCluster) {
                boost::add_edge(index[cationRealList[i]],index[anionRealList[j]],d,g);
                if (d > cutoffCluster -eps) edgesToBeChecked.push_back(std::make_pair(cationRealList[i],anionRealList[j]));
            }
        }
    nameMap = (boost::get(boost::vertex_name, g));
}

bool GcmcIonPair::checkClusterCriteria(Graph& g){
    std::vector<int> components(boost::num_vertices(g));
    int num = boost::connected_components(g,&components[0]);
    if (num == 1)
        return true;
    else
        return false;
}
     
void GcmcIonPair::generateMST(Graph&g, std::vector<Edge> &spanning_tree){
    spanning_tree.clear();
    boost::kruskal_minimum_spanning_tree(g,std::back_inserter(spanning_tree));
    std::cout << "MST edges " << std::endl;
    for (std::vector <Edge>::iterator ei = spanning_tree.begin();
        ei != spanning_tree.end(); ++ei) {
        int a = nameMap[boost::source(*ei, g)];
        int b = nameMap[boost::target(*ei, g)];
        if (a < b)
        std::cout << a << "," << b << std::endl;
        else
        std::cout << b << "," << a << std::endl;
    }
} 
     
void GcmcIonPair::updateForce(int type, double scaling){
    if (type == 0)
        simulation->nonbond->setParticleParameters(cationRandIndex,cationCharge*scaling,
            cationSigma,cationEpsilon*pow(scaling,2.0));
    else
        simulation->nonbond->setParticleParameters(anionRandIndex,anionCharge*scaling,
            anionSigma,anionEpsilon*pow(scaling,2.0));
    simulation->nonbond->updateParametersInContext(*(simulation->context));
}

void GcmcIonPair::addRandRefRestraint(int type){
    if (type == 0)
        simulation->custombond->setBondParameters(cationRandIndex/2*upperLimit+(anionRefIndex-1)/2,
            cationRandIndex,anionRefIndex,{1.0});
    else
        simulation->custombond->setBondParameters(cationRefIndex/2*upperLimit+(anionRandIndex-1)/2,
            cationRefIndex,anionRandIndex,{1.0});
    simulation->custombond->updateParametersInContext(*(simulation->context));
}

void GcmcIonPair::removeRandRefRestraint(int type){
    if (type == 0)
        simulation->custombond->setBondParameters(cationRandIndex/2*upperLimit+(anionRefIndex-1)/2,
            cationRandIndex,anionRefIndex,{0.0});
    else
        simulation->custombond->setBondParameters(cationRefIndex/2*upperLimit+(anionRandIndex-1)/2,
            cationRefIndex,anionRandIndex,{0.0});
    simulation->custombond->updateParametersInContext(*(simulation->context));
}

void GcmcIonPair::addMSTRestraints(Graph &g,std::vector<Edge> &spanning_tree){
    for (std::vector <Edge>::iterator ei = spanning_tree.begin();
        ei != spanning_tree.end(); ++ei) {
        int a = nameMap[boost::source(*ei, g)];
        int b = nameMap[boost::target(*ei, g)];
        if (a%2 == 0)
            simulation->custombond->setBondParameters(a/2*upperLimit+(b-1)/2,a,b,{1.0});
        else
            simulation->custombond->setBondParameters(b/2*upperLimit+(a-1)/2,b,a,{1.0});
    }
    simulation->custombond->updateParametersInContext(*(simulation->context));
}

void GcmcIonPair::removeMSTRestraints(Graph &g,std::vector<Edge> &spanning_tree){
    for (std::vector <Edge>::iterator ei = spanning_tree.begin();
        ei != spanning_tree.end(); ++ei) {
        int a = nameMap[boost::source(*ei, g)];
        int b = nameMap[boost::target(*ei, g)];
        if (a%2 == 0)
            simulation->custombond->setBondParameters(a/2*upperLimit+(b-1)/2,a,b,{0.0});
        else
            simulation->custombond->setBondParameters(b/2*upperLimit+(a-1)/2,b,a,{0.0});
    }
    simulation->custombond->updateParametersInContext(*(simulation->context));
}

void GcmcIonPair::setRandomPositions(int type){
    generateGraph(positions,g);
    generateMST(g,spanning_tree);
    if (type == 0){
        anionRefIndex = anionRealList[rand()%anionNumReal];
        std::pair<boost::graph_traits<Graph>::adjacency_iterator,
            boost::graph_traits<Graph>::adjacency_iterator> anionCoorList; 
        anionCoorList = boost::adjacent_vertices(index[anionRefIndex],g);
        cationNIn = std::distance(anionCoorList.first,anionCoorList.second);
        OpenMM::Vec3 vec;
        for (int i = 0; i < 3; i++)
            vec[i] = gauss(generator);
        double mag = sqrt(pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2));
        double length = pow(rand()/double(RAND_MAX),1.0/3.0); 
        for (int i = 0; i < 3; i++)
            vec[i] = vec[i]/mag * length;
        positions[cationRandIndex] = positions[anionRefIndex]+vec*cutoffCluster;
        simulation->context->setPositions(positions);
    } else {
        cationRefIndex = cationRealList[rand()%cationNumReal];
        std::pair<boost::graph_traits<Graph>::adjacency_iterator,
            boost::graph_traits<Graph>::adjacency_iterator> cationCoorList; 
        cationCoorList = boost::adjacent_vertices(index[cationRefIndex],g);
        anionNIn = std::distance(cationCoorList.first,cationCoorList.second); 
        OpenMM::Vec3 vec;
        for (int i = 0; i < 3; i++)
            vec[i] = gauss(generator);
        double mag = sqrt(pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2));
        double length = pow(rand()/double(RAND_MAX),1.0/3.0); 
        for (int i = 0; i < 3; i++)
            vec[i] = vec[i]/mag * length;
        positions[anionRandIndex] = positions[cationRefIndex]+vec*cutoffCluster;
        simulation->context->setPositions(positions);
    }
}

bool GcmcIonPair::getRandomAtoms(int type){
    generateGraph(positions,g);
    if (type == 0){
        anionRefIndex = anionRealList[rand()%anionNumReal];
        std::pair<boost::graph_traits<Graph>::adjacency_iterator,
            boost::graph_traits<Graph>::adjacency_iterator> anionCoorList; 
        anionCoorList = boost::adjacent_vertices(index[anionRefIndex],g);
        cationNIn = std::distance(anionCoorList.first,anionCoorList.second);
        cationRandIndex = nameMap[*(anionCoorList.first + rand()%cationNIn)];
        cationRealList.erase(std::remove(cationRealList.begin(),cationRealList.end(),
            cationRandIndex),cationRealList.end());
        boost::clear_vertex(index[cationRandIndex],g);
        boost::remove_vertex(index[cationRandIndex],g);
        // reject this step if we break the cluster criterion by deleting the atom
        if (!checkClusterCriteria(g)){
            cationRealList.push_back(cationRandIndex);
            generateGraph(positions,g);
            generateMST(g,spanning_tree);
            addMSTRestraints(g,spanning_tree);
            return false;
        }
        generateMST(g,spanning_tree);
        return true;
    } else {
        cationRefIndex = cationRealList[rand()%cationNumReal];
        std::pair<boost::graph_traits<Graph>::adjacency_iterator,
            boost::graph_traits<Graph>::adjacency_iterator> cationCoorList; 
        cationCoorList = boost::adjacent_vertices(index[cationRefIndex],g);
        anionNIn = std::distance(cationCoorList.first,cationCoorList.second); 
        anionRandIndex = nameMap[*(cationCoorList.first + rand()%anionNIn)];
        anionRealList.erase(std::remove(anionRealList.begin(),anionRealList.end(),
            anionRandIndex),anionRealList.end());
        boost::clear_vertex(index[anionRandIndex],g);
        boost::remove_vertex(index[anionRandIndex],g);
        // reject this step if we break the cluster criterion by deleting the atom
        if (!checkClusterCriteria(g)){
            anionRealList.push_back(anionRandIndex);
            generateGraph(positions,g);
            generateMST(g,spanning_tree);
            addMSTRestraints(g,spanning_tree);
            return false;
        }
        generateMST(g,spanning_tree);
        return true;
    }
}

void GcmcIonPair::setRandomVelocities(int type){
    velocities = 
        simulation->context->getState(OpenMM::State::Velocities).getVelocities();
    simulation->context->setVelocitiesToTemperature(temp,rand());
    newVelocities = 
        simulation->context->getState(OpenMM::State::Velocities).getVelocities();
    if (type == 0)
        velocities[cationRandIndex] = newVelocities[cationRandIndex];
    else 
        velocities[anionRandIndex] = newVelocities[anionRandIndex];
    simulation->context->setVelocities(velocities);
}

void GcmcIonPair::insertion(){
    if (numReal == upperLimit){
        // reject this step when cluster size reaches upper limit 
        generateGraph(positions,g);
        generateMST(g,spanning_tree);
        addMSTRestraints(g,spanning_tree);
    } else {
        numInsr++;
        type = m/(M/2);
        int n = m%(M/2);
        if (type == 0){
            if (n == 0){
                int randNumber = random()%cationNumGhost;
                cationRandIndex = cationGhostList[randNumber];
                cationGhostList.erase(cationGhostList.begin()+randNumber);
                setRandomPositions(type);
                addMSTRestraints(g,spanning_tree);
                addRandRefRestraint(type);
            }
            double oldPotential = simulation->context->getState(OpenMM::State::Energy).
                getPotentialEnergy();
            updateForce(type, lamda[n+1]);
            double newPotential = simulation->context->getState(OpenMM::State::Energy).
                getPotentialEnergy();
            double dU = newPotential - oldPotential;
            double acc = exp(-beta*1000*(dU-(lamda[n+1]-lamda[n])*mu/2.0)+eta[numReal*M+m+1]
                - eta[numReal*M+m])*pow(vIn/wavelengthCube,lamda[n+1]-lamda[n])
                *pow(cationNIn+lamda[n],lamda[n])/
                pow(cationNIn+lamda[n+1],lamda[n+1]);
            acc = std::min(acc,1.0);
            if (rand()/double(RAND_MAX) < acc){
                sucInsr++;
                if (n == 0){
                    m = 1;
                    setRandomVelocities(type);
                } else if (n == M/2-1){
                    m = M/2;
                    cationRealList.push_back(cationRandIndex);
                    cationNumReal++;
                    cationNumGhost--;
                } else
                    m++;
            } else{
                updateForce(type,lamda[n]);
                if (n == 0){
                    cationGhostList.push_back(cationRandIndex);
                    removeRandRefRestraint(type);
                }
            }
        } else {
            if (n == 0){
                int randNumber = random()%anionNumGhost;
                anionRandIndex = anionGhostList[randNumber];
                anionGhostList.erase(anionGhostList.begin()+randNumber);
                setRandomPositions(type);
                addMSTRestraints(g,spanning_tree);
                addRandRefRestraint(type);
            }
            double oldPotential = simulation->context->getState(OpenMM::State::Energy).
                getPotentialEnergy();
            updateForce(type, lamda[n+1]);
            double newPotential = simulation->context->getState(OpenMM::State::Energy).
                getPotentialEnergy();
            double dU = newPotential - oldPotential;
            double acc = exp(-beta*1000*(dU-(lamda[n+1]-lamda[n])*mu/2.0)+eta[numReal*M+m+1]
                - eta[numReal*M+m])*pow(vIn/wavelengthCube,lamda[n+1]-lamda[n])
                *pow(anionNIn+lamda[n],lamda[n])/
                pow(anionNIn+lamda[n+1],lamda[n+1]);
            acc = std::min(acc,1.0);
            if (rand()/double(RAND_MAX) < acc){
                sucInsr++;
                if (n == 0){
                    m = M/2+1;
                    setRandomVelocities(type);
                } else if (n == M/2-1){
                    m = 0;
                    anionRealList.push_back(anionRandIndex);
                    anionNumReal++;
                    anionNumGhost--;
                    numReal++;
                    numGhost--;
                } else
                    m++;
            } else{
                updateForce(type,lamda[n]);
                if (n == 0){
                    anionGhostList.push_back(anionRandIndex);
                    removeRandRefRestraint(type);
                }
            }
        }
    }
}

void GcmcIonPair::deletion(){
    if (numReal == lowerLimit && m == 0){
        // reject this step when cluster size reaches lower limit 
        generateGraph(positions,g);
        generateMST(g,spanning_tree);
        addMSTRestraints(g,spanning_tree);
    } else{
        numDel++;
        type = m/(M/2);
        int n = m%(M/2);
        if (n == 0)
            type = 1 - type;
        double oldPotential = 0.0;
        if (type == 0){
            if (n==0){
                if (!getRandomAtoms(type))
                    return;
                addMSTRestraints(g,spanning_tree);
                addRandRefRestraint(type);
                oldPotential = simulation->context->
                    getState(OpenMM::State::Energy).getPotentialEnergy();
                updateForce(type,lamda[M/2-1]);
            } else{
                oldPotential = simulation->context->
                    getState(OpenMM::State::Energy).getPotentialEnergy();
                updateForce(type,lamda[n-1]);
            }
            double newPotential = simulation->context->getState(OpenMM::State::Energy).
                getPotentialEnergy();
            double dU = newPotential - oldPotential;
            double acc = 1.0;
            if (n == 0)
                acc = exp(-beta*1000*(dU-(lamda[M/2-1]-1.0)*mu/2)+eta[numReal*M+m-1]-eta[numReal*M+m])
                    *pow(vIn/wavelengthCube,lamda[M/2-1]-1.0)
                    *pow(cationNIn,1.0)/
                    pow(cationNIn-1.0+lamda[M/2-1],lamda[M/2-1]);
            else
                acc = exp(-beta*1000*(dU-(lamda[n-1]-lamda[n])*mu/2)+eta[numReal*M+m-1]-eta[numReal*M+m])
                    *pow(vIn/wavelengthCube,lamda[n-1]-lamda[n])
                    *pow(cationNIn+lamda[n],lamda[n])/
                    pow(cationNIn+lamda[n-1],lamda[n-1]);
            acc = std::min(1.0,acc);
            if (rand()/double(RAND_MAX) < acc){
                sucDel++;
                if (n == 0){
                    m = M/2-1;
                    cationNumReal--;
                    cationNumGhost++;
                    cationNIn--;
                } else if (n == 1){
                    m = 0;
                    cationGhostList.push_back(cationRandIndex);
                    removeRandRefRestraint(type);
                } else
                    m--;
            } else{
                if (n == 0){
                    updateForce(type,1.0);
                    cationRealList.push_back(cationRandIndex);
                } else
                    updateForce(type,lamda[n]);
            }
        } else {
            if (n==0){
                if (!getRandomAtoms(type))
                    return;
                addMSTRestraints(g,spanning_tree);
                addRandRefRestraint(type);
                oldPotential = simulation->context->
                    getState(OpenMM::State::Energy).getPotentialEnergy();
                updateForce(type,lamda[M/2-1]);
            } else{
                oldPotential = simulation->context->
                    getState(OpenMM::State::Energy).getPotentialEnergy();
                updateForce(type,lamda[n-1]);
            }
            double newPotential = simulation->context->getState(OpenMM::State::Energy).
                getPotentialEnergy();
            double dU = newPotential - oldPotential;
            double acc = 1.0;
            if (n == 0)
                acc = exp(-beta*1000*(dU-(lamda[M/2-1]-1.0)*mu/2)+eta[numReal*M+m-1]-eta[numReal*M+m])
                    *pow(vIn/wavelengthCube,lamda[M/2-1]-1.0)
                    *pow(anionNIn,1.0)/
                    pow(anionNIn-1.0+lamda[M/2-1],lamda[M/2-1]);
            else
                acc = exp(-beta*1000*(dU-(lamda[n-1]-lamda[n])*mu/2)+eta[numReal*M+m-1]-eta[numReal*M+m])
                    *pow(vIn/wavelengthCube,lamda[n-1]-lamda[n])
                    *pow(anionNIn+lamda[n],lamda[n])/
                    pow(anionNIn+lamda[n-1],lamda[n-1]);
            acc = std::min(1.0,acc);
            if (rand()/double(RAND_MAX) < acc){
                sucDel++;
                if (n == 0){
                    m = M-1;
                    anionNumReal--;
                    anionNumGhost++;
                    anionNIn--;
                    numReal--;
                    numGhost++;
                } else if (n == 1){
                    m = M/2;
                    anionGhostList.push_back(anionRandIndex);
                    removeRandRefRestraint(type);
                } else
                    m--;
            } else{
                if (n == 0){
                    updateForce(type,1.0);
                    anionRealList.push_back(anionRandIndex);
                } else
                    updateForce(type,lamda[n]);
            }
        }
    }
}

void GcmcIonPair::printEnergyComponents(){
    double NB_kJmol = simulation->context->getState(OpenMM::State::Energy, false, 1).getPotentialEnergy();
    double cmmr_kJmol = simulation->context->getState(OpenMM::State::Energy, false, 2).getPotentialEnergy();
    double cbond_kJmol = simulation->context->getState(OpenMM::State::Energy, false, 4).getPotentialEnergy();
    double KE = simulation->context->getState(OpenMM::State::Energy+OpenMM::State::Velocities).getKineticEnergy();
    double PE = simulation->context->getState(OpenMM::State::Energy+OpenMM::State::Velocities).getPotentialEnergy();
    printf("Nonbonded interaction = %.4f kJ/mol, CMMotionRemover = %.4f kJ/mol, Half harmonic potential = %.4f kJ/mol\n", NB_kJmol, cmmr_kJmol, cbond_kJmol);
    printf("Potential energy = %.4f kJ/mol, Kinetic energy = %.4f kJ/mol, Total energy = %.4f kJ/mol\n", PE, KE, PE+KE);
}

void GcmcIonPair::printGraphEdges(Graph& g){
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        auto source = nameMap[boost::source(*ei, g)];
        auto target = nameMap[boost::target(*ei, g)];
        double d = getDistance(positions, source, target);
        std::cout << "(" << source << ", " << target << ") " << std::setprecision(7) << d << std::endl; // Print each edge descriptor
    }
}

void GcmcIonPair::printEdgeSet(std::set<std::pair<int, int>>& edgeSet){
    for (const auto& pair : edgeSet) {
        int a = pair.first;
        int b = pair.second;
        double d = getDistance(positions, a, b);
        std::cout << "(" << a << ", " << b << ") " << std::setprecision(7) << d << std::endl;
    }
}

void GcmcIonPair::printEdgeVector(std::vector<std::pair<int, int>>& edgeVec){
    for (int i=0; i<edgeVec.size(); ++i) {
        int a = edgeVec[i].first;
        int b = edgeVec[i].second;
        double d = getDistance(positions, a, b);
        std::cout << "(" << a << ", " << b << ") " << std::setprecision(7) << d << std::endl;
    }
}

void GcmcIonPair::MDstep(){
    simulation->integrator->step(1);
}
void GcmcIonPair::step(){
    //printEnergyComponents();
    // check whether constraints can be removed and update constraints
    positions = simulation->context->getState(OpenMM::State::Positions).getPositions();
    for (const auto& pair : constraints) {
        int a = pair.first; // atom indices
        int b = pair.second;
        double d = getDistance(positions, a, b);
        if (d < cutoffCluster - eps) {
            constraintsToBeRemoved.insert(std::make_pair(a, b));
            if (a%2 == 0)
                simulation->custombond->setBondParameters(a/2*upperLimit+(b-1)/2,a,b,{0.0});
            else
                simulation->custombond->setBondParameters(b/2*upperLimit+(a-1)/2,b,a,{0.0});
        }
    }
    for (const auto& pair : constraintsToBeRemoved) constraints.erase(pair);
    constraintsToBeRemoved.clear();

    //velocities = simulation->context->getState(OpenMM::State::Velocities).getVelocities();
    //for (int i=0; i < numReal; ++i){
    //// keeping the index of atoms unchanged, unlike how GCMC codes output PDB
    //    printf("%5d%-5.5s%5.5s%5d", (2*i+1) % 100000, "Na+", "Na", (2*i+1) % 100000);
    //    printf("%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n", positions[2*i][0], positions[2*i][1], positions[2*i][2], velocities[2*i][0], velocities[2*i][1], velocities[2*i][2]);

    //    printf("%5d%-5.5s%5.5s%5d", (2*i+2) % 100000, "Cl-", "Cl", (2*i+2) % 100000);
    //    printf("%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n", positions[2*i+1][0], positions[2*i+1][1], positions[2*i+1][2], velocities[2*i+1][0], velocities[2*i+1][1], velocities[2*i+1][2]);

    //}
    generateGraph2(positions, g, edgesToBeChecked);
    gConnectivity = checkClusterCriteria(g);
    //printf("Current graph edges:\n");
    //printGraphEdges(g);
    //printf("Edges to be checked:\n");
    //printEdgeVector(edgesToBeChecked);
    // make g_plus_constraints
    g_plus_constraints.clear();
    boost::copy_graph(g, g_plus_constraints);
    for (const auto& pair : constraints) boost::add_edge(index[pair.first], index[pair.second], g_plus_constraints); // vertex index
    //printf("Current graph+constraints edges:\n");
    //printGraphEdges(g_plus_constraints);

    // check g_plus_constraints connectivity
    if (!checkClusterCriteria(g_plus_constraints)) {
        double timeInPs = simulation->context->getState(OpenMM::State::Energy).getTime();
        printf("Time = %.3f ps\n", timeInPs);
        throw std::bad_exception();
    }
    else {
        for (const auto& pair : constraints) {
            std::vector<std::pair<int, int>>::iterator position = std::find(edgesToBeChecked.begin(), edgesToBeChecked.end(), pair);
            if (position != edgesToBeChecked.end()) edgesToBeChecked.erase(position);
        }
        std::shuffle(edgesToBeChecked.begin(), edgesToBeChecked.end(), gen);
        // check whether edges are important in graph
        for (int i=0; i<edgesToBeChecked.size(); ++i) { // a-b is longer than rc-eps
            int a = edgesToBeChecked[i].first;
            int b = edgesToBeChecked[i].second;
            boost::remove_edge(index[a], index[b], g_plus_constraints);
            if (!checkClusterCriteria(g_plus_constraints)) { // edge (a,b) is important and a-b is longer than rc-eps
                constraints.insert(std::make_pair(a, b));
                if (a%2 == 0)
                    simulation->custombond->setBondParameters(a/2*upperLimit+(b-1)/2,a,b,{1.0});
                else
                    simulation->custombond->setBondParameters(b/2*upperLimit+(a-1)/2,b,a,{1.0});
                double d = getDistance(positions, a, b);
                boost::add_edge(index[a], index[b], d, g_plus_constraints);
            }
        }
        edgesToBeChecked.clear();
        //printf("Current edge set with constraints:\n");
        //printEdgeSet(constraints);

        simulation->custombond->updateParametersInContext(*(simulation->context));
        simulation->integrator->step(2);

        //printf("One step finished.\n\n");
        //printf("\n");
    }
}

double GcmcIonPair::getInsrProbability(){
    return sucInsr/(double)numInsr;
}

double GcmcIonPair::getDelProbability(){
    return sucDel/(double)numDel;
}

double GcmcIonPair::getMDProbability(){
    return sucMD/(double)numMD;
}

void GcmcIonPair::rezeroStat(){
    numInsr = 0;
    numDel = 0;
    numMD = 0;
    sucInsr = 0;
    sucDel = 0;
    sucMD = 0;
    rezeroHistogram();
    rezeroEnergy();
}

void GcmcIonPair::rezeroHistogram(){
    for (int i = 0; i < upperLimit*M + 1; i++)
        histogram[i] = 0;
}

void GcmcIonPair::rezeroEnergy(){
    for (int i = 0; i < upperLimit*M + 1; i++)
        energy[i] = 0.0;
}

void GcmcIonPair::updateWangLaudauFactor(){
    double sumHistogram = 0;
    for (int i = lowerLimit*M; i < upperLimit*M+1; i++)
        sumHistogram += histogram[i];
    double aveHistogram = sumHistogram/(double)((upperLimit-lowerLimit)*M+1);
    int count = 0;
    for (int i = lowerLimit*M; i < upperLimit*M+1; i++)
        if (histogram[i] > 0.8*aveHistogram)
            count++;
    if (count == (upperLimit-lowerLimit)*M+1){
        f = 0.5*f;
        rezeroHistogram();
    }
}

double* GcmcIonPair::getEnergy(){
    for (int i = lowerLimit*M; i < upperLimit*M + 1; i++)
        if (histogram[i] > 0)
            energy[i] = energy[i]/histogram[i];
    return energy;
}

double* GcmcIonPair::getFreeEnergy(){
    for (int i = lowerLimit*M; i < upperLimit*M + 1; i++)
        freeEnergy[i] = -(log(histogram[i]/(1.0*histogram[lowerLimit*M]))-eta[i]+
            eta[lowerLimit*M])/beta;
    return freeEnergy;
}

