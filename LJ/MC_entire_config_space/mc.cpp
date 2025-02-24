#include <cstdio>
#include <random>
#include <cmath>

#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <time.h>

//Parameters that define the system
#define MAXPART 10

#define CLUSTER_DIST 1.5
#define CLUSTER_DIST2 (CLUSTER_DIST*CLUSTER_DIST)
#define T 0.7 //reduced temperature in the regime of liquid/vapor coexistence

//Parameters that define the length of the simulation
#define NCYCLE_EQUIL 1000000
#define NCYCLE_PROD  4000000

int npart = 10; // TODO: This is current hard-coded
double current_energy;

//initializes particles into a 1D line
void init_positions(double (&positions)[MAXPART][3], bool (&connectivity)[MAXPART][MAXPART]) {
	//initialize to all disconnected
	for (int ipart = 0; ipart < npart; ipart++)
		for (int jpart = 0; jpart < npart; jpart++)
			connectivity[ipart][jpart] = false;
		
	for (int ipart = 0; ipart < npart; ipart++) {
		positions[ipart][0] = 0.0;
		positions[ipart][1] = 0.0;
		positions[ipart][2] = ipart * 1.0;
		
		//also update the connectivity, initially forming a line
		if (ipart > 0) {
			connectivity[ipart-1][ipart] = true;
			connectivity[ipart][ipart-1] = true;
		}
	}
}

//write the current position out to an XYZ formatted file
void write_positions(FILE* posOutput, double (&positions)[MAXPART][3]) {
        for (int ipart = 0; ipart < npart; ipart++) {
                fprintf(posOutput, "C %f %f %f\n", positions[ipart][0], positions[ipart][1], positions[ipart][2]);
        }
        fprintf(posOutput, "\n");
}

void read_positions(const char* filename, double (&positions)[MAXPART][3], bool (&connectivity)[MAXPART][MAXPART]) {
	FILE *handle = fopen(filename, "r");
	char type[20];
	for (int ipart = 0; ipart < npart; ipart++) {
		fscanf(handle, "%s %lf %lf %lf\n", type, &positions[ipart][0], &positions[ipart][1], &positions[ipart][2]);
	}
	
	double dx,dy,dz,dr2;
	for (int ipart = 0; ipart < npart; ipart++) {
		for (int jpart = 0; jpart < ipart; jpart++) {
			dx = positions[ipart][0]-positions[jpart][0];
			dy = positions[ipart][1]-positions[jpart][1];
			dz = positions[ipart][2]-positions[jpart][2];
			
			dr2 = dx*dx + dy*dy + dz*dz;
			
			if (dr2 >= CLUSTER_DIST2) {
				connectivity[ipart][jpart] = false;
				connectivity[jpart][ipart] = false;
			} else {
				connectivity[ipart][jpart] = true;
				connectivity[jpart][ipart] = true;
			}
		}
	}
	fclose(handle);
}

//returns the energy in dimensionless units
double get_energy(double (&positions)[MAXPART][3], bool (&connectivity)[MAXPART][MAXPART]) {
	double dx,dy,dz,dr2,drm2,drm6,drm12,dr;
	double e = 0;
	for (int ipart = 0; ipart < npart; ipart++)
		for (int jpart = ipart + 1; jpart < npart; jpart++) {
			dx = positions[ipart][0]-positions[jpart][0];
			dy = positions[ipart][1]-positions[jpart][1];
			dz = positions[ipart][2]-positions[jpart][2];
			
			dr2 = (dx*dx + dy*dy + dz*dz);
			drm2 = 1.0/dr2;
			drm6 = drm2*drm2*drm2;
			drm12 = drm6*drm6;
			
			e += 4*(drm12-drm6);
		}
	
	return e;
}

double get_ave_dist(double (&positions)[MAXPART][3], bool (&connectivity)[MAXPART][MAXPART]) {
        double dx,dy,dz,dr2,drm2,dr;
        double d = 0;
        for (int ipart = 0; ipart < npart; ipart++)
                for (int jpart = ipart + 1; jpart < npart; jpart++) {
                        dx = positions[ipart][0]-positions[jpart][0];
                        dy = positions[ipart][1]-positions[jpart][1];
                        dz = positions[ipart][2]-positions[jpart][2];

                        dr2 = (dx*dx + dy*dy + dz*dz);
                        dr = sqrt(dr2);
                        d +=dr;
                }
        d /= 6;
        return d;
}

//returns true if the clusters satisfy the cluster criterion
bool check_connected(bool (&connectivity)[MAXPART][MAXPART]) {
	
	using namespace std;
	using namespace lemon;
	
	ListGraph graph; //graph of connectivity
	ListGraph::Node nodes[MAXPART]; //list of all nodes
	
	for (int ipart = 0; ipart < npart; ipart++) {
		nodes[ipart] = graph.addNode();
	}
	
	for (int ipart = 0; ipart < npart; ipart++)
		for (int jpart = 0; jpart < ipart; jpart++)
			if (connectivity[ipart][jpart])
				graph.addEdge(nodes[ipart], nodes[jpart]);
	
	return connected(graph);
}

//(potentially) executes a MC single particle displacement move
int ndisplacement_accepted = 0;
int ndisplacement_rejected = 0;
double do_displement(double (&positions)[MAXPART][3], bool (&connectivity)[MAXPART][MAXPART]) {
	
	bool prev_connectivity[MAXPART][MAXPART];
	std::copy(std::begin(connectivity), std::end(connectivity), std::begin(prev_connectivity));
		
	bool need_to_check_connectivity = false;
	const double max_displacmenet = 2*0.35; //(twice the) largest possible trial move
	
	//first pick the candidate particle
	int ipart = rand() % npart;
	//then do a trial displacement
	double trial_dx = (drand48()-0.5)*max_displacmenet;
	double trial_dy = (drand48()-0.5)*max_displacmenet;
	double trial_dz = (drand48()-0.5)*max_displacmenet;
	
	//calculate the CHANGE in energy
	double dx,dy,dz,dr2,drm2,drm6,drm12;
	double de = 0;
	for (int jpart = 0; jpart < npart; jpart++) {
		if (ipart == jpart)
			continue;
		
		dx = positions[ipart][0]-positions[jpart][0] + trial_dx;
		dy = positions[ipart][1]-positions[jpart][1] + trial_dy;
		dz = positions[ipart][2]-positions[jpart][2] + trial_dz;
		
		dr2 = dx*dx + dy*dy + dz*dz;
		//first check the cluster criterion
		if (dr2 >= CLUSTER_DIST2) {
			if (prev_connectivity[ipart][jpart]) {
				need_to_check_connectivity = true; //we just broke an existing connection!
				connectivity[ipart][jpart] = false;
				connectivity[jpart][ipart] = false;
			}
		} else {
			connectivity[ipart][jpart] = true;
			connectivity[jpart][ipart] = true;
		}
		
		drm2 = 1.0/dr2;
		drm6 = drm2*drm2*drm2;
		drm12 = drm6*drm6;
		//add on the new contribution to the energy
		de += 4*(drm12-drm6);
		
		dx = positions[ipart][0]-positions[jpart][0];
		dy = positions[ipart][1]-positions[jpart][1];
		dz = positions[ipart][2]-positions[jpart][2];
		
		dr2 = dx*dx + dy*dy + dz*dz;
		drm2 = 1.0/dr2;
		drm6 = drm2*drm2*drm2;
		drm12 = drm6*drm6;		
		//subtract off the old contribution to the energy
		de -= 4*(drm12-drm6);
		
	}
	
	//now accept or reject the move
	
	//immediately reject if this move will not obey the cluster criterion
	if (need_to_check_connectivity && !check_connected(connectivity)) {
		std::copy(std::begin(prev_connectivity), std::end(prev_connectivity), std::begin(connectivity));
		ndisplacement_rejected++;

		return current_energy;
	}

	//otherwise check the energy
	double uniform = drand48();
	if (de < 0 || uniform < exp(-de/T) ) {
		//accept move
		positions[ipart][0] += trial_dx;
		positions[ipart][1] += trial_dy;
		positions[ipart][2] += trial_dz;
		//update the energy
		current_energy += de;
		
		ndisplacement_accepted++;
		return current_energy;
	}
	//otherwise restore configuration
	std::copy(std::begin(prev_connectivity), std::end(prev_connectivity), std::begin(connectivity));
	ndisplacement_rejected++;
	
	return current_energy;
}

int main() {
        srand(time(0));
	double positions[MAXPART][3]; //all atomic positions
	bool connectivity[MAXPART][MAXPART]; //connectivity matrix; upper triangular only needed

	init_positions(positions, connectivity);
	//read_positions("input.xyz", positions, connectivity);
	current_energy = get_energy(positions, connectivity);
	printf("Initial energy: %f\n", current_energy);
	
	for (int icycle = 0; icycle < NCYCLE_EQUIL; icycle++)
		for (int imove = 0; imove < npart; imove++)
			do_displement(positions, connectivity);

        double avg_energy = 0;
        double pe = 0;
        //FILE *energy = fopen("energy_allgraph.txt", "w");
        //FILE *posF = fopen("cluster_allgraph.xyz", "w");

        for (int icycle = 0; icycle < NCYCLE_PROD; icycle++){
                for (int imove = 0; imove < npart; imove++){
                        pe = do_displement(positions, connectivity);
                        avg_energy += pe;
                        //fprintf(energy, "%.4f\n", pe);
                        //write_positions(posF, positions);
                }
        }

        avg_energy /= (NCYCLE_PROD*npart);
        //fclose(energy);
        //fclose(posF);

	printf("Average energy: %f\n", avg_energy);
	printf("Displacement percent accepted: %f\n", ndisplacement_accepted / (double) (ndisplacement_accepted + ndisplacement_rejected)); 
	
	return 0; //successful completion
}
