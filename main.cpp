//needed headers
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <numeric>
#include <string>
#include <getopt.h>
#include <sstream>


//Other connected classes
#include "Replica.h"
#include "searchparameters.h"
#include "searchresult.h"
#include "random.h"
#include "Conformation.h"
#include "const.h"

using namespace std;

void check_options( searchparameters &params, long int argc, char *argv[] );
void print_program_id();
void compareFolds( string, int );
void validateFold( string, string, int );

string fold = "";
int numSets = 1;
bool validate = false;

int main (int argc, char *argv[]) {
	searchparameters params;
	check_options( params, argc, argv );
	
	if( validate == true ) {
        validateFold( params.seq, fold, params.energy );
		return 0;
	}
    
	if( params.seedvalue != 0 ) {
		setSeed( params.seedvalue );
	} else {
		setSeed( time(0) );
	}
    
	ReplicaExchange simulation( params );
	vector<solution> bestSolutionTrace;
	
	cout << "Begin Simulation" << endl;
	simulation.runreplicaexchange(bestSolutionTrace);
	cout << "End Simulation" << endl << endl;
	
	return 0;
}

void validateFold( string sequence, string fold, int energy ) {
	Conformation *c = new Conformation(sequence,200,"",1);
    cout<<"Enter validation\n";
    
    c->loadConformation(fold);
	cout << c->getDirections() << endl;
	c->recalculateDirections();
	cout << c->getDirections() << endl;
	
	if( ! c->testLatticePositions() ) {
		cout << "Invalid fold" << endl << endl;
		c->printRelLatticePositions();
		return;
	}
	
	if( ! c->testFreeEnergy(energy) ) {
		cout << "Invalid energy" << endl << endl;
		c->printRelLatticePositions();
		return;
	}
	cout << "Valid fold and energy" << endl << endl;
}


void check_options(searchparameters &params, long int argc, char *argv[] )
{
    int opt;
    const string delims(", ");
    string line;
    string::size_type begIdx, endIdx;
    
    int option_index = 0;
    
    cin>>params.runs>>params.seq>>params.energy>>params.no_of_replicas>>params.mintemp>>params.maxtemp>>params.steps>>params.runtime>>params.seedvalue;
        
    
} 


