#ifndef CONFORMATION_H
#define CONFORMATION_H

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <sstream>
#include <map>
#include <assert.h>
#include "const.h"
#include "random.h"
#include "searchresult.h"
#include "searchparameters.h"
#include "Timer.h"
#include "stopwatch.h"

// default static initialization
#if DIM == 2
static const int VSF[2] = {1,0};
static const int VLF[2] = {0,1};
static const int VRF[2] = {0,-1};
static const int VBF[2] = {-1,0};

static const int VSB[2] = {-1,0};
static const int VLB[2] = {0,-1};
static const int VRB[2] = {0,1};
static const int VBB[2] = {1,0};

static const int BS[2][2] = {1,0, 0,1};
static const int BL[2][2] = {0,-1, 1,0};
static const int BR[2][2] = {0,1, -1,0};
static const int BB[2][2] = {-1,0,0,-1};

static const int DI[4] = {1,L,-1,-L};

#elif DIM == 3
static const int VSF[3] = {1,0,0};
static const int VLF[3] = {0,1,0};
static const int VRF[3] = {0,-1,0};
static const int VBF[3] = {-1,0,0};
static const int VUF[3] = {0,0,1};
static const int VDF[3] = {0,0,-1};

static const int VSB[3] = {-1,0,0};
static const int VLB[3] = {0,-1,0};
static const int VRB[3] = {0,1,0};
static const int VBB[3] = {1,0,0};
static const int VUB[3] = {0,0,-1};
static const int VDB[3] = {0,0,1};

static const int BS[3][3] = {1,0,0,   0,1,0,  0,0,1};
static const int BL[3][3] = {0,-1,0, 1,0,0,   0,0,1};
static const int BR[3][3] = {0,1,0,  -1,0,0,  0,0,1};
static const int BU[3][3] = {0,0,-1, 0,1,0, 1,0,0};
static const int BD[3][3] = {0,0,1, 0,1,0, -1,0,0};
static const int BB[3][3] = {-1,0,0, 0,-1,0, 0,0,-1};

static const int DI[6] = {1,L,LL,-1,-L,-LL};

#endif

using namespace std;

class Conformation {
private:
    
    static int OFFSET;
    static double probAcceptValues[MAX_NUM_REPLICAS][MAX_NUM_ENERGY_VALUES];
    static double Tvalues[MAX_NUM_REPLICAS];
    
    // local class properties
    int replicaID;
    int seqLength;
    int FreeEnergy;
    double pmWeight;
    double T;
    int Lattice[HashSize]; //Lattice hash
    int i0;	// initial starting position in lattice
    int ii[MAX_SEQ_LENGTH]; //hash keys into the lattice
    int FixedPoint; //folding point
    int Directions[MAX_SEQ_LENGTH]; //relative directions of aa
    
    int sequence[MAX_SEQ_LENGTH];
    int sequenceCode[MAX_SEQ_LENGTH];
    int bestRelDirections[MAX_SEQ_LENGTH];
    int bestFreeEnergy;
    
    int FreeEnergy_M;
    int Lattice_M[HashSize];
    int ii_M[MAX_SEQ_LENGTH];
    int move_startResidue;
    int move_endResidue;
    bool move_isPossible;
    //Mutation mutation;
    searchresult searchResult;
    float prevTime;
    long double prevMCStep;
    Timer t_clock;
    
    // used to track statistics
    long double acceptedMoves;
    long double rejectedMoves;
    long double impossibleMoves;
    long double numResidues;
    
    // private class methods
    int row(int);
    int col(int);
#if DIM == 3
    int height(int);
#endif
    void determineDirection(int i, int n1, int n2, int &Dir, int GBasisP[DIM][DIM], int fixedpoint);
#if DIM == 3
    int determineII(int TBasis [DIM][DIM], int Direction, bool forward,int newBasis[DIM][DIM]);
#endif
    void determineBasis( int n1, int n2, int Basis[DIM][DIM] );
    int findFreeNeighbours(int residue, int freeN[]);
    int findFreeNeighboursL(int lp, int freeN[DIM*2] );
    int findFreeNeighbours(int res1, int res2, int freeN[DIM*2]);
    int selectRandomFreeNeighbour(int freeN[DIM*2], int numFree);
    int countBonds(int n1, int lp, int plp, int nlp);
    void applyMutation();
    void resetMutation();
    void recenterOnLattice();
    
#if MOVESET == LOCAL || MOVESET == MIXED
    void tryEndMove(int residue);
    //#if DIM == 2
    //		void tryCrankshaftOrCornerMove(int residue);
    //#elif DIM == 3
    void tryCrankshaftMove(int residue);
    void tryCornerMove(int residue);
    //#endif
#endif
    
    
#if MOVESET == PULLMOVES || MOVESET == MIXED
    void tryEndPullMove(int residue);
    void tryPullMove(int residue);
#endif
    
    // private methods
    void printRelDirections();
    void printLatticePositions();
    
    int calcFreeEnergy();
    
	// public interface
public:
    // default constructor
    Conformation(const string &sequence, const double temperature, const string &experimentID, const int replicaNumber);
    
    static void initProbAcceptValues(const double tempValues[], int numReplicas, int absKnownEnergy) {
        // for each temperature value
        for(int i=0; i<numReplicas; i++) {
            double T = tempValues[i];
            Tvalues[i] = T;
            
            // for each possible energy delta
            for(int j=0; j<=absKnownEnergy; j++) {
                probAcceptValues[i][j] = exp(-j/(T*K_b));
            }
        }
        
    }
    
    static void swapOffset() {
        OFFSET = 1 - OFFSET;
    }
    
    bool testLatticePositions();
    bool testFreeEnergy(int);
    void printRelLatticePositions();
    
    // method to set the sequence
    void setSequence( const string &seq ) {
        seqLength = seq.length();
        for(int i=0; i<seqLength; i++) {
            if( seq[i] == 'W')
                sequence[i] = W;
            else
                sequence[i] = B;
        }
    }
    
    // method to get the sequence
    string getSequence() {
        ostringstream buf;
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] == W )
                buf << 'w';
            else
                buf << 'b';
        }
        return buf.str();
    }
    
    string getDirections() {
        char relDir[seqLength];
        for(int i=0; i<seqLength; i++) {
            switch(Directions[i]) {
				case straight: relDir[i] = 's'; 	break;
				case left:	 relDir[i] = 'l'; 		break;
				case right:	 relDir[i] = 'r'; 		break;
				case up:	 relDir[i] = 'u'; 		break;
				case down:	 relDir[i] = 'd'; 		break;
				case ground: relDir[i] = 'g'; 	break;
            }
        }
        string tmp = relDir;
        return tmp;
    }
    
    void recalculateDirections();
    
    // returns the current temperature of the system
    double getTemperature();
    
    // set the current temperature of the system
    void setTemperature(const double t);
    
    // returns the current free energy of the system
    int getFreeEnergy();
    
    // returns the replica id
    int getReplicaID();
    
    void setReplicaID(int a) {
        replicaID = a;
    }
    
    // runs the simulation for n steps
    searchresult& runSimulation(const long double n);
    
    // returns the search results
    searchresult& getSearchResult();
    
    // informs SearchResult to dump results to file
    void writeResults();
    
    // temporary test function
    void test2();
    
    // blank the lattice and all conformation information
    void reset();
    
    // reset the conformation to a straight conformation
    void makeStraight();
    
    // calculate the free energy of current conformation
    void computeFreeEnergy();
    
    // prints the conformation
    void printConformation();
    
    // returns the contact string
    string getContactString();
    
    // returns the symmetric contact string in case folding occured from the other end
    // in a symmetric sequence
    string getSymmetricContactString();
    
    // returns the relative h-h contact order
    double relativeWWContactOrder();
    
    // returns the absolute h-h contact order
    double absoluteWWContactOrder();
    
    // returns the solvent accessible area
    int solventAccessibility();
    
    // adds frequency of residue contacts to contact frequency matrix (parameter)
    void addToContactFrequencyMatrix(vector< vector< int > > &, int);
    
#if DIM == 2
    // quickly folds a near optimal solution
    void fastFold(double);
#endif
    // loads a conformation given by the parameters
    void loadConformation(const int [], const int=0);
    
    void loadConformation(const string &, const int=0);
    
    
};

// returns the relative row of a lattice position
inline
int Conformation:: row(int latticePos) {
#if DIM == 2
	return latticePos / L;
#elif DIM == 3
	return (latticePos & (LL-1)) / L;
#endif
}

// returns the relative column of a lattice position
inline
int Conformation:: col(int latticePos) {
	return latticePos & (L - 1);
}

#if DIM == 3
// returns the relative height of a lattice position
inline
int Conformation:: height(int latticePos) {
	return latticePos / LL;
}
#endif

#if DIM == 2
// given two residue n1 and n2, determines the basis of direction
// from n1 to n2 in 2D
inline
void Conformation:: determineBasis( int n1, int n2, int Basis[DIM][DIM] ) {
	int vect[DIM];
	const int *basis;
	
	if( n2>n1 ) {
		vect[0] = row(ii[n2]) - row(ii[n1]);
		vect[1] = col(ii[n2]) - col(ii[n1]);
		
		if( vect[0] == VSF[0]  && vect[1] == VSF[1] ) {
			basis = &BS[0][0];
		} else if( vect[0] == VLF[0]  && vect[1] == VLF[1] ) {
			basis = &BL[0][0];
		} else if( vect[0] == VRF[0]  && vect[1] == VRF[1] ) {
			basis = &BR[0][0];
		} else if( vect[0] == VBF[0]  && vect[1] == VBF[1] ) {
			basis = &BB[0][0];
		}
	} else {
		vect[0] = row(ii[n1]) - row(ii[n2]);
		vect[1] = col(ii[n1]) - col(ii[n2]);
		
		if( vect[0] == VSB[0]  && vect[1] == VSB[1] ) {
			basis = &BS[0][0];
		} else if( vect[0] == VLB[0]  && vect[1] == VLB[1] ) {
			basis = &BL[0][0];
		} else if( vect[0] == VRB[0]  && vect[1] == VRB[1] ) {
			basis = &BR[0][0];
		} else if( vect[0] == VBB[0]  && vect[1] == VBB[1] ) {
			basis = &BB[0][0];
		}
		
	}
	
	for(int row=0; row<DIM; row++) {
		for(int col=0; col<DIM; col++) {
			Basis[row][col] = *(basis + row*DIM + col);
			//cout << Basis[row][col] << " ";
		}
		//cout << endl;
	}
}

#endif

// given a free neighbour matrix, randomly selects one of the free positions
inline
int Conformation:: selectRandomFreeNeighbour(int freeN[DIM*2], int numFree) {
	double cum_p[DIM*2];
	int Q = random_number(1, numFree);
	
	int index = 0;
	int passes = 0;
	
	for(int i=0; i<DIM*2; i++) {
		if( freeN[i] != -1 ) {
			passes++;
		}
		if( passes == Q ) {
			index = i;
			break;
		}
	}
	
	return freeN[index];
}

// returns the free neighbours of a given residue
inline
int Conformation:: findFreeNeighbours(int residue, int freeN[DIM*2]) {
	return findFreeNeighboursL(ii[residue], freeN);
}

inline
int Conformation:: findFreeNeighboursL(int lp, int freeN[DIM*2] ) {
	int numFree = 0;
	
	for(int i=0; i<DIM*2; i++) {
		if( Lattice[lp + DI[i]] == LPempty ) {
			freeN[i] = lp + DI[i];
			numFree++;
		} else {
			freeN[i] = -1;
		}
	}
	return numFree;
}

// returns the free neighbours in common of two given residues
inline
int Conformation:: findFreeNeighbours(int res1, int res2, int freeN[DIM*2]) {
	int freeN1[DIM*2];
	int freeN2[DIM*2];
	int numFree=0;
	
	findFreeNeighbours(res1, freeN1);
	findFreeNeighbours(res2, freeN2);
	
	for(int i=0; i<DIM*2; i++) {
		freeN[i]=-1;
		if( freeN1[i] != -1 ) {
			for(int j=0; j<DIM*2; j++) {
				if( freeN1[i] == freeN2[j] ) {
					freeN[i] = freeN1[i];
					numFree++;
				}
			}
		}
	}
	return numFree;
}

// returns the number of hydrogen bonds between residue at
// lattice position lp and other residues.  plp and nlp are
// the lattice positions of the previous and next residue.
// these are required to ensure we don't count them as hydrogen
// bonds.  if plp or nlp = -1, we assume they don't exist
inline
int Conformation::countBonds(int n1, int lp, int plp, int nlp) {
	int bonds = 0;
	
	if( sequence[n1] < 0 ) return 0;
	
	for(int i=0; i < 2*DIM; i++) {
		int olp = lp + DI[i];
		if( Lattice[ olp ] >= 0 && olp != plp && olp != nlp && olp != ii[n1] ) bonds++;
	}
	return bonds;
}

//#if DIM == 2
// returns the new direction from n1 to n2 in 2D
inline
void Conformation:: determineDirection(int i, int n1, int n2,
                                       int &Dir, int GBasisP[DIM][DIM], int fixedpoint) {
	
 	int rvect[DIM];
	int co_i[DIM], co_j[DIM];
	
	if(i==fixedpoint-1 || i == fixedpoint-2){
		Dir=ground;
		GBasisP[0][0] = 1;
		GBasisP[0][1] = 0;
		
		GBasisP[1][0] = 0;
		GBasisP[1][1] = 1;
		
	}else {
		if(i==n2){
			
			co_i[0]=ii[i]/L;
			co_i[1]=(ii[i]&(L-1));
			co_j[0]=ii[i-1]/L;
			co_j[1]=(ii[i-1]&(L-1));
			
			for(int dim = 0; dim<DIM; dim++)
				rvect[dim]= co_i[dim]-co_j[dim];
			
			int vect[DIM];
            int multBasis [DIM][DIM];
			
            for (int row = 0; row<DIM; row++){
				int transp_row[DIM];
				
				for(int col = 0; col<DIM; col++){
					transp_row[col] = GBasisP[col][row];
				}
				
				vect[row] = 0;
				for(int col = 0; col<DIM; col++){
					vect[row]+=transp_row[col]*rvect[col];
				}
			}
            
			if(((vect[0] == VSF[0]) && (vect[1] == VSF[1]))){// && (vect[2] == VSF[2])) ){
				Dir = straight;
				
				for(int row = 0; row<DIM; row++){
					for(int col = 0; col<DIM; col++){
						multBasis[row][col] = BS[row][col];
					}
                }
			}
			else if((vect[0] == VLF[0] && vect[1] == VLF[1] )){//&& vect[2] == VLF[2])  ){
                Dir = left;
				
				for(int row = 0; row<DIM; row++)
					for(int col = 0; col<DIM; col++)
						multBasis[row][col] = BL[row][col];
			}
			else if((vect[0] == VRF[0] && vect[1] == VRF[1])){// && vect[2] == VRF[2])){
				Dir = right;
				for(int row = 0; row<DIM; row++)
					for(int col = 0; col<DIM; col++)
						multBasis[row][col] = BR[row][col];
			}
			
			int prevBasis[DIM][DIM];
            for(int row = 0; row<DIM; row++)
				for(int col = 0; col<DIM; col++)
					prevBasis[row][col] = GBasisP[row][col];
			for(int row=0; row<DIM; row++){
				for(int col=0; col<DIM; col++){
					GBasisP[row][col] = 0;
					for(int mf=0; mf<DIM; mf++)
						GBasisP[row][col]+=prevBasis[row][mf]*multBasis[mf][col];
				}
			}
		}
		
	} 
}

int Conformation::OFFSET=0;
double Conformation::probAcceptValues[MAX_NUM_REPLICAS][MAX_NUM_ENERGY_VALUES];
double Conformation::Tvalues[MAX_NUM_REPLICAS];

Conformation::Conformation(const string &seq, const double t, const string &expID, const int repID)
:searchResult(expID, repID)
{
	setSequence( seq );
	
	replicaID = repID;
	i0 = I0;
    
	for(int i=0; i<HashSize; i++) {
		Lattice[i] = LPempty;				Lattice_M[i] = LPempty;
	}
    int residueType;
	for(int i=0; i<seqLength; i++) {
		if( sequence[i] == B ) {
			Lattice[i0 + (i - (int)(seqLength/2))*L] = B; 	Lattice_M[i0 + (i - (int)(seqLength/2))*L] = B;
			sequenceCode[i] = B;
		} else {
			Lattice[i0 + (i - (int)(seqLength/2))*L] = i; 	Lattice_M[i0 + (i - (int)(seqLength/2))*L] = i;
			sequenceCode[i] = i;
		}
		ii[i] = i0 + (i - (int)(seqLength/2))*L; 	ii_M[i] = i0 + (i - (int)(seqLength/2))*L;
		Directions[i] = straight;
	}
	Directions[0] = ground;
	Directions[1] = ground;
	
	FixedPoint = 2;
	FreeEnergy = 0; FreeEnergy_M = 0;
	T = t;
	prevTime = 0;
	prevMCStep = 0;
	pmWeight = searchparameters::pmweight;
}
void Conformation::writeResults() {
	searchResult.write();
}
double Conformation::getTemperature() {
	return Tvalues[replicaID];
}
void Conformation::setTemperature(double t) {
	T = t;
}
int Conformation::getFreeEnergy() {
	return FreeEnergy;
}
int Conformation::getReplicaID() {
	return replicaID;
}

searchresult& Conformation::getSearchResult() {
	return searchResult;
}

searchresult& Conformation::runSimulation(long double n) {
	int mutationPoint;
    
	acceptedMoves = 0;
	rejectedMoves = 0;
	impossibleMoves = 0;
	numResidues = 0;
	long double i;
    
	for(i = 0; i < n; i++) {
		mutationPoint = random_number(0,seqLength-1);
        
#if MOVESET == LOCAL
		if( mutationPoint == 0 || mutationPoint == seqLength-1 ) {
			tryEndMove( mutationPoint );
        } else
        {
			double Q = ran01();
			if( Q >= 0.8 ) {
				tryCrankshaftMove( mutationPoint );
			} else {
				tryCornerMove( mutationPoint );
            }
            
		}
        
#elif MOVESET == PULLMOVES
		if( mutationPoint == 0 || mutationPoint == seqLength-1 ) {
			tryEndPullMove( mutationPoint );
		} else {
			tryPullMove( mutationPoint );
		}
        
#elif MOVESET == MIXED
        
		if( mutationPoint == 0 || mutationPoint == seqLength-1 ) {
			double Q = ran01();
			if( Q <= pmWeight ) {
				tryEndPullMove( mutationPoint );
			} else {
				tryEndMove( mutationPoint );
			}
		} else {
			double Q = ran01();
			if( Q <= pmWeight ) {
				tryPullMove( mutationPoint );
			} else {
				double Q = ran01();
				if( Q >= 0.8 ) {
					tryCrankshaftMove( mutationPoint );
				} else {
					tryCornerMove( mutationPoint );
				}
			}
		}
#endif
		
		if( move_isPossible ) {
			
			int energyDelta = FreeEnergy_M - FreeEnergy;
			
			if( energyDelta < 0 ) {
				applyMutation();
				acceptedMoves++;
			} else {
				
				double prob = probAcceptValues[replicaID][energyDelta];
				double Q = ran01();
                
				if( prob >= Q)
                {
					applyMutation();
					acceptedMoves++;
				} else {
					rejectedMoves++;
					resetMutation();
				}
			}
			
		} else {
			impossibleMoves++;
		}
		
		if( FreeEnergy < bestFreeEnergy ) {
			recalculateDirections();
			searchResult.addsol(prevMCStep + i, 0 /*prevTime + t_clock.elapsed()*/,
                                     FreeEnergy, Tvalues[replicaID], Directions, seqLength);
			bestFreeEnergy = FreeEnergy;
			for(int i=0; i<seqLength; i++)
				bestRelDirections[i] = Directions[i];
		}
        
	}
	prevMCStep += n - 1;
	
	searchResult.setma(acceptedMoves);
	searchResult.setmi(impossibleMoves);
	searchResult.setmr(rejectedMoves);
	searchResult.setresidue(numResidues);
	searchResult.setsteps(i-1);
    
	return searchResult;
}

void Conformation::recalculateDirections() {
	int GBasis[DIM][DIM];
	
	Directions[0] = ground;
    Directions[1] = ground;
	determineBasis(0, 1, GBasis);
	
	for(int i=2; i<seqLength; i++) {
		determineDirection(i, i-1, i, Directions[i], GBasis, 2);
	}
}

#if MOVESET == LOCAL || MOVESET == MIXED
void Conformation::tryEndMove(int residue) {
    int freeNeighbours[2*DIM];
	int bondedResidue;
	int plp, nlp;
	move_isPossible = false;
	
	if( residue == 0 ) {
		bondedResidue = 1;
		plp = -1;
		nlp = ii[residue + 1];
	}
	else if(residue == (seqLength-1)) {
		bondedResidue = seqLength-2;
		plp = ii[residue - 1];
		nlp = -1;
	}
	else
		return;
	
	int num = findFreeNeighbours(bondedResidue, freeNeighbours);
    
    if(num == 0)
		return;
    int lp = selectRandomFreeNeighbour(freeNeighbours, num);
	Lattice_M[ii_M[residue]] = LPempty;
	Lattice_M[lp] = sequenceCode[residue];
	ii_M[residue] = lp;
	int bondsBefore = countBonds(residue, ii[residue], plp, nlp);
	int bondsAfter = countBonds(residue, lp, plp, nlp);
    
	FreeEnergy_M += bondsBefore - bondsAfter;
	move_startResidue = residue;
	move_endResidue = residue;
	move_isPossible = true;
	
	return; ;
}

void Conformation::tryCrankshaftMove(int residue) {
#if DIM == 2
	static const int componentSum = (1 + L);
#elif DIM == 3
	static const int componentSum = (1 + L + LL);
#endif
    
	int newNiplus1 = -1, newNiplus2 = -1;
	bool canCrankshaft = false;
	
	move_isPossible = false;
    
	//cout << "Begin Crankshaft" << endl;
	
	// can't do a corner or crankshaft move on first or last residue
	if( residue == 0 || residue == seqLength-1 ) {
		return;
	}
	
	int abs_i_to_i_plus_3;
	
	if( residue + 3 < seqLength - 1 ) {
		abs_i_to_i_plus_3 = abs(ii[residue+3] - ii[residue]);
		
#if DIM == 2
		if( abs_i_to_i_plus_3 == 1 || abs_i_to_i_plus_3 == L )
#elif DIM == 3
            if( abs_i_to_i_plus_3 == 1 || abs_i_to_i_plus_3 == L || abs_i_to_i_plus_3 == LL  )
#endif
                canCrankshaft = true;
	}
	
	if( !canCrankshaft && residue - 3 >= 1 ) {
		residue = residue - 3;
		abs_i_to_i_plus_3 = abs(ii[residue+3] - ii[residue]);
        
#if DIM == 2
		if( abs_i_to_i_plus_3 == 1 || abs_i_to_i_plus_3 == L )
            
#endif
			canCrankshaft = true;
	}
	
	if( canCrankshaft ) {
		int abs_i_to_i_plus_1 = abs(ii[residue+1] - ii[residue]);
		int distance = componentSum - abs_i_to_i_plus_3 - abs_i_to_i_plus_1;
        
		// randomly select which crankshaft direction to try first
		int direction = random_number(1,2);
		if( direction == 2 ) distance = -distance;
        
		if( (Lattice[ii[residue] + distance] == LPempty) && (Lattice[ii[residue + 3] + distance] == LPempty) ) {
			newNiplus1 = ii[residue] + distance;
			newNiplus2 = ii[residue + 3] + distance;
		} else if( (Lattice[ii[residue] - distance] == LPempty) && (Lattice[ii[residue + 3] - distance] == LPempty) ) {
			newNiplus1 = ii[residue] - distance;
			newNiplus2 = ii[residue + 3] - distance;
		}
        
		if( newNiplus1 != -1 && newNiplus2 != -1 ) {
			// determine the energy differential if we move the residues here
			// residue + 1
			int bondsBefore = countBonds(residue + 1, ii[residue+1], ii[residue], ii[residue+2]);
			int bondsAfter = countBonds(residue + 1, newNiplus1, ii[residue], newNiplus2);
			FreeEnergy_M += bondsBefore - bondsAfter;
			// residue + 2
			bondsBefore = countBonds(residue + 2, ii[residue+2], ii[residue+1], ii[residue+3]);
			bondsAfter = countBonds(residue + 2, newNiplus2, newNiplus1, ii[residue+3]);
			FreeEnergy_M += bondsBefore - bondsAfter;
			
			int n1=residue+1;
			int n2=residue+2;
			
			Lattice_M[ii[n1]] = LPempty;
			Lattice_M[ii[n2]] = LPempty;
			
			ii_M[n1] = newNiplus1;
			ii_M[n2] = newNiplus2;
			
			Lattice_M[newNiplus1] = sequenceCode[n1];
			Lattice_M[newNiplus2] = sequenceCode[n2];
			
			if(n1 > n2) {
				int tmp = n1;
				n1 = n2;
				n2 = tmp;
			}
			move_startResidue = n1;
			move_endResidue = n2;
			move_isPossible = true;
		}
	}
	return;
}

void Conformation::tryCornerMove(int residue) {
	move_isPossible = false;
	
	if( (ii[residue] - ii[residue-1]) != (ii[residue+1] - ii[residue]) ) {
        int freeN[DIM*2];
		int numFree;
		
		if( (numFree=findFreeNeighbours( residue-1, residue+1, freeN )) > 0 ) {
			int new_n1 = selectRandomFreeNeighbour(freeN, numFree);
			// determine the energy differential if we move the residue here
			int bondsBefore = countBonds(residue, ii[residue], ii[residue-1], ii[residue+1]);
			int bondsAfter = countBonds(residue, new_n1, ii[residue-1], ii[residue+1]);
			FreeEnergy_M += bondsBefore - bondsAfter;
			
			Lattice_M[ii[residue]] = LPempty;
			ii_M[residue] = new_n1;
			Lattice_M[new_n1] = sequenceCode[residue];
			
			move_startResidue = residue;
			move_endResidue = residue;
			move_isPossible = true;
            
		}
	}
	return;
}

#endif

#if MOVESET == PULLMOVES || MOVESET == MIXED

void Conformation::tryEndPullMove(int residue) {
	int bondsBefore = 0, bondsAfter = 0;
	int newLP1 = -1, newLP2 = -1;
	int startResidue = residue, endResidue = residue;
	int offset = 1;
	move_isPossible = false;
	
	if( residue == 0 ) {
		offset = -1;
	}
	
	int freeNeighbours[2*DIM];
	int numFree = findFreeNeighbours(residue, freeNeighbours);
	if( numFree < 1 ) return;
	
	newLP1 = selectRandomFreeNeighbour(freeNeighbours, numFree);
	
	numFree = findFreeNeighboursL(newLP1, freeNeighbours);
	if( numFree < 1 ) return;
	
	newLP2 = selectRandomFreeNeighbour(freeNeighbours, numFree);
	
	if( newLP2 == -1 ) return;
	
	int lastMovePos, lastPastMovePos, secondLastPastMovePos;
	int residueJ;
    
	Lattice_M[ii_M[residue]] = LPempty;
	Lattice_M[newLP2] = sequenceCode[residue];
	ii_M[residue] = newLP2;
    Lattice_M[ii_M[residue - offset]] = LPempty;
	Lattice_M[newLP1] = sequenceCode[residue - offset];
	ii_M[residue - offset] = newLP1;
	endResidue = residue - offset;
	secondLastPastMovePos = ii[residue];
	lastPastMovePos = ii[residue - offset];
	lastMovePos = newLP1;
	residueJ = residue - offset - offset;
	
	while( residueJ >= 0 && residueJ <= seqLength - 1 &&
#if DIM == 2
          ( abs(ii[residueJ] - lastMovePos) != L && abs(ii[residueJ] - lastMovePos) != 1 ) ) {
#endif
		
		endResidue = residueJ;
		Lattice_M[ii_M[residueJ]] = LPempty;
		ii_M[residueJ] = secondLastPastMovePos;
		Lattice_M[secondLastPastMovePos] = sequenceCode[residueJ];
        
		lastMovePos = secondLastPastMovePos;
		secondLastPastMovePos = lastPastMovePos;
		lastPastMovePos = ii[residueJ];
		residueJ = residueJ - offset;
	}
	
	if( startResidue > endResidue ) {
		int tmp = startResidue;
		startResidue = endResidue;
		endResidue = tmp;
	}
    
	for( int res = startResidue; res <= endResidue; res += 1 )
	{
		
		if( sequenceCode[res] >= 0 ) {
			int lp = ii[res];
			int olp;
			int oRes;
			
			for( int i = 0; i < 2*DIM; i += 1 )
			{
				
				olp = lp + DI[i];
				oRes = Lattice[olp];
				if( oRes >= 0 ) {
                    
					if( abs(oRes - res) > 1 ) {
						bondsBefore++;
                        
						if( oRes > res && oRes <= endResidue ) {
							bondsBefore--;
						}
					}
		    	}
			}
			
			for( int i = 0; i < 2*DIM; i += 1 )
			{
				lp = ii_M[res];
				olp = lp + DI[i];
				oRes = Lattice_M[olp];
				if( oRes >= 0 ) {
                    
					if( abs(oRes - res) > 1 ) {
						bondsAfter++;
                        
						if( oRes > res && oRes <= endResidue ) {
							bondsAfter--;
						}
					}
		    	}
			}
			
		}
	}
	
	FreeEnergy_M += bondsBefore - bondsAfter;
	move_startResidue = startResidue;
	move_endResidue = endResidue;
	move_isPossible = true;
}

void Conformation::tryPullMove(int residue) {
	int bondsBefore = 0, bondsAfter = 0;
	int offset = 1, Lpos = -1, Cpos = -1;
	int startResidue = residue, endResidue = residue;
	move_isPossible = false;
	
#if DIM == 2
	int residue2, offset2;
	int direction = random_number(1,2);
	if( direction == 1 ) offset = -offset;
	
	for(int i = 0; i < 2; i++) {
		residue2 = residue + offset;
		int d = abs(ii[residue] - ii[residue2]);
		if( d == L ) {
			offset2 = 1;
		} else {
			offset2 = L;
		}
        direction = random_number(1,2);
		if( direction == 1 ) {
			offset2 = -offset2;
		}
		if( Lattice[ii[residue2] + offset2] == LPempty
           && (Lattice[ii[residue] + offset2] == LPempty
               || ii[residue] + offset2 == ii[residue - offset] ) ) {
               
               Lpos = ii[residue2] + offset2;
               Cpos = ii[residue] + offset2;
               
           } else if( Lattice[ii[residue2] - offset2] == LPempty
                     && (Lattice[ii[residue] - offset2] == LPempty
                         || ii[residue] - offset2 == ii[residue - offset] )) {
                         
                         Lpos = ii[residue2] - offset2;
                         Cpos = ii[residue] - offset2;
                     }
		
		if( Lpos != -1 ) break;
        offset = -offset;
	}
    
#elif DIM == 3
	int residue2, offset2, offset2b;
	// randomly select which pull direction to try first
	int direction = random_number(1,2);
	if( direction == 1 ) offset = -offset;
	
	for(int i = 0; i < 2; i++) {
		residue2 = residue + offset;
        
		int d = abs(ii[residue] - ii[residue2]);
        
		if( d == L ) {
			if( random_number(1,2) == 1 ) {
				offset2 = 1; offset2b = LL;
			} else {
				offset2 = LL; offset2b = 1;
			}
		} else if( d == LL ) {
			if( random_number(1,2) == 1 ) {
				offset2 = 1; offset2b = L;
			} else {
				offset2 = L; offset2b = 1;
			}
		} else {
			// d = 1
			if( random_number(1,2) == 1 ) {
				offset2 = L; offset2b = LL;
			} else {
				offset2 = LL; offset2b = L;
			}
		}
		
		for(int j=0; j<2; j++) {
            
			// randomly select which direction to check first for L
			direction = random_number(1,2);
			if( direction == 1 ) {
				offset2 = -offset2;
			}
            
			if( Lattice[ii[residue2] + offset2] == LPempty
               && (Lattice[ii[residue] + offset2] == LPempty
                   || ii[residue] + offset2 == ii[residue - offset] ) ) {
                   
                   Lpos = ii[residue2] + offset2;
                   Cpos = ii[residue] + offset2;
                   
               } else if( Lattice[ii[residue2] - offset2] == LPempty
                         && (Lattice[ii[residue] - offset2] == LPempty
                             || ii[residue] - offset2 == ii[residue - offset] )) {
                             
                             Lpos = ii[residue2] - offset2;
                             Cpos = ii[residue] - offset2;
                         }
            
			if( Lpos != -1 ) break;
			offset2 = offset2b;
		}
		
		if( Lpos != -1 ) break;
        offset = -offset;
	}
    
#endif
	
	if( Lpos == -1 ) {
		return;
	}
	
	int lastMovePos, lastPastMovePos, secondLastPastMovePos;
	int residueJ;
	
	Lattice_M[ii_M[residue]] = LPempty;
	Lattice_M[Lpos] = sequenceCode[residue];
	ii_M[residue] = Lpos;
    if( ii[residue - offset] != Cpos  ) {
		
		Lattice_M[ii_M[residue - offset]] = LPempty;
		ii_M[residue - offset] = Cpos;
		Lattice_M[Cpos] = sequenceCode[residue - offset];
		endResidue -= offset;
		
		secondLastPastMovePos = ii[residue];
		lastPastMovePos = ii[residue - offset];
		lastMovePos = Cpos;
		residueJ = residue - offset - offset;
		
		while( residueJ >= 0 && residueJ <= seqLength - 1 &&
#if DIM == 2
              ( abs(ii[residueJ] - lastMovePos) != L && abs(ii[residueJ] - lastMovePos) != 1 ) ) {
#elif DIM == 3
            ( abs(ii[residueJ] - lastMovePos) != L && abs(ii[residueJ] - lastMovePos) != 1 ) && abs(ii[residueJ] - lastMovePos) != LL ) {
#endif
                
                endResidue -= offset;
                Lattice_M[ii_M[residueJ]] = LPempty;
                ii_M[residueJ] = secondLastPastMovePos;
                Lattice_M[secondLastPastMovePos] = sequenceCode[residueJ];
                
                lastMovePos = secondLastPastMovePos;
                secondLastPastMovePos = lastPastMovePos;
                lastPastMovePos = ii[residueJ];
                residueJ = residueJ - offset;
            }
            
        }
        
        if( startResidue > endResidue ) {
            int tmp = startResidue;
            startResidue = endResidue;
            endResidue = tmp;
        }
		
        for( int res = startResidue; res <= endResidue; res += 1 )
        {
            
            if( sequenceCode[res] >= 0 ) {
                int lp = ii[res];
                int olp;
                int oRes;
                
                for( int i = 0; i < 2*DIM; i += 1 )
                {
                    
                    olp = lp + DI[i];
                    oRes = Lattice[olp];
                    if( oRes >= 0 ) {
                        
                        if( abs(oRes - res) > 1 ) {
                            bondsBefore++;
                            
                            if( oRes > res && oRes <= endResidue ) {
                                bondsBefore--;
                            }
                        }
                    }
                }
                
                for( int i = 0; i < 2*DIM; i += 1 )
                {
                    lp = ii_M[res];
                    olp = lp + DI[i];
                    oRes = Lattice_M[olp];
                    if( oRes >= 0 ) {
                        
                        if( abs(oRes - res) > 1 ) {
                            bondsAfter++;
                            
                            if( oRes > res && oRes <= endResidue ) {
                                bondsAfter--;
                            }
                        }
                    }
                }
                
            }
        }
        
        FreeEnergy_M += bondsBefore - bondsAfter;
        move_startResidue = startResidue;
        move_endResidue = endResidue;
        move_isPossible = true;
    }
    
#endif
    
    void Conformation::applyMutation() {
        int residue;
        int lp;
        int bondedResidue;
        int boundaryHit = false;
        if( !move_isPossible ) return;
        for(residue=move_startResidue; residue <= move_endResidue; residue += 1) {
            Lattice[ii[residue]] = LPempty;
        }
        
        for(residue=move_startResidue; residue <= move_endResidue; residue += 1) {
            lp = ii_M[residue];
            int x = row(lp);
            int y = col(lp);
            
#if DIM == 2
            if( x<=3 || y<=3 || x>=L-3 || y>=L-3) {
                boundaryHit = true;
            }
#endif
            Lattice[lp] = sequenceCode[residue];
            ii[residue] = ii_M[residue];
        }
        
        
        FreeEnergy = FreeEnergy_M;
        
        if( boundaryHit ) {
            recenterOnLattice();
        }
        numResidues += move_endResidue - move_startResidue + 1;
        
    }
    
    void Conformation::resetMutation() {
        int residue;
        
        for(residue=move_startResidue; residue <= move_endResidue; residue += 1) {
            Lattice_M[ii_M[residue]] = LPempty;
        }
        
        for(residue=move_startResidue; residue <= move_endResidue; residue += 1) {
            ii_M[residue] = ii[residue];
            Lattice_M[ii_M[residue]] = sequenceCode[residue];
            FreeEnergy_M = FreeEnergy;
        }
    }
    void Conformation:: recenterOnLattice() {
        
        for(int i = 0; i<seqLength; i++) {
            Lattice[ii[i]]=LPempty;
        }
        
        int minRow = L + 1, minCol = L + 1, minHeight = L + 1, maxRow = 0, maxCol = 0, maxHeight = 0;
        
        for(int i=0; i<seqLength; i++) {
            int lrow = row(ii_M[i]);
            int lcol = col(ii_M[i]);
#if DIM == 3
            int lheight = height(ii_M[i]);
#endif
            
            if( lrow < minRow ) minRow = lrow;
            if( lrow > maxRow ) maxRow = lrow;
            if( lcol < minCol ) minCol = lcol;
            if( lcol > maxCol ) maxCol = lcol;
#if DIM == 3
            if( lheight < minHeight ) minHeight = lheight;
            if( lheight > maxHeight ) maxHeight = lheight;
#endif
        }
        
        // step 3) Determine the centroid of the embedding
        int r = minRow + (int)((maxRow - minRow)/2);
        int c = minCol + (int)((maxCol - minCol)/2);
#if DIM == 3
        int h = minHeight + (int)((maxHeight - minHeight)/2);
#endif
        
        int C = c + r * L;
#if DIM == 3
        C += h * LL;
#endif
        
        int delta = i0 - C;
        
        
        
        // step 4) copy residues from Lattice_M to Lattice, translating by delta
        for(int i=0; i<seqLength; i++) {
            Lattice[ii_M[i]+delta]=Lattice_M[ii[i]];
            ii[i]=ii_M[i]+delta;
        }
        
        // step 5) erase Lattice_M
        for(int i=0; i<seqLength; i++) {
            Lattice_M[ii_M[i]] = LPempty;
        }
        
        // step 6) copy Lattice to Lattice_M and ii to ii_M
        for(int i=0; i<seqLength; i++) {
            Lattice_M[ii[i]] = Lattice[ii[i]];
            ii_M[i] = ii[i];
        }
        
#ifdef TRACE_VERBOSE
        cout << endl << "Done translation move" << endl;
        printLatticePositions();
        cout << endl;
#endif
    }
    
    void Conformation:: printRelDirections() {
        for(int i=0; i<seqLength; i++)
            cout << Directions[i] << ",";
        cout << endl;
    }
    
    void Conformation:: printLatticePositions() {
        cout << "i0: (" << row(i0) << "," << col(i0) << ")" << endl;
        for(int i=0; i<seqLength; i++) {
            cout << "(" << row(ii[i]) << "," << col(ii[i]);
            cout << ")" << ",";
        }
        cout << endl;
    }
    
    void Conformation:: printRelLatticePositions() {
#if DIM == 2
        cout << "(0,0)" << ",";
        cout << endl;
#endif
        for(int i=1; i<seqLength; i++) {
            cout << "(" << row(ii[i])-row(ii[0]) << ","
            << col(ii[i])-col(ii[0])
            << ")" << ",";
            cout << endl;
        }
        cout << endl;
    }
    
    bool Conformation:: testLatticePositions() {
        bool flag;
        
        flag=false;
        for(int j=0; j<2*DIM; j++) {
            if(ii[0] == ii[1]+DI[j] ) {
                flag=true;
                break;
            }
        }
        
        if( flag == false ) {
            cerr << "1. Something's wrong" << endl;
            return false;
            //throw "Something's wrong";
        }
        
        for(int i=1; i<seqLength; i++) {
            flag = false;
            for(int j=0; j<2*DIM; j++) {
                if( ii[i] == ii[i-1] + DI[j] ) {
                    flag=true;
                    break;
                }
            }
            if( flag==false ) {
                cerr << "2. Something's wrong" << endl;
                return false;
                //throw "Something's wrong";
            }
        }
        
        // now make sure no lattice points are shared by two residues
        for(int i = 0; i<seqLength-1; i++) {
            for(int j=i+1; j<seqLength; j++) {
                if(ii[i] == ii[j]) {
                    cerr << "3. Something's wrong" << endl;
                    return false;
                }
            }
        }
        
        return true;
    }
    
    int Conformation:: calcFreeEnergy() {
        int freeE = 0;
        
        freeE += countBonds(0,ii[0],-1,ii[1]);
        
        for(int i = 1; i<seqLength-1; i++) {
            freeE += countBonds(i,ii[i],ii[i-1],ii[i+1]);
        }
        
        freeE += countBonds(seqLength-1,ii[seqLength-1],ii[seqLength-2],-1);
        return -freeE/2;
    }
    
    bool Conformation:: testFreeEnergy(int E) {
        int realEnergy = calcFreeEnergy();
        
        if( realEnergy != E ) {
            cout << "Actual Free Energy: " << realEnergy << endl;
            cout << "Reported Free Energy: " << E << endl;
        }
        return realEnergy == E;
    }
    
    void Conformation::reset() {
        i0 = I0;
        for(int i=0; i<seqLength; i++) {
            Lattice[ii[i]] = LPempty;
            Lattice_M[ii_M[i]] = LPempty;
        }
        
        for(int i=0; i<seqLength; i++) {
            Directions[i] = straight;
			bestRelDirections[i] = straight;
        }
        
        FixedPoint = 2;
        FreeEnergy = 0;
		searchResult.resetminsol();
		bestFreeEnergy = 0;
		
    }
    
    void Conformation::makeStraight() {
        reset();
        for(int i=0; i<seqLength; i++) {
            Lattice[i0 + i*L] = sequence[i];
			Lattice_M[i0 + i*L] = Lattice[i0 + i*L];
            ii[i] = i0 + i*L;
			ii_M[i] = ii[i];
            Directions[i] = straight;
        }
        Directions[0] = ground;
        Directions[1] = ground;
    }
    
    void Conformation::computeFreeEnergy() {
        FreeEnergy = calcFreeEnergy();
    }
    
    string Conformation::getContactString() {
        ostringstream cs;
        bool firstContact = true;
        
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] == B) continue;
            
            for(int j=i+2; j<seqLength; j++) {
                if( sequence[j] == B) continue;
                
                int distance = abs(ii[i] - ii[j]);
                
                if( distance == 1 || distance == L ) {
                    if(!firstContact) cs << ",";
                    cs << "(" << (i+1) << "," << (j+1) << ")";
                    firstContact = false;
                }
            }
        }
        
        return cs.str();
    }
    
    string Conformation::getSymmetricContactString() {
        ostringstream cs;
        bool firstContact = true;
        
        for(int i=seqLength-1; i>=0; i--) {
            if( sequence[i] == B) continue;
            
            for(int j=i-2; j>=0; j--) {
                if( sequence[j] == B) continue;
                
                int distance = abs(ii[i] - ii[j]);
                
#if DIM == 2
                if( distance == 1 || distance == L ) {
                    if(!firstContact) cs << ",";
                    cs << "(" << (seqLength-i) << "," << (seqLength-j) << ")";
                    firstContact = false;
                }
#endif
            }
        }
        
        return cs.str();
    }
    
    void Conformation::addToContactFrequencyMatrix(vector< vector< int > > &freqMatrix, int count) {
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] == B) continue;
            
            for(int j=i+2; j<seqLength; j++) {
                if( sequence[j] == B) continue;
                
                int distance = abs(ii[i] - ii[j]);
                
                if( distance == 1 || distance == L ) {
                    freqMatrix[i][j] += count;
                    freqMatrix[j][i] += count;
                }
            }
        }
    }

    double Conformation::relativeWWContactOrder() {
        int sumOfDistance = 0;
        int n = 0;
        int l = 0;
        
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] == B) continue;
            n++;
            
            for(int j=i+2; j<seqLength; j++) {
                if( sequence[j] == B) continue;
                
                int distance = abs(ii[i] - ii[j]);
                
#if DIM == 2
                if( distance == 1 || distance == L ) {
                    l++;
                    sumOfDistance += abs(i-j);
                }
                
#endif
            }
        }
        if( n == 0 || l == 0) return 0.0;
        
        return (1.0 / (double)(l * n)) * sumOfDistance;
    }
    
    double Conformation::absoluteWWContactOrder() {
        int n = 0;
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] != W ) n++;
        }
        
        return n * relativeWWContactOrder();
    }
    
    int Conformation::solventAccessibility() {
        int saCnt = 0;
        
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] == B) continue;
            
            for(int j=0; j<2*DIM; j++) {
                if( Lattice[ii[i] + DI[j]] == LPempty ) {
                    saCnt++;
                }
            }
        }
        
        return saCnt;
    }
    
    
    void Conformation::loadConformation(const string &relDirections, const int freeEnergy) {
        int relDir[MAX_SEQ_LENGTH];
        int len = relDirections.size();
        for(int i=0; i<len; i++) {
            switch(relDirections[i]) {
                case 's': 	relDir[i] = straight; 	break;
                case 'l':	relDir[i] = left; 		break;
                case 'r':	relDir[i] = right; 		break;
                case 'u':	relDir[i] = up; 		break;
                case 'd':	relDir[i] = down; 		break;
                case 'g':	relDir[i] = ground; 	break;
            }
        }
        loadConformation(relDir, freeEnergy);
    }
    
#if DIM == 2
    void Conformation::loadConformation(const int relDirections[], const int freeEnergy) {
        reset();
        Directions[0] = relDirections[0];
        Directions[1] = relDirections[1];
        Lattice[i0] = sequenceCode[0]; ii[0] = i0;
        Lattice[i0+L] = sequenceCode[1]; ii[1] = i0+L;
        
        int lp, dv;
        
        for(int i=2; i<seqLength; i++) {
            dv = ii[i-1] - ii[i-2];
            if( relDirections[i] == straight ) {
                lp = ii[i-1] + dv; 
            } else if( relDirections[i] == left ) {
                
                if( dv == L ) {
                    lp = ii[i-1] + 1;
                } else if ( dv == -L ) {
                    lp = ii[i-1] - 1;
                } else if ( dv == 1 ) {
                    lp = ii[i-1] - L;
                } else if( dv == -1 ) {
                    lp = ii[i-1] + L;
                } 					
            } else if( relDirections[i] == right ) {
                
                if( dv == L ) {
                    lp = ii[i-1] - 1;
                } else if ( dv == -L ) {
                    lp = ii[i-1] + 1;
                } else if ( dv == 1 ) {
                    lp = ii[i-1] + L;
                } else if( dv == -1 ) {
                    lp = ii[i-1] - L;
                }
            }
            Lattice[lp] = sequenceCode[i]; ii[i] = lp;
            Directions[i] = relDirections[i];
        }
        
        computeFreeEnergy();
        searchResult.setminsol(Directions, FreeEnergy, seqLength);
        bestFreeEnergy = FreeEnergy;
        for(int i=0; i<seqLength; i++) {
            bestRelDirections[i] = Directions[i];
            ii_M[i] = ii[i];
            Lattice_M[ii_M[i]] = Lattice[ii[i]];
        }
    }
    
    void Conformation::fastFold(double timeout) {
        cout << "Reseeding worst conformation - FreeEnergy: " << FreeEnergy << " -> ";
        makeStraight();
        stopwatch timer;
        timer.start();
        
        double oldTemp = getTemperature();
        setTemperature(190);
        
        while( timer.read() < timeout ) {
            runSimulation(1000);
        }
        
        loadConformation(bestRelDirections, bestFreeEnergy);
        setTemperature(oldTemp);
        cout << FreeEnergy << endl;
    }
    
    void Conformation:: printConformation() {
        
        for(int i=0; i<HashSize; i++) {
            Lattice_M[i] = 0;
        }
        
        int code;
        for(int i=0; i<seqLength; i++) {
            if( sequence[i] == W ) {
                code = 2*(i+1) - 1;
            } else {
                code = 2*(i+1);
            }
            Lattice_M[ ii[i] ] = code;
        }
        
        int col=0, minCol = L, maxCol = 0;
        int row=-1, minRow = L, maxRow = 0;
        for(int i=0; i<HashSize; i++) {
            if( i % L == 0 ) {
                col=0;
                row++;
            }
            if( Lattice[i] != LPempty ) {
                if( row < minRow ) minRow = row;
                if( row > maxRow ) maxRow = row;
                if( col < minCol ) minCol = col;
                if( col > maxCol ) maxCol = col;
            }
            col++;
        }
        
        cout << endl;
        
        for(int r=minRow; r<=maxRow; r++) {
            for(int c=minCol; c<=maxCol; c++) {
                if( r == 0 ) continue;
                int lp = r * L + c;
                if( Lattice_M[lp] == 0 || Lattice_M[lp-L] == 0) {
                    cout << "   ";
                } else {
                    int curCode = Lattice_M[lp];
                    int topCode = Lattice_M[lp-L];
                    int diffCode = curCode - topCode;
                    if( (curCode % 2) == 0 ) {
                        
                        if( diffCode == 2 || diffCode == 3 || diffCode == -1 || diffCode == -2 ) {
                            cout << "|  ";
                        } else {
                            cout << "   ";
                        }
                    } else {
                        if( diffCode == 1 || diffCode == 2 || diffCode == -2 || diffCode == -3 ) {
                            cout << "|  ";
                        } else {
                            
                            if( searchparameters::showwwbonds && (topCode % 2) == 1 ) {
                                cout << ":  ";
                            } else {
                                cout << "   ";
                            }
                        }
                    }
                }
            }
            cout << endl;
			
            for(int c=minCol; c<maxCol; c++) {
                int lp = r * L + c;
                int curCode = Lattice_M[lp];
                
                if( curCode == 0 ) {
                    cout << "   ";
                    continue;
                }
                
                int rightCode = Lattice_M[lp+1];
                
                if( rightCode == 0 ) {
                    if( (curCode % 2) == 0 ) {
                        if( curCode == 2 || curCode == 2 * seqLength )
                            cout << "B  ";
                        else
                            cout << "b  ";
                    } else {
                        if( curCode == 1 || curCode == 2 * seqLength - 1)
                            cout << "W  ";
                        else
                            cout << "w  ";
                    }
                    continue;
                }
                
                int diffCode = curCode - rightCode;
                
                if( (curCode % 2) == 0 ) {
                    if( curCode == 2 || curCode == 2 * seqLength )
                        cout << "B";
                    else
                        cout << "b";
                    
                    if( diffCode == -1 || diffCode == -2 || diffCode == 2 || diffCode == 3 ) {
                        cout << "--";
                    } else {
                        cout << "  ";
                    }
                } else {
                    if( curCode == 1 || curCode == 2 * seqLength - 1)
                        cout << "W";
                    else
                        cout << "w";
                    
                    if( diffCode == 1 || diffCode == 2 || diffCode == -2 || diffCode == -3 ) {
                        cout << "--";
                    } else {
                        if( searchparameters::showwwbonds && (rightCode % 2) == 1 ) {
                            cout << "..";
                        } else {
                            cout << "  ";
                        }
                    }
                }
                
            }
			
            int lp = r * L + maxCol;
            int curCode = Lattice_M[lp];
			
            if( curCode == 0 ) {
                cout << " ";
            } else if( (curCode % 2) == 0 ) {
                if( curCode == 1 || curCode == 2 * seqLength - 1)
                    cout << "B";
                else
                    cout << "b";
            } else {
                if( curCode == 1 || curCode == 2 * seqLength - 1)
                    cout << "W";
                else
                    cout << "w";
            }
			
            cout << endl;
        }
        
        cout << endl;
    }
#endif
#endif


