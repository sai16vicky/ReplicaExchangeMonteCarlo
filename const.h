#ifndef CONST_H
#define CONST_H
#define TRACE_SQT_GENERAL_IMMEDIATE		// general sqt trace, dumped to file on the
#define TRACE_ARI_MOVES // Accept, Reject, Impossible moves
#define TRACE_VERBOSE	// Verbose output messages
#define LOCAL		1
#define PULLMOVES   2
#define MIXED		3
#ifndef DIM
#define DIM		2
#endif
#ifndef MOVESET
#define MOVESET	LOCAL
#endif
#if DIM == 3
#define L			256
#define LL			65536
#define HashSize	16777216
//	#define L			128
//	#define LL			16384
//	#define HashSize	2097152
#define MinBoundary LL
#define MaxBoundary HashSize - LL
#define I0			(L/2*LL + L/2*L + L/2)
#elif DIM == 2
#define L			2048
#define LL			4096
#define HashSize	4194304
#define MinBoundary LL
#define MaxBoundary HashSize - LL
#define I0			(L/2*L + L/2)
#endif
#if defined LOWMEM
#define L		128
#define LL		512
#define HashSize 	16384
#define MinBoundary 	LL
#define MaxBoundary	HashSize - LL
#define I0		(L/2*L + L/2)
#MAX_NUM		128
#endif

static const double K_b = 0.0019872;
static const double K_b2 = 0.0001679010;
static const int MAX_SEQ_LENGTH = 1024;
static const int MAX_NUM_REPLICAS = 400;
static const int MAX_NUM_ENERGY_VALUES = 1000;
#define straight		1
#define right			2
#define left 			4
#define up				8
#define down			16
#define ground			32
#define LPempty			-1

#define W 				1
#define B 				-2

#endif
