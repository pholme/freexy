// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <stdint.h>
#include <float.h>
#include <time.h>
#include <assert.h>

#define NSAVE 100
#define NTHERM 10000 // NSWEEP for initial thermalization
#define NSWEEP 100 // number of sweeps (N updates) per x-mc trial
#define NSPINXTRA 10 // how many spin sweeps per links swaps / 2
#define NT 100 // number of temperature levels
#ifdef FREER
#define TMIN 0.1 // lowest temperature
#define TMAX 10.0 // highest temperature
#else
#define TMIN 0.1 // lowest temperature
#define TMAX 100.0 // highest temperature
#endif
#define TRATIO exp(log(TMAX/TMIN)/(NT-1)) // makes the lowest temperature TMIN
#define NANG 0x10000 // angles coded as a 16-bit integer (to facilitate tabulizing)
#define ANGMAX 0xFFFF // largest angle
#define ANGMIN 0 // smallest trial angle (for infinite FERMIFAC)
#define FERMIMID (0.6*NT) // inflection point for Fermi function for trial angles
#ifdef FREER
#define FERMIFAC 5.0 // setting the shape of the Fermi function for trial angles
#else
#define FERMIFAC 10.0 // setting the shape of the Fermi function for trial angles
#endif
#define TRAVERSE_FAC(x) ((x) ? 0.25 : 0.5) // for the convergence criterion, longer first time for thermalization
#define NONE 0xFFFF

// result vector
#define MAG1 0
#define MAG2 1 // for Binder's cumulant
#define MAG4 2
#define NRG1 3 // energy
#define NRG2 4 // energy square for CV
#define NRES 5
#define DIA1 0
#define DIA2 1
#define NC1 2 // number of components
#define NC2 3 // number of components
#define SLARGEST 4 // size of largest component
#define S2LARGEST 5 // size of second largest component
#define M1 6 // number of links
#define M2 7 // square number of links
#ifdef FREER
#define NIRES 8
#else
#define NIRES 6
#endif

#define T(x) (TMIN*pow(TRATIO,(x)))
#define E(x) g.energy[(uint16_t) (x)]
#define SQ(x) ((x) * (x))

typedef struct GLOBAL {
	// NETWORK SPECS
	unsigned int n, mmax, m0; // number of nodes, max number of links, initial number of links
	double k; // avg. (for freer, initial) degree
	double beta[NT]; // temperature sequence
	unsigned int nwk_with_t[NT]; // for book-keeping of the replicas
	uint16_t *node; // all nodes in random order
	uint16_t *save; // scratch space
	uint16_t *cluster_members;
	// to speed up exp and log
	double energy[NANG], mlog[NANG];
	long double cos[NANG], sin[NANG];
	uint16_t maxang[NT];
	// for results
	unsigned int navg[NT];
	unsigned long ires[NT][NIRES];
	long double res[NT][NRES];
} GLOBAL;

typedef struct NODE {
	uint16_t s;
	uint16_t *nb; // the list of neighbors, the diagonal stores the degrees
	uint16_t *a; // the adjacency matrix's row at this node. this stores the index of the node in nb (this index starts from 1, so adjacency matrix equal 0 means no link, as usual)
	uint16_t dist, cluster;
} NODE;

typedef struct LINK {
	uint16_t left, right;
} LINK; // always left < right

typedef struct NWK { // info about the different copies of the network
	unsigned int t, visited[NT]; // current temperature, set of temperatures visited
	uint32_t m;
	double energy; // the true energy is this energy - m (the number of links (to increase precision))
	NODE *node;
	LINK *link;
} NWK;

extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded (uint32_t);
extern uint16_t pcg_16 ();
extern uint16_t pcg_16_bounded (uint16_t);
extern uint32_t coin_flip ();

extern void add_link (NWK *, LINK *);
extern void update_spin (NWK *);
extern void update_link_swap (NWK *);
#ifdef FREER
extern void update_link_add_delete (NWK *);
#endif

extern void measure (NWK *);


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
