#include "freexy.h"

GLOBAL g;
uint64_t pcg_state; // for RNG
NWK nwk[NT];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine initializes one network for every temperature

void setup_nwk () {
	unsigned int i, it;
	uint16_t v, u;
	LINK *ladd;
	NWK *nw;

	for (nw = nwk, it = 0; it < NT; nw++, it++) {
#ifdef FREER
		nw->link = calloc(g.mmax, sizeof(LINK));
		for (i = 0, u = 1; u < g.n; u++) 
			for (v = 0; v < u; v++, i++) {
				nw->link[i].left = v;
				nw->link[i].right = u;
			}
#else
		nw->link = calloc(g.m0, sizeof(LINK));
#endif
		nw->node = calloc(g.n, sizeof(NODE));

		for (v = 0; v < g.n; v++) {
			nw->node[v].a = calloc(g.n, sizeof(uint16_t));
			for (u = 0; u < v; u++)
				nw->node[v].a[u] = nw->node[u].a[v] = NONE;
			nw->node[v].nb = malloc(g.n * sizeof(uint16_t));
			nw->node[v].s = pcg_16();
		}

		nw->t = it;
		nw->m = 0;

		for (i = 0; i < g.m0; i++) {
			ladd = nw->link + i;
#ifndef FREER
			do {
				ladd->left = pcg_16_bounded(g.n);
				ladd->right = pcg_16_bounded(g.n);
			} while ((ladd->left >= ladd->right) || (nw->node[ladd->left].a[ladd->right] != NONE));
#endif
			nw->energy += E(nw->node[ladd->left].s - nw->node[ladd->right].s);
			add_link(nw, ladd);
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// saves the network and spin states to file

void save (uint16_t it) {
	unsigned int i;
	FILE *fp;
	char fname[100];

#ifdef FREER
	sprintf(fname, "outer/%02u/%u.txt", it, g.n);
#else
	sprintf(fname, "out/%02u/%3.1lf/%u.txt", it, g.k, g.n);
#endif

	fp = fopen(fname, "a");
	if (!fp) {
		fprintf(stderr, "can't open save file (%s)\n", fname);
		exit(2);
	}


	for (i = 0; i < NRES; i++) g.res[it][i] /= g.navg[it];

	fprintf(fp, "%.30Lg", g.res[it][NRG1]); // energy
	fprintf(fp, " %.30Lg", g.res[it][MAG1]); // magnetization
	fprintf(fp, " %.30Lg", 1 - g.res[it][MAG4] / (3 * SQ(g.res[it][MAG2]))); // binder
	fprintf(fp, " %.15Lg", (g.res[it][NRG2] - SQ(g.res[it][NRG1])) * SQ(g.beta[it])); // specific heat
	fprintf(fp, " %.15Lg", (g.res[it][MAG2] - SQ(g.res[it][MAG1])) * g.beta[it]); // susceptibility
	fprintf(fp, " %g", g.ires[it][DIA1] / (double) g.navg[it]); // diameter
	fprintf(fp, " %g", g.ires[it][NC1] / (double) g.navg[it]); // number of components
	fprintf(fp, " %g", g.ires[it][SLARGEST] / (double) g.navg[it]); // size of largest component
	fprintf(fp, " %g", g.ires[it][S2LARGEST] / (double) g.navg[it]); // size of 2nd largest component

#ifdef FREER
	fprintf(fp, " %g", g.ires[it][M1] / (double) g.navg[it]); // number of links
#endif

	fprintf(fp, "\n");

	fclose(fp);

	g.navg[it] = 0;
	for (i = 0; i < NRES; i++) g.res[it][i] = 0.0;
	for (i = 0; i < NIRES; i++) g.ires[it][i] = 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// testing for a swap of temperatures according to the exchange MC criterion

void exchange () {
	int i, j0, j1;
	static int i0; // used to swap between odd / even exchange pairs
	double de;

	for (i = i0; i < NT - 1; i += 2) {
		j0 = g.nwk_with_t[i];
		j1 = g.nwk_with_t[i + 1];
		de = (g.beta[i] - g.beta[i + 1]) *
			(nwk[g.nwk_with_t[i + 1]].energy - nwk[g.nwk_with_t[i]].energy);

		if ((de < 0.0) || (de < g.mlog[pcg_16()])) {
			nwk[j0].visited[nwk[j0].t = i + 1] = 1;
			nwk[j1].visited[nwk[j1].t = i] = 1;

			g.nwk_with_t[i] = j1;
			g.nwk_with_t[i + 1] = j0;
		}
	}

	i0 = (i0) ? 0 : 1; 
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// checking that all networks have traversed a fraction TRAVERSE_FAC of all
// the temperatures (which is our criteria for saving the state)

int not_enough () {
	unsigned int i, it, s;
	static unsigned int not_first;

	for (it = 0; it < NT; it++) {
		for (s = i = 0; i < NT; i++) if (nwk[it].visited[i]) s++;
		if (s < NT * TRAVERSE_FAC(not_first)) return 1;
	}

	not_first = 1;

	for (it = 0; it < NT; it++) for (i = 0; i < NT; i++) nwk[it].visited[i] = 0;

	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// To get trial angles giving intermediate acceptance ratios (aiming for 50%
// acceptance for average nodes), this is all very sketchy and could be improved.
// It does not affect the quality of the output, but it with better tuned
// angles one could decrease the number of updates to get reasonably
// uncorrelated samples.

uint16_t get_maxang (unsigned int it) {

	return (uint16_t) rint(ANGMAX - (ANGMAX - ANGMIN) /
		(1.0 + exp(FERMIFAC * ((int) it - FERMIMID) / (double) NT)));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// main function handling i/o
// input limits:
// 2 < n < 2^15
// 1 < k < n
// 0 < seed < 2^64

int main (int argc, char *argv[]) {
	unsigned int i, j, k, it;
	double x;

#ifdef FREER
	if (argc != 3) {
		fprintf(stderr, "usage: ./freerxy [n] [seed]\n");
		return 1;
	}
	pcg_state = strtoull(argv[2], NULL, 10); // reading RNG seed / state
#else
	if (argc != 4) {
		fprintf(stderr, "usage: ./freexy [n] [k] [seed]\n");
		return 1;
	}
	g.k = atof(argv[2]);
	pcg_state = strtoull(argv[3], NULL, 10); // reading RNG seed / state
#endif

	// initialization - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	g.n = atoi(argv[1]);
	g.mmax = (g.n * (g.n - 1)) / 2;
	g.save = malloc(2 * g.mmax * sizeof(uint16_t));
	g.node = malloc(g.n * sizeof(uint16_t));
	g.cluster_members = malloc(g.n * sizeof(uint16_t));
	for (i = 0; i < g.n; i++) g.node[i] = i;

#ifdef FREER
	g.m0 = g.mmax / 2;
#else
	g.m0 = rint(g.n * g.k / 2.0);
#endif

	for (it = 0; it < NT; it++) {
		g.beta[it] = 1.0 / T(it);
		g.nwk_with_t[it] = it;
		g.maxang[it] = get_maxang(it);
	}

	for (i = 0; i < NANG; i++) { // tabulizing for faster calculation
		x = sin((M_PI * i) / NANG); // using 1-cos rather than -cos (and the formula 1 - cos x = 2 sin^2 x/2) to avoid subtracting values very close to 1
		g.energy[i] = 2.0 * x * x;
#ifdef FREER
		g.energy[i] /= g.n; // to make energy extensive
#endif
		g.sin[i] = sinl((2 * 3.14159265358979323846264338327950288419716939937510L * i) / NANG);
		g.cos[i] = cosl((2 * 3.14159265358979323846264338327950288419716939937510L * i) / NANG);
		g.mlog[i] = (i) ? -log(i / (double) NANG) : DBL_MAX;
	}

	setup_nwk();

	// run  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	for (it = 0; it < NT; it++) for (i = 0; i < NTHERM; i++) {
		for (j = 0; j < NSPINXTRA; j++) update_spin(nwk + it);
#ifdef FREER
		update_link_add_delete(nwk + it);
#endif
		update_link_swap(nwk + it);
		for (j = 0; j < NSPINXTRA; j++) update_spin(nwk + it);
	}

	for (i = 0; i < NSAVE; i++) {
		do {
			for (it = 0; it < NT; it++) {
				for (j = 0; j < NSWEEP; j++) {
#ifdef FREER
					update_link_add_delete(nwk + it);
#endif
					update_link_swap(nwk + it);
					for (k = 0; k < NSPINXTRA; k++) update_spin(nwk + it);
				}
				measure(nwk + it);
			}
			exchange();
		} while (not_enough());

		for (it = 0; it < NT; it++) save(it);
	}

	// cleaning  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	for (it = 0; it < NT; it++) {
		for (i = 0; i < g.n; i++) {
			free(nwk[it].node[i].a);
			free(nwk[it].node[i].nb);
		}
		free(nwk[it].node);
		free(nwk[it].link);
	}
	free(g.save); free(g.node); free(g.cluster_members);

	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
