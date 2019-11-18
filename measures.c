#include "freexy.h"

extern GLOBAL g;
extern NWK nwk[NT];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void measure (NWK *nw) {
	uint16_t me, you, at;
	unsigned int i, nc = 0, slargest = 0, s2largest = 0, dia = 0, first;
	unsigned int write_at, read_at;
	long double smag = 0.0, d, ssin = 0.0, scos = 0.0;
	NODE *n = nw->node;

	for (i = 0; i < g.n; i++) n[i].cluster = n[i].dist = NONE;

	for (me = 0; me < g.n; me++) {
		if (n[me].cluster == NONE) {
			ssin = scos = 0.0;
			first = 1;
			n[me].cluster = nc++;
		} else first = 0;

		write_at = 1;
		read_at = 0;
		g.cluster_members[0] = me;
		n[me].dist = 0;

		while (read_at < write_at) {
			at = g.cluster_members[read_at++];
			if (first) {
				n[at].cluster = n[me].cluster;
				ssin += g.sin[n[at].s];
				scos += g.cos[n[at].s];
			}
			for (i = 0; i < n[at].a[at]; i++) {
				you = n[at].nb[i];
				if (n[you].dist == NONE) {
					n[you].dist = n[at].dist + 1;
					if (n[you].dist > dia) dia = n[you].dist;
					g.cluster_members[write_at++] = you;
				}
			}
		}
		for (i = 0; i < write_at; i++) n[g.cluster_members[i]].dist = NONE;

		if (first) {
			smag += sqrtl(SQ(scos) + SQ(ssin));

			if (write_at > slargest) {
				s2largest = slargest;
				slargest = write_at;
			} else if (write_at > s2largest) s2largest = write_at;
		}
	}

	smag /= g.n;

	g.navg[nw->t]++;
	g.res[nw->t][MAG1] += smag;
	d = SQ(smag);
	g.res[nw->t][MAG2] += d;
	g.res[nw->t][MAG4] += SQ(d);
	g.ires[nw->t][SLARGEST] += slargest;
	g.ires[nw->t][S2LARGEST] += s2largest;
	d = nw->energy - nw->m;
	g.res[nw->t][NRG1] += d;
	g.res[nw->t][NRG2] += SQ(d);
#ifdef FREER
	g.ires[nw->t][M1] += nw->m;
	g.ires[nw->t][M2] += SQ((double) nw->m);
#endif
	g.ires[nw->t][DIA1] += dia;
	g.ires[nw->t][DIA2] += SQ(dia);
	g.ires[nw->t][NC1] += nc;
	g.ires[nw->t][NC2] += SQ(nc);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
