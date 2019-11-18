#include "freexy.h"

extern GLOBAL g;
extern NWK nwk[NT];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// does the book keeping for the adjacency lists and matrices when deleting a link
// i is the index of n in the nw->node struct, j is the index at the other
// end of the link

void adjacency_delete (NODE *n, uint16_t i, uint16_t j) {
	NODE *ni = n + i;
	uint16_t k = ni->nb[--(ni->a[i])]; // getting the index of the node to repl j by, decreasing degree

	ni->a[k] = ni->a[j]; // update k's position in i's nb list
	ni->nb[ni->a[j]] = k; // overwrite j
	ni->a[j] = NONE; // reset the adjacency matrix
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// does the book keeping for the adjacency lists and matrices when adding a link
// i is the index of n in the nw->node struct, j is the index at the other
// end of the link

void adjacency_add (NODE *n, uint16_t i, uint16_t j) {
	NODE *ni = n + i;

	ni->a[j] = ni->a[i]; // because j will be added to the end of i's adjacency list (n->a[i] is i's degree)
	ni->nb[(ni->a[i])++] = j; // write j and increase the degree
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// help routine to swap links in the link list

void lswap (LINK *l0, LINK *l1) {
	uint16_t tmp;

	tmp = l0->left;
	l0->left = l1->left;
	l1->left = tmp;

	tmp = l0->right;
	l0->right = l1->right;
	l1->right = tmp;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void add_link (NWK *nw, LINK *ladd) {

	adjacency_add(nw->node, ladd->left, ladd->right);
	adjacency_add(nw->node, ladd->right, ladd->left);

#ifdef FREER
	lswap(ladd, nw->link + nw->m);
#endif

	nw->m++;
}

#ifdef FREER //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void delete_link (NWK *nw, LINK *ldel) {
	LINK *lrepl;

	nw->m--;
	lrepl = nw->link + nw->m;

	// delete link in adjacency lists
	adjacency_delete(nw->node, ldel->left, ldel->right);
	adjacency_delete(nw->node, ldel->right, ldel->left);

	// delete link in link struct
	lswap(lrepl, ldel);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// try to update links by adding a random one

void update_link_add (NWK *nw, LINK *ladd) {
	double de, dde;

	dde = E(nw->node[ladd->left].s - nw->node[ladd->right].s);
	de = dde - 1.0;

	// metropolis condition
	if ((de < 0.0) || (g.beta[nw->t] * de < g.mlog[pcg_16()])) {
		// updating energy
		nw->energy += dde;

		add_link(nw, ladd);
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// try to update links by deleting a random one

void update_link_delete (NWK *nw, LINK *ldel) {
	double de, dde;

	// getting energy difference
	dde = - E(nw->node[ldel->left].s - nw->node[ldel->right].s);
	de = 1.0 + dde;

	// metropolis condition
	if ((de < 0.0) || (g.beta[nw->t] * de < g.mlog[pcg_16()])) {
		// updating energy
		nw->energy += dde;

		delete_link(nw, ldel);
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void update_link_add_delete (NWK *nw) {
	int i, ilim = nw->m;

	for (i = ilim; i < (int) g.mmax; i++) if (coin_flip())
		update_link_add(nw, nw->link + i);

	for (i = ilim - 1; i >= 0; i--) if (coin_flip())
		update_link_delete(nw, nw->link + i);
}

#endif // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void swap_link (NWK *nw, LINK *ldel, LINK *ladd) {

	adjacency_delete(nw->node, ldel->left, ldel->right);
	adjacency_delete(nw->node, ldel->right, ldel->left);

	adjacency_add(nw->node, ladd->left, ladd->right);
	adjacency_add(nw->node, ladd->right, ladd->left);

	lswap(ladd, ldel);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// sweep all links and try to update them by swapping

void update_link_swap (NWK *nw) {
	unsigned int i;
	double de;
	LINK *ladd, *ldel;
#ifdef FREER
	if (nw->m == 0) return;
#else
	LINK ltmp;
#endif

	ldel = nw->link + pcg_32_bounded(nw->m);

	for (i = 0; i < nw->m; i++) {
#ifdef FREER
		if ((nw->m < 1) || (nw->m >= g.mmax)) return;
		ladd = nw->link + nw->m + pcg_32_bounded(g.mmax - nw->m);
#else
		ladd = &ltmp;
		do {
			ladd->left = pcg_16_bounded(g.n);
			ladd->right = pcg_16_bounded(g.n);
		} while ((ladd->left >= ladd->right) || (nw->node[ladd->left].a[ladd->right] != NONE));
#endif

		// getting energy difference
		de = E(nw->node[ladd->left].s - nw->node[ladd->right].s) - E(nw->node[ldel->left].s - nw->node[ldel->right].s);

		// metropolis condition
		if ((de < 0.0) || (g.beta[nw->t] * de < g.mlog[pcg_16()])) {
			// updating energy
			nw->energy += de;

			swap_link(nw, ldel, ladd);
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// makes a sweep through all nodes and attempts to update them

void update_spin (NWK *nw) {
	unsigned int i, j;
	uint16_t v, su, snew;
	double de;
	NODE *n;

	// Fischer-Yates shuffling
	for (i = g.n; i > 0; ) {
		j = pcg_16_bounded(i--);
		v = g.node[i];
		g.node[i] = g.node[j];
		g.node[j] = v;
	}

	for (i = 0; i < g.n; i++) {
		v = g.node[i]; // go thru all nodes in random order
		n = nw->node + v;
		snew = n->s + pcg_16_bounded(g.maxang[nw->t]) - (g.maxang[nw->t] / 2); // get trial angle

		// calculate energy difference
		for (j = 0, de = 0.0; j < n->a[v]; j++) {
			su = nw->node[n->nb[j]].s;
			de += E(snew - su) - E(n->s - su);
		}

		// metropolis condition
		if ((de < 0.0) || (g.beta[nw->t] * de < g.mlog[pcg_16()])) {
			nw->energy += de;
			n->s = snew; 
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
