// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for temporal network SIR by Petter Holme (2018)

// this file contains the random number generator, derived from the PCG
// RNG v0.94 http://www.pcg-random.org under the Apache License 2.0
// http://www.apache.org/licenses/LICENSE-2.0

// 32-bit Output, 64-bit State: PCG-XSH-RS

#include <stdint.h>

extern uint64_t pcg_state; // for RNG seed/state

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t pcg_32 () {
	uint64_t state = pcg_state;
	uint32_t value, rot;

	pcg_state = pcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
	value = ((state >> 18u) ^ state) >> 27u;
	rot = state >> 59u;
	return (value >> rot) | (value << ((- rot) & 31));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t pcg_32_bounded (uint32_t bound) {
	uint32_t threshold = -bound % bound, r;

	for ( ; ; ) {
		r = pcg_32();
		if (r >= threshold) return r % bound;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint16_t pcg_16 () {
	static uint32_t exist, rmem;

	if (exist) {
		exist = 0;
		return (uint16_t) (rmem >> 16);
	}
	exist = 1;
	rmem = pcg_32();
	return (uint16_t) rmem;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint16_t pcg_16_bounded (uint16_t bound) {
	uint16_t threshold = -bound % bound, r;

	for ( ; ; ) {
		r = pcg_16();
		if (r >= threshold) return r % bound;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;

uint32_t coin_flip () {
	static uint32_t mask, rmem;

	if (mask == 0) {
		rmem = pcg_32();
		mask = 1;
	}

	mask <<= 1;

	return rmem & mask;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
