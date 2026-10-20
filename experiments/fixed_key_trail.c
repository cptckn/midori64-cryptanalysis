/* Fixed-key probability of the four-round Midori-64 characteristic
 *   2000000000200002 -> ... -> 2220000022022022  (DDT-product 2^-32)
 * by conditional sampling.
 *
 * A pair follows the characteristic iff the LEFT text's value at every active
 * S-box lies in the conforming set of that transition (differences through the
 * linear layers are then forced). We sample the internal state x0 = P xor WK
 * directly (WK only relabels plaintexts), force the three round-1 active cells
 * into triples that already satisfy rounds 1 and 2 (enumerated exactly), and
 * estimate rounds 3-4 by sampling.
 *
 *   gcc -O3 -march=native -o fixed_key_trail fixed_key_trail.c -lm
 *   ./fixed_key_trail key  <32-hex master key> <samples> <seed>
 *   ./fixed_key_trail multi <number of random keys> <samples per key> <seed>
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int S[16] = {0xC,0xA,0xD,0x3,0xE,0xB,0xF,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6};
static const int PERM[16] = {0,10,5,15,14,4,11,1,9,3,12,6,7,13,2,8};
static const int BETA[3][16] = {{0,0,0,1,0,1,0,1,1,0,1,1,0,0,1,1},
                                {0,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0},
                                {1,0,1,0,0,1,0,0,0,0,1,1,0,1,0,1}};
static const int R4_CELLS[9] = {0,1,2,8,9,11,12,14,15};
static int c24[16], c42[16];          /* conforming inputs for 2->4 and 4->2 */
static uint64_t rng_state;

static uint64_t rnd(void) {           /* xorshift64 */
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}

/* x = MixColumn(ShuffleCell(y)) xor rk */
static void linear(const int *y, const int *rk, int *x) {
    int z[16];
    for (int i = 0; i < 16; i++) z[i] = y[PERM[i]];
    for (int c = 0; c < 4; c++) {
        int t = z[4*c] ^ z[4*c+1] ^ z[4*c+2] ^ z[4*c+3];
        for (int r = 0; r < 4; r++) x[4*c+r] = t ^ z[4*c+r] ^ rk[4*c+r];
    }
}

typedef struct { int good_triples; long long samples, pass_r3, hits; } result_t;

static result_t run(const char *key, long long n) {
    int k0[16], k1[16], rk[3][16], tri[64][3];
    result_t res = {0, n, 0, 0};
    for (int i = 0; i < 16; i++) {
        char a[2] = {key[i], 0}, b[2] = {key[16 + i], 0};
        k0[i] = (int)strtol(a, 0, 16); k1[i] = (int)strtol(b, 0, 16);
    }
    for (int r = 0; r < 3; r++)
        for (int i = 0; i < 16; i++) rk[r][i] = (r % 2 ? k1[i] : k0[i]) ^ BETA[r][i];
    /* round 1: cells 0,10,15 must be in c24; round 2: x1[2] = S(x0[0])^S(x0[10])^S(x0[15])^RK0[2] in c42 */
    for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) for (int c = 0; c < 16; c++)
        if (c24[a] && c24[b] && c24[c] && c42[S[a] ^ S[b] ^ S[c] ^ rk[0][2]]) {
            tri[res.good_triples][0] = a; tri[res.good_triples][1] = b; tri[res.good_triples][2] = c;
            res.good_triples++;
        }
    if (res.good_triples == 0) return res;   /* impossible already at round 2 */
    for (long long s = 0; s < n; s++) {
        int x0[16], y0[16], x1[16], y1[16], x2[16], y2[16], x3[16];
        uint64_t r = rnd();
        for (int i = 0; i < 16; i++) { x0[i] = r & 15; r >>= 4; }
        int t = (int)(rnd() % (uint64_t)res.good_triples);
        x0[0] = tri[t][0]; x0[10] = tri[t][1]; x0[15] = tri[t][2];
        for (int i = 0; i < 16; i++) y0[i] = S[x0[i]];
        linear(y0, rk[0], x1);
        for (int i = 0; i < 16; i++) y1[i] = S[x1[i]];
        linear(y1, rk[1], x2);
        if (!(c24[x2[12]] && c24[x2[13]] && c24[x2[15]])) continue;
        res.pass_r3++;
        for (int i = 0; i < 16; i++) y2[i] = S[x2[i]];
        linear(y2, rk[2], x3);
        int ok = 1;
        for (int k = 0; k < 9 && ok; k++) ok = c42[x3[R4_CELLS[k]]];
        res.hits += ok;
    }
    return res;
}

static double trail_log2(result_t r) {    /* log2 of the estimated fixed-key trail probability */
    if (r.good_triples == 0 || r.hits == 0) return -INFINITY;
    return log2(r.good_triples / 4096.0) + log2((double)r.hits / r.samples);
}

int main(int argc, char **argv) {
    for (int v = 0; v < 16; v++) { c24[v] = ((S[v] ^ S[v ^ 2]) == 4); c42[v] = ((S[v] ^ S[v ^ 4]) == 2); }
    if (argc == 5 && !strcmp(argv[1], "key")) {
        rng_state = strtoull(argv[4], 0, 10) | 1;
        result_t r = run(argv[2], atoll(argv[3]));
        if (r.good_triples == 0) { printf("%s: impossible (round-2 condition fails for every pair)\n", argv[2]); return 0; }
        printf("%s: P(R1,R2)=2^%.2f  P(R3|.)=2^%.2f  P(R4|.)=2^%.2f  hits=%lld/%lld  P_trail ~ 2^%.2f\n",
               argv[2], log2(r.good_triples / 4096.0), log2((double)r.pass_r3 / r.samples),
               r.pass_r3 ? log2((double)r.hits / r.pass_r3) : -INFINITY, r.hits, r.samples, trail_log2(r));
        return 0;
    }
    if (argc == 5 && !strcmp(argv[1], "multi")) {
        int nkeys = atoi(argv[2]); long long n = atoll(argv[3]);
        rng_state = strtoull(argv[4], 0, 10) | 1;
        int imp_r2 = 0, no_hit = 0, seen = 0; char key[33];
        for (int k = 0; k < nkeys; k++) {
            snprintf(key, sizeof key, "%016llx%016llx", (unsigned long long)rnd(), (unsigned long long)rnd());
            result_t r = run(key, n);
            if (r.good_triples == 0) { imp_r2++; continue; }
            if (r.hits == 0) { no_hit++; continue; }
            seen++;
            printf("  key %s  P_trail ~ 2^%.2f\n", key, trail_log2(r));
        }
        printf("%d keys: impossible at round 2: %d | no hit in %lld samples: %d | trail observed: %d (fraction 2^%.2f)\n",
               nkeys, imp_r2, n, no_hit, seen, seen ? log2((double)seen / nkeys) : -INFINITY);
        return 0;
    }
    fprintf(stderr, "usage: %s key <32hex> <samples> <seed> | %s multi <nkeys> <samples> <seed>\n", argv[0], argv[0]);
    return 1;
}
