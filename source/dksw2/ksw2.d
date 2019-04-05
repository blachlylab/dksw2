module dksw2.ksw2;

import core.stdc.stdint;
import core.stdc.stdio;
import core.stdc.string;
import core.stdc.stdlib;
import std.bitmanip:bitfields;

import dksw2.kalloc;

extern(C):

enum KSW_NEG_INF =-0x40000000;

enum KSW_EZ_SCORE_ONLY  =  0x01; // don't record alignment path/cigar
enum KSW_EZ_RIGHT       =  0x02; // right-align gaps
enum KSW_EZ_GENERIC_SC  =  0x04; // without this flag: match/mismatch only; last symbol is a wildcard
enum KSW_EZ_APPROX_MAX  =  0x08; // approximate max; this is faster with sse
enum KSW_EZ_APPROX_DROP =  0x10; // approximate Z-drop; faster with sse
enum KSW_EZ_EXTZ_ONLY   =  0x40; // only perform extension
enum KSW_EZ_REV_CIGAR   =  0x80; // reverse CIGAR in the output
enum KSW_EZ_SPLICE_FOR  =  0x100;
enum KSW_EZ_SPLICE_REV  =  0x200;
enum KSW_EZ_SPLICE_FLANK = 0x400;


struct ksw_extz_t{
    mixin(bitfields!(uint32_t,"max",31,uint32_t,"zdropped",1));
	// uint32_t max:31, zdropped:1;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be KSW_NEG_INF
	int m_cigar, n_cigar;
	int reach_end;
	uint32_t *cigar;
}

/**
 * NW-like extension
 *
 * @param km        memory pool, when used with kalloc
 * @param qlen      query length
 * @param query     query sequence with 0 <= query[i] < m
 * @param tlen      target length
 * @param target    target sequence with 0 <= target[i] < m
 * @param m         number of residue types
 * @param mat       m*m scoring mattrix in one-dimension array
 * @param gapo      gap open penalty; a gap of length l cost "-(gapo+l*gape)"
 * @param gape      gap extension penalty
 * @param w         band width (<0 to disable)
 * @param zdrop     off-diagonal drop-off to stop extension (positive; <0 to disable)
 * @param flag      flag (see KSW_EZ_* macros)
 * @param ez        (out) scores and cigar
 */
void ksw_extz(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const(int8_t) *mat,
			  int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);

void ksw_extz2_sse(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const(int8_t) *mat,
				   int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

void ksw_extd(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const(int8_t) *mat,
			  int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int flag, ksw_extz_t *ez);

void ksw_extd2_sse(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const int8_t *mat,
				   int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

void ksw_exts2_sse(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const int8_t *mat,
				   int8_t gapo, int8_t gape, int8_t gapo2, int8_t noncan, int zdrop, int flag, ksw_extz_t *ez);

void ksw_extf2_sse(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t mch, int8_t mis, int8_t e, int w, int xdrop, ksw_extz_t *ez);

/**
 * Global alignment
 *
 * (first 10 parameters identical to ksw_extz_sse())
 * @param m_cigar   (modified) max CIGAR length; feed 0 if cigar==0
 * @param n_cigar   (out) number of CIGAR elements
 * @param cigar     (out) BAM-encoded CIGAR; caller need to deallocate with kfree(km, )
 *
 * @return          score of the alignment
 */
int ksw_gg(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const(int8_t) *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
int ksw_gg2(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const(int8_t) *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
int ksw_gg2_sse(void *km, int qlen, const(uint8_t) *query, int tlen, const(uint8_t) *target, int8_t m, const(int8_t) *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

void *ksw_ll_qinit(void *km, int size, int qlen, const(uint8_t) *query, int m, const(int8_t) *mat);
int ksw_ll_i16(void *q, int tlen, const(uint8_t) *target, int gapo, int gape, int *qe, int *te);


pragma(inline,true)
uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = cast(uint32_t*)krealloc(km, cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

// In the backtrack matrix, value p[] has the following structure:
//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
pragma(inline,true)
void ksw_backtrack(void *km, int is_rot, int is_rev, int min_intron_len, const(uint8_t) *p, const(int) *off, const(int) *off_end, int n_col, int i0, int j0,
								 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{ // p[] - lower 3 bits: which type gets the max; bit
	int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
	uint32_t *cigar = *cigar_;
    uint32_t tmp;
	while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
		int force_state = -1;
		if (is_rot) {
			r = i + j;
			if (i < off[r]) force_state = 2;
			if (off_end && i > off_end[r]) force_state = 1;
			tmp = force_state < 0? p[cast(size_t)r * n_col + i - off[r]] : 0;
		} else {
			if (j < off[i]) force_state = 2;
			if (off_end && j > off_end[i]) force_state = 1;
			tmp = force_state < 0? p[cast(size_t)i * n_col + j - off[i]] : 0;
		}
		if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
		else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
		if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
		if (force_state >= 0) state = force_state;
		if (state == 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
		else if (state == 1 || (state == 3 && min_intron_len <= 0)) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
		else if (state == 3 && min_intron_len > 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 3, 1), --i; // intron
		else cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
	}
	if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len? 3 : 2, i + 1); // first deletion
	if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
	if (!is_rev)
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
	*m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

pragma(inline,true)
void ksw_reset_extz(ksw_extz_t *ez)
{
	ez.max_q = ez.max_t = ez.mqe_t = ez.mte_q = -1;
	ez.max = 0, ez.score = ez.mqe = ez.mte = KSW_NEG_INF;
	ez.n_cigar = 0, ez.zdropped = 0, ez.reach_end = 0;
}

pragma(inline,true)
int ksw_apply_zdrop(ksw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e)
{
	int r, t;
	if (is_rot) r = a, t = b;
	else r = a + b, t = a;
	if (H > cast(int32_t)ez.max) {
		ez.max = H, ez.max_t = t, ez.max_q = r - t;
	} else if (t >= ez.max_t && r - t >= ez.max_q) {
		int tl = t - ez.max_t, ql = (r - t) - ez.max_q, l;
		l = tl > ql? tl - ql : ql - tl;
		if (zdrop >= 0 && ez.max - H > zdrop + l * e) {
			ez.zdropped = 1;
			return 1;
		}
	}
	return 0;
}

ubyte[256] seq_nt4_table = [
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
];

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

//this is more for demonstration than actually being used
static void global_aln(const(char) *algo, void *km, const(char) *qseq_, const(char) *tseq_, int8_t m, const(int8_t) *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
					   int w, int zdrop, int flag, ksw_extz_t *ez)
{
	int i, qlen, tlen;
	uint8_t *qseq, tseq;
	ez.max_q = ez.max_t = ez.mqe_t = ez.mte_q = -1;
	ez.max = 0, ez.mqe = ez.mte = KSW_NEG_INF;
	ez.n_cigar = 0;
	qlen = cast(int)strlen(qseq_);
	tlen = cast(int)strlen(tseq_);
	qseq = cast(uint8_t*)calloc(qlen + 33, 1); // 32 for gaba
	tseq = cast(uint8_t*)calloc(tlen + 33, 1);
	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[cast(uint8_t)qseq_[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[cast(uint8_t)tseq_[i]];
	if (strcmp(algo, "gg") == 0) {
		if (flag & KSW_EZ_SCORE_ONLY) ez.score = ksw_gg(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, null, null, null);
		else ez.score = ksw_gg(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
	} else if (strcmp(algo, "gg2") == 0) {
		if (flag & KSW_EZ_SCORE_ONLY) ez.score = ksw_gg2(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, null, null, null);
		else ez.score = ksw_gg2(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
	}
	else if (strcmp(algo, "gg2_sse") == 0)     ez.score = ksw_gg2_sse(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
	else if (strcmp(algo, "extz") == 0)        ksw_extz(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, zdrop, flag, ez);
	else if (strcmp(algo, "extz2_sse") == 0)   ksw_extz2_sse(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, w, zdrop, 0, flag, ez);
	else if (strcmp(algo, "extd") == 0)        ksw_extd(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	else if (strcmp(algo, "extd2_sse") == 0)   ksw_extd2_sse(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, 0, flag, ez);
	else if (strcmp(algo, "extf2_sse") == 0)   ksw_extf2_sse(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, mat[0], mat[1], e, w, zdrop, ez);
	else if (strcmp(algo, "exts2_sse") == 0) {
        // I don't understand this code's purpose
		// int8_t[25] mat;
		// ksw_gen_simple_mat(5, mat, 1, 2);
		ksw_exts2_sse(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, 2, 1, 32, 4, zdrop, flag|KSW_EZ_SPLICE_FOR, ez);
	}
	else if (strcmp(algo, "test") == 0) ksw_extd2_sse(km, qlen, cast(uint8_t*)qseq, tlen, cast(uint8_t*)tseq, m, mat, 4, 2, 24, 1, 751, 400, 0, 8, ez);
	else {
		fprintf(stderr, "ERROR: can't find algorithm '%s'\n", algo);
		exit(1);
	}
	free(qseq); free(tseq);
}

void print_aln(const(char) *tname, const(char) *qname, ksw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez.score);
	printf("\t%d\t%d\t%d", ez.max, ez.max_t, ez.max_q);
	if (ez.n_cigar > 0) {
		int i;
		putchar('\t');
		for (i = 0; i < ez.n_cigar; ++i)
			printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
	}
	putchar('\n');
}

