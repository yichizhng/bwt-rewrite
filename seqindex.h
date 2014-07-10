#ifndef _SEQINDEX_H
#define _SEQINDEX_H

// The function to build the sequence index are here, as are the functions
// relating to the actual FM-index, as well as the struct definition thereof


long long **seq_index(unsigned char *, long long, long long, const unsigned char *);

long long seq_rank(unsigned char *, long long **, long long, long long, char, const unsigned char *);

unsigned char *lookup_table();

typedef struct _fmi {
	unsigned char *bwt;
	long long *idxs;
	long long **rank_index;
	unsigned char* lookup;
	long long endloc;
	long long C[5];
	long long len;
} fm_index;

// As the name suggests; deallocates all memory allocated for fmi, including
// fmi itself
void destroy_fmi(fm_index *fmi);

// Creates a FM-index from a given sequence using SACA-K
// (allocating memory dynamically)
fm_index *make_fmi(const unsigned char *str, unsigned long long len);

// Calculates the rank of a given symbol at a given index (i.e. the number
// of times the symbol has appeared up to that polong long) using the FM-index
// (Roughly constant time; this depends on implementation)
long long rank(const fm_index *fmi, unsigned char c, long long idx);

// Calculates the LF column mapping using the FM-index (constant time)
long long lf(const fm_index *fmi, long long idx);

// Searches for an exact pattern over the indexed sequence; returns
// the "first" (in this context, this means the rotation which appears
// first lexicographically) match, printing a warning message if there
// are multiple matches.
// pattern should be given uncompressed but in 0-3 form.
// Linear time in len * complexity of rank()
long long reverse_search(const fm_index *fmi, const unsigned char *pattern, long long len);

// Calculates SA[idx] from the FM-index
long long unc_sa(const fm_index *fmi, long long idx);

// Same as reverse_search
long long locate(const fm_index *fmi, const unsigned char *pattern, long long len);

// Same as locate, but returns the indices all matches (in the BWT; this means
// that you need to retrieve their indices using unc_sa) via sp and ep
void loc_search(const fm_index *fmi, const unsigned char *pattern, long long len, long long *sp, long long *ep);

// Finds the maximum mappable suffix of a given pattern; returns the number
// of bases matched, storing matches in sp and ep as per loc_search
long long mms(const fm_index *fmi, const unsigned char *pattern, long long len, long long *sp, long long *ep);

// Prlong longs part of a compressed sequence in human-readable form
void printseq(const unsigned char *seq, long long startidx, long long len);

#endif /* _SEQINDEX_H */
