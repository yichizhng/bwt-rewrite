// Tries aligning reads from a file against an index and sequence read from
// file, assuming that they are not spliced reads
// This, of course, requires that we put another function together.

// usage: single_align seqfile indexfile readfile

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include "seqindex.h"
#include "csacak.h"
#include "fileio.h"
#include "rdtscll.h"
#include "time.h"
#include "smw.h"
#include "stack.h"

static inline unsigned char getbase(const unsigned char *str, int idx) {
  if (idx<0) idx=0;
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// Continues a MMS search
int mms_continue(const fm_index *fmi, const unsigned char *pattern, int len, int *sp, int *ep) {
  int start, end, i;
  start = *sp;
  end = *ep;
  for (i = len-1; i >= 0; --i) {
    if (end <= start) {
      break;
    }
    *sp = start;
    *ep = end;
    start = fmi->C[(ptrdiff_t)pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[(ptrdiff_t)pattern[i]] + rank(fmi, pattern[i], end);
  }
  if (end <= start) // Didn't finish matching
    return len - i - 2;
  else { // Finished matching
    *sp = start;
    *ep = end;
    return len - i - 1;
  }
}

// Tries continuing a mms search with mismatch; returns upon finding any continuation with at least 6 matching nts
// Last argument is the difference between the return value and the number of nts on the genome matched (from -3 to 3).
int mms_mismatch(const fm_index *fmi, const unsigned char *seq, const unsigned char *pattern, int len, long long int *sp, long long int *ep, int *genomeskips) {
  // If there are too many matches, don't even bother
  //  if (*ep - *sp > 10)
  //    return -1;
  if (len < 2) { // nothing to do, really
    int loc = unc_sa(fmi, *sp);
    unsigned char sub_c = getbase(seq, loc-1);
    *sp = fmi->C[(ptrdiff_t)sub_c] + rank(fmi, sub_c, *sp);
    *ep = *sp + 1;
    *genomeskips = 0;
    return 1;
  }
  int best_align = 0;
  int best_pos = -1;
  for (long long int i = *sp; i < *ep; ++i) {
    // Reads the start and end from sp and ep instead of using the last
    // character of the sequence. It assumes that we have a mismatch at that
    // point (mms returns if that happens or it finished)
    // and tries the following things to try aligning it
    
    // 1) Assume that there was a substitution at that point. Use LF() to skip
    // to the next nt and decrement len, then try aligning
    {
      int loc = unc_sa(fmi, i);
      char sub_c = getbase(seq, loc-1);
      int sub_idx = fmi->C[(ptrdiff_t)sub_c] + rank(fmi, sub_c, i), ins_idx = sub_idx;
      int sub_end = sub_idx + 1, sub_align;
      sub_align = mms_continue(fmi, pattern, len-1, &sub_idx, &sub_end) + 1;
      best_align = sub_align;
      best_pos = sub_idx;
      if (sub_align > 6 || sub_align == len) {
	*genomeskips = 0;
	break;
      }

      // 1.5) Assume that there was an insertion (on the genome) at that point of up to three nts
      // Use LF() to skip one, two, and three nts and _don't_ decrement len, then try aligning for each of those
      int bleh = ins_idx;

      int ins_end = ins_idx + 1, ins_align;
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*genomeskips = 1;
	break;
      }

      // two!
      sub_c = getbase(seq, loc-2);
      ins_idx = fmi->C[(ptrdiff_t)sub_c] + rank(fmi, sub_c, bleh);
      int blah = ins_idx;
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*genomeskips = 2;
	break;
      }

      // three!
      sub_c = getbase(seq, loc-3);
      ins_idx = fmi->C[(ptrdiff_t)sub_c] + rank(fmi, sub_c, blah);
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*genomeskips = 3;
	break;
      }
    }
    
    // 2) Assume that there was a deletion (on the genome) at that point.
    // Ignore up to three nts and start aligning again
    {
      // This one is a lot simpler because we don't actually need to
      // figure out the character
      int del_idx = i, del_end = del_idx + 1, del_align;
      del_align = mms_continue(fmi, pattern, len-1, &del_idx, &del_end) + 1;
      if (del_align > 6 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*genomeskips = -1;
	break;
      }

      del_idx = i;
      del_end = del_idx + 1;
      del_align = mms_continue(fmi, pattern, len-2, &del_idx, &del_end) + 2;
      if (del_align > 7 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*genomeskips = -2;
	break;
      }

      del_idx = i;
      del_end = del_idx + 1;
      del_align = mms_continue(fmi, pattern, len-3, &del_idx, &del_end) + 3;
      if (del_align > 8 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*genomeskips = -3;
	break;
      }
    }
  }
  *sp = best_pos;
  *ep = best_pos + 1;
  return best_align;
}

// Pass in the required anchor length. No mismatch will be allowed.
unsigned long long align_read_anchored(const fm_index *fmi, const unsigned char *seq, const unsigned char *pattern, int len, int anchor_len, stack *s) {
  int score;
  int indels;
  const int olen = len;
  long long curgap = 0;
  long long curpos = -1;
  long long endpos;
  int anchlen;
  while ((len > anchor_len) && ((olen - len) < 2 * anchor_len)) {
    score = -1;
    while ((len > anchor_len)  && ((olen - len) < 2 * anchor_len)) {
      int seglen = mms(fmi, pattern, len, &curpos, &endpos);
      if (seglen < anchor_len || endpos - curpos > 1) {
	len -= 3;
	continue;
      }
      else {
	len -= seglen;
	anchlen = seglen;
	score = (int) (0.6 * (1 + olen));
	indels = 5; // Should be adjustable?
	curpos = unc_sa(fmi, curpos);
	//fprintf(stderr, "%d %d %d\n", anchlen, olen, len);

	// And use N-W to align the "tail" of the read
	int buflen = indels + (olen - (len + seglen));
	if (buflen + curpos + seglen > fmi->len)
	  buflen = fmi->len - curpos - seglen;
	unsigned char *buf = malloc(buflen);
	for (int i = 0; i < buflen; ++i)
	  buf[i] = getbase(seq, curpos + seglen + i);
	nw_fast(pattern + len + seglen, olen - (len + seglen),
		buf, buflen, s, &score, &indels);
	// TODO: it is probably better to instead reverse the tail of
	// the buffer and align it from there

	// We can ignore the return value (we don't really care where the
	// end of the read ends up; we can calculate that from the CIGAR)
	free(buf);
	// Then push this anchor onto it
	stack_push(s, 'M', seglen);
	break;
      }
    }
    
    if (score < 0)
      continue;

    // In the second loop we try to extend our anchor backwards
    while ((len > anchor_len) && (score >= 0) && (indels >= 0)) {
      for (curgap = 1; curgap + anchor_len < len && curgap < anchor_len && curgap < len; ++curgap) {
	long long int start, end;
	int seglen = mms(fmi, pattern, len-curgap, &start, &end);
	if (seglen < anchor_len)
	  continue;
	int matched = 0;
	for (long long i = start; i < end; ++i) {
	  long long cpos = unc_sa(fmi, i);
	  if (abs(cpos + seglen - curpos) - curgap <= 3) {
	    matched = 1;
	    // Align the stuff in between. In this case we don't need to
	    // copy pattern to a new buffer, but we do still need to copy
	    // the genome
	    int buflen = curpos - (cpos + seglen);
	    // There's a semi-theoretical problem that this might actually
	    // be negative, but that's easy to resolve
	    if (buflen < 0) {
	      stack_push(s, 'I', -buflen);
	      indels += buflen;
	    }
	    else {
	      unsigned char *buf = malloc(buflen);
	      for (int j = 0; j < buflen; ++j)
		buf[j] = getbase(seq, cpos + seglen + j);
	      // And compare
	      sw_fast(pattern + (len - curgap), curgap, buf, buflen, s, &score, &indels);
	      free(buf);
	    }
	    stack_push(s, 'M', seglen);
	    curpos = cpos;
	    len -= seglen + curgap;
	    curgap = 0;
	    break;
	  }
	}
	if (matched)
	  break;
	else
	  continue;
      }
      if (curgap) {
	break;
      }
    }
    if ((score >= 0) && (indels >= 0)) {
      if (len < 0) {
	// I don't even know when this happens
	return 0;
      }
      // Set up matrix for N-W alignment
      int buflen = len + indels;
      if (buflen > curpos)
	buflen = curpos;
      unsigned char *buf = malloc(buflen);
      for (int i = 0; i < buflen; ++i)
	buf[i] = getbase(seq, curpos - 1 - i);
      unsigned char *buf2 = malloc(len);
      for (int i = 0; i < len; ++i)
	buf2[i] = pattern[len-1-i];
      int x = nw_fast(buf2, len, buf, buflen, s, &score, &indels);
      free(buf);
      free(buf2);
      //printf("%d %d\t", x, len);
      printf("\t%d\n", indels);
      if ((score >= 0) && (indels >= 0))
	return curpos - x;
      //return 0; // Give up early to save some time
    }
    // reset the stack
    len-= anchlen;
    s->size = 0;
  }
  if ((score < 0) || (indels < 0)) {
    return 0;
  }

  int buflen = len + indels;
  if (buflen > curpos)
    buflen = curpos;
  unsigned char *buf = malloc(buflen);
  for (int i = 0; i < buflen; ++i)
    buf[i] = getbase(seq, curpos - 1 - i);
  unsigned char *buf2 = malloc(len);
  for (int i = 0; i < len; ++i)
    buf2[i] = pattern[len-1-i];
  int x = nw_fast(buf2, len, buf, buflen, s, &score, &indels);
  free(buf);
  free(buf2);
  if ((score >= 0) && (indels >= 0))
    return curpos - x;
  return 0;
}

int align_read(const fm_index *fmi, const unsigned char *seq, const unsigned char *pattern, int len, int thresh) {
  int starts[10], lens[10], nsegments;
  int penalty;
  int nmisses = len/10;
  int olen = len;
  for (nsegments = 0; nsegments < 10; nsegments++) {
    if (len < 10)
      break;
    long long int start, end;
    int seglen = mms(fmi, pattern, len, &start, &end);
    if (seglen < thresh) {
      int mlen = mms_mismatch(fmi, seq, pattern, len - seglen, &start, &end, &penalty);
      if (mlen + seglen > 2 * thresh) {
	len -= seglen + mlen + 3;
	starts[nsegments] = start;
	lens[nsegments] = seglen + mlen;
	continue;
      }
      if (!nmisses--)
	return 0;
      len -= 3;
      nsegments--;
      if (nsegments > -1) {
	starts[nsegments] -= 3;
	lens[nsegments] += 3;
      }
      continue;
    }
    if ((len - seglen == 0) || ((len - seglen > 10) && end - start == 1)) {
      starts[nsegments] = start;
      lens[nsegments] = seglen;
      len -= seglen + 3;
      continue;
    }
    // Otherwise try continuing the search
    int mlen = mms_mismatch(fmi, seq, pattern, len - seglen, &start, &end, &penalty);
    len -= seglen + mlen + 3;
    starts[nsegments] = start;
    lens[nsegments] = seglen + mlen;
  }
  int totlen = lens[0];
  if (nsegments == 10)
    return 0; // Too many segments
  else {
    // For each segment check whether it's within 6 nts of the next
    
    for (int i = 0; i < nsegments - 1; ++i) {
      if (abs(unc_sa(fmi, starts[i+1]) + lens[i+1] - unc_sa(fmi, starts[i])) < 7) {
	totlen += lens[i+1];
	continue;
      } 
      else
	return 0; // Gapped
    }
  }
  if (3 * totlen > 2 * olen)
    return unc_sa(fmi, starts[nsegments-1]) - len;
  return 0;
}

// Reminder to self: buf length (i.e. maximum read length) is currently
// hardcoded; change to a larger value (to align longer reads) or make it
// dynamic

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s seqfile indexfile readfile\n", argv[0]);
    exit(-1);
  }
  char *seqfile, *indexfile, *readfile;
  unsigned char *seq, *buf = malloc(256 * 256), *revbuf = malloc(256 * 256), c;
  fm_index *fmi;
  int len;
  int i;
  FILE *sfp, *ifp, *rfp;
  seqfile = argv[1];
  indexfile = argv[2];
  readfile = argv[3];
  sfp = fopen(seqfile, "rb");
  if (sfp == 0) {
    fprintf(stderr, "Could not open sequence\n");
    exit(-1);
  }
  fseek(sfp, 0L, SEEK_END);
  len = ftell(sfp);
  rewind(sfp);
  seq = malloc(len/4+1);
  for (i = 0; i < len/4 + 1; ++i) {
    switch(fgetc(sfp)) {
    case 'C': c = 64; break;
    case 'G': c = 128; break;
    case 'T': c = 192; break;
    default: c = 0;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 16; break;
    case 'G': c ^= 32; break;
    case 'T': c ^= 48;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 4; break;
    case 'G': c ^= 8; break;
    case 'T': c ^= 12;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 1; break;
    case 'G': c ^= 2; break;
    case 'T': c ^= 3;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < (len&3); ++i) {
    switch(fgetc(sfp)) {
    case 'C': c ^= 64 >> (2 * i); break;
    case 'G': c ^= 128 >> (2 * i); break;
    case 'T': c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  fclose(sfp);
  
  // Open index file
  ifp = fopen(indexfile, "rb");
  if (ifp == 0) {
    fprintf(stderr, "Could not open index file");
    exit(-1);
  }
  fmi = read_index(ifp);
  fclose(ifp);

  // And now we go read the index file
  rfp = fopen(readfile, "r");
  if (rfp == 0) {
    fprintf(stderr, "Could not open reads file");
    exit(-1);
  }
  // Read one line ("read") and try aligning it
  
  int naligned = 0;
  int nread = 0;
  while (!feof(rfp)) {
    if (!fgets((char *)buf, 256*256-1, rfp))
      break;
    nread++;
    if (buf[strlen((char *)buf)-1] == '\n')
      buf[strlen((char *)buf)-1] = 0;
    int len = strlen((char *)buf);
    for (int i = 0; i < len; ++i) {
      // Replace with "compressed" characters
      switch(buf[i]) {
      case 'A':
	buf[i] = 0;
	revbuf[len-i-1] = 3;
	break;
      case 'C':
	buf[i] = 1;
	revbuf[len-i-1] = 2;
	break;
      case 'T':
	buf[i] = 3;
	revbuf[len-i-1] = 0;
	break;
      case 'G':
	buf[i] = 2;
	revbuf[len-i-1] = 1;
	break;
      default: // 'N'
	buf[i] = 5;
	revbuf[len-i-1] = 5;
	break;
      }
    }

    stack *s = stack_make();
    //    int thresh = (int) (-1.2 * (1+len));

    //    int pos = align_read(fmi, seq, buf, len, 10);
    int pos = align_read_anchored(fmi, seq, buf, len, 12, s);
    if (pos) {
      naligned++;
      printf("%d\n", pos + 1);
      stack_print_destroy(s);
    }
    else {
      stack_destroy(s);
      s = stack_make();
      //      pos = align_read(fmi, seq, revbuf, len, 10);
      pos = align_read_anchored(fmi, seq, revbuf, len, 12, s);
      if (pos) {
	naligned++;
	printf("%d\n", pos + 1);
	stack_print_destroy(s);
      }
      else {
	printf("0\n");
	stack_destroy(s);
      }
    }
  }
  fclose(rfp);
  fprintf(stderr, "%d of %d reads aligned\n", naligned, nread);
  
  free(buf);
  free(revbuf);
  destroy_fmi(fmi);
  free(seq);
  return 0;
}
