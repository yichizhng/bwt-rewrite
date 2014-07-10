// Slowish implementation of Smith-Waterman algorithm for local alignments
// There are optimizations involving SIMD instructions but they don't really
// seem that worthwhile
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "rdtscll.h"
#include "stack.h"

static inline int max(int a, int b, int c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

// Optimization of needleman-wunsch by using a 1d output array; is faster by a
// pretty hilarious factor (>40x) due to the lack of another level of
// indirection and/or cache optimizations and/or lack of malloc() calls

// str1 should be the read (allowed characters are 0-3 and 5), str2 the
// genome (allowed characters 0-3), both unpacked form. 'N' on the read will
// be treated as if it matches all characters.

// Returns the position on str2 that the last character of str1 was aligned
// to. Outputs some CIGARs to the given stack. This function can be used
// to align both the head and tail; the head should be passed in backwards
// (i.e. str1[0] and str2[0] are the characters that the MMS failed on)
int nw_fast(const unsigned char *str1, int len1, const unsigned char *str2, int len2, stack *s, int *score, int *indels) {
  //  fprintf(stderr, "%d %d\n", len1, len2);
  if (len1 == 0) { // happens more often than you'd think
    return 0; // Nothing at all to do
  }
  int *values, i, j;
  int mx = -1000000, maxloc = 0;
  values = malloc((len1 + 1) * (len2 + 1) * sizeof(int));
  char *pointers = malloc((len1 + 1) * (len2 + 1));
  stack *flips = stack_make();
  // "Zero" first row
  values[0] = 0;
  for (j = 1; j <= len2; ++j) {
    pointers[j] = 2;
    values[j] = -5-3*j;
  }
  for (i = 1; i <= len1; ++i) {
    // Zero first column
    pointers[i * (len2 + 1)] = 1;
    values[i * (len2 + 1)] = -5-3*i;
    for (j = 1; j <= len2; ++j) {
      int skip1, skip2;
      skip2 = (pointers[i * (len2 + 1) + j - 1] == 2) ? 0 : -5;
      skip1 = (pointers[(i - 1) * (len2 + 1) + j] == 1) ? 0 : -5;
      // Update cell appropriately
      values[i * (len2 + 1) + j] =
	max(values[(i-1) * (len2 + 1) + j - 1] + (((str1[i-1]==5)||(str1[i-1] == str2[j-1]))?0:-6),
	    values[i * (len2 + 1) + j - 1] - 3 + skip2,
	    values[(i-1) * (len2 + 1) + j] - 3 + skip1);
      if (values[i * (len2 + 1) + j] == values[i * (len2 + 1) + j - 1] - 3 + skip2) {
	pointers[i * (len2 + 1) + j] = 2; // skip on 2
      }
      else if (values[i * (len2 + 1) + j] == values[(i-1) * (len2 + 1) + j] - 3 + skip1) {
	pointers[i * (len2 + 1) + j] = 1; // skip on 1
      }
      else
	pointers[i * (len2 + 1) + j] = 0; // no skip
    }
  }
  for (j = 0; j <= len2; ++j) {
    if (values[len1 * (len2 + 1) + j] > mx) {
      mx = values[len1 * (len2 + 1) + j];
      maxloc = j;
    }
  }
  i = len1;
  j = maxloc;
  *score += values[i * (len2 + 1) + j];
  // In theory something bad might have happened (some insertions on the end)
  // if str1 wasn't long enough. I'm going to ignore that possibility :)
  // But it's something to think about implementing
  while (i && j) {
    // Figure out which way the correct path goes; we can infer this directly
    // from the score, in fact, but it's just easier to keep track of that
    // from a separate pointer array ("pointer"; in this case it's just
    // storing which way we're going)
    
    char dir = pointers[i * (len2 + 1) + j];
    switch(dir) {
    case 1:
      i--;
      stack_push(flips, 'I', 1);
      (*indels)--;
      break;
    case 2:
      j--;
      stack_push(flips, 'D', 1);
      (*indels)--;
      break;
    default: // 0
      i--;
      j--;
      stack_push(flips, 'M', 1);
      break;
    }
  }
  while (i) {
    i--;
    stack_push(flips, 'I', 1);
    (*indels)--;
  }
  while (j) {
    j--;
    stack_push(flips, 'D', 1);
    (*indels)--;
  }
  stack_flip (flips, s);
  free(pointers);
  free(values);
  return maxloc - 1;
}

// The same thing, except this time we always backtrack from the end of both
// strings (so obviously we don't need to return the position on str2 that
// we aligned to), outputting to the stack
// str1 still refers to the pattern and str2 to the genome, for consistency
void sw_fast(const unsigned char *str1, int len1, const unsigned char *str2, int len2, stack *s, int *score, int *indels) {
  int *values, i, j;
  char *pointers;
  values = malloc((len1 + 1) * (len2 + 1) * sizeof(int));
  pointers = malloc((len1 + 1) * (len2 + 1));
  // "Zero" first row
  values[0] = 0;
  for (j = 1; j <= len2; ++j) {
    values[j] = -5-3*j;
    pointers[j] = 2;
  }
  for (i = 1; i <= len1; ++i) {
    // Zero first column
    values[i * (len2 + 1)] = -5-3*i;
    pointers[i * (len2 + 1)] = 1;
    for (j = 1; j <= len2; ++j) {
      int skip1, skip2;
      skip2 = (pointers[i * (len2 + 1) + j - 1] == 2) ? 0 : -5;
      skip1 = (pointers[(i - 1) * (len2 + 1) + j] == 1) ? 0 : -5;
      // Update cell appropriately
      values[i * (len2 + 1) + j] =
	max(values[(i-1) * (len2 + 1) + j - 1] + (((str1[i-1]==5)||(str1[i-1] == str2[j-1]))?0:-6),
	    values[i * (len2 + 1) + j - 1] - 3 + skip2,
	    values[(i-1) * (len2 + 1) + j] - 3 + skip1);
      if (values[i * (len2 + 1) + j] == values[i * (len2 + 1) + j - 1] - 3 + skip2)
	pointers[i * (len2 + 1) + j] = 2; // skip on 2
      else if (values[i * (len2 + 1) + j] == values[(i-1) * (len2 + 1) + j] - 3 + skip1)
	pointers[i * (len2 + 1) + j] = 1; // skip on 1
      else
	pointers[i * (len2 + 1) + j] = 0; // no skip
    }
  }
  // Now we need to iterate back from the end of the match (len1, len2)
  // and figure out the best (or A best) path. Because these segments are
  // passed in forward order, we do not need to allocate another stack to
  // do a flip (we are pushing them to the stack in back to front order, which
  // is correct)
  i = len1;
  j = len2;
  *score += values[i * (len2 + 1) + j];
  while (i && j) {
    // Figure out which way the correct path goes; we can infer this directly
    // from the score, in fact, but it's just easier to keep track of that
    // from a separate pointer array ("pointer"; in this case it's just
    // storing which way we're going)
    
    // Since, as before, str1 refers to the pattern and str2 the genome,
    // a skip on 1 is to be called a deletion (the base is present on the genome
    // but not the read) and a skip on 2 an insertion. We do not output
    // characters other than M, I, and D. Score can be calculated directly
    // and trivially from the CIGAR if necessary.
    char dir = pointers[i * (len2 + 1) + j];
    switch(dir) {
    case 1:
      i--;
      stack_push(s, 'I', 1);
      (*indels)--;
      break;
    case 2:
      j--;
      stack_push(s, 'D', 1);
      (*indels)--;
      break;
    default: // 0
      i--;
      j--;
      stack_push(s, 'M', 1);
      break;
    }
  }
  while (i) {
    i--;
    stack_push(s, 'I', 1);
    (*indels)--;
  }
  while (j) {
    j--;
    stack_push(s, 'D', 1);
    (*indels)--;
  }
  // verify contents of stack
  //for (int k = 0; k < s->size; ++k) {
  //  fprintf(stderr, "%c\n", s->chars[k]);
  //}

  free(values);
  free(pointers);
  return;
}

// Note that this implementation takes the full O(m*n) memory; it is possible
// to do with less (especially if we only want the optimal
// alignment), but much more annoying
// The slow and obnoxious implementation, returns a 2D array with results
int** smw(const char *str1, int len1, const char *str2, int len2) {
  // Allocate a 2-D array for values
  int **values, i, j;
  values = malloc((len1+1) * sizeof(void *));
  for (i = 0; i <= len1; ++i) {
    values[i] = malloc((len2+1) * sizeof(int));
    values[i][0] = -i;
  }
  for (j = 1; j <= len2; ++j) {
    values[0][j] = -j;
  }
  // The first row and column of the matrix are 0's, and have already
  // been assigned in the previous loop
  
  // For now we use +2 for a match and -1 for indel or mismatch
  
  for (i = 1; i <= len1; ++i) {
    for (j = 1; j <= len2; ++j) {
      if (str1[i-1] == str2[j-1]) {
	values[i][j] = 2 + values[i-1][j-1];
	continue;
	// Taking a match is always better than taking an indel if
	// we use linear gap penalty (instead of affine)
      }
      // Find the maximum of the three possibilities and subtract 1
      values[i][j] = max(values[i][j-1], values[i-1][j], values[i-1][j-1]) - 1;
      //if (values[i][j] < 0) // uncomment to use smith-waterman
      //values[i][j] = 0; // keep commented to use needleman-wunsch
      // smith-waterman is harder to backtrack, if you're wondering
      // what the problem is, so it's better to just pass appropriate length
      // strings in to get local alignments that way
    }
  }
  // In order to retrieve the correct path we need to iterate backward from
  // (i, j), but we can leave that for a different function
  return values;
}
