#ifndef _SMW_H
#define _SMW_H

#include "stack.h"

int **smw(const char*, int, const char*, int);

int nw_fast(const unsigned char *str1, int len1, const unsigned char *str2, int len2, stack *s, int *score, int *indels);

void sw_fast(const unsigned char *str1, int len1, const unsigned char *str2, int len2, stack *s, int *score, int *indels);

#endif /* _SMW_H */
