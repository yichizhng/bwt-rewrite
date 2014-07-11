#ifndef _RDTSCLL_H
#define _RDTSCLL_H


  #define rdtscll(val)				\
  __asm__ __volatile__ ("rdtsc" : "=A" (val));

/*/

#define rdtscll(val) \
  __asm__ __volatile__ ("rdtsc\n\tshl $32, %%rdx\n\tor %%rdx, %%rax" :	\
			"=A" (val) : : "rdx");


#endif /* _RDTSCLL_H */
