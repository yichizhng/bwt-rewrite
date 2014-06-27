// Author: Ge Nong,  Email: issng@mail.sysu.edu.cn
// Department of Computer Science, Sun Yat-sen University, 
// Guangzhou, China
// Date: December 24, 2012
//
// This is the demo source code for the algorithm SACA-K presented in this article:
// G. Nong, Practical Linear-Time O(1)-Workspace Suffix Sorting for Constant Alphabets, 
// ACM Transactions on Information Systems, Scheduled to Appear in July 2013.
// A draft for this article can be retrieved from http://code.google.com/p/ge-nong/.

#include <stdlib.h>

// set only the highest bit as 1, i.e. 1000...
const long long EMPTY = 0x8000000000000000;

static __inline__ unsigned char getbase(unsigned char *str, unsigned long long idx) {
  return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// Minor adaptations to make this compile in C (which has no bool type)
// Compile with -std=gnu99, or fix the for loop initial declaration yourself.
// This version is adapted to work with the compressed representation of
// nucleotide sequences (4nts / byte).
typedef char bool;
#define true ((char) 1)
#define false ((char) 0)
// -- Yichi Zhang

// get getbase(s,i) at a certain level
#define chr(i) ((level==0)?(getbase((unsigned char *)s,i)):((long long *)s)[i])

void getBuckets(unsigned char *s, 
		unsigned long long *bkt, unsigned long long n,
		unsigned long long K, bool end) { 
  unsigned long long i, sum=0;
  
  // clear all buckets .
  for(i=0; i<K; i++) bkt[i]=0; 
  
  // compute the size of each bucket .
  for(i=0; i<n; i++) bkt[getbase(s,i)]++; 
  
  for(i=0; i<K; i++) { 
    sum+=bkt[i]; 
    bkt[i]=end ? sum-1 : sum-bkt[i]; 
  }
}

void putSuffix0(unsigned long long *SA, 
		unsigned char *s, unsigned long long *bkt, 
		unsigned long long n, unsigned long long K, long long n1) {
  unsigned long long i, j;
  
  // find the end of each bucket.
  getBuckets(s, bkt, n, K, true);

  // put the suffixes long longo their buckets.
  for(i=n1-1; i>0; i--) {
    j=SA[i]; SA[i]=0;
    SA[bkt[getbase(s,j)]--]=j;
  }
  SA[0]=n-1; // set the single sentinel suffix.
}

void induceSAl0(unsigned long long *SA,
		unsigned char *s, unsigned long long *bkt,
		unsigned long long n, unsigned long long K, bool suffix) {
  unsigned long long i, j;

  // find the head of each bucket.
  getBuckets(s, bkt, n, K, false);

  bkt[0]++; // skip the virtual sentinel.
  for(i=0; i<n; i++)
    if(SA[i]>0) {
      j=SA[i]-1;
      if(getbase(s,j)>=getbase(s,j+1)) {
        SA[bkt[getbase(s,j)]]=j;
        bkt[getbase(s,j)]++;
        if(!suffix && i>0) SA[i]=0;
      }
    }
}

void induceSAs0(unsigned long long *SA,
		unsigned char *s, unsigned long long *bkt,
		unsigned long long n, unsigned long long K, bool suffix) {
  unsigned long long i, j;

  // find the end of each bucket.
  getBuckets(s, bkt, n, K, true);

  for(i=n-1; i>0; i--)
    if(SA[i]>0) {
      j=SA[i]-1;
      if(getbase(s,j)<=getbase(s,j+1) && bkt[getbase(s,j)]<i) {
        SA[bkt[getbase(s,j)]]=j;
        bkt[getbase(s,j)]--;
        if(!suffix) SA[i]=0;
      }
    }
}

void putSubstr0(unsigned long long *SA,
		unsigned char *s, unsigned long long *bkt,
		unsigned long long n, unsigned long long K) {
  unsigned long long i, cur_t, succ_t;

  // find the end of each bucket.
  getBuckets(s, bkt, n, K, true);

  // set each item in SA as empty.
  for(i=0; i<n; i++) SA[i]=0;

  succ_t=0; // getbase(s,n-2) must be L-type.
  for(i=n-2; i>0; i--) {
    cur_t=(getbase(s,i-1)<getbase(s,i) ||
           (getbase(s,i-1)==getbase(s,i) && succ_t==1)
          )?1:0;
    if(cur_t==0 && succ_t==1) SA[bkt[getbase(s,i)]--]=i;
    succ_t=cur_t;
  }

  // set the single sentinel LMS-substring.
  SA[0]=n-1;
}

void putSuffix1(long long *SA, long long *s, long long n1) {
  long long i, j, pos, cur, pre=-1;
  
  for(i=n1-1; i>0; i--) {
    j=SA[i]; SA[i]=EMPTY;
    cur=s[j];
    if(cur!=pre) {
      pre=cur; pos=cur;
    }
    SA[pos--]=j;
  }
}

void induceSAl1(long long *SA, long long *s, 
		long long n, bool suffix) {
  long long h, i, j, step=1;
  
  for(i=0; i<n; i+=step) {
    step=1; j=SA[i]-1;
    if(SA[i]<=0) continue;
    long long c=s[j], c1=s[j+1];
    bool isL=c>=c1;
    if(!isL) continue;

    // getbase(s,j) is L-type.

    long long d=SA[c];
    if(d>=0) {
      // SA[c] is borrowed by the left
      //   neighbor bucket.
      // shift-left the items in the
      //   left neighbor bucket.
      long long foo, bar;
      foo=SA[c];
      for(h=c-1; SA[h]>=0||SA[h]==EMPTY; h--)
      { bar=SA[h]; SA[h]=foo; foo=bar; }
      SA[h]=foo;
      if(h<i) step=0;

      d=EMPTY;
    }

    if(d==EMPTY) { // SA[c] is empty.
      if(c<n-1 && SA[c+1]==EMPTY) {
        SA[c]=-1; // init the counter.
        SA[c+1]=j;
      }
      else        
        SA[c]=j; // a size-1 bucket.
    }
    else { // SA[c] is reused as a counter.
        long long pos=c-d+1;
        if(pos>n-1 || SA[pos]!=EMPTY) {
          // we are running long longo the right
          //   neighbor bucket.
          // shift-left one step the items
          //   of bucket(SA, S, j).
          for(h=0; h<-d; h++)
            SA[c+h]=SA[c+h+1];
          pos--;
          if(c<i) step=0;
        }
        else
          SA[c]--;

        SA[pos]=j;
    }

    long long c2;
    bool isL1=(j+1<n-1) && (c1>(c2=s[j+2]) || (c1==c2 && c1<i));  // is s[SA[i]] L-type?
    if((!suffix || !isL1) && i>0) {
      long long i1=(step==0)?i-1:i;
      SA[i1]=EMPTY;
    }
  }

  // scan to shift-left the items in each bucket 
  //   with its head being reused as a counter.
  for(i=1; i<n; i++) {
    j=SA[i];
    if(j<0 && j!=EMPTY) { // is SA[i] a counter?
      for(h=0; h<-j; h++)
        SA[i+h]=SA[i+h+1];
      SA[i+h]=EMPTY;
    }
  }
}

void induceSAs1(long long *SA, long long *s, 
		long long n, bool suffix) {
  long long h, i, j, step=1;
  
  for(i=n-1; i>0; i-=step) {
    step=1; j=SA[i]-1;
    if(SA[i]<=0) continue;
    long long c=s[j], c1=s[j+1];
    bool isS=(c<c1) || (c==c1 && c>i);
    if(!isS) continue;

    // getbase(s,j) is S-type

    long long d=SA[c];
    if(d>=0) {
      // SA[c] is borrowed by the right
      //   neighbor bucket.
      // shift-right the items in the
      //   right neighbor bucket.
      long long foo, bar;
      foo=SA[c];
      for(h=c+1; SA[h]>=0||SA[h]==EMPTY; h++)
      { bar=SA[h]; SA[h]=foo; foo=bar; }
      SA[h]=foo;
      if(h>i) step=0;

      d=EMPTY;
    }

    if(d==EMPTY) { // SA[c] is empty.
      if(SA[c-1]==EMPTY) {
        SA[c]=-1; // init the counter.
        SA[c-1]=j;
      }
      else
        SA[c]=j; // a size-1 bucket.
    }
    else { // SA[c] is reused as a counter.
        long long pos=c+d-1;
        if(SA[pos]!=EMPTY) {
          // we are running long longo the left
          //   neighbor bucket.
          // shift-right one step the items 
          //   of bucket(SA, S, j).
          for(h=0; h<-d; h++)
            SA[c-h]=SA[c-h-1];
          pos++;
          if(c>i) step=0;
        }
        else
          SA[c]--;

        SA[pos]=j;
    }

    if(!suffix) {
      long long i1=(step==0)?i+1:i;
      SA[i1]=EMPTY;
    }
  }

  // scan to shift-right the items in each bucket
  //   with its head being reused as a counter.
  if(!suffix)
    for(i=n-1; i>0; i--) {
      j=SA[i];
      if(j<0 && j!=EMPTY) { // is SA[i] a counter?
        for(h=0; h<-j; h++)
          SA[i-h]=SA[i-h-1];
        SA[i-h]=EMPTY;
      }
    }
}

void putSubstr1(long long *SA, long long *s, long long n) {
  long long h, i, j;

  for(i=0; i<n; i++) SA[i]=EMPTY;

  long long c, c1, t, t1;
  c1=s[n-2];
  t1=0; 
  for(i=n-2; i>0; i--) {
    c=c1; t=t1; 
    c1=s[i-1];
    t1=c1<c || (c1==c && t);
    if(t && !t1) {
      if(SA[c]>=0) {
        // SA[c] is borrowed by the right
        //   neighbor bucket.
        // shift-right the items in the
        //   right neighbor bucket.
        long long foo, bar;
        foo=SA[c];
        for(h=c+1; SA[h]>=0; h++)
        { bar=SA[h]; SA[h]=foo; foo=bar; }
        SA[h]=foo;

        SA[c]=EMPTY;
      }

      long long d=SA[c];
      if(d==EMPTY) { // SA[c] is empty.
        if(SA[c-1]==EMPTY) {
          SA[c]=-1; // init the counter.
          SA[c-1]=i;
        }
        else
          SA[c]=i; // a size-1 bucket.
      }
      else { // SA[c] is reused as a counter
          long long pos=c+d-1;
          if(SA[pos]!=EMPTY) {
            // we are running long longo the left
            //   neighbor bucket.
            // shift-right one step the items 
            //   of bucket(SA, S, i).
            for(h=0; h<-d; h++)
              SA[c-h]=SA[c-h-1];
            pos++;
          }
          else
            SA[c]--;

          SA[pos]=i;
      }
    }
  }

  // scan to shift-right the items in each bucket
  //   with its head being reused as a counter.
  for(i=n-1; i>0; i--) {
    j=SA[i];
    if(j<0 && j!=EMPTY) { // is SA[i] a counter?
      for(h=0; h<-j; h++)
        SA[i-h]=SA[i-h-1];
      SA[i-h]=EMPTY;
    }
  }

  // put the single sentinel LMS-substring.
  SA[0]=n-1;
}

unsigned long long getLengthOfLMS(unsigned char *s, 
			    unsigned long long n, long long level, unsigned long long x) {
  if(x==n-1) return 1;  
  
  unsigned long long dist, i=1;  
  while(1) {
    if(chr(x+i)<chr(x+i-1)) break;
    i++;
  }  
  while(1) {
    if(x+i>n-1 || chr(x+i)>chr(x+i-1)) break;
    if(x+i==n-1 || chr(x+i)<chr(x+i-1)) dist=i;
    i++;
  }
  
  return dist+1;
}

unsigned long long nameSubstr(unsigned long long *SA, 
			unsigned char *s, unsigned long long *s1, unsigned long long n, 
			unsigned long long m, unsigned long long n1, long long level) {
  unsigned long long i, j, cur_t, succ_t;

  // init the name array buffer
  for(i=n1; i<n; i++) SA[i]=EMPTY;

  // scan to compute the long longerim s1
  unsigned long long name, name_ctr=0;
  unsigned long long pre_pos, pre_len=0;
  for(i=0; i<n1; i++) {
    bool diff=false;
    unsigned long long len, pos=SA[i];

    len=getLengthOfLMS(s, n, level, pos);
    if(len!=pre_len) diff=true;
    else
      for(unsigned long long d=0; d<len; d++)
        if(pos+d==n-1 || pre_pos+d==n-1 ||
           chr(pos+d)!=chr(pre_pos+d)) {
          diff=true; break;
        }

    if(diff) {
      name=i; name_ctr++;
      SA[name]=1; // a new name.
      pre_pos=pos; pre_len=len;
    }
    else
      SA[name]++; // count this name.

    SA[n1+pos/2]=name;
  }

  // compact the long longerim s1 sparsely stored 
  //   in SA[n1, n-1] long longo SA[m-n1, m-1].
  for(i=n-1, j=m-1; i>=n1; i--)
    if(SA[i]!=EMPTY) SA[j--]=SA[i];

  // rename each S-type character of the
  //   long longerim s1 as the end of its bucket
  //   to produce the final s1.
  succ_t=1;
  for(i=n1-1; i>0; i--) {
    long long ch=s1[i], ch1=s1[i-1];
    cur_t=(ch1< ch || (ch1==ch && succ_t==1))?1:0;
    if(cur_t==1) {
      s1[i-1]+=SA[s1[i-1]]-1;
    }
    succ_t=cur_t;
  }

  return name_ctr;
}

void getSAlms(unsigned long long *SA, 
  unsigned char *s, 
  unsigned long long *s1, unsigned long long n, 
  unsigned long long n1, long long level ) {
  unsigned long long i, j, cur_t, succ_t;

  j=n1-1; s1[j--]=n-1;
  succ_t=0; // getbase(s,n-2) must be L-type
  for(i=n-2; i>0; i--) {
    cur_t=(chr(i-1)<chr(i) ||
          (chr(i-1)==chr(i) && succ_t==1))?1:0;
    if(cur_t==0 && succ_t==1) s1[j--]=i;
    succ_t=cur_t;
  }

  for(i=0; i<n1; i++) SA[i]=s1[SA[i]];
  
  // init SA[n1..n-1]
  for(i=n1; i<n; i++) SA[i]=level?EMPTY:0; 
}


void SACA_K(unsigned char *s, unsigned long long *SA,
	    unsigned long long n, unsigned long long K,
	    unsigned long long m, long long level) {
  unsigned long long i;
  unsigned long long *bkt=NULL;

  // stage 1: reduce the problem by at least 1/2.

  if(level==0) {
    bkt=(unsigned long long *)malloc(sizeof(long long)*K);
    putSubstr0(SA, s, bkt, n, K);
    induceSAl0(SA, s, bkt, n, K, false);
    induceSAs0(SA, s, bkt, n, K, false);
  }
  else {
    putSubstr1((long long *)SA, (long long *)s,(long long)n);
    induceSAl1((long long *)SA, (long long *)s, n ,false);
    induceSAs1((long long *)SA, (long long *)s, n, false);
  }

  // now, all the LMS-substrings are sorted and 
  //   stored sparsely in SA.

  // compact all the sorted substrings long longo
  //   the first n1 items of SA.
  // 2*n1 must be not larger than n.
  unsigned long long n1=0;
  for(i=0; i<n; i++) 
    if((!level&&SA[i]>0) || (level&&((long long *)SA)[i]>0))
      SA[n1++]=SA[i];

  unsigned long long *SA1=SA, *s1=SA+m-n1;
  unsigned long long name_ctr;
  name_ctr=nameSubstr(SA,s,s1,n,m,n1,level);

  // stage 2: solve the reduced problem.

  // recurse if names are not yet unique.
  if(name_ctr<n1)
    SACA_K((unsigned char *)s1, SA1, 
          n1, 0, m-n1, level+1);
  else // get the suffix array of s1 directly.
    for(i=0; i<n1; i++) SA1[s1[i]]=i;

  // stage 3: induce SA(S) from SA(S1).

  getSAlms(SA, s, s1, n, n1, level);
  if(level==0) {
    putSuffix0(SA, s, bkt, n, K, n1);
    induceSAl0(SA, s, bkt, n, K, true);
    induceSAs0(SA, s, bkt, n, K, true);
    free(bkt);
  }
  else {
    putSuffix1((long long *)SA, (long long *)s, n1);
    induceSAl1((long long *)SA, (long long *)s, n, true);
    induceSAs1((long long *)SA, (long long *)s, n, true);
  }
}

// This function is a drop-in replacement for histsort() (except
// that it expects a 0 bp after the sequence :]), which, while
// slower for len<10^9, also uses a lot less memory
unsigned long long *csuff_arr(const unsigned char *seq, unsigned long long len) {
  // seq is assumed to be given in compressed form form and be
  // null-terminated (having an long longernal zero byte is fine)
  // Testing, ahoy!
  unsigned long long *SA = malloc((len+1) * sizeof(long long));
  SACA_K((unsigned char *)seq, SA, len+1, 4 /* Not 256*/, len+1, 0);
  return SA;
}
