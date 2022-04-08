//////////////////////////////////////////////////////////////////////
// CountBySNP.cpp
// (c) 2010-2019 Wei-Min Chen
//
// This file is distributed as part of the KING source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile KING.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Oct 3, 2019

#include "IntArray.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void ComputeAlleleFrequency64Bit(IntArray &subset, IntArray &AACounts, IntArray &AaCounts, IntArray &missingCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount)
{
   const int BLOCKSIZE=255;
   const int CACHESIZE=16; // total cache size: 2^(6 + 4 - 2 + 8) = 2^16 = 64KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;  // 512
   unsigned long long int word, genobit[2][CACHESIZE_BITPAR];
   unsigned char *pchar;
   unsigned char byte;

   int base[8];
   for(int i = 0; i < 8; i++)
      base[i] = 1 << i;
   char revbase[256];
   for(int i = 0; i < 8; i++)
      revbase[base[i]] = i;
   char rightmost[256];
   for(int i = 0; i < 256; i++)
      rightmost[i] = revbase[i&(-i)];
   int thread = 0;
   int **genobitCount[3];
   for(int k = 0; k < 3; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
   int mymarkerCount = localmarkerCount;
   int mylongCount = (localmarkerCount-1)/64+1;
   int idCount = subset.Length();
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, byte, thread)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [(mylongCount<<6)];
      for(int m = 0; m < mymarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [mylongCount<<6];
      for(int m = 0; m < mymarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < idCount; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > idCount-BLOCKSIZE) ? idCount: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < mylongCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= mylongCount-CACHESIZE) ? mylongCount: blockb+CACHESIZE;
         for(int j = 0; j < CACHESIZE_BITPAR; j++)
            genobit[0][j] = genobit[1][j] = 0;
         for(int i = bi; i < iMax; i++){
             int id = subset[i];
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               for(int j = 0; j < 2; j++){
                  word = localLG[j][id][b];  // jth bit word
                  genobit[j][bbb] += (word&0x0101010101010101);
                  for(int k = 1; k < 8; k++)
                     genobit[j][bbb|k] += ((word>>k)&0x0101010101010101);
               }  // end of bit
               word = ~(localLG[0][id][b] | localLG[1][id][b]);  // Miss
               pchar = (unsigned char*)&word;
               for(int k = 0; k < 8; k++)
                  for(byte = pchar[k]; byte; byte &= (byte-1))
                     genobitCount[2][thread][(b<<6)|(k<<3)|rightmost[byte]] ++;
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 2; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   AACounts.Dimension(mymarkerCount);
   AACounts.Zero();
   AaCounts.Dimension(mymarkerCount);
   AaCounts.Zero();
   missingCounts.Dimension(mymarkerCount);
   missingCounts.Zero();
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < mymarkerCount; m++){
            AACounts[m] += genobitCount[0][i][m] + genobitCount[1][i][m] + genobitCount[2][i][m];
            AaCounts[m] += genobitCount[1][i][m];
            missingCounts[m] += genobitCount[2][i][m];
         }
         for(int k = 0; k < 3; k++)
            delete []genobitCount[k][i];
      }
   for(int m = 0; m < mymarkerCount; m++){
      AACounts[m] -= idCount;
      AaCounts[m] -= AACounts[m];
   }
   for(int k = 0; k < 3; k++)
      delete []genobitCount[k];
}

void ComputeMZBySNP64Bit(IntArray &L0, IntArray &nonmissingMZCounts, IntArray &ibs1MZCounts, IntArray &HetHetMZCounts, IntArray &ibs0MZCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount)
{
   int locallongCount = (localmarkerCount-1)/64+1;
   const int BLOCKSIZE=255;
   const int CACHESIZE=16; // total cache size: 2^(6 + 4 - 2 + 8 + 1) = 2^17 = 128KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;  // 128
   unsigned long long int word, genobit[4][CACHESIZE_BITPAR];
   unsigned char *pchar;
   int id1, id2;
   int thread = 0;
   int L0Count = (L0.Length()>>1);
   HetHetMZCounts.Dimension(localmarkerCount);
   ibs0MZCounts.Dimension(localmarkerCount);
   ibs1MZCounts.Dimension(localmarkerCount);
   nonmissingMZCounts.Dimension(localmarkerCount);
   HetHetMZCounts.Zero();
   ibs0MZCounts.Zero();
   ibs1MZCounts.Zero();
   nonmissingMZCounts.Zero();
   int **genobitCount[4];
   for(int k = 0; k < 4; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, thread, id1, id2)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 4; k++){
      genobitCount[k][thread] = new int [(locallongCount<<6)];
      for(int m = 0; m < localmarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 4; k++){
      genobitCount[k][thread] = new int [locallongCount<<6];
      for(int m = 0; m < localmarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < L0Count; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > L0Count-BLOCKSIZE) ? L0Count: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < locallongCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= locallongCount-CACHESIZE) ? locallongCount: blockb+CACHESIZE;
         for(int  k = 0;  k < 4; k++)
            for(int j = 0; j < CACHESIZE_BITPAR; j++)
               genobit[k][j] = 0;
         for(int i = bi; i < iMax; i++){
            id1 = L0[i*2];
            id2 = L0[i*2+1];
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               word = (localLG[0][id1][b] | localLG[1][id1][b]) & (localLG[0][id2][b] | localLG[1][id2][b]); // nonmissing
               genobit[0][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[0][bbb|k] += ((word>>k)&0x0101010101010101);
               word &= (localLG[0][id1][b] ^ localLG[0][id2][b]);  // IBS1
               genobit[1][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[1][bbb|k] += ((word>>k)&0x0101010101010101);
               word = (~localLG[0][id1][b]) & localLG[1][id1][b] & (~localLG[0][id2][b]) & localLG[1][id2][b];   // HetHet
               genobit[2][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[2][bbb|k] += ((word>>k)&0x0101010101010101);
               word = localLG[0][id1][b] & localLG[0][id2][b] & (localLG[1][id1][b]^localLG[1][id2][b]);   // IBS0
               genobit[3][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[3][bbb|k] += ((word>>k)&0x0101010101010101);
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 3; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }  // end of bit
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < localmarkerCount; m++){
            nonmissingMZCounts[m] += genobitCount[0][i][m];
            ibs1MZCounts[m] += genobitCount[1][i][m];
            HetHetMZCounts[m] += genobitCount[2][i][m];
            ibs0MZCounts[m] += genobitCount[3][i][m];
         }
         for(int k = 0; k < 4; k++)
            delete []genobitCount[k][i];
      }
   for(int k = 0; k < 4; k++)
      delete []genobitCount[k];
}

void ComputeTrioBySNP64Bit(IntArray &Ltrio, IntArray &MItrioCounts, IntArray &HetInOffspringCounts, IntArray &nonmissingtrioCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount)
{
   int locallongCount = (localmarkerCount-1)/64+1;
   const int BLOCKSIZE=255;
   const int CACHESIZE=8; // total cache size: 2^(6 + 4 - 2 + 8 + 1) = 2^17 = 128KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;
   unsigned long long int word, genobit[3][CACHESIZE_BITPAR];
   unsigned char *pchar;
   int id1, id2, id3;
   int thread = 0;
   int LtrioCount = Ltrio.Length()/3;
   MItrioCounts.Dimension(localmarkerCount);
   HetInOffspringCounts.Dimension(localmarkerCount);
   nonmissingtrioCounts.Dimension(localmarkerCount);
   MItrioCounts.Zero();
   HetInOffspringCounts.Zero();
   nonmissingtrioCounts.Zero();
   int **genobitCount[3];
   for(int k = 0; k < 3; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, thread, id1, id2, id3)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [(locallongCount<<6)];
      for(int m = 0; m < localmarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [locallongCount<<6];
      for(int m = 0; m < localmarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < LtrioCount; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > LtrioCount-BLOCKSIZE) ? LtrioCount: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < locallongCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= locallongCount-CACHESIZE) ? locallongCount: blockb+CACHESIZE;
         for(int  k = 0;  k < 3; k++)
            for(int j = 0; j < CACHESIZE_BITPAR; j++)
               genobit[k][j] = 0;
         for(int i = bi; i < iMax; i++){
            id1 = Ltrio[i*3+1];
            id2 = Ltrio[i*3+2];
            id3 = Ltrio[i*3];   // id3 is the child
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               word = localLG[0][id1][b] & localLG[0][id2][b] & (~(localLG[1][id1][b] ^ localLG[1][id2][b]))
               & (~localLG[0][id3][b]) & localLG[1][id3][b];   // AA x AA -> Aa, or aa x aa -> Aa
               genobit[0][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[0][bbb|k] += ((word>>k)&0x0101010101010101);
               word = (localLG[0][id1][b] | localLG[1][id1][b]) & (localLG[0][id2][b] | localLG[1][id2][b])
                  & (localLG[0][id3][b] | localLG[1][id3][b]);  // nonmissing
               genobit[1][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[1][bbb|k] += ((word>>k)&0x0101010101010101);
               word &= ((~localLG[0][id3][b]) & localLG[1][id3][b]); // HetInOffspring
               genobit[2][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[2][bbb|k] += ((word>>k)&0x0101010101010101);
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 3; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }  // end of bit
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < localmarkerCount; m++){
            MItrioCounts[m] += genobitCount[0][i][m];
            nonmissingtrioCounts[m] += genobitCount[1][i][m];
            HetInOffspringCounts[m] += genobitCount[2][i][m];
         }
         for(int k = 0; k < 3; k++)
            delete []genobitCount[k][i];
      }
   for(int k = 0; k < 3; k++)
      delete []genobitCount[k];
}

void ComputePOBySNP64Bit(IntArray &Lpo, IntArray &HomHomCounts, IntArray &ibs0Counts, IntArray &nonmissingCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount)
{
   int locallongCount = (localmarkerCount-1)/64+1;
   const int BLOCKSIZE=255;
   const int CACHESIZE=16; // total cache size: 2^(6 + 4 - 2 + 8 + 1) = 2^17 = 128KB
   const int CACHESIZE_BITPAR=CACHESIZE<<3;  // 128
   unsigned long long int word, genobit[3][CACHESIZE_BITPAR];
   unsigned char *pchar;
   int id1, id2;
   int thread = 0;
   int LpoCount = (Lpo.Length()>>1);
   HomHomCounts.Dimension(localmarkerCount);
   ibs0Counts.Dimension(localmarkerCount);
   nonmissingCounts.Dimension(localmarkerCount);
   HomHomCounts.Zero();
   ibs0Counts.Zero();
   nonmissingCounts.Zero();    
   int **genobitCount[3];
   for(int k = 0; k < 3; k++){
      genobitCount[k] = new int *[defaultMaxCoreCount];
      for(int i = 0; i < defaultMaxCoreCount; i++)
         genobitCount[k][i] = NULL;
   }
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(word, genobit, pchar, thread, id1, id2)
{
   thread = omp_get_thread_num();
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [(locallongCount<<6)];
      for(int m = 0; m < localmarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
   #pragma omp for
#else
   for(int k = 0; k < 3; k++){
      genobitCount[k][thread] = new int [locallongCount<<6];
      for(int m = 0; m < localmarkerCount; m++)
         genobitCount[k][thread][m] = 0;
   }
#endif
   for(int bi = 0; bi < LpoCount; bi+=BLOCKSIZE){   // blocked individuals
      int iMax = (bi > LpoCount-BLOCKSIZE) ? LpoCount: bi+BLOCKSIZE;
      for(int blockb = 0; blockb < locallongCount; blockb += CACHESIZE){  // blocked SNP words
         int bMax = (blockb >= locallongCount-CACHESIZE) ? locallongCount: blockb+CACHESIZE;
         for(int  k = 0;  k < 3; k++)
            for(int j = 0; j < CACHESIZE_BITPAR; j++)
               genobit[k][j] = 0;
         for(int i = bi; i < iMax; i++){
            id1 = Lpo[i*2];
            id2 = Lpo[i*2+1];
            for(int b = blockb; b < bMax; b++){
               int bbb = ((b-blockb)<<3);
               word = localLG[0][id1][b] & localLG[0][id2][b];   // HomHom
               genobit[0][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[0][bbb|k] += ((word>>k)&0x0101010101010101);
               word &= (localLG[1][id1][b]^localLG[1][id2][b]);   // IBS0
               genobit[1][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[1][bbb|k] += ((word>>k)&0x0101010101010101);
               word = (localLG[0][id1][b] | localLG[1][id1][b]) & (localLG[0][id2][b] | localLG[1][id2][b]); // nonmissing
               genobit[2][bbb] += (word&0x0101010101010101);
               for(int k = 1; k < 8; k++)
                  genobit[2][bbb|k] += ((word>>k)&0x0101010101010101);
            }  // end of blocked marker
         }  // end of blocked individual
         for(int b = blockb; b < bMax; b++){
            int bbb = ((b-blockb)<<3);
            for(int k = 0; k < 8; k++){   // shift
               int pos = (b<<6) | k;
               for(int bit = 0; bit < 3; bit++){
                  pchar = (unsigned char *)&genobit[bit][bbb|k];
                  for(int j = 0; j < 8; j++) // bit parallel
                     genobitCount[bit][thread][pos | (j<<3)] += pchar[j];
               }  // end of bit
            }  // end of 8 shifts
         }  // end of b
      }  // end of blockb
   }  // end of bi
#ifdef _OPENMP
}
#endif
   for(int i = 0; i < defaultMaxCoreCount; i++)
      if(genobitCount[0][i]){  // thread effectively used
         for(int m = 0; m < localmarkerCount; m++){
            HomHomCounts[m] += genobitCount[0][i][m];
            ibs0Counts[m] += genobitCount[1][i][m];
            nonmissingCounts[m] += genobitCount[2][i][m];
         }
         for(int k = 0; k < 3; k++)
            delete []genobitCount[k][i];
      }
   for(int k = 0; k < 3; k++)
      delete []genobitCount[k];
}


