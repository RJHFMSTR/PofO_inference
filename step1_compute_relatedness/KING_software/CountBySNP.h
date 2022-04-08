#ifndef __CountBySNP_h__
#define __CountBySNP_h__                      
#include "IntArray.h"
void ComputeMZBySNP64Bit(IntArray &L0, IntArray &nonmissingMZCounts, IntArray &ibs1MZCounts, IntArray &HetHetMZCounts, IntArray &ibs0MZCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount);
void ComputeTrioBySNP64Bit(IntArray &Ltrio, IntArray &MItrioCounts, IntArray &HetInOffspringCounts, IntArray &nonmissingtrioCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount);
void ComputePOBySNP64Bit(IntArray &Lpo, IntArray &HomHomCounts, IntArray &ibs0Counts, IntArray &nonmissingCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount);
void ComputeAlleleFrequency64Bit(IntArray &subset, IntArray &AACounts, IntArray &AaCounts, IntArray &missingCounts, unsigned long long int **localLG[2], int localmarkerCount, int defaultMaxCoreCount);
#endif
