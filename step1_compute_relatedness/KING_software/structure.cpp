//////////////////////////////////////////////////////////////////////
// structure.cpp
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
// Oct 11, 2019

#include <math.h>
#include <stdint.h>
#include <memory.h>
#include "analysis.h"
#include "Kinship.h"
#include "KinshipX.h"
#include "MathStats.h"
#include "MathSVD.h"
#include "QuickIndex.h"
#include "rplot.h"
#include "CountBySNP.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef WITH_LAPACK
extern "C" void dsyevr_(char*, char*, char*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, int*, double*, int*, int*, int*, int*);
extern "C" double dlamch_(char*);
#endif

#ifdef WITH_OPENBLAS
extern "C" void openblas_set_num_threads(int);
#endif

void Engine::WeightedInnerProduct64Bit(IntArray & subset, IntArray & AACounts, IntArray & AaCounts, IntArray & missingCounts, Matrix & IP)
{
    const double MISSINGRATE = 0.01;
    const double MINMAF = 0.01;
    int subsetCount = subset.Length();
    IP.Dimension(subsetCount, subsetCount); IP.Zero();
    const int WEIGHTSIZE = 32;
    const int DICTSIZE = 256;
    float weights[WEIGHTSIZE][DICTSIZE];
    unsigned long long int IBS[3];
    double freq[64];
    const int BLOCKSIZE = 16;
    const int CACHESIZE = 4;
    double localIP[BLOCKSIZE][BLOCKSIZE];
    IntArray loopIndex[2]; loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
    for (int i = 0; i < subsetCount; i += BLOCKSIZE)
        for (int j = i; j < subsetCount; j += BLOCKSIZE) {
            loopIndex[0].Push(i);
            loopIndex[1].Push(j);
        }
    int loopIndexLength = loopIndex[0].Length();
    double varweight[8];
    for (int w = 0; w < longCount; w += CACHESIZE) {
        for (int i = 0; i < WEIGHTSIZE; i++)
            for (int j = 0; j < DICTSIZE; j++)
                    weights[i][j] = 0.0;
        int wMax = ((w > longCount - CACHESIZE) ? longCount : w + CACHESIZE);
        for (int ww = w; ww < wMax; ww++) {
            for (int m = 0; m < 64; m++) {
                int pos = (ww << 6) | m;
                double p = ((missingCounts[pos] < subsetCount * MISSINGRATE) ?
                    (AACounts[pos] + AaCounts[pos] * 0.5) / (subsetCount - missingCounts[pos]) : 0);
                freq[m] = ((p < MINMAF || p > 1-MINMAF) ? 0 : p);
            }
            for (int block = 0; block < 8; block++) {
                for (int b = 0; b < 8; b++) {
                    double p = freq[(block<<3)|b];
                    varweight[b] = (p < MINMAF? 0: 0.25 / (p * (1 - p)));
                }
                int p2 = ((ww - w) << 3) | block;
                for (int byte = 0; byte < 256; byte++)
                    for (int bit = 0; bit < 8; bit++) 
                        if (byte & base[bit])
                            weights[p2][byte] += varweight[bit];
            }   // End of loop block
        }   // End of loop ww
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount) private(localIP, IBS)
#endif
        for (int k = 0; k < loopIndexLength; k++) {
            int i = loopIndex[0][k];
            int iMax = (i < subsetCount - BLOCKSIZE ? i + BLOCKSIZE : subsetCount);
            int j = loopIndex[1][k];
            int jMax = (j < subsetCount - BLOCKSIZE ? j + BLOCKSIZE : subsetCount);
            int jMin = j;
            for (int ii = 0; ii < BLOCKSIZE; ii++)
                for (int jj = 0; jj < BLOCKSIZE; jj++)
                    localIP[ii][jj] = 0;
            for (int ww = w; ww < wMax; ww++) {
                int offset_m = ((ww - w) << 3);
                for (int i1 = i; i1 < iMax; i1++) {
                    char ii = i1 - i;
                    int id1 = subset[i1];
                    for (int i2 = j; i2 < jMax; i2++) {
                        char jj = i2 - j;
                        int id2 = subset[i2];
                        IBS[0] = (LG[0][id1][ww] ^ LG[0][id2][ww]) & (LG[0][id1][ww] | LG[1][id1][ww]) & (LG[0][id2][ww] | LG[1][id2][ww]);// AA x Aa or Aa x aa
                        IBS[1] = LG[0][id1][ww] & LG[0][id2][ww] & (LG[1][id1][ww] ^ LG[1][id2][ww]);// AA x aa
                        double sum = 0.0;
                        for (int b = 0; b < 8; b++) {
                            int b8 = (b << 3);
                            sum += weights[offset_m | b][(IBS[0] >> b8) & 0xFF];
                            sum += weights[offset_m | b][(IBS[1] >> b8) & 0xFF] * 4;
                        }
                        localIP[ii][jj] += sum;
                    }   // End of loop i2
                }   // End of loop i1
            }   // End of loop ww
            for (int i1 = i; i1 < iMax; i1++) {
                char ii = i1 - i;
                if (i == j) jMin = i1+1;
                for (int i2 = jMin; i2 < jMax; i2++) {
                    char jj = i2 - j;
                    IP[i1][i2] -= localIP[ii][jj];
                }   // End of loop i2
            }   // End of loop i1
        }   // End of loop k
    }   // End of loop w
    
    int snpCount = 0;
    float weights0[8][DICTSIZE][3];
    double numerator[8][3];
    Vector IP2(subsetCount); IP2.Zero();
    for (int w = 0; w < longCount; w ++) {
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < DICTSIZE; j++)
                for (int k = 0; k < 3; k++)
                    weights0[i][j][k] = 0.0;
        for (int m = 0; m < 64; m++) {
            int pos = (w << 6) | m;
            double p = ((missingCounts[pos] < subsetCount*MISSINGRATE) ?
                (AACounts[pos] + AaCounts[pos] * 0.5) / (subsetCount - missingCounts[pos]) : 0);
            freq[m] = ((p < MINMAF || p > 1-MINMAF) ? 0 : p);
            if (freq[m] > MINMAF) snpCount++;
        }
        for (int block = 0; block < 8; block++) {
            for (int b = 0; b < 8; b++) {
                double p = freq[(block << 3) | b];
                if (p < MINMAF) varweight[b] = 0;
                else {
                    double q = 1 - p;
                    varweight[b] = 0.25 / (p * q);
                    numerator[b][0] = p * p * 4;
                    numerator[b][1] = (1 - p * 2)*(1 - p * 2);
                    numerator[b][2] = q * q * 4;
                }
            }
            for (int byte = 0; byte < 256; byte++)
                for (int bit = 0; bit < 8; bit++)
                    if (byte & base[bit])
                        for (int t = 0; t < 3; t++)
                            weights0[block][byte][t] += varweight[bit] * numerator[bit][t];
        }   // End of loop block
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount) private(IBS)
#endif
        for (int i = 0; i < subsetCount; i++) {
            int id = subset[i];
            IBS[0] = LG[0][id][w] & (~LG[1][id][w]);// aa
            IBS[1] = (~LG[0][id][w]) & LG[1][id][w];// Aa
            IBS[2] = LG[0][id][w] & LG[1][id][w];// AA
            double sum = 0.0;
            for (int b = 0; b < 8; b++) {
                int b8 = (b << 3);
                for (int t = 0; t < 3; t++)
                    sum += weights0[b][(IBS[t] >> b8) & 0xFF][t];
            }
            IP2[i] += sum;
        }   // End of loop i
    }   // End of loop w

    printf("  %d SNPs are used in PCA.\n", snpCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
    for (int i = 0; i < subsetCount; i++){
        IP[i][i] += IP2[i];
        for (int j = 0; j < subsetCount; j++) 
            IP[i][j] += IP2[i];
    }
    for (int i = 0; i < subsetCount; i++)
        for (int j = 0; j < i; j++) {
            IP[j][i] += IP[i][j];
            IP[i][j] = IP[j][i];
        }
}

void TopEigenForMDS(Matrix & D, int nPC, Vector & Values, Matrix & Vectors, bool distFlag, int coreCount)
{
    uintptr_t dimN = D.cols;
    int dimPC = dimN > nPC ? nPC : dimN;
    Vector Dmean(dimN);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int i = 0; i < dimN; i++) {
        double sum = 0.0;
        for (int j = 0; j < dimN; j++)
            sum += D[i][j];
        Dmean[i] = sum / dimN;
    }   // End of OMP i loop
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int i = 0; i < dimN; i++) {
        D[i].Subtract(Dmean);   // (I-11'/N) * D
        D[i].Add(-D[i].Sum() / dimN);    // D * (I-11'/N)
        if (distFlag) D[i].Multiply(-0.5);
    }   // End of OMP i loop
    Values.Dimension(dimPC);
    Vectors.Dimension(dimPC, dimN);
#ifdef WITH_LAPACK
    printf("  LAPACK is being used...\n"); fflush(stdout);
#ifdef WITH_OPENBLAS
    openblas_set_num_threads(coreCount);
#endif
    char JOBZ = 'V';
    char RANGE = 'I';
    char UPLO = 'U';
    int N = dimN;
    int info;
    double *A = new double[dimN*dimN];
    int LDA = N;
    for (uintptr_t i = 0; i < dimN; i++)    // transpose is actually used for symmetry
        memcpy(&(A[i*dimN]), &(D[i][0]), dimN * sizeof(double));//assign ith row to ith col
    double VL, VU;
    int M = dimPC;
    int IL = N - dimPC + 1;
    int IU = N;
    char CMACH = 'S';
    double ABSTOL = dlamch_(&CMACH);
    double *W = new double[N];
    double *Z = new double[M*N]; // N rows x M columns
    int LDZ = N;
    int *ISUPPZ = new int[M * 2];
    int INFO;
    double optim_lwork;
    int optim_liwork;
    int LWORK = -1;
    int LIWORK = -1;
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, &optim_lwork, &LWORK, &optim_liwork, &LIWORK, &INFO);
    LWORK = (int)optim_lwork;
    LIWORK = optim_liwork;
    double *WORK = new double[LWORK];
    int *IWORK = new int[LIWORK];
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    if (INFO != 0) error("SVD failed with INFO=%d.", INFO);
    delete[]ISUPPZ;
    delete[]IWORK;
    delete[]WORK;
    delete[]A;
    for (int j = 0; j < dimPC; j++) {
        Values[j] = W[dimPC - 1 - j];   //Next: assign (M-j)th col of Z to jth row of Vectors
        memcpy(&(Vectors[j][0]), &(Z[dimN*(dimPC - 1 - j)]), dimN * sizeof(double));
    }// for(int i = 0; i < N; i++) for(int j = 0; j < dimPC; j++) Vectors[i][j] = Z[(dimPC - 1 - j)*N + i];
    delete[]W;
    delete[]Z;
#else
    printf("  Please re-compile KING with LAPACK library.\n");
    SVD svd;
    int N = dimN;
    svd.Decompose(D);
    printf("done\n");
    if (svd.n == 0) return;
    QuickIndex idx;
    idx.Index(svd.w);
    for (int j = 0; j < dimPC; j++)
        Values[j] = svd.w[idx[N - 1 - j]];
    printf("\n");
    for (int i = 0; i < N; i++)
        for (int j = 0; j < dimPC; j++)
            Vectors[j][i] = svd.v[i][idx[N - 1 - j]];
#endif
}

void TopEigenForPCA(Matrix & D, int nPC, Vector & Values, Matrix & Vectors, int coreCount)
{
    uintptr_t dimN = D.cols;
    int dimPC = dimN > nPC ? nPC : dimN;
    Values.Dimension(dimPC);
    Vectors.Dimension(dimPC, dimN);
#ifdef WITH_LAPACK
    printf("  LAPACK is being used...\n"); fflush(stdout);
#ifdef WITH_OPENBLAS
    openblas_set_num_threads(coreCount);
#endif
    char JOBZ = 'V';
    char RANGE = 'I';
    char UPLO = 'U';
    int N = dimN;
    int info;
    double *A = new double[dimN*dimN];
    int LDA = N;
    for (uintptr_t i = 0; i < dimN; i++)    // transpose is actually used for symmetry
        memcpy(&(A[i*dimN]), &(D[i][0]), dimN * sizeof(double));//assign ith row to ith col
    double VL, VU;
    int M = dimPC;
    int IL = N - dimPC + 1;
    int IU = N;
    char CMACH = 'S';
    double ABSTOL = dlamch_(&CMACH);
    double *W = new double[N];
    double *Z = new double[M*N]; // N rows x M columns
    int LDZ = N;
    int *ISUPPZ = new int[M * 2];
    int INFO;
    double optim_lwork;
    int optim_liwork;
    int LWORK = -1;
    int LIWORK = -1;
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, &optim_lwork, &LWORK, &optim_liwork, &LIWORK, &INFO);
    LWORK = (int)optim_lwork;
    LIWORK = optim_liwork;
    double *WORK = new double[LWORK];
    int *IWORK = new int[LIWORK];
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    if (INFO != 0) error("SVD failed with INFO=%d.", INFO);
    delete[]ISUPPZ;
    delete[]IWORK;
    delete[]WORK;
    delete[]A;
    for (int j = 0; j < dimPC; j++) {
        Values[j] = W[dimPC - 1 - j];   //Next: assign (M-j)th col of Z to jth row of Vectors
        memcpy(&(Vectors[j][0]), &(Z[dimN*(dimPC - 1 - j)]), dimN * sizeof(double));
    }// for(int i = 0; i < N; i++) for(int j = 0; j < dimPC; j++) Vectors[i][j] = Z[(dimPC - 1 - j)*N + i];
    delete[]W;
    delete[]Z;
#else
    printf("  Please re-compile KING with LAPACK library.\n");
    SVD svd;
    int N = dimN;
    svd.Decompose(D);
    printf("done\n");
    if (svd.n == 0) return;
    QuickIndex idx;
    idx.Index(svd.w);
    for (int j = 0; j < dimPC; j++)
        Values[j] = svd.w[idx[N - 1 - j]];
    printf("\n");
    for (int i = 0; i < N; i++)
        for (int j = 0; j < dimPC; j++)
            Vectors[j][i] = svd.v[i][idx[N - 1 - j]];
#endif
}

void Engine::pca_projection_Internal64bit(IntArray &refList, IntArray &projList, Matrix &EV)
{
    int dimM = markerCount;
    if (dimM < 10) error("There are only %d autosomal SNPs", dimM);
    int dimN = refList.Length();
    if (dimN < 2) {
        printf("The number of reference samples for PC calculation is < 2.\n");
        return;
    }
    else if (dimN >= 65536) {
        printf("The number of reference samples for PC calculation is %d >= 65536.\n", dimN);
        return;
    }
    int projCount = projList.Length();
    int validCount = dimN + projCount;
    printf("Preparing matrix (%d x %d) for PCA...\n", dimN, dimN);
    Matrix IP;
    IntArray AACounts, AaCounts, missingCounts;
    ComputeAlleleFrequency64Bit(refList, AACounts, AaCounts, missingCounts, LG, longCount*64, defaultMaxCoreCount);
    WeightedInnerProduct64Bit(refList, AACounts, AaCounts, missingCounts, IP);
    printf("SVD starts at %s", currentTime()); fflush(stdout);
    int dimPC = dimN > nPC ? nPC : dimN;
    if (dimPC > 100) {
        dimPC = 100;
        printf("Only up to 100 PCs is allowed.\n");
    }
    Vector Values;
    Matrix Vectors;
    TopEigenForPCA(IP, dimPC, Values, Vectors, defaultMaxCoreCount);
    EV.Dimension(dimPC, validCount); EV.Zero();
    for (int k = 0; k < dimPC; k++)
        for (int i = 0; i < dimN; i++)
            EV[k][i] = Vectors[k][i];
    printf("Largest %d eigenvalues:", dimPC);
    for (int j = 0; j < dimPC; j++) printf(" %.2lf", sqrt(Values[j]));
    printf("\nPreparing scores at %s", currentTime()); fflush(stdout);
    Vector VectorSum(dimPC);
#ifdef _OPENMP
    int coreCount = defaultMaxCoreCount > dimPC ? dimPC : defaultMaxCoreCount;
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int k = 0; k < dimPC; k++) {
        Vectors[k].Multiply(1.0 / Values[k]);
        VectorSum[k] = Vectors[k].Sum();
    }
    Matrix scoreWT(longCount*64, dimPC);
    double weight[64][100];
    unsigned long long int G[2];
    double localPC[2][100];
    unsigned char MaskArray[256][9];
    for (int i = 0; i < 256; i++) {
        MaskArray[i][0] = 0;
        for (int k = 0; k < 8; k++)
            if (i & (1 << k)) {
                MaskArray[i][0]++;
                MaskArray[i][MaskArray[i][0]] = k;
            }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount) private(weight, G, localPC)
#endif
    for (int w = 0; w < longCount; w++) {
        for (int m = 0; m < 64; m++) 
            for (int k = 0; k < dimPC; k++)
                weight[m][k] = 0.0;
        for (int i = 0; i < dimN; i++) {
            for (int k = 0; k < dimPC; k++) {
                localPC[0][k] = Vectors[k][i];
                localPC[1][k] = Vectors[k][i] * 2;
            }
            int id = refList[i];
            G[0] = (~LG[0][id][w]) & LG[1][id][w];  // Aa
            G[1] = LG[0][id][w] & LG[1][id][w];     // AA
            for (int g = 0; g < 2; g++)
                for (int b = 0; b < 8; b++) {
                    int mbase = (b << 3);
                    unsigned char byte = ((G[g] >> mbase) & 0xFF);                    
                    int count = MaskArray[byte][0];
                    for (int m = 0; m < count; m++) {
                        int pos = mbase + MaskArray[byte][m + 1];
                        for (int k = 0; k < dimPC; k++)
                            weight[pos][k] += localPC[g][k];
                    }
                }
        }
        for (int m = 0; m < 64; m++) {
            int pos = (w << 6) | m;
            double p = (missingCounts[pos] < dimN * 0.01 ? (AACounts[pos] + AaCounts[pos] * 0.5) / (dimN - missingCounts[pos]) : 0);
            if (p < 0.001 || p > 0.999) 
                for (int k = 0; k < dimPC; k++)
                    scoreWT[pos][k] = 0;
            else {
                double mean = p * 2;
                double var = p * (1 - p) * 2;
                for (int k = 0; k < dimPC; k++)
                    scoreWT[pos][k] = (weight[m][k] - mean * VectorSum[k]) / var;
            }
        }
    }
    printf("Projecting %d samples starts at %s", validCount - dimN, currentTime()); fflush(stdout);
    float *weights[8][256];
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 256; j++)
            weights[i][j] = new float[dimPC];
    Vector mean(dimPC);
    for (int k = 0; k < dimPC; k++) mean[k] = 0.0;
    for (int w = 0; w < longCount; w++)
        for (int m = 0; m < 64; m++) {
            int pos = (w << 6) | m;
            double p = (missingCounts[pos] < dimN * 0.01 ? (AACounts[pos] + AaCounts[pos] * 0.5) / (dimN - missingCounts[pos]) : 0);
            if (p < 0.001 || p > 0.999) p = 0.0;
            for (int k = 0; k < dimPC; k++)
                mean[k] += p * 2 * scoreWT[w * 64 + m][k];
        }
    for (int k = 0; k < dimPC; k++)
        for (int i = 0; i < projCount; i++)
            EV[k][dimN + i] -= mean[k];
    for (int w = 0; w < longCount; w++) {
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 256; j++) {
                for (int k = 0; k < dimPC; k++)
                    weights[i][j][k] = 0;
                for (int bit = 0; bit < 8; bit++)
                    if (j & (1 << bit))
                        for (int k = 0; k < dimPC; k++)
                            weights[i][j][k] += scoreWT[w * 64 + i * 8 + bit][k];
            }
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount) private(G)
#endif
        for (int i = 0; i < projCount; i++) {
            int id = projList[i];
            G[0] = (~LG[0][id][w]) & LG[1][id][w];  // Aa
            G[1] = LG[0][id][w] & LG[1][id][w]; // AA
            for (int k = 0; k < dimPC; k++) {
                double sum = 0.0;
                for (int b = 0; b < 8; b++) {
                    int b8 = (b << 3);
                    sum += weights[b][(G[0] >> b8) & 0xFF][k];
                    sum += weights[b][(G[1] >> b8) & 0xFF][k] * 2;
                }
                EV[k][dimN + i] += sum;
            }
        }
    }
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 256; j++)
            delete[]weights[i][j];
}



void Engine::mds_projection_Internal64bit(IntArray &refList, IntArray &projList, Matrix &EV)
{
    int dimM = markerCount;
    if (dimM < 10) error("There are only %d autosomal SNPs", dimM);
    int dimN = refList.Length();
    if (dimN < 2) {
        printf("The number of reference samples for PC calculation is < 2.\n");
        return;
    }else if (dimN >= 65536) {
        printf("The number of reference samples for PC calculation is %d >= 65536.\n", dimN);
        return;
    }
    int validCount = dimN + projList.Length();
    int *missingInOnePersonCount = new int[idCount];
    for (int i = 0; i < idCount; i++) missingInOnePersonCount[i] = markerCount - ped[phenoid[i]].ngeno;
    IntArray *counts = new IntArray[dimN];
    for (int i = 0; i < dimN; i++)
        counts[i].Dimension(dimN);
    printf("Preparing matrix (%d x %d) for MDS...\n", dimN, dimN);
    ComputeInnerProduct64Bit(refList, counts);
    printf("SVD starts at %s", currentTime()); fflush(stdout);
    Matrix D(dimN, dimN);
    Vector Dmean(dimN);
#ifdef _OPENMP
    #pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
    for (int i = 0; i < dimN; i++) {
        for (int j = 0; j < i; j++)
            D[i][j] = (double)counts[j][i] / counts[i][j];
        D[i][i] = (double)counts[i][i] / (markerCount - missingInOnePersonCount[refList[i]]);
        for (int j = i + 1; j < dimN; j++)
            D[i][j] = (double)counts[i][j] / counts[j][i];
        double sum = 0.0;
        for (int j = 0; j < dimN; j++)
            sum += D[i][j];
        Dmean[i] = sum / dimN;
    }   // End of OMP i loop
    delete[]counts;
    int dimPC = dimN > nPC ? nPC : dimN;
    Vector Values;
    Matrix Vectors;
    TopEigenForMDS(D, dimPC, Values, Vectors, false, defaultMaxCoreCount);
    printf("Largest %d eigenvalues:", dimPC);
    for (int j = 0; j < dimPC; j++) printf(" %.2lf", Values[j]);
    printf("\nProjecting %d samples starts at %s", validCount - dimN, currentTime()); fflush(stdout);
    Vector subtractV(dimPC);
    Matrix rightMatrix = Vectors;
    EV.Dimension(dimPC, validCount);
#ifdef _OPENMP
    int coreCount = defaultMaxCoreCount > dimPC ? dimPC : defaultMaxCoreCount;
    #pragma omp parallel for num_threads(coreCount)
#endif
    for (int j = 0; j < dimPC; j++) {
        Vector tempV = rightMatrix[j];
        tempV.Multiply(1.0 / Values[j]);
        tempV.Add(-tempV.Sum() / dimN);
        subtractV[j] = Dmean.InnerProduct(tempV);
        rightMatrix[j] = tempV;
        memcpy(&(EV[j][0]), &(Vectors[j][0]), dimN * sizeof(double));
    }   // End of OMP dimPC loop
    const int BLOCKSIZE1 = 32;
    const int BLOCKSIZE2 = 32;
    int CACHESIZE = 128;   // cache size: BLOCKSIZE*CACHESIZE*32 = 2^17 = 128KB
    int ip[BLOCKSIZE1][BLOCKSIZE2], localMiss[BLOCKSIZE1][BLOCKSIZE1];
    unsigned long long int word0, word, word1, word2;
#ifdef _OPENMP
    #pragma omp parallel for num_threads(defaultMaxCoreCount) private(ip, localMiss, word, word0, word1, word2)
#endif
    for (int ii = dimN; ii < validCount; ii += BLOCKSIZE1) {
        int iMax = ii < validCount - BLOCKSIZE1 ? ii + BLOCKSIZE1 : validCount;
        for (int k = 0; k < dimPC; k++)
            for (int i = ii; i < iMax; i++)
                EV[k][i] = -subtractV[k];
        for (int jj = 0; jj < dimN; jj += BLOCKSIZE2) {
            int jMax = jj < dimN - BLOCKSIZE2 ? jj + BLOCKSIZE2 : dimN;
            for (int i = 0; i < BLOCKSIZE1; i++)
                for (int j = 0; j < BLOCKSIZE2; j++)
                    ip[i][j] = localMiss[i][j] = 0;
            for (int mm = 0; mm < longCount; mm += CACHESIZE) {
                int mMax = (mm > longCount - CACHESIZE) ? longCount : mm + CACHESIZE;
                for (int i = ii; i < iMax; i++) {
                    int id1 = projList[i - dimN];
                    for (int j = jj; j < jMax; j++) {
                        int id2 = refList[j];
                        word1 = word2 = 0;
                        int missmissCount = 0;
                        for (int m = mm; m < mMax; m++) {
                            for (word = ~(LG[0][id1][m] | LG[0][id2][m] | LG[1][id1][m] | LG[1][id2][m]);
                                word; word &= (word - 1), missmissCount++);  // MissMiss
                            word0 = LG[0][id1][m] & LG[0][id2][m];          // HomHom
                            word = word0 - ((word0 >> 1) & 0x5555555555555555);
                            word = (word & 0x3333333333333333) + ((word >> 2) & 0x3333333333333333);
                            word = (word + (word >> 4)) & 0x0F0F0F0F0F0F0F0F;
                            word1 += (word + (word >> 8)) & 0x00FF00FF00FF00FF;    //HomHom
                            word0 &= (LG[1][id1][m] ^ LG[1][id2][m]);       // IBS0
                            word = word0 - ((word0 >> 1) & 0x5555555555555555);  // IBS0
                            word = (word & 0x3333333333333333) + ((word >> 2) & 0x3333333333333333);
                            word = (word + (word >> 4)) & 0x0F0F0F0F0F0F0F0F;
                            word2 += (word + (word >> 8)) & 0x00FF00FF00FF00FF;    // IBS0
                        }  // ip: Inner product between g[i] and G[j]
                        word1 = (word1 + (word1 >> 16)) & 0x0000FFFF0000FFFF;
                        word2 = (word2 + (word2 >> 16)) & 0x0000FFFF0000FFFF;
                        ip[i - ii][j - jj] += ((word1 + (word1 >> 32)) & 0xFFFFFFFF) - (((word2 + (word2 >> 32)) & 0xFFFFFFFF) << 1);
                        localMiss[i - ii][j - jj] += missmissCount;
                    }  // End of j loop
                }  // End of i loop
            }  // End of mm loop
            for (int k = 0; k < dimPC; k++)
                for (int i = ii; i < iMax; i++) {
                    double sum = 0.0;
                    for (int j = jj; j < jMax; j++)
                        sum += ip[i - ii][j - jj] * rightMatrix[k][j] / (markerCount - missingInOnePersonCount[projList[i - dimN]] - missingInOnePersonCount[refList[j]] + localMiss[i - ii][j - jj]);
                    EV[k][i] += sum;
                }
        }// End of jj loop
    }   // End of OMP ii loop
    delete[]missingInOnePersonCount;
}

void Engine::mds_projection(bool flagMDS)
{
    if (Bit64 != 64) { printf("Only available for 64-bit system.\n"); return; }
    printf("\nOptions in effect:\n");
    if(flagMDS) printf("\t--mds\n");
    else printf("\t--pca\n");
    if (projectStart)
        printf("\t--projection %d\n", projectStart);
    else
        printf("\t--projection\n");
    if (Bit64Flag)
        printf("\t--sysbit 64\n");
    if (CoreCount)
        printf("\t--cpus %d\n", CoreCount);
    if (rplotFlag)
        printf("\t--rplot\n");
    if (prefix != "king")
        printf("\t--prefix %s\n", (const char*)prefix);
    printf("\n");

    String MDSorPCA = flagMDS ? "MDS" : "PCA";
    printf("%s projection starts at %s", (const char*)MDSorPCA, currentTime());
    int *missingInOnePersonCount = new int[idCount];
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
    for (int i = 0; i < idCount; i++) {
        int sum = 0;
        for (int m = 0; m < longCount; m++)   // not all non-missing
            for (unsigned long long int word = ~(LG[0][i][m] | LG[1][i][m]); word; word &= (word - 1), sum++);
        missingInOnePersonCount[i] = sum;
    }  // parallel among individuals ends
    for (int i = 0; i < ped.count; i++) ped[i].ngeno = 0;
    for (int i = 0; i < idCount; i++) ped[phenoid[i]].ngeno = markerCount - missingInOnePersonCount[i];
    int maxCount = projectStart ? projectStart : idCount;
    IntArray ID_AFF(0), ID_UN(0);
    bool KGflag = false;
    if (!popreffile.IsEmpty()) {
        FILE *fp = fopen((const char*)popreffile, "rt");
        if (fp == NULL) printf("Cannot open %s to read. Option --popref %s is ignored.\n",
            (const char*)popreffile, (const char*)popreffile);
        else {
            String line;
            line.ReadLine(fp);
            StringArray tokens;
            tokens.Clear();
            tokens.AddTokens(line);
            if (tokens.Length() < 3 || tokens[0] != "FID") error("Format: FID IID Population");
            StringArray FID(0), PID(0), POP(0);
            while (!feof(fp)) {
                line.ReadLine(fp);
                tokens.Clear();
                tokens.AddTokens(line);
                if (tokens.Length() < 3) continue;
                PID.Push(tokens[1]);
                POP.Push(tokens[2]);
            }
            fclose(fp);
            StringIntHash HashPop;
            HashPop.Clear();
            int count = POP.Length();
            printf("  %d samples are in the reference file %s\n", count, (const char*)popreffile);
            int popCount = 0;
            StringArray allPops(0);
            for (int i = 0; i < count; i++) {
                if (HashPop.Integer(POP[i]) == -1) {
                    HashPop.SetInteger(POP[i], popCount++);
                    allPops.Push(POP[i]);
                }
            }
            printf("  %d populations in the Reference:\n    ", popCount);
            for (int i = 0; i < popCount; i++)
                printf(" %s", (const char*)allPops[i]);
            printf("\n");
            StringIntHash HashRef;
            HashRef.Clear();
            for (int i = 0; i < count; i++)
                HashRef.SetInteger(PID[i], HashPop.Integer(POP[i]));
            String reffile = prefix;
            reffile.Add("_popref.txt");
            if (popreffile == reffile) {
                for (int i = 0; i < idCount; i++) {
                    int id = phenoid[i];
                    if (ped[id].ngeno >= MINSNPCOUNT) {
                        if (HashRef.Integer(ped[id].pid) > -1) ID_UN.Push(i);
                        else ID_AFF.Push(i);
                    }
                }
            }
            else {
                fp = fopen((const char*)reffile, "wt");
                if (fp == NULL) error("Cannot open %s to write.", (const char*)reffile);
                fprintf(fp, "FID IID Population\n");
                for (int i = 0; i < idCount; i++) {
                    int id = phenoid[i];
                    if (ped[id].ngeno >= MINSNPCOUNT) {
                        int pop = HashRef.Integer(ped[id].pid);
                        if (pop > -1) {
                            ID_UN.Push(i);
                            fprintf(fp, "%s %s %s\n",
                                (const char*)ped[id].famid, (const char*)ped[id].pid,
                                (const char*)allPops[pop]);
                        }
                        else ID_AFF.Push(i);
                    }
                }
                fclose(fp);
                printf("  Population reference saved in %s\n", (const char*)reffile);
            }
            KGflag = true;
            printf("  %d samples are in reference for MDS analysis.\n", ID_UN.Length());
        }
    }
    if (!KGflag) {   // Ref file not provided
        int unknownCount = 0;
        for (int i = 0; i < maxCount; i++) {
            int id = phenoid[i];
            if (ped[id].ngeno >= MINSNPCOUNT && ped[id].affections[0] == 0)
                unknownCount++;
        }
        if (unknownCount > 400) {
            StringIntHash HashKG;
            MakeHashKG(HashKG);
            int overlapCounts[5];
            for (int i = 0; i < 5; i++) overlapCounts[i] = 0;
            for (int i = 0; i < maxCount; i++) {
                int id = phenoid[i];
                if (ped[id].famid != ped[id].pid) continue;
                int pop = HashKG.Integer(ped[id].pid);
                if (ped[id].ngeno >= MINSNPCOUNT && ped[id].affections[0] == 0 && pop > -1 && pop < 6)
                    overlapCounts[pop - 1] ++;
            }
            int overlapCount = 0;
            for (int i = 0; i < 5; i++) overlapCount += overlapCounts[i];
            if (overlapCount > 2000 || (overlapCount > 400 && overlapCounts[3] == overlapCount)) KGflag = true;
            if (KGflag) {
                if (overlapCount > 2000) {
                    printf("%d 1000 Genomes samples are detected and used as reference.\n", overlapCount);
                    String reffile = prefix;
                    reffile.Add("_popref.txt");
                    FILE *fp = fopen((const char*)reffile, "wt");
                    if (fp == NULL) error("Cannot open %s to write.", (const char*)reffile);
                    fprintf(fp, "FID IID Population\n");
                    for (int i = 0; i < idCount; i++) {
                        int id = phenoid[i];
                        if (ped[id].ngeno >= MINSNPCOUNT) {
                            int pop = HashKG.Integer(ped[id].pid);
                            if (ped[id].famid == ped[id].pid && ped[id].affections[0] == 0 && pop > -1 && pop < 6) {
                                ID_UN.Push(i);
                                fprintf(fp, "%s %s ", (const char*)ped[id].famid, (const char*)ped[id].pid);
                                switch (pop) {
                                case 1: fprintf(fp, "AFR\n"); break;
                                case 2: fprintf(fp, "AMR\n"); break;
                                case 3: fprintf(fp, "EAS\n"); break;
                                case 4: fprintf(fp, "EUR\n"); break;
                                case 5: fprintf(fp, "SAS\n"); break;
                                }
                            }
                            else if (pop != 6)
                                ID_AFF.Push(i);
                        }
                    }
                    fclose(fp);
                }
                else {   // EUR only inference
                    StringIntHash HashEUR;
                    MakeHashEUR(HashEUR);
                    String reffile = prefix;
                    reffile.Add("_popref.txt");
                    FILE *fp = fopen((const char*)reffile, "wt");
                    if (fp == NULL) error("Cannot open %s to write.", (const char*)reffile);
                    fprintf(fp, "FID IID Population\n");
                    for (int i = 0; i < idCount; i++) {
                        int id = phenoid[i];
                        if (ped[id].ngeno >= MINSNPCOUNT) {
                            int pop = HashEUR.Integer(ped[id].pid);
                            if (ped[id].famid == ped[id].pid && ped[id].affections[0] == 0 && pop > -1 && pop < 4) {
                                ID_UN.Push(i);
                                fprintf(fp, "%s %s ", (const char*)ped[id].famid, (const char*)ped[id].pid);
                                switch (pop) {
                                case 1: fprintf(fp, "NEUR\n"); break;
                                case 2: fprintf(fp, "SEUR\n"); break;
                                case 3: fprintf(fp, "FIN\n"); break;
                                }
                            }
                            else if (pop != 4)
                                ID_AFF.Push(i);
                        }
                    }
                    fclose(fp);
                    printf("%d 1000 Genomes European samples are detected and used as reference.\n", ID_UN.Length());
                }
            }
        }
    }
    if (!KGflag) {   // No reference
        if (projectStart) {
            for (int i = 0; i < projectStart; i++)
                if (ped[phenoid[i]].ngeno >= MINSNPCOUNT)
                    ID_UN.Push(i);
            for (int i = projectStart; i < idCount; i++)
                if (ped[phenoid[i]].ngeno >= MINSNPCOUNT)
                    ID_AFF.Push(i);
            printf("The first %d samples are used as reference.\n", projectStart);
        }
        else {
            for (int i = 0; i < idCount; i++) {
                int id = phenoid[i];
                if (ped[id].ngeno >= MINSNPCOUNT) {
                    if (ped[id].affections[0] != 2)
                        ID_UN.Push(i);
                    else
                        ID_AFF.Push(i);
                }
            }
            printf("%d samples with unknown or unaffected status are used as reference\n", ID_UN.Length());
        }
    }
    int dimN = ID_UN.Length();
    int validCount = dimN + ID_AFF.Length();
    int dimPC = dimN > nPC ? nPC : dimN;
    if (dimPC > 100) {
        dimPC = 100;
        printf("Only up to 100 PCs is allowed.\n");
    }
    Matrix EV;
    if(flagMDS) mds_projection_Internal64bit(ID_UN, ID_AFF, EV);
    else pca_projection_Internal64bit(ID_UN, ID_AFF, EV);
    if (EV.cols != validCount) return;
    String pedfile = prefix;
    pedfile.Add("pc.txt");
    FILE *fp = fopen((const char*)pedfile, "wt");
    if(fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
    fprintf(fp, "FID IID FA MO SEX AFF");
    for(int j = 0; j < dimPC; j++)
        fprintf(fp, " PC%d", j+1);
    fprintf(fp, "\n");
    for(int i = 0; i < validCount; i++){
        int id = phenoid[i < dimN? ID_UN[i]: ID_AFF[i-dimN]];
        fprintf(fp, "%s %s %s %s %d %s",
            (const char*)ped[id].famid, (const char*)ped[id].pid,
            (const char*)ped[id].fatid, (const char*)ped[id].motid,
            ped[id].sex, i < dimN? "1": "2");
        for(int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("%s projection ends at %s", (const char*)MDSorPCA, currentTime());
    printf("%d principal components saved in file %s\n", dimPC, (const char*)pedfile);
}

void Engine::mds_family64Bit(bool flagMDS)
{
    printf("\nOptions in effect:\n");
    if(flagMDS) printf("\t--mds\n");
    else printf("\t--pca\n");
    if (Bit64Flag)
        printf("\t--sysbit 64\n");
    if (CoreCount)
        printf("\t--cpus %d\n", CoreCount);
    if (rplotFlag)
        printf("\t--rplot\n");
    if (prefix != "king")
        printf("\t--prefix %s\n", (const char*)prefix);
    printf("\n");
    String MDSorPCA = flagMDS ? "MDS" : "PCA";
    printf("%s for families starts at %s", (const char*)MDSorPCA, currentTime());
    printf("Genotypes stored in %d words for each of %d individuals.\n", longCount, idCount);
    Kinship kin;
    int *missingInOnePersonCount = new int[idCount];
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
    for (int i = 0; i < idCount; i++) {
        int sum = 0;
        for (int m = 0; m < longCount; m++)   // not all non-missing
            for (unsigned long long int word = ~(LG[0][i][m] | LG[1][i][m]); word; word &= (word - 1), sum++);
        missingInOnePersonCount[i] = sum;
    }  // parallel OMP among individuals ends
    for (int i = 0; i < ped.count; i++) ped[i].ngeno = 0;
    for (int i = 0; i < idCount; i++) ped[phenoid[i]].ngeno = markerCount - missingInOnePersonCount[i];
    IntArray ID(idCount); ID.Zero();
    IntArray ID_UN(0), ID_unrelated;
    for (int f = 0; f < ped.familyCount; f++) {
        ID_unrelated.Dimension(0);
        for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            if (ped[i].ngeno >= MINSNPCOUNT) {
                if (ped[i].isFounder())
                    ID_unrelated.Push(i);
            }
        kin.Setup(*ped.families[f]);
        for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++) {
            if (ped[i].ngeno < MINSNPCOUNT) continue;
            bool isUnrelated = true;
            for (int j = 0; j < ID_unrelated.Length(); j++)
                if (kin(ped[i], ped[ID_unrelated[j]]) > 0) {
                    isUnrelated = false;
                    break;
                }
            if (isUnrelated) ID_unrelated.Push(i);
        }
        for (int i = 0; i < ID_unrelated.Length(); i++) {
            ID_UN.Push(geno[ID_unrelated[i]]);
            ID[geno[ID_unrelated[i]]] = 1;
        }
    }
    IntArray ID_AFF(0);
    for (int i = 0; i < idCount; i++)
        if (ID[i] == 0 && ped[phenoid[i]].ngeno >= MINSNPCOUNT) ID_AFF.Push(i);
    int dimN = ID_UN.Length();
    int validCount = dimN + ID_AFF.Length();
    int dimPC = dimN > nPC ? nPC : dimN;
    if (dimPC > 100) {
        dimPC = 100;
        printf("Only up to 100 PCs is allowed.\n");
    }
    Matrix EV;
    if(flagMDS) mds_projection_Internal64bit(ID_UN, ID_AFF, EV);
    else pca_projection_Internal64bit(ID_UN, ID_AFF, EV);
    if (EV.cols != validCount) return;
    String pedfile = prefix;
    pedfile.Add("pc.txt");
    FILE *fp = fopen((const char*)pedfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
    fprintf(fp, "FID IID FA MO SEX AFF");
    for (int j = 0; j < dimPC; j++)
        fprintf(fp, " PC%d", j + 1);
    fprintf(fp, "\n");
    for (int i = 0; i < validCount; i++) {
        int id = phenoid[i < dimN ? ID_UN[i] : ID_AFF[i - dimN]];
        fprintf(fp, "%s %s %s %s %d %s",
            (const char*)ped[id].famid, (const char*)ped[id].pid,
            (const char*)ped[id].fatid, (const char*)ped[id].motid,
            ped[id].sex, i < dimN ? "1" : "2");
        for (int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", EV[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("%s for families ends at %s", (const char*)MDSorPCA, currentTime());
    printf("%d principal components saved in file %s\n", dimPC, (const char*)pedfile);
}

void Engine::mds()
{
    if (xflag) error("To be implemented");
    printf("\nOptions in effect:\n");
    printf("\t--mds\n");
    if (Bit64Flag)
        printf("\t--sysbit 64\n");
    if (CoreCount)
        printf("\t--cpus %d\n", CoreCount);
    if (rplotFlag)
        printf("\t--rplot\n");
    if (prefix != "king")
        printf("\t--prefix %s\n", (const char*)prefix);
    printf("\n");
    printf("MDS starts at %s", currentTime());
    printf("Genotypes stored in %d words for each of %d individuals.\n",
        Bit64 == 64 ? longCount : shortCount, idCount);
    IntArray ID(0);
    if (Bit64 == 64)
        for (int i = 0; i < idCount; i++) {
            int count = 0;
            for (int m = 0; m < longCount; m++)
                count += popcount(LG[0][i][m] | LG[1][i][m]);
            if (count >= MINSNPCOUNT) ID.Push(i);
        }
    else
        for (int i = 0; i < idCount; i++) {
            int count = 0;
            for (int m = 0; m < shortCount; m++) {
                int k = GG[0][i][m] | GG[1][i][m];
                count += oneCount[k & 255] + oneCount[(k >> 8) & 255];
            }
            if (count >= MINSNPCOUNT) ID.Push(i);
        }
    int dimN = ID.Length();
    if (dimN < 1) {
        printf("The number of individuals is < 1.\n");
        return;
    }
    else if (dimN > 46340) {
        printf("The number of individuals is > 46340.\n");
        return;
    }
    printf("Preparing matrix (%d x %d) for MDS...\n", dimN, dimN);
    Matrix D;
    if (Bit64 == 64)
        ComputeDistanceMatrix64Bit(ID, D); // computationally intensive here
    else
        ComputeDistanceMatrix(D);// computationally intensive here

    printf("SVD starts at %s", currentTime()); fflush(stdout);
    int dimPC = dimN > nPC ? nPC : dimN;
    Vector Values;
    Matrix Vectors;
    TopEigenForMDS(D, dimPC, Values, Vectors, true, defaultMaxCoreCount);

    String pedfile = prefix;
    pedfile.Add("pc.txt");
    FILE *fp = fopen((const char*)pedfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
    fprintf(fp, "FID IID FA MO SEX AFF");
    for (int j = 0; j < dimPC; j++)
        fprintf(fp, " PC%d", j + 1);
    fprintf(fp, "\n");
    for (int i = 0; i < ID.Length(); i++) {
        int id = ID[i];
        fprintf(fp, "%s %s %s %s %d",
            (const char*)ped[phenoid[id]].famid, (const char*)ped[phenoid[id]].pid,
            (const char*)ped[phenoid[id]].fatid, (const char*)ped[phenoid[id]].motid,
            ped[phenoid[id]].sex);
        fprintf(fp, " 1");
        for (int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", Vectors[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("MDS ends at %s", currentTime());
    printf("Largest %d eigenvalues:", dimPC);
    for (int j = 0; j < dimPC; j++) printf(" %.2lf", Values[j]);
    printf("\n");
    printf("%d principal components saved in file %s\n", dimPC, (const char*)pedfile);
}

void Engine::pca64Bit()
{
    printf("\nOptions in effect:\n");
    printf("\t--pca\n");
    if (Bit64Flag)
        printf("\t--sysbit 64\n");
    if (CoreCount)
        printf("\t--cpus %d\n", CoreCount);
    if (rplotFlag)
        printf("\t--rplot\n");
    if (prefix != "king")
        printf("\t--prefix %s\n", (const char*)prefix);
    printf("\n");
    printf("PCA starts at %s", currentTime());
    printf("Genotypes stored in %d words for each of %d individuals.\n", longCount, idCount);
    int *missingInOnePersonCount = new int[idCount];
#ifdef _OPENMP
#pragma omp parallel for num_threads(defaultMaxCoreCount)
#endif
    for (int i = 0; i < idCount; i++) {
        int sum = 0;
        for (int m = 0; m < longCount; m++)// not all non-missing
            for (unsigned long long int word = ~(LG[0][i][m] | LG[1][i][m]); word; word &= (word - 1), sum++);
        missingInOnePersonCount[i] = sum;
    }  // parallel among individuals ends
    IntArray ID(0);
    for (int i = 0; i < idCount; i++)   // Call Rate > 80%
        if (missingInOnePersonCount[i] < markerCount * 0.2) ID.Push(i);
    delete[]missingInOnePersonCount;
    int dimN = ID.Length();
    int dimM = markerCount;
    Matrix IP;
    printf("Preparing matrix (%d x %d) for PCA...\n", dimN, dimN);
    IntArray AACounts, AaCounts, missingCounts;
    ComputeAlleleFrequency64Bit(ID, AACounts, AaCounts, missingCounts, LG, longCount * 64, defaultMaxCoreCount);
    WeightedInnerProduct64Bit(ID, AACounts, AaCounts, missingCounts, IP);
    printf("SVD starts at %s", currentTime()); fflush(stdout);
    int dimPC = dimN > nPC ? nPC : dimN;
    Vector Values;
    Matrix Vectors;
    TopEigenForPCA(IP, dimPC, Values, Vectors, defaultMaxCoreCount);

    String pedfile = prefix;
    pedfile.Add("pc.txt");
    FILE *fp = fopen((const char*)pedfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)pedfile);
    fprintf(fp, "FID IID FA MO SEX AFF");
    for (int j = 0; j < dimPC; j++)
        fprintf(fp, " PC%d", j + 1);
    fprintf(fp, "\n");
    for (int i = 0; i < ID.Length(); i++) {
        int id = phenoid[ID[i]];
        fprintf(fp, "%s %s %s %s %d",
            (const char*)ped[id].famid, (const char*)ped[id].pid,
            (const char*)ped[id].fatid, (const char*)ped[id].motid,
            ped[id].sex);
        fprintf(fp, " 1");
        for (int j = 0; j < dimPC; j++)
            fprintf(fp, " %.4lf", Vectors[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("PCA ends at %s", currentTime());
    printf("Largest %d eigenvalues:", dimPC);
    for (int j = 0; j < dimPC; j++) printf(" %.2lf", sqrt(Values[j]));
    printf("\n");
    printf("%d principal components saved in file %s\n", dimPC, (const char*)pedfile);
}

