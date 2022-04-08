#ifndef __analysis_h__
#define __analysis_h__
#include "KingCore.h"

inline unsigned char popcount(unsigned long long int word)
{
   word = word - ((word>>1)&0x5555555555555555);
   word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
   word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
   word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
   word = (word+(word>>16)) & 0x0000FFFF0000FFFF;
   return((word+(word>>32)) & 0xFF);
}

inline unsigned char popcount(unsigned short int word)
{
   word = word - ((word>>1)&0x5555);
   word = (word&0x3333) + ((word>>2)&0x3333);
   word = (word+(word>>4)) & 0x0F0F;
   return((word+(word>>8)) & 0xFF);
}

class Engine:public KingEngine{
   protected:
      void printRelationship(int *beforeCount, int *afterCount);
      int *pAACount, *pAaCount, *paaCount;
      int *pxAACount, *pxAaCount, *pxaaCount;
      int *pyAACount, *pyAaCount, *pyaaCount;
      int *pmtAACount, *pmtAaCount, *pmtaaCount;
      IntArray covariatePC;
      int cAge;
      int missingBase;
      bool uniqueIID;
      Matrix *pedKin;
      void SemifamilyKinship(void);
      double EmpP(double, Vector &);
      void ComputeDistanceMatrix(Matrix & Dist);
      void ComputeDistanceMatrix64Bit(IntArray & subset, Matrix & Dist);
      void ComputeInnerProduct64Bit(IntArray & subset, IntArray * counts);
      double ComputeAUC(Vector & risks, IntArray & diseaseStatus, int printFlag);
      void mds_projection_Internal64bit(IntArray &refList, IntArray &projList, Matrix &EV);
      void pca_projection_Internal64bit(IntArray &refList, IntArray &projList, Matrix &EV);

   public:
      Engine(Pedigree &ped);
      ~Engine();

      int xflag;
// Relationship Inference
      int faster;
      int slower;
      String popreffile;
      double kinFilter;
      int relativedegree;
      int mincons;
      bool homogeneity;
      bool adjustFamily;
      bool rplotFlag;
      void internalKING(int);
      void IBD2SegInOnePair64Bit(int id1, int id2, IntArray &chrSeg, double totalLength, double &prop, double &maxLength);
      bool PreSegment(bool printFlag=true);
      void ROH();
	   void ROH64BitMain(unsigned long long int **localLG[2], int localmarkerCount, IntArray &localchrSeg, double localtotalLength, IntArray &localbp, Vector &rohprops, Vector &maxLengths, IntArray *allsegments);
      void ROHOnly(IntArray & idList, int segment, IntArray & rohStorage, IntArray & rohIndex, bool LengthOnly=false);
      void PopulationROH();
      void ComputeIBDSegment64Bit();
      void ComputeIBDSegment64BitWithFilter(IntArray & subsetRef, IntArray & subsetProj);
      void ComputeIBDSegmentMain64Bit(unsigned long long int **localLG[2], int localmarkerCount, IntArray & localchrSeg, double localtotalLength, IntArray &localbp);

      void ComputeBigDataDuplicate();
      void ComputeBigDataDuplicate64Bit();
      void ComputeBigDataSecondDegree();
      void ComputeBigDataDistant();
      void ComputeBigDataPO();
      void IntegratedRelationshipInference();
      int ScreenDuplicates(IntArray & rpList);
      long long int ScreenDuplicates64Bit(IntArray rpList[]);
      void ScreenCloseRelativesInSubset(IntArray rpList[]);
      void ScreenCloseRelativesInSubset64Bit(IntArray rpList[]);

      void ComputeShortFastHomoKinship();
      void ComputeShortExtendedIBS();
      void ComputeExtendedIBS64Bit();
      void ComputeShortRobustKinship();
      void ComputeLongRobustKinship64Bit();
      void ComputeLongRobustKinship64BitWithFilter(IntArray ids[], bool WriteFlag=true);
      void ComputeShortRobustXKinship();
      void ComputeLongRobustXKinship64Bit();

      void DuplicateInSubset64Bit(IntArray & pairList, IntArray & HetHetCounts,
         IntArray & DiffHomCounts, IntArray & HomHomCounts, IntArray & notMissingHetCounts);
      void KinshipInSubset(IntArray & pairList, IntArray & HetHetCounts,
         IntArray & IBS0Counts, IntArray & het1Counts, IntArray & het2Counts, IntArray & HomHomCounts, IntArray & IBS1Counts);
      void KinshipInSubset64Bit(IntArray & pairList, IntArray & HetHetCounts,
         IntArray & IBS0Counts, IntArray & het1Counts, IntArray & het2Counts, IntArray & HomHomCounts, IntArray & IBS1Counts);
      void IBD2SegInSubset(IntArray & pairList, Vector & ibd2props, Vector & maxLengths);
      void IBD2SegInSubset64Bit(IntArray & pairList, Vector & ibd2props, Vector & maxLengths);
      void IBDSegInSubset64Bit(IntArray & pairList, Vector & ibdprops, Vector & maxLengths, Vector & ibd2props, Vector & maxLengths2, unsigned long long int **localLG[2], int localmarkerCount, IntArray & localchrSeg, double localtotalLength, IntArray &localbp, IntArray *ibd1segs=NULL, IntArray *ibd2segs=NULL);
      void IBDSegOnly(IntArray & pairList, int segment, IntArray & ibdsegStorage1, IntArray & ibdsegIndex1, IntArray & ibdsegStorage2, IntArray & ibdsegIndex2, bool LengthOnly=false, int MINSEGLENGTH=2500000, int MINCCOUNT=200);
      void NPL();
      void HEreg();
      void IBDmapping(int nperm);
      void HomozygosityMapping();
      void HomozygosityMappingMH(const char *popName="Reference");
      void HomozygosityMappingForQT();
      void HomozygosityMappingForQTMH(const char *popName="Reference");
      void IBDGDT();
      void IBDMI();
      void IBDVC();
      void AUCmapping();
      void AncestryInference();
      void PopulationDistance();
      void PopulationIBD();
      void AUCpredicting(IntArray &allchr, IntArray &allpos);
      void LocalH2();

// IBD Analysis
      IntArray chrSeg, chrSegX;
      double totalLength, totalLengthX;
      String segmessage;

// Population Structure Analysis
      int nPC;
      void pca();
      void pca64Bit();
      void pca_family();
      void pca_projection();
      void mds();
      void mds_family64Bit(bool flagMDS);
      void mds_projection(bool flagMDS);
      void WeightedInnerProduct64Bit(IntArray & subset, IntArray & AACounts, IntArray & AaCounts, IntArray & missingCounts, Matrix & IP);
          
      bool semifamilyFlag;
      bool projectFlag;
      int projectStart;
      bool unrelatedExtraction;

// Pedigree Reconstruction
      char specialChar;
      bool rebuild(int id_added=1);
      int BuildOneFamily(int f, IntArray, double, String &);
      int ClusterFamily(int pedrebuildFlag, int degree);
      void rebuild_semifamily();
      bool SplitPedigree();

// QC
      void MakeFamilyForMI(void);
      void countGenotype(void);
      void OutputIndividualInfo(void);
      void QC_By_SNP(void);
      void QC_By_SNP64Bit(void);
      void QC_By_Sample(void);
      void QC_By_Sample64Bit(void);
      void QC_WipeMI(void);

// autoQC & autoplots
      char *idMask;
      unsigned short int *gMask, *xMask, *yMask, *mtMask;
      void Monomorphic_SNP(IntArray & removeList);
      void xHeterozygosity_SNP(IntArray & removeList, double xHeterozygosity);
      void CallRate_SNP(double rateFilter, IntArray & removeList);
      void CallRate_xSNP(double rateFilter, IntArray & removeList);
      void CallRate_ySNP(double rateFilter, IntArray & removeList, bool lessthanFlag=true);
      void CallRate_Sample(double rateFilter, IntArray & removeList);
      void autoQC(double samplecallrate, double snpcallrate/*, double xHeterozygosity*/);
// Tools
      int WriteFile;
      void WriteMerlin();
      void WritePlink();

// Binary Trait Association
      void ResidualQuantitative(int trait, int * vid, Vector & adjResid);
      void ResidualDisease(int * vid, Vector & adjResid);
      double ComputeDaviesP(double Q, Matrix & kernel);
      char KernelShort;
      void InvNorm(int tr);
      void SKAT();
      void SKAT_Batch_WeightedLinear(const char* skatfile);
      void SKAT_WeightedLinear();
      void SKAT_family_quantitative(int trait, StringArray &, IntArray *, Vector *, FILE *fp);
      void SKAT_family_quantitatives(const char* skatfile);
      void VCR_family_quantitative(int trait, StringArray &, IntArray *, Vector *, FILE *fp);
      void VCR_family_quantitatives(const char* vcrfile);
      void AdmixtureMapping(double, double);
      void MakeTrioForTDT(void);
      void TDT();
      void CountTDT(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts);
      void CountTDT64Bit(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts, unsigned long long int **localLG[2], int localmarkerCount);
      void CountTDTinX(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts);
      //      void CountTDTinX64Bit(IntArray &trioList, IntArray &ncounts, IntArray & ntcounts);
      void rareTDT();
      void permuteTDT();
      int permuteCount;
      double OnePermuteTDT(IntArray &);
      void permuteRareTDT();
      double OnePermuteRareTDT(IntArray &, int);
      int TDT_T[2];
      int TDT_PeakPos;
      IntArray rmarkers;
      double rareMAF;
      double chisqFilter;
      void HetSharing();
      void ExtractUnrelated(IntArray &unrelatedList, bool printFlag=false);
      void Haplotyping(bool SaveFlag=false);
      IntArray *hap[2];

// Quantitative Trait Association
      static double MACHEPS;
      static double MACHEPS_SQRT;
      static double cbrent;
      double minimize(double a, double b, double eps, double &funcx, int &numiter, int maxIter, bool quiet);

      int normalization;
      bool HeritFlag;
      bool slowFlag;
      void ROADTRIPS();

      int FixedEff;
      bool noiterFlag;
      bool svdoutFlag;
      String svdinfile;
      IntArray ID;
      double **UT;
      double **VR;
      double **D;
      double *varVR;
      void ReadSVD();
      void ComputeSimilarity();
      void ComputeScoreSVD(bool quiet);
      void ComputeSVD();
      void ComputeMTSVD();

      IntArray ID0;
      IntArray validpheno;
      int **missingpheno;
      String dosagefile;
      String dfamfile;
      String dmapfile;
      Vector lambda0, tau0;
      void PreScan_FixedEff();
      void PreScan_MTSCORE();
      void PreScan_RSCORE();
      void GenomeScan();
      void GenomeScan64Bit();
      void ScoreScan(FILE *fp);
      void GenomeScanWithPermutation(FILE *fp);
      void DosageScan(FILE *fp);
      Vector *chisqs;
      void PreVC();
      void PostVC();
      void VC();
      void LMMSCORE();
      void FASTASSOC();
      void GRAMMAR();
      void RSCORE();
      void LMM();
      double fLL(double x);
      Matrix UX;
      Vector UY;
      Vector EV;
      int currentT;
      void polygenic();
      void PrintPolygenic();
      void IBDGRM();
      void AllIBDSegments();

      void IBDMDS();
      void IBDMDS_Internal(IntArray & subset, Vector & EigenValue, Matrix & EigenVector);
      void IBDMDS_Projection();

      Matrix means;
      Vector variances, heritabilities;
      IntArray traits, covariates, SampleSize, diseases;
      StringArray traitList;
      StringArray covariateList;
      int unrelatedFlag;
      bool effectFlag;
      bool CheckCovariates(Person & p);
      void PrintGeneticRiskScoreSNPMajor(const char* weightfile, bool noflipFlag=false);
      double prevalence;
};
#endif

//      int *quality;
//      Vector BetaSum;
//      Matrix BetaSquareSum;
//      Matrix freqBeta;
 //     int diagnosis(void);
 //     int qualityT;
 //     void ComputeBigDataKinshipAdjustPC();

      //void mds_kin_family();
      //void mds_kin();
      //void mds_moving();
//     void TopEigenForMDS(Matrix & D, int nPC, Vector & Values, Matrix & Vectors, bool distFlag=true);
//      void pca_family_check();
//      void OneWindowProjection(int mv_chr, double mv_start, double mv_stop, Vector & ancestry);
//      void SlidingWindows();
//      Matrix localAncestry;
//     void mds_family_projection();
//      void mds_family();
//      void IBDLength_3Pairs(IntArray &trio, int order, Vector &lengths);
//      void IBDLength_trioMI(IntArray &trio, int order, Vector &lengths);
//      void AssocFamHistory();
