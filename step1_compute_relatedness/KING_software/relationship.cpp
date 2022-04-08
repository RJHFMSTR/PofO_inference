//////////////////////////////////////////////////////////////////////
// relationship.cpp
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
#include "analysis.h"
#include "Kinship.h"
#include "KinshipX.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

Engine::Engine(Pedigree &pedigree) :KingEngine(pedigree)
{
    traits.Dimension(0);
    homogeneity = false;
    adjustFamily = false;
    mafFilterFlag = false;
    detailFlag = false;
    pAACount = pAaCount = paaCount = pxAACount = pxAaCount = pxaaCount = NULL;
    pmtAACount = pmtAaCount = pmtaaCount = pyAACount = pyAaCount = pyaaCount = NULL;
    //   quality = NULL;
    relativedegree = 0;
    pedKin = NULL;
    semifamilyFlag = false;
    xflag = 0;
    chrList.Dimension(0);
    autoflipFlag = false;
    uniqueIID = false;
    rareMAF = _NAN_;
    hap[0] = hap[1] = NULL;
    chisqFilter = _NAN_;
    permuteCount = 0;
    projectFlag = false;
    projectStart = 0;
    unrelatedExtraction = false;
    normalization = 0;
    HeritFlag = false;
    slowFlag = false;
    svdinfile = "";
    svdoutFlag = false;
    noiterFlag = false;
    FixedEff = 1;
    dosagefile = dfamfile = dmapfile = "";
    D = NULL;
    effectFlag = false;
    bigdataFlag = false;
    faster = 0;
    slower = 0;
    specialChar = '_';
    autosomeOnly = true;
    genotypeOnly = true;
    CoreCount = 0;
    SaveFormat = "";
    idMask = NULL;
    gMask = xMask = yMask = mtMask = NULL;
    kinFilter = _NAN_;
    //   BetaSum.Dimension(0);
    lessmemFlag = false;
    prevalence = _NAN_;
    mincons = 0;
    rplotFlag = false;
    nPC = 10;
    popreffile.Clear();
}

void Engine::ComputeLongRobustXKinship64Bit()
{
   if(xmarkerCount < MINSNPCOUNT){
      printf("No sufficient X-chromosome SNPs (%d) available for the robust analysis.\n", xmarkerCount);
      return;
   }
   printf("X-chromosome genotypes stored in %d 64-bit words for each of %d individuals.\n",
      xlongCount, idCount);

   pxAACount = new int [idCount];
   pxAaCount = new int [idCount];
   pxaaCount = new int [idCount];
   for(int i = 0; i < idCount; i++)
      pxAACount[i] = pxAaCount[i] = pxaaCount[i] = 0;
   for(int k = 0; k < idCount; k++)
      for(int m = 0; m < xlongCount; m++){
         pxAACount[k] += popcount(XLG[0][k][m] & XLG[1][k][m]);
         pxAaCount[k] += popcount((~XLG[0][k][m]) & XLG[1][k][m]);
         pxaaCount[k] += popcount(XLG[0][k][m] & (~XLG[1][k][m]));
      }
   Vector hets, hetf;
   hetf.Dimension(ped.familyCount);
   hetf.Set(_NAN_);
   int type;
   String sextype[4];
   sextype[0] = "MM"; sextype[1] = "FM"; sextype[2] = "MF"; sextype[3] = "FF";
   for(int f = 0; f < ped.familyCount; f++){
      hets.Dimension(0);
      for(int i = 0; i < id[f].Length(); i++){
         if(ped[id[f][i]].sex == 1) continue;
         int k = geno[id[f][i]];
         double het = pxAaCount[k]*1.0/(pxAACount[k]+pxAaCount[k]+pxaaCount[k]);
         if(het > 0.05) hets.Push(het);
      }
      if(hets.Length()==1)
         hetf[f] = hets[0];
      else if(hets.Length()==2)
         hetf[f] = hets[0] < hets[1]? hets[0]: hets[1];
      else if(hets.Length()>2){
         hets.Sort();
         hetf[f] = hets[(hets.Length()-1)/2];
      }
   }
   String outfile;
   outfile.Copy(prefix);
   outfile.Add("X.kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tSex\tN_SNP\tPhiX\tHet\tIBS0\tKinshipX\n");
   double kinship;
   int id1, id2;
   KinshipX kinX;
   int HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi;
   double ibs0, het;

   for(int f = 0; f < ped.familyCount; f++){
      kinX.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            if(ped[id[f][i]].sex == 1 && ped[id[f][j]].sex == 1 && hetf[f]==_NAN_) // male-male
               continue;
            if(ped[id[f][i]].sex == 0 || ped[id[f][j]].sex == 0)
               continue;
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            HetHetCount = IBS0Count = het1Count = het2Count = notMissingCount = 0;
            type = 0;
            if(ped[id[f][i]].sex == 2) type |= 1;
            if(ped[id[f][j]].sex == 2) type |= 2;
            for(int m = 0; m < xlongCount; m++){
               IBS0Count += popcount(XLG[0][id1][m] & XLG[0][id2][m] & (XLG[1][id1][m] ^ XLG[1][id2][m]));
               notMissingCount += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (XLG[0][id2][m] | XLG[1][id2][m]));
            }
            if(type==3)
               for(int m = 0; m < xlongCount; m++){
                  HetHetCount += popcount((~XLG[0][id1][m]) & (XLG[1][id1][m]) & (~XLG[0][id2][m]) & XLG[1][id2][m]);
                  het1Count += popcount((XLG[0][id2][m] | XLG[1][id2][m]) & (~XLG[0][id1][m]) & XLG[1][id1][m]);
                  het2Count += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (~XLG[0][id2][m]) & XLG[1][id2][m]);
               }
            else{
               if(ped[id[f][i]].sex == 2) // id1 is female
                  for(int m = 0; m < xlongCount; m++)
                     het1Count += popcount((XLG[0][id2][m] | XLG[1][id2][m]) & (~XLG[0][id1][m]) & XLG[1][id1][m]);
               else if(ped[id[f][j]].sex == 2)// id2 is female
                  for(int m = 0; m < xlongCount; m++)
                     het1Count += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (~XLG[0][id2][m]) & XLG[1][id2][m]);
            }
            if(type == 3){
               if(het1Count+het2Count == 0) continue;
               kinship = (HetHetCount - IBS0Count*2.0) / (het1Count+het2Count);
               het = (het1Count+het2Count)*0.5/notMissingCount;
            }else if(type == 0){
               kinship = 0.75 - IBS0Count * 1.0 / notMissingCount / hetf[f];
               het = hetf[f];
            }else{
               if(het1Count == 0) continue;
               kinship = 0.5 - IBS0Count * 1.0 / het1Count;
               het = het1Count*1.0/notMissingCount;
            }
            phi = kinX(ped[id[f][i]], ped[id[f][j]]);
            ibs0 = IBS0Count*1.0/notMissingCount;
            fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.3lf\t%.4lf\t%.4lf\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, (const char*)sextype[type], notMissingCount, phi,
               het, ibs0, kinship);
         }
   }
   fclose(fp);

   printf("Within-family kinship data saved in file %s\n", (const char*)outfile);

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      return;
   }

   printf("Relationship inference across families starts at %s", currentTime());
   char buffer[1024];
   StringArray buffers;
   buffers.Dimension(ped.familyCount);
   for(int i = 0; i < buffers.Length(); i++)
      buffers[i].Clear();
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(HetHetCount, IBS0Count, notMissingCount, het1Count, het2Count, \
      type, id1, id2, kinship, phi, het, ibs0, buffer)
#endif
   for(int f1 = 0; f1 < ped.familyCount; f1++)
      for(int i = 0; i < id[f1].Length(); i++){
         if(ped[id[f1][i]].sex == 0) continue;
         id1 = geno[id[f1][i]];
         for(int f2 = f1+1; f2 < ped.familyCount; f2++)
         for(int j = 0; j < id[f2].Length(); j++){
            if(ped[id[f1][i]].sex == 1 && ped[id[f2][j]].sex == 1 && hetf[f1] == _NAN_ && hetf[f2] == _NAN_) // male-male
               continue;
            if(ped[id[f2][j]].sex == 0) continue;
            id2 = geno[id[f2][j]];
            HetHetCount = IBS0Count = het1Count = het2Count = notMissingCount = 0;
            type = 0;
            if(ped[id[f1][i]].sex == 2) type |= 1;
            if(ped[id[f2][j]].sex == 2) type |= 2;
            for(int m = 0; m < xlongCount; m++)
               IBS0Count += popcount(XLG[0][id1][m] & XLG[0][id2][m] & (XLG[1][id1][m] ^ XLG[1][id2][m]));
            if(type == 3){
               for(int m = 0; m < xlongCount; m++)
                  HetHetCount += popcount(~(XLG[0][id1][m]  | XLG[0][id2][m]) & (XLG[1][id1][m]) & XLG[1][id2][m]);
               if(relativedegree && HetHetCount < IBS0Count*2) continue;
               for(int m = 0; m < xlongCount; m++){
                  het1Count += popcount((XLG[0][id2][m] | XLG[1][id2][m]) & (~XLG[0][id1][m]) & XLG[1][id1][m]);
                  het2Count += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (~XLG[0][id2][m]) & XLG[1][id2][m]);
               }
            }else{
               if(ped[id[f1][i]].sex == 2) // id1 is female
                  for(int m = 0; m < xlongCount; m++)
                     het1Count += popcount((XLG[0][id2][m] | XLG[1][id2][m]) & (~XLG[0][id1][m]) & XLG[1][id1][m]);
               else if(ped[id[f2][j]].sex == 2) // id2 is female
                  for(int m = 0; m < xlongCount; m++)
                     het1Count += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (~XLG[0][id2][m]) & XLG[1][id2][m]);
            }
            if(type == 3){
               if(het1Count+het2Count == 0) continue;
               kinship = (HetHetCount - IBS0Count*2.0) / (het1Count+het2Count);
               if(relativedegree && kinship < 0) continue;
               for(int m = 0; m < xlongCount; m++)
                  notMissingCount += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (XLG[0][id2][m] | XLG[1][id2][m]));
               het = (het1Count+het2Count)*0.5/notMissingCount;
            }else if(type == 0){ // male-male
               if(hetf[f1]!=_NAN_ && hetf[f2]!=_NAN_)
                  het = hetf[f1] < hetf[f2]? hetf[f1]: hetf[f2];
               else if(hetf[f1]!=_NAN_)
                  het = hetf[f1];
               else
                  het = hetf[f2];
               for(int m = 0; m < xlongCount; m++)
                  notMissingCount += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (XLG[0][id2][m] | XLG[1][id2][m]));
               kinship = 0.75 - IBS0Count * 1.0 / notMissingCount / het;
               if(relativedegree && kinship < 0) continue;
            }else{ // male-female
               if(het1Count == 0) continue;
               kinship = 0.5 - IBS0Count * 1.0 / het1Count;
               if(relativedegree && kinship < 0) continue;
               for(int m = 0; m < xlongCount; m++)
                  notMissingCount += popcount((XLG[0][id1][m] | XLG[1][id1][m]) & (XLG[0][id2][m] | XLG[1][id2][m]));
               het = het1Count * 1.0 / notMissingCount;
            }
            sprintf(buffer, "%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\n",
               (const char*)ped[id[f1][i]].famid, (const char*)ped[id[f1][i]].pid,
               (const char*)ped[id[f2][j]].famid, (const char*)ped[id[f2][j]].pid,
               (const char*)sextype[type], notMissingCount, het, IBS0Count*1.0/notMissingCount, kinship);
            buffers[f1].Add(buffer);
         }
   }

   outfile.Copy(prefix);
   outfile.Add("X.kin0");
   fp = fopen(outfile, "wt");
   fprintf(fp, "FID1\tID1\tFID2\tID2\tSex\tN_SNP\tHet\tIBS0\tKinshipX\n");
   for(int b = 0; b < buffers.Length(); b++)
      buffers[b].Write(fp);
   fclose(fp);
   printf("                                         ends at %s", currentTime());
   printf("Between-family kinship data saved in file %s\n", (const char*)outfile);
}

void Engine::ComputeLongRobustKinship64BitWithFilter(IntArray ids[], bool WriteFlag)
{
   int m1, m2, m3;
   unsigned long long int word, word1, word2;
   const int cutoffMissingCount=markerCount-MINSNPCOUNT;
   int id1, id2;
   String outfile;
   FILE *fp;
   int notMissingCount;
   double kinshipcutoff = 0.3535534;
   if(relativedegree)
      kinshipcutoff /= (1<<relativedegree);  // 1: /2; 2: /4; 3: /8;
   int *missingInOnePersonCount = new int[idCount];
   int *hetInOnePersonCount = new int[idCount];
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) private(m1, m3, word)
#endif
   for(int i = 0; i < idCount; i++){
      for(m3 = m1 = 0; m1 < longCount; m1++)   // not all non-missing
         for(word = ~(LG[0][i][m1] | LG[1][i][m1]); word; word &= (word-1), m3++);
      missingInOnePersonCount[i] = m3;
      for(m3 = m1 = 0; m1 < longCount; m1++)   // Het
         m3 += popcount(~LG[0][i][m1] & LG[1][i][m1]);
      hetInOnePersonCount[i] = m3;
   }  // parallel among individuals ends

   if(WriteFlag){
      printf("Autosome genotypes stored in %d", longCount);
      printf(" words for each of %d individuals.\n", idCount);
      printf("\nOptions in effect:\n");
      printf("\t--kinship\n");
      if(relativedegree)
         printf("\t--degree %d\n", relativedegree);
      if(CoreCount)
         printf("\t--cpus %d\n", CoreCount);
      if(Bit64Flag)
         printf("\t--sysbit 64\n");
      if(lessmemFlag)
         printf("\t--lessmem\n");
      if(prefix!="king")
         printf("\t--prefix %s\n", (const char*)prefix);
      printf("\n");
      outfile.Copy(prefix);
      outfile.Add(".kin");
      fp = fopen(outfile, "wt");
      fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n");
      double kinship, inflation;
      Kinship kin;
      int HetHetCount, IBS0Count, het1Count;
      double phi, pi0, errorFlag;
      int beforeCount[6], afterCount[6];
      int degree; double ibs0;
      for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
      for(int f = 0; f < ped.familyCount; f++){
         kin.Setup(*ped.families[f]);
         for(int i = 0; i < id[f].Length(); i++)
            for(int j = i+1; j < id[f].Length(); j++){
               id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
               for(HetHetCount=IBS0Count=0, m1 = 0; m1 < longCount; m1++){
                  HetHetCount += popcount((~LG[0][id1][m1]) & (LG[1][id1][m1]) & (~LG[0][id2][m1]) & LG[1][id2][m1]);   // HetHet
                  for(word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
                     word; word &= (word-1), IBS0Count++);
               }
               if(HetHetCount==0 && IBS0Count==0) continue;
               het1Count = hetInOnePersonCount[id1] + hetInOnePersonCount[id2];
               notMissingCount = longCount*64 - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
               for(m1 = 0; m1 < longCount; m1++){
                  word1 = ~(LG[0][id1][m1] | LG[1][id1][m1]);
                  word2 = ~(LG[0][id2][m1] | LG[1][id2][m1]);
                  for(word = word1 & word2; word; word &= (word-1), notMissingCount++);
                  for(word = (~LG[0][id1][m1]) & LG[1][id1][m1] & word2; word; word &= (word-1), het1Count--);
                  for(word = word1 & (~LG[0][id2][m1]) & LG[1][id2][m1]; word; word &= (word-1), het1Count--);
               }
               kinship = (HetHetCount - IBS0Count*2.0) / het1Count;
               phi = kin(ped[id[f][i]], ped[id[f][j]]);
               pi0 = 0.0;
               if(phi < 0.2)
                  pi0 = 1-4*phi;
               else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
                  pi0 = 0.25;
               errorFlag = 0;
               inflation = (phi > 0)? kinship / phi: -1;
               if(phi > 0.03){  // up to 4th-degree relative
                  if(inflation > 2 || inflation < 0.5)
                     errorFlag = 1;
                  else if(inflation > 1.4142 || inflation < 0.70711)
                     errorFlag = 0.5;
               }else if(phi < 0.005){  // unrelated pair
                  if(kinship > 0.0442) errorFlag = 1;
                  else if(kinship > 0.0221) errorFlag = 0.5;
               }else{   // distant relatives
                  if(kinship < -0.0221) errorFlag = 1;
                  else errorFlag = 0.5;
               }
               ibs0 = IBS0Count*1.0/notMissingCount;
               degree = 4;
               if(phi > 0.0442)
                  degree = int(-log(phi)/log(2.0) - 0.5);
               if(degree < 4){
                  beforeCount[degree] ++;
                  if(degree == 1)
                     if(pi0==0) beforeCount[5] ++;
               }else
                  beforeCount[4] ++;
               degree = 4;
               if(kinship > 0.0442)
                  degree = int(-log(kinship)/log(2.0) - 0.5);
               if(degree < 4){
                  afterCount[degree] ++;
                  if(degree == 1){
                     if(errorrateCutoff == _NAN_){
                        if(ibs0 < 0.008)
                           afterCount[5] ++;
                     }else{
                        if(ibs0 < errorrateCutoff)
                           afterCount[5]++;
                     }
                  }
               }else
                  afterCount[4] ++;

               fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.3lf\t%.4lf\t%.4lf\t%G\n",
                  (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
                  (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
                  HetHetCount*1.0/notMissingCount, ibs0,
                  kinship, errorFlag);
            }
      }
      fclose(fp);
      bool pedigreeFlag = false;
      for(int i = 0; i < 6; i++)
         if(beforeCount[i]) pedigreeFlag = true;
      if(pedigreeFlag){
         printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
         printRelationship(beforeCount, afterCount);
      }else
         printf("Each family consists of one individual.\n");

      if(ped.familyCount < 2) {
         if(ped.familyCount==1 && ped.families[0]->famid=="0")
            warning("All individuals with family ID 0 are considered as relatives.\n");
         if(xmarkerCount >= MINSNPCOUNT){
            printf("\nX-chromosome analysis...\n");
            if(Bit64==64)
               ComputeLongRobustXKinship64Bit();
            else
               ComputeShortRobustXKinship();
         }
         printf("There is only one family.\n");
         return;
      }
      printf("Relationship inference across families starts at %s", currentTime());
   }

   IntArray subsetRef(0), subsetProj(0);
   bool projValid = projectFlag? true: false;
   if(projValid && projectStart==0 && ped.affectionCount==0){
      printf("Either affection status needs to be assigned for --projection analysis,\n");
      printf("  or a projection start position N can be specified as --projection N.\n");
      return;
   }
   if(projectStart){ // --kinship --proj N
      for(int i = 0; i < projectStart; i++) subsetRef.Push(i);
      for(int i = projectStart; i < idCount; i++) subsetProj.Push(i);
   }else if(projectFlag){  // --kinship --proj
      for(int i = 0; i < idCount; i++){
         int id = phenoid[i];
         if(ped[id].affections[0]==1) subsetRef.Push(i);
         else if(ped[id].affections[0]==2) subsetProj.Push(i);
      }
   }else{   // standard --kinship
      for(int i = 0; i < idCount; i++) subsetRef.Push(i);
      subsetProj = subsetRef;
   }
   int refCount = subsetRef.Length();
   int projCount = subsetProj.Length();
   if(refCount==0 || projCount==0){
      printf("No valid samples for inference.\n");
      return;
   }
   if(projectStart)
      printf("Only pairs between the first %d and the last %d individuals are inferred.\n", projectStart, idCount-projectStart);
   else if(projValid)
      printf("Only pairs between %d unaffected and %d affected individuals are inferred.\n",
         refCount, projCount);

   int thread = 0;
   int BLOCKSIZE=32;
   char SHIFTSIZE=5; // takes up to 2.1 million samples
   int CACHESIZE=128;   // cache size: BLOCKSIZE*CACHESIZE*32 = 2^17 = 128KB
   if(idCount > (1<<20)){  // sample size over a million
      int count = (idCount >> 20);
      int F = 0;
      for(; count; F++) count >>= 1;
      BLOCKSIZE <<= F;
      SHIFTSIZE += F;
      CACHESIZE >>= F;
   }
   subsetRef.Dimension(((refCount+BLOCKSIZE-1)>>SHIFTSIZE)<<SHIFTSIZE);
   subsetProj.Dimension(((projCount+BLOCKSIZE-1)>>SHIFTSIZE)<<SHIFTSIZE);
   unsigned int blockCount = (idCount-1)/BLOCKSIZE+1;
   unsigned int blockCountRef = (refCount-1)/BLOCKSIZE+1;
   unsigned int blockCountProj = (projCount-1)/BLOCKSIZE+1;
   unsigned int loopIndexLength = projValid? blockCountRef * blockCountProj:
      ((blockCountRef & 1) ? blockCountRef * ((blockCountRef+1)/2): (blockCountRef/2) * (blockCountRef+1));
   unsigned int *loopIndex = new unsigned int [loopIndexLength];
   int index = 0;
   for(unsigned int i = 0; i < refCount; i += BLOCKSIZE){
      unsigned int iShifted = (i<<(16-SHIFTSIZE));
      for(int j = (projValid? 0: i); j < projCount; j += BLOCKSIZE)
         loopIndex[index++] = iShifted | (j>>SHIFTSIZE);
   }
#ifdef _OPENMP
   if(WriteFlag) printf("%d CPU cores are used.\n", defaultMaxCoreCount);
#endif
   IntArray *ibuffer = new IntArray[defaultMaxCoreCount];
   Vector *rbuffer = new Vector[defaultMaxCoreCount];
   for(int c = 0; c < defaultMaxCoreCount; c++){
      ibuffer[c].Dimension(0);
      rbuffer[c].Dimension(0);
   }
   int firstCount = loopIndexLength / defaultMaxCoreCount;
   int onepercent = firstCount? (firstCount+99) / 100: 1;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(notMissingCount, m1, m2, m3, thread, word, word1, word2)
{
   thread = omp_get_thread_num();
#endif
   int **localD[5];
   for(int i = 0; i < 5; i++){
      localD[i] = new int * [BLOCKSIZE];
      for(int j = 0; j < BLOCKSIZE; j++)
         localD[i][j] = new int [BLOCKSIZE];
   }//   int localD[5][BLOCKSIZE][BLOCKSIZE];
   int **threshold = new int * [BLOCKSIZE];
   bool **stopFlag = new bool * [BLOCKSIZE];
   for(int i = 0; i < BLOCKSIZE; i++){
      threshold[i] = new int [BLOCKSIZE];
      stopFlag[i] = new bool [BLOCKSIZE];
   }//   int threshold[BLOCKSIZE][BLOCKSIZE];
   int *idBlock[2];
   for(int i = 0; i < 2; i++)
      idBlock[i] = new int [BLOCKSIZE];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(unsigned int k = 0; k < loopIndexLength; k++){
      if(k < firstCount && (k % onepercent) == 0) {
         printf("%d%%\r", k/onepercent);
         fflush(stdout);
      }
      int i = (loopIndex[k]>>16)<<SHIFTSIZE;
      int iMax = i<refCount-BLOCKSIZE? i+BLOCKSIZE: refCount;
      int j = (loopIndex[k]&0xFFFF)<<SHIFTSIZE;
      int jMax = j<projCount-BLOCKSIZE? j+BLOCKSIZE: projCount;
      int jMin = j;
      for(int c = 0; c < 5; c++)
         for(int ii = 0; ii < BLOCKSIZE; ii++)
            for(int jj = 0; jj < BLOCKSIZE; jj++)
               localD[c][ii][jj] = 0;
      for(int ii = 0; ii < BLOCKSIZE; ii++)
         for(int jj = 0; jj < BLOCKSIZE; jj++)
            stopFlag[ii][jj] = false;
      for(int s = 0; s < BLOCKSIZE; s++){
         idBlock[0][s] = subsetRef[i+s];
         idBlock[1][s] = subsetProj[j+s];
      }
      for(int i1 = i; i1 < iMax; i1++){
         char ii = i1 - i;
         int id1 = idBlock[0][ii];
         if(!projValid && i == j) jMin = i1 + 1;
         for(int i2 = jMin; i2 < jMax; i2++){
            char jj = i2 - j;
            int id2 = idBlock[1][jj];
            threshold[ii][jj] = int((2.0 - kinshipcutoff*4.0) *
               ((hetInOnePersonCount[id1] < hetInOnePersonCount[id2])?
               hetInOnePersonCount[id1]: hetInOnePersonCount[id2])) +
               missingInOnePersonCount[id1] + missingInOnePersonCount[id2];
            if( (missingInOnePersonCount[id1] >= cutoffMissingCount) ||
               (missingInOnePersonCount[id2] >= cutoffMissingCount) ||
               (ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) )
               stopFlag[ii][jj] = true;
         }
      }
      for(int mm = 0; mm < longCount; mm += CACHESIZE){
         int mMax = (mm > longCount-CACHESIZE) ? longCount: mm+CACHESIZE;
         for(int i1 = i; i1 < iMax; i1++){
            char ii = i1 - i;
            int id1 = idBlock[0][ii];
            if(!projValid && i == j) jMin = i1 + 1;
            for(int i2 = jMin; i2 < jMax; i2++){
               char jj = i2 - j;
               int id2 = idBlock[1][jj];
               if(stopFlag[ii][jj]) continue;
               for(word2=0, m1 = mm; m1 < mMax; m1++){
                  word = LG[0][id1][m1]^LG[0][id2][m1]; // HetHom
                  word = word - ((word>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               } // Lower bound of HetHom Count, minus HetMiss and MissMiss
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               localD[0][ii][jj] += ((word2+(word2>>32)) & 0xFFFFFFFF);//HetHomCount

               for(word2=0, m1 = mm; m1 < mMax; m1++){   // IBS0
                  word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
                  word = word - ((word>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               }
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               localD[4][ii][jj] += ((word2+(word2>>32)) & 0xFFFFFFFF);//IBS0Count

               if(localD[0][ii][jj] + (localD[4][ii][jj]<<2) > threshold[ii][jj])
                  stopFlag[ii][jj] = true;
            }  // end of id2
         }  // end of id1
      }// end of mm
      for(int i1 = i; i1 < iMax; i1++){
         char ii = i1 - i;
         int id1 = idBlock[0][ii];
         if(!projValid && i == j) jMin = i1 + 1;
         for(int i2 = jMin; i2 < jMax; i2++){
            char jj = i2 - j;
            if(stopFlag[ii][jj]) continue;
            int id2 = idBlock[1][jj];
            for(notMissingCount = m2 = m3 = 0, m1 = 0; m1 < longCount; m1++){
               word1 = ~(LG[0][id1][m1] | LG[1][id1][m1]);
               word2 = ~(LG[0][id2][m1] | LG[1][id2][m1]);
               for(word = word1 & word2; word; word &= (word-1), notMissingCount++);
               for(word = (~LG[0][id1][m1]) & LG[1][id1][m1] & word2; word; word &= (word-1), m2++);
               for(word = word1 & (~LG[0][id2][m1]) & LG[1][id2][m1]; word; word &= (word-1), m3++);
            }
            localD[0][ii][jj] += m2 + m3 + notMissingCount*2   // HetHomCount
               -(missingInOnePersonCount[id1] + missingInOnePersonCount[id2]);
            localD[1][ii][jj] = hetInOnePersonCount[id1]-m2;   // het1Count
            localD[2][ii][jj] = hetInOnePersonCount[id2]-m3;   // het2Count
            localD[3][ii][jj] = (longCount<<6) - missingInOnePersonCount[id1] - missingInOnePersonCount[id2] + notMissingCount;  // notMissingCount
            int het1Count = localD[1][ii][jj];
            int het2Count = localD[2][ii][jj];
            m3 = het1Count < het2Count? het1Count: het2Count;
            double kinship = (m3 == 0? 0.0: 0.5 - (localD[0][ii][jj]*0.25+localD[4][ii][jj]) / m3);
            notMissingCount = localD[3][ii][jj];
            if(kinship >= kinshipcutoff){ // save to ibuffer and rbuffer
               ibuffer[thread].Push(id1);
               ibuffer[thread].Push(id2);
               ibuffer[thread].Push(notMissingCount);
               rbuffer[thread].Push((het1Count+het2Count-localD[0][ii][jj])*0.5/notMissingCount);
               rbuffer[thread].Push(localD[4][ii][jj]*1.0/notMissingCount);
               rbuffer[thread].Push(kinship);
            }
         }  // end of id2 loop
      }  // end of id1 loop
   }  // end of pair loop
   for(int i = 0; i < 5; i++){
      for(int j = 0; j < BLOCKSIZE; j++)
         delete []localD[i][j];
      delete []localD[i];
   }
   for(int i = 0; i < BLOCKSIZE; i++){
      delete []stopFlag[i];
      delete []threshold[i];
   }
   delete []stopFlag;
   delete []threshold;
#ifdef _OPENMP
}  // extra bracket for omp
#endif
   delete []loopIndex;
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   if(WriteFlag){
      int pbuffer = 0;
      char buffer[0x20000];
      outfile.Copy(prefix);
      outfile.Add(".kin0");
      fp = fopen(outfile, "wb");
      pbuffer += sprintf(&buffer[pbuffer], "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n");
      long long int totalpairCount=0;
      for(int c = 0; c < defaultMaxCoreCount; c++){
         int count = ibuffer[c].Length()/3;
         for(int i = 0; i < count; i++){
            id1 = phenoid[ibuffer[c][i*3]];
            id2 = phenoid[ibuffer[c][i*3+1]];
            pbuffer += sprintf(&buffer[pbuffer], "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\n",
               (const char*)ped[id1].famid, (const char*)ped[id1].pid,
               (const char*)ped[id2].famid, (const char*)ped[id2].pid,
               ibuffer[c][i*3+2],
               rbuffer[c][i*3], rbuffer[c][i*3+1], rbuffer[c][i*3+2]);
            if(pbuffer > 0xFFFF){
               fwrite(buffer, 1, pbuffer, fp);
               pbuffer = 0;
            }
         }
         totalpairCount += count;
      }
      if(pbuffer) fwrite(buffer, 1, pbuffer, fp);
      fclose(fp);
      printf("                                         ends at %s", currentTime());
      printf("Between-family kinship data (up to degree %d, %lli pairs in total) saved in file %s\n",
         relativedegree, totalpairCount, (const char*)outfile);
   }else{
      for(int c = 0; c < defaultMaxCoreCount; c++){
         ids[c].Dimension(0);
         int count = ibuffer[c].Length()/3;
         for(int i = 0; i < count; i++){
            ids[c].Push(ibuffer[c][i*3]);
            ids[c].Push(ibuffer[c][i*3+1]);
         }
      }
   }
   delete []ibuffer;
   delete []rbuffer;
}


Engine::~Engine()
{
   if(pAACount) delete pAACount;
   if(pAaCount) delete pAaCount;
   if(paaCount) delete paaCount;
   if(pxAACount) delete pxAACount;
   if(pxAaCount) delete pxAaCount;
   if(pxaaCount) delete pxaaCount;
//   if(quality) delete quality;
   if(pedKin) delete []pedKin;
   if(hap[0]) {delete [](hap[0]); delete [](hap[1]);}
}

void Engine::ComputeExtendedIBS64Bit()
{
   printf("Autosome genotypes stored in %d", longCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--ibs\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   bool IBDvalidFlag = PreSegment(/*chrSeg, totalLength, segmessage*/);
   if(!IBDvalidFlag)
      printf("%s\n", (const char*)segmessage);

   int *missingInOnePersonCount = new int[idCount];
   int *hetInOnePersonCount = new int[idCount];
   int m1, m2, m3;
   unsigned long long int word, word1, word2;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(m1, m3, word)
#endif
   for(int i = 0; i < idCount; i++){
      for(m3 = m1 = 0; m1 < longCount; m1++)   // not all non-missing
         for(word = ~(LG[0][i][m1] | LG[1][i][m1]); word; word &= (word-1), m3++);
      missingInOnePersonCount[i] = m3;
      for(m3 = m1 = 0; m1 < longCount; m1++)   // Het
         m3 += popcount(~LG[0][i][m1] & LG[1][i][m1]);
      hetInOnePersonCount[i] = m3;
   }  // parallel among individuals ends
   const int cutoffMissingCount=markerCount-MINSNPCOUNT;
   IntArray TooFewSNPsList(0);
   for(int i = 0; i < idCount; i++)
      if(missingInOnePersonCount[i]>=cutoffMissingCount)
         TooFewSNPsList.Push(i);
   if(TooFewSNPsList.Length()){
      int count = TooFewSNPsList.Length();
      printf("The following %d samples are excluded from the kinship analysis (M<%d):\n",
         count, MINSNPCOUNT);
      for(int i = 0; i < count; i++)
         printf("\t(%s %s)", (const char*)ped[phenoid[i]].famid, (const char*)ped[phenoid[i]].pid);
      printf("\n\n");
   }
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".ibs");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tZ0\tPhi\tN_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\t");
   fprintf(fp, "IBS\tDist\tHetConc\tHet2|1\tHet1|2\tHomConc\tKinship");
   if(IBDvalidFlag)
      fprintf(fp, "\tMaxIBD2\tPr_IBD2");
   fprintf(fp, "\n");
   double kinship, inflation;
   int id1, id2;
   Kinship kin;
   int HetHomCount, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            for(HetHetCount=IBS0Count=0, m1 = 0; m1 < longCount; m1++){
               HetHetCount += popcount((~LG[0][id1][m1]) & (LG[1][id1][m1]) & (~LG[0][id2][m1]) & LG[1][id2][m1]);   // HetHet
               for(word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
                  word; word &= (word-1), IBS0Count++);
            }
            if(HetHetCount==0 && IBS0Count==0) continue;
            het1Count = hetInOnePersonCount[id1];
            het2Count = hetInOnePersonCount[id2];
            notMissingCount = longCount*64 - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            for(m1 = 0; m1 < longCount; m1++){
               word1 = ~(LG[0][id1][m1] | LG[1][id1][m1]);  // missingness in ID1
               word2 = ~(LG[0][id2][m1] | LG[1][id2][m1]);  // missingness in ID2
               for(word = word1 & word2; word; word &= (word-1), notMissingCount++);
               for(word = (~LG[0][id1][m1]) & LG[1][id1][m1] & word2; word; word &= (word-1), het1Count--);
               for(word = word1 & (~LG[0][id2][m1]) & LG[1][id2][m1]; word; word &= (word-1), het2Count--);
            }
            kinship = (HetHetCount - IBS0Count*2.0) / (het1Count+het2Count);

            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            errorFlag = 0;
            inflation = (phi > 0)? kinship / phi: -1;
            if(phi > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phi < 0.005){  // unrelated pair
               if(kinship > 0.0442) errorFlag = 1;
               else if(kinship > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinship < -0.0221) errorFlag = 1;
               else errorFlag = 0.5;
            }
            ibs0 = IBS0Count*1.0/notMissingCount;
            degree = 4;
            if(phi > 0.0442)
               degree = int(-log(phi)/log(2.0) - 0.5);
            if(degree < 4){
               beforeCount[degree] ++;
               if(degree == 1)
                  if(pi0==0) beforeCount[5] ++;
            }else
               beforeCount[4] ++;
            degree = 4;
            if(kinship > 0.0442)
               degree = int(-log(kinship)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1){
                  if(errorrateCutoff == _NAN_){
                    if(ibs0 < 0.008)
                        afterCount[5] ++;
                  }else{
                     if(ibs0 < errorrateCutoff)
                        afterCount[5]++;
                  }
               }
            }else
               afterCount[4] ++;
            fprintf(fp, "%s\t%s\t%s\t%.3lf\t%.4lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf",
               (const char*)ped[id[f][i]].famid,
               (const char*)ped[id[f][i]].pid, (const char*)ped[id[f][j]].pid,
               pi0, phi, notMissingCount, IBS0Count,  // N, N_IBS0
               het1Count+het2Count-HetHetCount*2,     // N_IBS1
               notMissingCount-IBS0Count-het1Count-het2Count+HetHetCount*2,   //N_IBS2
               HetHetCount, notMissingCount-het1Count-het2Count+HetHetCount,   // N_HetHet, N_HomHom
               het1Count, het2Count,                  // N_Het1, N_Het2
               2+(HetHetCount*2.0-IBS0Count*2.0-het1Count-het2Count)/notMissingCount, //IBS
               (IBS0Count*4.0 + het1Count + het2Count - HetHetCount*2.0) / notMissingCount,//Dist
               HetHetCount*1.0/(het1Count+het2Count-HetHetCount), // CHet
               HetHetCount*1.0/het1Count, HetHetCount*1.0/het2Count, // Pr(Het2|Het1), Pr(Het1|Het2)
               1-IBS0Count*1.0/(notMissingCount-het1Count-het2Count+HetHetCount), // CHom
               kinship);
            if(IBDvalidFlag){
               double prop;
               double maxLength;
               IBD2SegInOnePair64Bit(id1, id2, chrSeg, totalLength, prop, maxLength);
               fprintf(fp, "\t%.3lf\t%.4lf", maxLength, prop);
            }
            fprintf(fp, "\n");
         }
   }
   fclose(fp);
   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      printf("Within-family IBS data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      printf("There is only one family.\n");
      return;
   }

   printf("IBS and relationship inference across families starts at %s", currentTime());
   int thread = 0;
   FILE **fps;
   outfile.Copy(prefix);
   outfile.Add(".ibs0");

   const int BLOCKSIZE=8;
   const int CACHESIZE=256;
   int localD[5][BLOCKSIZE][BLOCKSIZE];
   IntArray loopIndex[2];
   loopIndex[0].Dimension(0); loopIndex[1].Dimension(0);
   for(int i = 0; i < idCount; i += BLOCKSIZE)
      for(int j = i; j < idCount; j += BLOCKSIZE){
         loopIndex[0].Push(i);
         loopIndex[1].Push(j);
      }
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   fps = new FILE *[defaultMaxCoreCount];
   StringArray outfiles(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.kin0");
      fps[c] = fopen(outfiles[c], "wt");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wt");
   fprintf(fps[0], "FID1\tID1\tFID2\tID2\tN_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\t");
   fprintf(fps[0], "IBS\tDist\tHetConc\tHet2|1\tHet1|2\tHomConc\tKinship");
   if(IBDvalidFlag)
      fprintf(fps[0], "\tMaxIBD2\tPr_IBD2");
   fprintf(fps[0],"\n");
   int loopIndexLength = loopIndex[0].Length();
   int firstCount = loopIndexLength / defaultMaxCoreCount;
   int onepercent = firstCount? (firstCount+99) / 100: 1;
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHomCount, IBS0Count, notMissingCount, het1Count, het2Count, \
      id1, id2, m1, m2, m3, thread, word, word1, word2, localD)
{
   thread = omp_get_thread_num();
   #pragma omp for
#endif
   for(int k = 0; k < loopIndexLength; k++){
      if(k < firstCount && (k % onepercent) == 0) {
         printf("%d%%\r", k/onepercent);
         fflush(stdout);
      }
      int i = loopIndex[0][k];
      int iMax = i<idCount-BLOCKSIZE? i+BLOCKSIZE: idCount;
      int j = loopIndex[1][k];
      int jMax = j<idCount-BLOCKSIZE? j+BLOCKSIZE: idCount;
      int jMin = j;
      for(int c = 0; c < 5; c++)
         for(int ii = 0; ii < BLOCKSIZE; ii++)
            for(int jj = 0; jj < BLOCKSIZE; jj++)
               localD[c][ii][jj] = 0;
         for(int mm = 0; mm < longCount; mm += CACHESIZE){
            int mMax = (mm > longCount-CACHESIZE) ? longCount: mm+CACHESIZE;
            for(id1 = i; id1 < iMax; id1++){
               char ii = id1 - i;
               if(i == j) jMin = id1 + 1;
               for(id2 = jMin; id2 < jMax; id2++){
                  char jj = id2 - j;
                  for(word2=0, m1 = mm; m1 < mMax; m1++){
                     word = LG[0][id1][m1]^LG[0][id2][m1]; // HetHom
                     word = word - ((word>>1)&0x5555555555555555);
                     word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                     word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                     word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
                  } // Lower bound of HetHom Count, minus HetMiss and MissMiss
                  word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
                  HetHomCount = (word2+(word2>>32)) & 0xFFFFFFFF;
                  for(word2=0, m1 = mm; m1 < mMax; m1++){   // IBS0
                     word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
                     word = word - ((word>>1)&0x5555555555555555);
                     word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                     word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                     word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
                  }
                  word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
                  IBS0Count = (word2+(word2>>32)) & 0xFFFFFFFF;

                  for(notMissingCount = m2 = m3 = 0, m1 = mm; m1 < mMax; m1++){
                     word1 = ~(LG[0][id1][m1] | LG[1][id1][m1]);
                     word2 = ~(LG[0][id2][m1] | LG[1][id2][m1]);
                     for(word = word1 & word2; word; word &= (word-1), notMissingCount++);
                     for(word = (~LG[0][id1][m1]) & LG[1][id1][m1] & word2; word; word &= (word-1), m2++);
                     for(word = word1 & (~LG[0][id2][m1]) & LG[1][id2][m1]; word; word &= (word-1), m3++);
                  }
                  localD[1][ii][jj] -= m2;   // het1Count
                  HetHomCount += m2;
                  localD[2][ii][jj] -= m3;   // het2Count
                  HetHomCount += m3;
                  localD[0][ii][jj] += HetHomCount + notMissingCount*2; // HetHomCount
                  localD[3][ii][jj] += notMissingCount;  // notMissingCount
                  localD[4][ii][jj] += IBS0Count;  // IBS0Count
               }  // end of id2
            }  // end of id1
         }// end of mm
      for(id1 = i; id1 < iMax; id1++){
         if(missingInOnePersonCount[id1]>=cutoffMissingCount) continue;
         char ii = id1 - i;
         if(i == j) jMin = id1 + 1;
         for(id2 = jMin; id2 < jMax; id2++){
            if(missingInOnePersonCount[id2]>=cutoffMissingCount) continue;
            if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
            char jj = id2 - j;
            localD[0][ii][jj] -= (missingInOnePersonCount[id1] + missingInOnePersonCount[id2]);
            localD[1][ii][jj] += hetInOnePersonCount[id1];  // het1Count
            localD[2][ii][jj] += hetInOnePersonCount[id2];  // het2Count
            localD[3][ii][jj] += longCount*64 - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];   // notMissingCount
            het1Count = localD[1][ii][jj];
            het2Count = localD[2][ii][jj];
            m3 = het1Count < het2Count? het1Count: het2Count;
            double kinship = (m3 == 0? 0.0: 0.5 - (localD[0][ii][jj]*0.25+localD[4][ii][jj]) / m3);
            notMissingCount = localD[3][ii][jj];
            fprintf(fps[thread], "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount, localD[4][ii][jj], localD[0][ii][jj],//N,N_IBS0,N_IBS1
               notMissingCount - localD[0][ii][jj] - localD[4][ii][jj],//N_IBS2
               (het1Count+het2Count-localD[0][ii][jj])/2, // N_HetHet
               notMissingCount-(het1Count+het2Count+localD[0][ii][jj])/2,  //N_HomHom
               het1Count, het2Count,                     // N_Het1, N_Het2
               2 - (localD[0][ii][jj] + localD[4][ii][jj]*2.0) / notMissingCount,//IBS
               (localD[0][ii][jj] + localD[4][ii][jj]*4.0) / notMissingCount,//Dist
               (het1Count+het2Count-localD[0][ii][jj])*1.0/(het1Count+het2Count+localD[0][ii][jj]),   // CHet
               (het1Count+het2Count-localD[0][ii][jj])*0.5 / het1Count, // Pr(Het2|Het1)
               (het1Count+het2Count-localD[0][ii][jj])*0.5 / het2Count, // Pr(Het1|Het2)
               1-localD[4][ii][jj]/(notMissingCount-(het1Count+het2Count+localD[0][ii][jj])*0.5),  // CHom
               kinship);
            if(IBDvalidFlag){
               if(kinship > 0.0884){
                  double prop;
                  double maxLength;
                  IBD2SegInOnePair64Bit(id1, id2, chrSeg, totalLength, prop, maxLength);
                  fprintf(fps[thread], "\t%.3lf\t%.4lf", maxLength, prop);
               }else
                  fprintf(fps[thread], "\t-9\t-9");
            }
            fprintf(fps[thread], "\n");
         }
      }
   }
#ifdef _OPENMP
}  // extra bracket for omp
#endif
#ifdef _OPENMP
   for(int c = 0; c < defaultMaxCoreCount; c++)
      fclose(fps[c]);
   fps[0] = fopen(outfile, "at");
   char buffer[1024];
   for(int c = 1; c < defaultMaxCoreCount; c++){
      fps[c] = fopen(outfiles[c], "rt");
      while(fgets(buffer, 1024, fps[c]) != NULL){
         fputs(buffer, fps[0]);
      }
      fclose(fps[c]);
      remove(outfiles[c]);
   }
#endif
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   fclose(fps[0]);
   delete []fps;
   printf("                                         ends at %s", currentTime());
   printf("Between-family IBS data saved in file %s\n", (const char*)outfile);
}

void Engine::ComputeLongRobustKinship64Bit()
{
   printf("Autosome genotypes stored in %d", longCount);
   printf(" words for each of %d individuals.\n", idCount);

   printf("\nOptions in effect:\n");
   printf("\t--kinship\n");
   if(projectStart) printf("\t--projection %d\n", projectStart);
   else if(projectFlag) printf("\t--projection\n");
   if(CoreCount)
      printf("\t--cpus %d\n", CoreCount);
   if(Bit64Flag)
      printf("\t--sysbit 64\n");
   if(lessmemFlag)
      printf("\t--lessmem\n");
   if(prefix!="king")
      printf("\t--prefix %s\n", (const char*)prefix);
   printf("\n");

   double kinshipcutoff = 0.3535534;
   int *missingInOnePersonCount = new int[idCount];
   int *hetInOnePersonCount = new int[idCount];
   int m1, m2, m3;
   unsigned long long int word, word1, word2;
#ifdef _OPENMP
   #pragma omp parallel for num_threads(defaultMaxCoreCount) \
      private(m1, m3, word)
#endif
   for(int i = 0; i < idCount; i++){
      for(m3 = m1 = 0; m1 < longCount; m1++)   // not all non-missing
         for(word = ~(LG[0][i][m1] | LG[1][i][m1]); word; word &= (word-1), m3++);
      missingInOnePersonCount[i] = m3;
      for(m3 = m1 = 0; m1 < longCount; m1++)   // Het
         m3 += popcount(~LG[0][i][m1] & LG[1][i][m1]);
      hetInOnePersonCount[i] = m3;
   }  // parallel among individuals ends     
   const int cutoffMissingCount=markerCount-MINSNPCOUNT;
   IntArray TooFewSNPsList(0);
   for(int i = 0; i < idCount; i++)
      if(missingInOnePersonCount[i]>=cutoffMissingCount)
         TooFewSNPsList.Push(i);
   if(TooFewSNPsList.Length()){
      int count = TooFewSNPsList.Length();
      printf("The following %d samples are excluded from the kinship analysis (M<%d):\n",
         count, MINSNPCOUNT);
      for(int i = 0; i < count; i++)
         printf("\t(%s %s)", (const char*)ped[phenoid[i]].famid, (const char*)ped[phenoid[i]].pid);
      printf("\n\n");
   }
   String outfile;
   outfile.Copy(prefix);
   outfile.Add(".kin");
   FILE *fp = fopen(outfile, "wt");
   fprintf(fp, "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n");
   double kinship, inflation;
   int id1, id2;
   Kinship kin;
   int HetHomCount, HetHetCount, IBS0Count, het1Count, het2Count, notMissingCount;
   double phi, pi0, errorFlag;
   int beforeCount[6], afterCount[6];
   int degree; double ibs0;
   for(int i = 0; i < 6; i++) beforeCount[i] = afterCount[i] = 0;
   for(int f = 0; f < ped.familyCount; f++){
      kin.Setup(*ped.families[f]);
      for(int i = 0; i < id[f].Length(); i++)
         for(int j = i+1; j < id[f].Length(); j++){
            id1 = geno[id[f][i]]; id2 = geno[id[f][j]];
            for(HetHetCount=IBS0Count=0, m1 = 0; m1 < longCount; m1++){
               HetHetCount += popcount((~LG[0][id1][m1]) & (LG[1][id1][m1]) & (~LG[0][id2][m1]) & LG[1][id2][m1]);   // HetHet
               for(word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
                  word; word &= (word-1), IBS0Count++);
            }
            if(HetHetCount==0 && IBS0Count==0) continue;
            het1Count = hetInOnePersonCount[id1] + hetInOnePersonCount[id2];
            notMissingCount = longCount*64 - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];
            for(m1 = 0; m1 < longCount; m1++){
               word1 = ~(LG[0][id1][m1] | LG[1][id1][m1]);
               word2 = ~(LG[0][id2][m1] | LG[1][id2][m1]);
               for(word = word1 & word2; word; word &= (word-1), notMissingCount++);
               for(word = (~LG[0][id1][m1]) & LG[1][id1][m1] & word2; word; word &= (word-1), het1Count--);
               for(word = word1 & (~LG[0][id2][m1]) & LG[1][id2][m1]; word; word &= (word-1), het1Count--);
            }
            kinship = (HetHetCount - IBS0Count*2.0) / het1Count;

            phi = kin(ped[id[f][i]], ped[id[f][j]]);
            pi0 = 0.0;
            if(phi < 0.2)
               pi0 = 1-4*phi;
            else if(phi < 0.3 && ped[id[f][i]].isSib(ped[id[f][j]]))
               pi0 = 0.25;
            errorFlag = 0;
            inflation = (phi > 0)? kinship / phi: -1;
            if(phi > 0.03){  // up to 4th-degree relative
               if(inflation > 2 || inflation < 0.5)
                  errorFlag = 1;
               else if(inflation > 1.4142 || inflation < 0.70711)
                  errorFlag = 0.5;
            }else if(phi < 0.005){  // unrelated pair
               if(kinship > 0.0442) errorFlag = 1;
               else if(kinship > 0.0221) errorFlag = 0.5;
            }else{   // distant relatives
               if(kinship < -0.0221) errorFlag = 1;
               else errorFlag = 0.5;
            }
            ibs0 = IBS0Count*1.0/notMissingCount;
            degree = 4;
            if(phi > 0.0442)
               degree = int(-log(phi)/log(2.0) - 0.5);
            if(degree < 4){
               beforeCount[degree] ++;
               if(degree == 1)
                  if(pi0==0) beforeCount[5] ++;
            }else
               beforeCount[4] ++;
            degree = 4;
            if(kinship > 0.0442)
               degree = int(-log(kinship)/log(2.0) - 0.5);
            if(degree < 4){
               afterCount[degree] ++;
               if(degree == 1){
                  if(errorrateCutoff == _NAN_){
                    if(ibs0 < 0.008)
                        afterCount[5] ++;
                  }else{
                     if(ibs0 < errorrateCutoff)
                        afterCount[5]++;
                  }
               }
            }else
               afterCount[4] ++;
            fprintf(fp, "%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.3lf\t%.4lf\t%.4lf\t%G\n",
               (const char*)ped[id[f][i]].famid, (const char*)ped[id[f][i]].pid,
               (const char*)ped[id[f][j]].pid, notMissingCount, pi0, phi,
               HetHetCount*1.0/notMissingCount, ibs0,
               kinship, errorFlag);
         }
   }
   fclose(fp);
   bool pedigreeFlag = false;
   for(int i = 0; i < 6; i++)
      if(beforeCount[i]) pedigreeFlag = true;
   if(pedigreeFlag){
      printf("Within-family kinship data saved in file %s\n", (const char*)outfile);
      printRelationship(beforeCount, afterCount);
   }else
      printf("Each family consists of one individual.\n");

   if(ped.familyCount < 2) {
      if(ped.familyCount==1 && ped.families[0]->famid=="0")
         warning("All individuals with family ID 0 are considered as relatives.\n");
      if(xmarkerCount >= MINSNPCOUNT){
         printf("\nX-chromosome analysis...\n");
         if(Bit64==64)
            ComputeLongRobustXKinship64Bit();
         else
            ComputeShortRobustXKinship();
      }
      printf("There is only one family.\n");
      return;
   }

   printf("Relationship inference across families starts at %s", currentTime());
   IntArray subsetRef(0), subsetProj(0);
   bool projValid = projectFlag? true: false;
   if(projValid && projectStart==0 && ped.affectionCount==0){
      printf("Either affection status needs to be assigned for --projection analysis,\n");
      printf("  or the count of the first part (N) can be specified as --projection N.\n");
      return;
   }
   if(projectStart){ // --kinship --proj N
      for(int i = 0; i < projectStart; i++) subsetRef.Push(i);
      for(int i = projectStart; i < idCount; i++) subsetProj.Push(i);
   }else if(projectFlag){  // --kinship --proj
      for(int i = 0; i < idCount; i++){
         int id = phenoid[i];
         if(ped[id].affections[0]==1) subsetRef.Push(i);
         else if(ped[id].affections[0]==2) subsetProj.Push(i);
      }
   }else{   // standard --kinship
      for(int i = 0; i < idCount; i++) subsetRef.Push(i);
      subsetProj = subsetRef;
   }
   int refCount = subsetRef.Length();
   int projCount = subsetProj.Length();
   if(refCount==0 || projCount==0){
      printf("No valid samples for inference.\n");
      return;
   }
   if(projectStart)
      printf("Only pairs between the first %d and the last %d individuals are inferred.\n", projectStart, idCount-projectStart);
   else if(projValid)
      printf("Only pairs between %d unaffected and %d affected individuals are inferred.\n",
         refCount, projCount);
   int thread = 0;
   FILE **fps;
   outfile.Copy(prefix);
   outfile.Add(".kin0");
   int BLOCKSIZE=32;
   char SHIFTSIZE=5; // takes up to 2.1 million samples
   int CACHESIZE=128;   // cache size: BLOCKSIZE*CACHESIZE*32 = 2^17 = 128KB
   if(idCount > (1<<20)){  // sample size over a million
      int count = (idCount >> 20);
      int F = 0;
      for(; count; F++) count >>= 1;
      BLOCKSIZE <<= F;
      SHIFTSIZE += F;
      CACHESIZE >>= F;
   }
   subsetRef.Dimension(((refCount+BLOCKSIZE-1)>>SHIFTSIZE)<<SHIFTSIZE);
   subsetProj.Dimension(((projCount+BLOCKSIZE-1)>>SHIFTSIZE)<<SHIFTSIZE);
   unsigned int blockCount = (idCount-1)/BLOCKSIZE+1;
   unsigned int blockCountRef = (refCount-1)/BLOCKSIZE+1;
   unsigned int blockCountProj = (projCount-1)/BLOCKSIZE+1;
   unsigned int loopIndexLength = projValid? blockCountRef * blockCountProj:
      ((blockCountRef & 1) ? blockCountRef * ((blockCountRef+1)/2): (blockCountRef/2) * (blockCountRef+1));
   unsigned int *loopIndex = new unsigned int [loopIndexLength];
   int index = 0;
   for(unsigned int i = 0; i < refCount; i += BLOCKSIZE){
      unsigned int iShifted = (i<<(16-SHIFTSIZE));
      for(int j = (projValid?0:i); j < projCount; j += BLOCKSIZE)
         loopIndex[index++] = iShifted | (j>>SHIFTSIZE);
   }
   int firstCount = loopIndexLength / defaultMaxCoreCount;
   int onepercent = firstCount? (firstCount+99) / 100: 1;
   StringArray outfiles;
   fps = new FILE *[defaultMaxCoreCount];
   int pbuffer = 0;
   char buffer[0x20000];
#ifdef _OPENMP
   printf("%d CPU cores are used.\n", defaultMaxCoreCount);
   outfiles.Dimension(defaultMaxCoreCount);
   for(int c = 1; c < defaultMaxCoreCount; c++){
      outfiles[c].Copy(prefix);
      outfiles[c] += (c+1);
      outfiles[c].Add("$$$.kin0");
      fps[c] = fopen(outfiles[c], "wb");
   }
#else
   fps = new FILE *[1];
#endif
   fps[0] = fopen(outfile, "wb");
#ifdef _OPENMP
   #pragma omp parallel num_threads(defaultMaxCoreCount) \
      private(HetHomCount, IBS0Count, notMissingCount, het1Count, het2Count, \
      m1, m2, m3, thread, word, word1, word2, buffer, pbuffer)
{
   pbuffer = 0;
   thread = omp_get_thread_num();
   if(thread==0)
      pbuffer += sprintf(&buffer[pbuffer], "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n");
#endif
   int **localD[5];
   for(int i = 0; i < 5; i++){
      localD[i] = new int * [BLOCKSIZE];
      for(int j = 0; j < BLOCKSIZE; j++)
         localD[i][j] = new int [BLOCKSIZE];
   }//   int localD[5][BLOCKSIZE][BLOCKSIZE];
   int *idBlock[2];
   for(int i = 0; i < 2; i++)
      idBlock[i] = new int [BLOCKSIZE];
#ifdef _OPENMP
   #pragma omp for
#endif
   for(unsigned int k = 0; k < loopIndexLength; k++){
      if(k < firstCount && (k % onepercent) == 0) {
         printf("%d%%\r", k/onepercent);
         fflush(stdout);
      }
      int i = ((loopIndex[k]>>16)<<SHIFTSIZE);
      int iMax = i<refCount-BLOCKSIZE? i+BLOCKSIZE: refCount;
      int j = ((loopIndex[k]&0xFFFF)<<SHIFTSIZE);
      int jMax = j<projCount-BLOCKSIZE? j+BLOCKSIZE: projCount;
      int jMin = j;
      for(int c = 0; c < 5; c++)
         for(int ii = 0; ii < BLOCKSIZE; ii++)
            for(int jj = 0; jj < BLOCKSIZE; jj++)
               localD[c][ii][jj] = 0;
      for(int s = 0; s < BLOCKSIZE; s++){
         idBlock[0][s] = subsetRef[i+s];
         idBlock[1][s] = subsetProj[j+s];
      }
      for(int mm = 0; mm < longCount; mm += CACHESIZE){
         int mMax = (mm > longCount-CACHESIZE) ? longCount: mm+CACHESIZE;
         for(int i1 = i; i1 < iMax; i1++){
            char ii = i1 - i;
            int id1 = idBlock[0][ii];
            if(!projValid && i == j) jMin = i1 + 1;
            for(int i2 = jMin; i2 < jMax; i2++){
               char jj = i2 - j;
               int id2 = idBlock[1][jj];
               for(word2=0, m1 = mm; m1 < mMax; m1++){
                  word = LG[0][id1][m1]^LG[0][id2][m1]; // HetHom
                  word = word - ((word>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               } // Lower bound of HetHom Count, minus HetMiss and MissMiss
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               HetHomCount = (word2+(word2>>32)) & 0xFFFFFFFF;
               for(word2=0, m1 = mm; m1 < mMax; m1++){   // IBS0
                  word = LG[0][id1][m1] & LG[0][id2][m1] & (LG[1][id1][m1] ^ LG[1][id2][m1]);
                  word = word - ((word>>1)&0x5555555555555555);
                  word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
                  word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
                  word2 += (word+(word>>8)) & 0x00FF00FF00FF00FF;
               }
               word2 = (word2+(word2>>16)) & 0x0000FFFF0000FFFF;
               IBS0Count = (word2+(word2>>32)) & 0xFFFFFFFF;

               for(notMissingCount = m2 = m3 = 0, m1 = mm; m1 < mMax; m1++){
                  word1 = ~(LG[0][id1][m1] | LG[1][id1][m1]);
                  word2 = ~(LG[0][id2][m1] | LG[1][id2][m1]);
                  for(word = word1 & word2; word; word &= (word-1), notMissingCount++);
                  for(word = (~LG[0][id1][m1]) & LG[1][id1][m1] & word2; word; word &= (word-1), m2++);
                     for(word = word1 & (~LG[0][id2][m1]) & LG[1][id2][m1]; word; word &= (word-1), m3++);
               }
               localD[1][ii][jj] -= m2;   // het1Count
               HetHomCount += m2;
               localD[2][ii][jj] -= m3;   // het2Count
               HetHomCount += m3;
               localD[0][ii][jj] += HetHomCount + notMissingCount*2; // HetHomCount
               localD[3][ii][jj] += notMissingCount;  // notMissingCount
               localD[4][ii][jj] += IBS0Count;  // IBS0Count
            }  // end of id2
         }  // end of id1
      }// end of mm
      for(int i1 = i; i1 < iMax; i1++){
         char ii = i1 - i;
         int id1 = idBlock[0][ii];
         if(missingInOnePersonCount[id1]>=cutoffMissingCount) continue;
         if(!projValid && i == j) jMin = i1 + 1;
         for(int i2 = jMin; i2 < jMax; i2++){
            char jj = i2 - j;
            int id2 = idBlock[1][jj];
            if(missingInOnePersonCount[id2]>=cutoffMissingCount) continue;
            if(ped[phenoid[id1]].famid == ped[phenoid[id2]].famid) continue;
            localD[0][ii][jj] -= (missingInOnePersonCount[id1] + missingInOnePersonCount[id2]);
            localD[1][ii][jj] += hetInOnePersonCount[id1];  // het1Count
            localD[2][ii][jj] += hetInOnePersonCount[id2];  // het2Count
            localD[3][ii][jj] += (longCount<<6) - missingInOnePersonCount[id1] - missingInOnePersonCount[id2];   // notMissingCount
            het1Count = localD[1][ii][jj];
            het2Count = localD[2][ii][jj];
            m3 = het1Count < het2Count? het1Count: het2Count;
            double kinship = (m3 == 0? 0.0: 0.5 - (localD[0][ii][jj]*0.25+localD[4][ii][jj]) / m3);
            notMissingCount = localD[3][ii][jj];
            pbuffer += sprintf(&buffer[pbuffer], "%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\n",
               (const char*)ped[phenoid[id1]].famid, (const char*)ped[phenoid[id1]].pid,
               (const char*)ped[phenoid[id2]].famid, (const char*)ped[phenoid[id2]].pid,
               notMissingCount,
               (het1Count+het2Count-localD[0][ii][jj])*0.5/notMissingCount,
               localD[4][ii][jj]*1.0/notMissingCount,
               kinship);
            if(pbuffer > 0xFFFF){
               fwrite(buffer, 1, pbuffer, fps[thread]);
               pbuffer=0;
            }
         }
      }
   }
   for(int i = 0; i < 5; i++){
      for(int j = 0; j < BLOCKSIZE; j++)
         delete []localD[i][j];
      delete []localD[i];
   }
   if(pbuffer>0)
      fwrite(buffer, 1, pbuffer, fps[thread]);
   fclose(fps[thread]);
#ifdef _OPENMP
}  // extra bracket for omp
#endif
   delete []loopIndex;
   int totalpairCount=0;
#ifdef _OPENMP
   FILE *fp0 = fopen(outfile, "ab");
   for(int c = 1; c < defaultMaxCoreCount; c++){
      FILE *fp = fopen(outfiles[c], "rb");
      int count = fread(buffer, 1, 0x20000, fp);
      for(; count == 0x20000; count = fread(buffer, 1, 0x20000, fp))
         fwrite(buffer, 1, 0x20000, fp0);
      if(count)
         fwrite(buffer, 1, count, fp0);
      fclose(fp);
      remove(outfiles[c]);
   }
   fclose(fp0);
#endif
   delete []missingInOnePersonCount;
   delete []hetInOnePersonCount;
   delete []fps;
   printf("                                         ends at %s", currentTime());
   printf("Between-family kinship data saved in file %s\n", (const char*)outfile);
   printf("Note --kinship --degree <n> can filter & speed up the kinship computing.\n");
   if(xmarkerCount >= MINSNPCOUNT){
      printf("\nX-chromosome analysis...\n");
      if(Bit64==64)
         ComputeLongRobustXKinship64Bit();
      else
         ComputeShortRobustXKinship();
   }
}

void Engine::WritePlink()
{
   if(geno.Length()==0)  // binary genotype does not exist
      BuildShortBinary();
   WritePlinkBinary(prefix);
}

void Engine::printRelationship(int *beforeCount, int *afterCount)
{
   int beforeTotal=0, afterTotal=0;
   for(int i = 0; i < 4; i++){
      if(beforeCount) beforeTotal += beforeCount[i];
      afterTotal += afterCount[i];
   }
   if(beforeTotal || afterTotal){
      printf("\n");
      printf("Relationship summary (total relatives: %d by pedigree, %d by inference)",
         beforeTotal, afterTotal);
      if(beforeCount){
         printf("\n  Source\tMZ\tPO\tFS\t2nd\t3rd\tOTHER\n");
         printf("  ===========================================================\n");
         printf("  Pedigree\t%d\t%d\t%d\t%d\t%d\t%d\n",
         beforeCount[0], beforeCount[5], beforeCount[1]-beforeCount[5], beforeCount[2], beforeCount[3], beforeCount[4]);
         printf("  Inference\t%d\t%d\t%d\t%d\t%d\t%d\n\n",
            afterCount[0], afterCount[5], afterCount[1]-afterCount[5], afterCount[2], afterCount[3], afterCount[4]);
      }else{
         printf("\n        \tMZ\tPO\tFS\t2nd\n");
         printf("  =====================================================\n");
         printf("  Inference\t%d\t%d\t%d\t%d\n\n",
            afterCount[0], afterCount[5], afterCount[1]-afterCount[5], afterCount[2]);
      }
   }
}


