//////////////////////////////////////////////////////////////////////
// FASTASSOC.cpp
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
// July 1, 2019


#include "analysis.h"
#include <math.h>
#include "MathStats.h"
#include "QuickIndex.h"
#include "MathCholesky.h"

void Engine::FASTASSOC()
{
   if(detailFlag) countGenotype();
   if(geno.Length()==0) {
      individualInfo = true;
      if(shortFlag) BuildShortBinary();
      else BuildBinary();
   }
   PreVC();

   int N_Trait = traits.Length();
   int N_Covariate = covariates.Length();
   int N_ID = ID.Length();
   Matrix Y(N_ID, N_Trait);
   Matrix X(N_ID, 1+N_Covariate);
   for(int i = 0; i < N_ID; i++)
      X[i][0] = 1.0;
   for(int t = 0; t < N_Trait; t++)
      for(int i = 0; i < N_ID; i++)
         Y[i][t] = ped[ID[i]].traits[traits[t]];
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_Covariate; j++)
         X[i][j+1] = ped[ID[i]].covariates[covariates[j]];

   Vector meanY(N_Trait);
   meanY.Set(0.0);
   int n;
   for(int t = 0; t < N_Trait; t++){
      n = 0;
      for(int i = 0; i < N_ID; i++)
         if(Y[i][t] != _NAN_) {
            meanY[t] += Y[i][t];
            n ++;
         }
      meanY[t] /= n;
      for(int i = 0; i < N_ID; i++)
         if(Y[i][t] == _NAN_)
            Y[i][t] = meanY[t];
   }

   Matrix UYs(N_ID, N_Trait); // U'Y += UT[i][k] * Y[k]
   UYs.Zero();
   UX.Dimension(N_ID, N_Covariate+1);
   UX.Zero();
   Vector tV(N_ID);

   Cholesky chol;
   QuickIndex idx;
   for(int t = 0; t < N_Trait; t++){
      if(normalization){
         for(int i = 0; i < N_ID; i++)
            tV[i] = Y[i][t];
         idx.Index(tV);
         for(int i = 0; i < N_ID;){
            int start = i, end = i + 1;
            while(end < N_ID && tV[idx[end]] == tV[idx[start]]) end++;
            end --;
            double q = ninv((start + (end-start)/2.0 + 0.5)/N_ID);
            for(int j = start; j <= end; j++)
               Y[idx[j]][t] = q;
            i = end + 1;
         }
      }
      for(int i = 0; i < N_ID; i++)
         for(int k = 0; k < N_ID; k++)
            UYs[i][t] += UT[i][k] * Y[k][t];
   }
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_Covariate+1; j++)
         for(int k = 0; k < N_ID; k++)
            UX[i][j] += UT[i][k] * X[k][j];

   UY.Dimension(N_ID);
   double a, b, LB, T, x, temp, funcx;
   int numiter;
   int maxIter=10000;
   bool quiet=true;
   T = 0.001;
   double R1R;

   Vector lambda0(N_Trait);
   lambda0.Zero();
   Matrix score0(N_Trait, N_Covariate+1);
   score0.Zero();
   Matrix *Info0 = new Matrix[N_Trait];
   printf("\nPolygenic parameter estimates\n");
   printf("%-15s %7s %7s %7s %9s",
      "TraitName", "Herit", "Lambda", "Tau", "Mu");
   for(int i = 0; i < N_Covariate; i++)
      printf(" %9s", (const char*)ped.covariateNames[covariates[i]]);
   printf("\n");

   Vector beta0(N_Covariate+1);
   Matrix Resid0(N_ID, N_Trait);
   Vector tau0(N_Trait);
   double loglik0;

   for(int t = 0; t < N_Trait; t++){
      currentT = t;
      for(int i = 0; i < N_ID; i++)
         UY[i] = UYs[i][currentT];
      loglik0 = 1E10;
      for(int i = 0; i < 100; i++){
         a = exp(i*0.1 - 5);
         b = exp(i*0.1 - 4.9);
         x = minimize(a, b, T, funcx, numiter, maxIter, quiet);
         if(funcx < loglik0) {
            loglik0 = funcx;
            lambda0[currentT] = 1/x;
         }
      }
      if(lambda0[currentT]<0.0067){  // border
         lambda0[currentT] = 0;
      }
      Info0[currentT].Dimension(N_Covariate+1, N_Covariate+1);
      Info0[currentT].Zero();
      for(int i = 0; i < N_Covariate+1; i++)
         for(int j = 0; j < N_Covariate+1; j++)
            for(int k = 0; k < N_ID; k++)
               Info0[currentT][i][j] += UX[k][i] * UX[k][j] / (lambda0[currentT] * EV[k] + 1);
      chol.Decompose(Info0[currentT]);
      for(int i = 0; i < N_Covariate+1; i++)
         for(int k = 0; k < N_ID; k++)
            score0[currentT][i] += UX[k][i] * UY[k] / (lambda0[currentT] * EV[k] + 1);
      chol.BackSubst(score0[currentT]);
      beta0 = chol.x;
      R1R = 0;
      for(int i = 0; i < N_ID; i++){
         Resid0[i][currentT] = UY[i];
         for(int j = 0; j < N_Covariate+1; j++)
            Resid0[i][currentT] -= UX[i][j] * beta0[j];
         R1R += Resid0[i][currentT] * Resid0[i][currentT] / (lambda0[currentT] * EV[i]+1);
      }
      tau0[currentT] = N_ID / R1R;
      printf("%-15s %7.4lf %7.4lf %7.3lf",
         (const char*)ped.traitNames[traits[currentT]], lambda0[currentT]/(1+lambda0[currentT]),
         lambda0[currentT], tau0[currentT]);
      for(int i = 0; i < N_Covariate+1; i++)
         printf(" %9.3lf", beta0[i]);
      printf("\n");
   }
   if(HeritFlag) {
      if(Info0) delete []Info0;
      for(int i = 0; i < N_ID; i++)
         delete UT[i];
      delete []UT;
      return;
   }
   printf("Variance component genome scan starts at %s", currentTime());

   int id, AA, Aa, missing, lowb, highb;
   double pvalue, h2;
   double **UG;
   UG = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UG[i] = new double [16];

   int MaskArray[256][9];
   for(int i = 0; i < 256; i++){
      MaskArray[i][0] = 0;
      for(int k = 0; k < 8; k++)
         if(i & (1<<k)){
            MaskArray[i][0]++;
            MaskArray[i][MaskArray[i][0]] = k;
         }
   }
   double **UT2 = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT2[i] = new double [N_ID];
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         UT2[i][j] = UT[i][j] * 2;

   Vector U1(N_ID);
   U1.Zero();
   for(int i = 0; i < N_ID; i++)
      for(int j = 0; j < N_ID; j++)
         U1[i] += UT[i][j];

   double freq[16];
   int AACount[16], AaCount[16], missingCount[16];
   int pos;
   String lmmfile = prefix;
   lmmfile.Add("score.txt");
   FILE *fp = fopen(lmmfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)lmmfile);
   fprintf(fp, "SNP");
   if(N_Trait > 1)
      fprintf(fp, "\tTraitName");
   if(chromosomes.Length()) fprintf(fp, "\tChr");
   if(bp.Length()) fprintf(fp, "\tPos");
   fprintf(fp, "\tLabelA\tLabela\tFreqA\tN\tBeta\tSE\tChisq\tH2\tPvalue\n");

   int k, mm;
   unsigned short int word;
   char revbase[65536];
   for(int i = 0; i < 16; i++)
      revbase[shortbase[i]] = i;
   char rightmost[65536];
   for(int i = 0; i < 65536; i++)
      rightmost[i] = revbase[i&(-i)];
   
   double statistic, score_local, var_local;
   for(int b = 0; b < shortCount; b++){ // loop over short integer
      for(int i = 0; i < 16; i++){
         AACount[i] = AaCount[i] = missingCount[i] = 0;
         freq[i] = 0.0;
      }     
      for(int i = 0; i < N_ID; i++){
         id = geno[ID[i]];
         AA = GG[0][id][b] & GG[1][id][b];
         for(word = AA; word; word &= (word-1))
            AACount[rightmost[word]] ++;
         Aa = (~GG[0][id][b]) & GG[1][id][b];
         for(word = Aa; word; word &= (word-1))
            AaCount[rightmost[word]] ++;
         missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
         for(word = missing; word; word &= (word-1))
            missingCount[rightmost[word]] ++;
      }
      for(int j = 0; j < 16; j++)
         if(missingCount[j] < N_ID)  // frequency of allele A
            freq[j] = (AACount[j] + AaCount[j]*0.5) / (N_ID - missingCount[j]);

      for(int i = 0; i < N_ID; i++)
         for(int m = 0; m < 16; m++)
            UG[i][m] = 0;
      for(int i = 0; i < N_ID; i++){ // O(N^2): most time-consuming loop in genome scan
         id = geno[ID[i]];
         AA = GG[0][id][b] & GG[1][id][b];
         Aa = (~GG[0][id][b]) & GG[1][id][b];
         missing = (~GG[0][id][b]) & (~GG[1][id][b]) & 65535;
         if(Aa & 255){
            lowb = Aa & 255;
            for(int m = 0; m < MaskArray[lowb][0]; m++){
               mm = MaskArray[lowb][m+1];
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT[k][i];
            }
         }
         if(Aa & 65280){
            highb = Aa >> 8;
            for(int m = 0; m < MaskArray[highb][0]; m++){
               mm = MaskArray[highb][m+1]+8;
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT[k][i];
            }
         }
         if(AA & 255){
            lowb = AA & 255;
            for(int m = 0; m < MaskArray[lowb][0]; m++){
               mm = MaskArray[lowb][m+1];
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i];
            }
         }
         if(AA & 65280){
            highb = AA >> 8;
            for(int m = 0; m < MaskArray[highb][0]; m++){
               mm = MaskArray[highb][m+1]+8;
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i];
            }
         }

         if(missing & 255){
            lowb = missing & 255;
            for(int m = 0; m < MaskArray[lowb][0]; m++){
               mm = MaskArray[lowb][m+1];
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i] * freq[mm];
            }
         }
         if(missing & 65280){
            highb = missing >> 8;
            for(int m = 0; m < MaskArray[highb][0]; m++){
               mm = MaskArray[highb][m+1]+8;
               for(k = 0; k < N_ID; k++)
                  UG[k][mm] += UT2[k][i] * freq[mm];
            }
         }
      }

      double temp, tempG;
      for(int m = 0; m < 16; m++){
         pos = b*16 + m;
         if(pos >= markerCount) break;
         // Average freq is used to avoid EV=0
         for(int t = 0; t < N_Trait; t++){
            statistic = 0;
            score_local = var_local = 0;
            for(int i = 0; i < N_ID; i++){
               temp = lambda0[t] * EV[i] + 1;
               tempG = UG[i][m] - 2*freq[m]*U1.data[i];
               score_local += tempG * Resid0.data[i]->data[t]/ temp;
               var_local += tempG * tempG / temp;
            }
            if(var_local > 0){
               statistic = tau0[t] * score_local * score_local / var_local;
               pvalue = ndist(sqrt(statistic))*2;
               chisqs[t].Push(statistic);
               h2 = score_local * score_local / var_local / var_local * 2 * freq[m] * (1-freq[m]);
               h2 = h2/((1+lambda0[t])/tau0[t]+h2);
            }else{
               statistic = _NAN_;
               pvalue = 1;
            }
            if(snpName.Length())
               fprintf(fp, "%s", (const char*)snpName[pos]);
            else
               fprintf(fp, "SNP%d", pos+1);
            if(N_Trait>1)
               fprintf(fp, "\t%s", (const char*)ped.traitNames[t]);
            if(chromosomes.Length())
               fprintf(fp, "\t%d", chromosomes[pos]);
            if(bp.Length())
               fprintf(fp, "\t%.6lf", bp[pos]*0.000001);
            fprintf(fp, "\t%s\t%s\t%.3lf\t%d",
               (const char*)alleleLabel[0][pos], (const char*)alleleLabel[1][pos],
               freq[m], N_ID-missingCount[m]);
            if(statistic != _NAN_)
               fprintf(fp, "\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3G\n",
                  score_local / var_local, sqrt(1/var_local/tau0[t]), statistic, h2, pvalue);
            else
               fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\n", "NA", "NA", "NA", "NA", "NA");
         }
      }
   }
   fclose(fp);
   printf("Variance component genome scan ends at %s", currentTime());
   printf("Variance component genome scan results saved in file %s\n\n", (const char*)lmmfile);

   PostVC();
   if(Info0) delete []Info0;
   for(int i = 0; i < N_ID; i++)
      delete UT[i];
   delete []UT;
   for(int i = 0; i < N_ID; i++)
      delete UG[i];
   delete []UG;
   for(int i = 0; i < N_ID; i++)
      delete UT2[i];
   delete []UT2;
}

         /*
         // est allele freq
         temp = 0;
         freq[m] = 0;
         for(int i = 0; i < N_ID; i++){
            freq[m] += U1[i] * UG[i][m] / EV[i];
            temp += U1[i] * U1[i] / EV[i];
         }
         freq[m] /= (temp*2);
           */

   /*
   if(N_Trait == 0) error("Quantitative traits are not found for variance component analysis.");
   IntArray validFlag(ped.count);
   validFlag.Set(1);
   for(int i = 0; i < ped.count; i++){
      bool allmissing = true;
      for(int t = 0; t < N_Trait; t++)
         if(ped[i].isPhenotyped(traits[t])) {
            allmissing = false;
            break;
         }
      if(allmissing) validFlag[i] = 0;
      bool somemissing = false;
      for(int j = 0; j < N_Covariate; j++)
         if(ped[i].covariates[covariates[j]] == _NAN_){
            somemissing = true;
            break;
         }
      if(somemissing) validFlag[i] = 0;
      if(geno[i]<0 || ped[i].ngeno < MINSNPCOUNT) validFlag[i] = 0;
   }
   ID.Dimension(0);
   for(int i = 0; i < ped.count; i++)
      if(validFlag[i]) ID.Push(i);

   printf("Genotypes stored in %d integers for each of %d individuals used in analysis.\n",
      shortCount, N_ID);
   EV.Dimension(N_ID);
   UT = new double * [N_ID];
   for(int i = 0; i < N_ID; i++)
      UT[i] = new double [N_ID];

   if(svdinfile!="")
      ReadSVD();
   else
      ComputeSVD();
   if(normalization)
      printf("Inverse normal transformation is applied to phenotypes.\n");
   */
   /*
   if(N_Trait>1){
      printf("\nTraitName\tN_SNP\tSmallestP\tGC_lambda\n");
      for(int t = 0; t < N_Trait; t++){
         idx.Index(chisq[t]);
         printf("%s\t%d\t%.3G\t%.3lf\n",
         (const char*)ped.traitNames[traits[t]],
         chisq[t].Length(),
         chidist(chisq[t].Max(),1),
         chisq[t][idx[chisq[t].Length()/2]]/.456);
      }
   }else{
      idx.Index(chisq[0]);
      printf("\nGC lambda for the genome scan of trait %s in %d autosome SNPS is %.3lf\n",
         (const char*)ped.traitNames[traits[0]],
         chisq[0].Length(), chisq[0][idx[chisq[0].Length()/2]]/.456);
   }
   printf("\n");
   for(int i = 0; i < N_ID; i++)
      delete UT[i];
   delete []UT;
   */

