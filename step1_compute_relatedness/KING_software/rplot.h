#ifndef __rplot_h__
#define __rplot_h__

#include "IntArray.h"
#include "MathVector.h"
#include "StringHash.h"

void MakeHashKG(StringIntHash & HashKG);
void MakeHashEUR(StringIntHash & HashEUR);
int plotAncestry(const char *prefix, const char *rpath=NULL);
void plotROH(const char *prefix, const char *rpath = NULL);
void plotMIerror(const char *prefix, const char *rpath = NULL);
void plotUniqueFamily(const char *prefix, int degree, const char *analysis, const char *rpath = NULL);
void plotDuplicate(const char *prefix, const char *rpath = NULL);
void plotBuild(const char *prefix, const char *rpath = NULL);
void plotSplitped(const char *prefix, const char *rpath = NULL);
void plotCluster(const char *prefix, const char *rpath = NULL);
void plotGenderError(const char *prefix, IntArray & plotx, Vector & ploty, IntArray & plotz, double xHeterozygosity, int gendererrorCount, const char *rpath = NULL);
void plotRelationship(const char *prefix, const char *rpath = NULL);
void plotIBDSeg(const char *prefix, const char *rpath = NULL);
void plotPopStructure(const char *prefix, int projectFlag, const char *rpath = NULL);

// not released yet
void plotAUCmapping(const char *prefix, int SEXCHR, const char *rpath = NULL);
void plotNPL(const char *prefix, int SEXCHR, const char *rpath = NULL);
void plotHEreg(const char *prefix, int SEXCHR, const char *rpath = NULL);
void plotIBDmapping(const char *prefix, int SEXCHR, const char *rpath = NULL);
void plotROHmapping(const char *prefix, const char *stratName, int SEXCHR, const char *rpath = NULL);
void plotROHforQT(const char *prefix, int SEXCHR, const char *rpath = NULL);
void plotPopROH(const char *prefix, int SEXCHR, const char *rpath = NULL);
void plotPopDist(const char *prefix, const char *rpath = NULL);

#endif
