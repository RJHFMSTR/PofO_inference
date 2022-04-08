#!/bin/bash


# 1. build group of relatives, trios, duos and sibs
Rscript src/step0_grouping.R


# 2. build the final set of target-groups pairs to use for the association testing: Duos and trios target with the direct parents, other target with the surrogate parents.
Rscript src/step1_final_call_set.R


