/////////////////////////////////////////////////////////////////////////////
//
//                          ROSIA
//
// This package contains the source code which implements the
// rotation-search-based star identification algorithm proposed in
//
// @article{chng2022rosia,
//  title={ROSIA: Rotation-Search-Based Star Identification Algorithm},
//  author={Chng, Chee-Kheng and Bustos, Alvaro Parra and McCarthy, Benjamin and Chin, Tat-Jun},
//  journal={arXiv preprint arXiv:2210.00429},
//  year={2022}
//}
//
// Copyright (c) 2023 Chee-Kheng Chng
// Australian Institute for Machine Learning (AIML)
// School of Computer Science, The University of Adelaide, Australia
// Please acknowledge the authors by citing the above paper in any academic
// publications that have made use of this package or part of it.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef REG_BINARYTREE_H
#define REG_BINARYTREE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"


namespace reg
{
namespace binarytree {


/* Definition for interval.*/
typedef struct interval
{
    double lw;
    double up;
} interval;

/* Definitions for binary search tree.*/
typedef struct payload
{
    double val;
    int order;
} payload;

typedef struct treeNode
{
    payload data;
    struct treeNode *left;
    struct treeNode *right;
} treeNode;

treeNode *Insert(treeNode*, payload, treeNode*);
void free_Binarytree(treeNode *node);
int queryLower(treeNode*,double,treeNode*);
int queryUpper(treeNode*,double,treeNode*);
double queryMiddle(treeNode *,double,treeNode *);
void PrintInorder(treeNode*);
int size_Binarytree(treeNode*);
int count_pointers_Binarytree(treeNode*);

} // End namespace binarytree
} // End namespace reg

#endif
