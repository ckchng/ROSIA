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

#include <stdlib.h>
#include <string.h>
#include "reg_common.h"


#ifndef REG_SORT_H
#define REG_SORT_H


namespace reg
{
namespace util
{

//TODO: use a template
struct value_index
{
    double value;
    unsigned int index;
};

struct value_index_uint
{
    unsigned int value;
    unsigned int index;
};


int compare_value_index(const void * a, const void * b);
int compare_value_index_uint(const void * a, const void * b);

unsigned int* sort_index(double *a, unsigned int len);
unsigned int* sort_index(unsigned int *a, unsigned int len);

void sorted_by_index(double *a, unsigned int* idx, unsigned int len);
void sorted_by_index(int *a, unsigned int* idx, unsigned int len);


/*
 * Sort an array according to indexes in idx. An aux array is given.
 */
void sorted_by_index2(double *a, unsigned int* idx, unsigned int len, double *tmp);


} // End namespace util
} // End namespace reg

#endif
