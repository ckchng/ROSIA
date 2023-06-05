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

#ifndef REG_DATA_INDEXATION_
#define REG_DATA_INDEXATION_

#include <vector>

namespace reg {
namespace search {

template<class SSR>
class DataIndexation
{
public:
    virtual ~DataIndexation(){}
    virtual int evalUpperBound(SSR ssr, int lwbnd) const = 0;
    virtual int evalUpperBound(SSR ssr, int lwbnd,
                               int *matchList, int matchListSize,
                               std::vector<bool> &matches) const = 0;

    virtual int evalLowerBound(SSR ssr) const = 0;
    virtual int evalLowerBound(SSR ssr, int *matchList, int matchListSize) const = 0;
    virtual int evalLowerBound(SSR ssr, int *matchList, int matchListSize,
                               std::vector<bool> &matches) const = 0;

    virtual int sweep() = 0;
    virtual int size() const = 0;

};


inline double dist_sq( Vector3 &a1, double *a2)
{
  double dist_sq = 0, diff;
  for (int i=0; i<3;i++)
  {
    diff = (a1[i] - a2[i]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}


} // End namespace sarch
} // End namespace reg

#endif

