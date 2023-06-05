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

#ifndef REG_MCIRC_INDEXATION_
#define REG_MCIRC_INDEXATION_

#include "reg_common.h"
#include "data_indexation.h"
#include "state.h"
#include "reg_rtree.h"

namespace reg {
namespace search {

template<class SSR>
class MCircIndexation : public DataIndexation<SSR>
{
public:
    MCircIndexation( const Matrix3X &X, const Matrix3X &Y,
                     const reg::Matrix1X &X_mag,
                     const reg::Matrix1X &Y_mag,
                     const reg::Matrix1X &SD1,
                     const reg::Matrix1X &SD2,
                     const reg::Matrix1X &CD1,
                     const reg::Matrix1X &CD2, double dist_th, double mag_th );
    ~MCircIndexation();

    int sweep();
    int sweep(int *matchList, int matchListSize, std::vector<bool> &matches);

    int size() const;

    int evalUpperBound(SSR ssr, int lwbnd) const;

    int evalUpperBound(SSR ssr, int lwbnd,
                       int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    int evalLowerBound(SSR ssr) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize) const;

    int evalLowerBound(SSR ssr, int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;



private:
    const Matrix3X &M_in;
    const Matrix3X &B_in;
    const Matrix1X &M_mag_in;
    const Matrix1X &B_mag_in;
    const Matrix1X &SD1_in;
    const Matrix1X &SD2_in;
    const Matrix1X &CD1_in;
    const Matrix1X &CD2_in;
    const double dist_th;
    const double mag_th;

    Matrix3X M;

    RTree **tree_array;

    int _size;
};

} // End namespace sarch
} // End namespace reg

#include "mcirc_indexation.hpp"


#endif
