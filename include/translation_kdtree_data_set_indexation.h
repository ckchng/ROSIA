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


#ifndef REG_TRANSLATION_KDT_DATA_INDEXATION_
#define REG_TRANSLATION_KDT_DATA_INDEXATION_

#include "reg_common.h"
#include "data_indexation.h"
#include "state.h"
#include "reg_rtree.h"
#include <nanoflann.hpp>


using namespace nanoflann;

namespace reg {
namespace search {
   

template<class SSR>
class TranslationKDTreeDataSetIndexation : public DataIndexation<SSR>
{
public:

    TranslationKDTreeDataSetIndexation(const Matrix3X &M, const Matrix3X &B, double th);

    ~TranslationKDTreeDataSetIndexation();

    int sweep();
    int size() const;

    int evalUpperBound(SSR ssr, int lwbnd) const;
    int evalUpperBound(SSR ssr, int lwbnd,
                       int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    int evalLowerBound(SSR ssr) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    void setM( Matrix4X *M);

private:

    const Matrix3X &M_in;
    const Matrix3X &B_in;
    const double th;

    Matrix4X *M;
    //pcl::KdTreeFLANN<Point> treeB;
    //flann::Index<flann::L2<double> > treeB;
    
    typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud > ,
    PointCloud, 3  > KdTree;

    
    PointCloud *cloud;
    KdTree *treeB;
    
    int _size;
};



} // End namespace sarch
} // End namespace reg

#include "translation_kdtree_data_set_indexation.hpp"

#endif
