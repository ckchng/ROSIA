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

#ifndef REG_NEAREST_NEIGHBOUR_SEARCH_
#define REG_NEAREST_NEIGHBOUR_SEARCH_

#include "reg_common.h"
#include "reg_dtransf.h"
#include <pcl/kdtree/kdtree_flann.h>


namespace reg
{
namespace search
{

typedef pcl::PointXYZ Point;
typedef pcl::PointCloud<Point>::Ptr Cloud_ptr;


class NearestNeighbourSearch
{
public:
    //virtual double distance(Point p)=0;
    virtual double distance(const Eigen::Vector3f &p)  =0 ;
    virtual double distance(const Point &p)  =0;
    virtual ~NearestNeighbourSearch(){}
};


class KDTreeNearestNeighbourSearch: public NearestNeighbourSearch
{
private:
    std::vector<int> idx;
    std::vector<float> sqrdist;
    pcl::KdTreeFLANN<Point>::Ptr kdtree_ptr;

public:
    KDTreeNearestNeighbourSearch(const Cloud_ptr &B);
    ~KDTreeNearestNeighbourSearch();

    double distance(const Eigen::Vector3f &p) ;
    double distance(const Point &p) ;
};


class DTransfNearestNeighbourSearch: public NearestNeighbourSearch
{
private:
    reg::DTransf dtransf;

public:
    DTransfNearestNeighbourSearch(const Cloud_ptr B,
                                  int resolution, double borderRate=0);
    ~DTransfNearestNeighbourSearch();

    double distance(const Point &p) ;
    double distance(const Eigen::Vector3f &p) ;
};



} // End namespace geometry
} // End namespace reg
#endif
