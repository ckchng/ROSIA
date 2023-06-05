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

#include "translation_kdtree_data_set_indexation.h"
#include "util_sort.h"
#include "reg_common.h"
#include "geometry.h"

using namespace reg::search;

template <class SSR>
TranslationKDTreeDataSetIndexation<SSR>::TranslationKDTreeDataSetIndexation(
const reg::Matrix3X &M, const reg::Matrix3X &B, double th ):
    M_in(M), B_in(B), th(th), _size(0)
{

    cloud = new PointCloud(B_in);

    treeB = new KdTree(3, *cloud);

}

template <class SSR>
TranslationKDTreeDataSetIndexation<SSR>::~TranslationKDTreeDataSetIndexation()
{
    delete treeB;
    
}

template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::size() const
{
    return _size;
}

template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::sweep()
{
    treeB->buildIndex();

    _size=1;
   return _size;
}


template <class SSR>
void TranslationKDTreeDataSetIndexation<SSR>::setM(reg::Matrix4X *M)
{
    this->M = M; //.topLeftCorner(3, M_in.cols());
    _size=1;
}


template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::evalUpperBound(
        SSR ssr, int lwbnd) const
{
    int bnd = 0;
    double srad = th+ssrUncertainty(ssr);
    srad *= srad; //nanoflann uses squared distances
    
    std::vector<std::pair<size_t,double>> ret_matches;
    nanoflann::SearchParams params;
    params.sorted = false;
    
    const Vector3 tr = ssrCentre(ssr);
    const int msize = M->cols();
    
    for (int i=0; i<msize && (bnd+msize-i > lwbnd); i++)
    {
        double p[3] = {tr.x+M->x[4*i], tr.y+M->x[4*i+1], tr.z+M->x[4*i+2] };
        //if (treeB.radiusSearch(p, srad, nbIdxB, nbSqrDistB, 1))
        if ( treeB->radiusSearch(p, srad, ret_matches, params) )
        {
            bnd++;
        }
    }
    return bnd;
}

template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::evalUpperBound(
        SSR ssr, int lwbnd,
        int *matchList, int matchListSize,
        std::vector<bool> &matches) const
{
//    mxAssert(lwbnd>=0 && matchListSize>=0 && matchListSize<=M->cols(), "invalid input");

    int bnd=0;
    double srad = th+ssrUncertainty(ssr);
    srad *= srad; //nanoflann uses squared distances
    
    std::vector<std::pair<size_t,double> > ret_matches;
    nanoflann::SearchParams params;
    params.sorted = false;
    
    const Vector3 tr = ssrCentre(ssr);
    
    
    for (size_t i=0; i<matchListSize && (bnd+matchListSize-i > lwbnd); i++)
    {
        //Eigen::Vector3d ep = tr+M.col(matchList[i]);
        const int &col = matchList[i];
        double p[3] = {tr.x+M->x[4*col], tr.y+M->x[4*col+1], tr.z+M->x[4*col+2] };

        //if (treeB.radiusSearch(p, srad, nbIdxB, nbSqrDistB, 1))
        if ( treeB->radiusSearch(p, srad, ret_matches, params) )
        {
            bnd++;
            matches[i]=true;
        }
        else
        {
            matches[i]=false;
        }
    }
    return bnd;
}

template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::evalLowerBound(SSR ssr) const
{
    int qual=0;

    const Vector3 tr = ssrCentre(ssr);
    const size_t msize = M->cols();
    
    std::vector<std::pair<size_t,double> > ret_matches;
    
    nanoflann::SearchParams params;
    params.sorted = false;
    double srad = th*th;
    
    for (size_t i=0; i<msize; i++)
    {
        double p[3] = {tr.x+M->x[4*i], tr.y+M->x[4*i+1], tr.z+M->x[4*i+2] };

      //  if (treeB.radiusSearch(p, th, nbIdxB, nbSqrDistB, 1))
        if(treeB->radiusSearch(p, srad, ret_matches, params))
        {
            qual++;
        }
    }

    //std::cout<<"xxxx lower bound "<<qual<<std::endl;
    return qual;
}

template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::evalLowerBound(
        SSR ssr, int *matchList, int matchListSize) const
{
//    mxAssert(matchListSize>=0 && matchListSize<=M->cols(), "invalid input");
    int qual=0;
    
    nanoflann::SearchParams params;
    params.sorted = false;
    double srad =th*th;
    std::vector<std::pair<size_t,double> > ret_matches;


    const Vector3 tr = ssrCentre(ssr);

    std::vector<std::pair<size_t,double> > indices_dists;
    RadiusResultSet<double,size_t> resultSet(th,indices_dists);

    for (int i=0; i<matchListSize; i++)
    {
        const int &col = matchList[i];
        double p[3] = {tr.x+M->x[4*col], tr.y+M->x[4*col+1], tr.z+M->x[4*col+2] };
        
//        Eigen::Vector3d ep = tr+M.col(matchList[i]);
//        pcl::PointXYZ p(ep(0), ep(1), ep(2));
//        if (treeB.radiusSearch(p, th, nbIdxB, nbSqrDistB, 1))
        if(treeB->radiusSearch(p, srad, ret_matches, params))
        {
            qual++;
        }
    }

    return qual;
}


template <class SSR>
int TranslationKDTreeDataSetIndexation<SSR>::evalLowerBound(
        SSR ssr, int *matchList, int matchListSize,
        std::vector<bool> &matches) const
{
//    mxAssert(matchListSize>=0 && matchListSize<=M->cols(), "invalid parameters");
    int qual=0;
    
    nanoflann::SearchParams params;
    params.sorted = false;
    double srad =th*th;
    std::vector<std::pair<size_t,double> > ret_matches;
    
    const Vector3 tr = ssrCentre(ssr);

    for (size_t i=0; i<matchListSize; i++)
    {
        const int &col = matchList[i];
        double p[3] = {tr.x+M->x[4*col], tr.y+M->x[4*col+1], tr.z+M->x[4*col+2] };
        
//        Eigen::Vector3d ep = tr+M.col(matchList[i]);
//        pcl::PointXYZ p(ep(0), ep(1), ep(2));
//        if (treeB.radiusSearch(p, th, nbIdxB, nbSqrDistB, 1))
        if(treeB->radiusSearch(p, srad, ret_matches, params))
        {
            qual++;
            matches[i]=1;
        }
        else
        {
            matches[i]=0;
        }
    }

    return qual;
}


