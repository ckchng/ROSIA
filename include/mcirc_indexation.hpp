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

#include "mcirc_indexation.h"
#include "util_sort.h"
#include "reg_binary_tree.h"
#include "reg_common.h"
#include "geometry.h"

using namespace reg::search;

template <class SSR>
MCircIndexation<SSR>::MCircIndexation(const reg::Matrix3X &X,
                                      const reg::Matrix3X &Y,
                                      const reg::Matrix1X &X_mag,
                                      const reg::Matrix1X &Y_mag,
                                      const reg::Matrix1X &SD1,
                                      const reg::Matrix1X &SD2,
                                      const reg::Matrix1X &CD1,
                                      const reg::Matrix1X &CD2,
                                      double dist_th,
                                      double mag_th):
    M_in(X), B_in(Y), M_mag_in(X_mag), B_mag_in(Y_mag), SD1_in(SD1), SD2_in(SD2), CD1_in(CD1), CD2_in(CD2), dist_th(dist_th), mag_th(mag_th), _size(0)
{}

template <class SSR>
MCircIndexation<SSR>::~MCircIndexation()
{
    if(_size>0)
    {
        for (size_t i=0; i<_size; i++)
        {
            delete tree_array[i];
        }
        delete []tree_array;
    }
}

template <class SSR>
int MCircIndexation<SSR>::size() const
{
    return _size;
}

template <class SSR>
int MCircIndexation<SSR>::sweep()
{
    int i, j;
    int lw, up;

    //-----------------------------
    //        Sort
    //-----------------------------
   
    RTree **rtreeArrayTmp;

    const size_t msize = M_in.n;
    const size_t bsize = B_in.n;


    rtreeArrayTmp = new RTree*[msize]();

    for(i=0; i<msize; i++)
    {
        rtreeArrayTmp[i] = NULL;
    }
    
    Vector M_len(msize);

    for(i=0; i<msize; i++)
    {
        M_len(i) = M_mag_in(i);

    }



    //TODO: Move sort to matrix class
    unsigned int *minx;
    minx = reg::util::sort_index(M_len.x, msize);
    Matrix3X M_tmp = M_in; //copy
    M_tmp.sort(minx);
    free(minx);

    //-----------------------------
    //        Sweep
    //-----------------------------

    RTree *rtree_i;
    binarytree::payload newdata;
    binarytree::treeNode *status;
    
   
    status = NULL;
    for (i=0;i<msize;i++)
    {
        // Create plane sweep status.
        newdata.val = M_len(i);
        newdata.order = i;
        status = Insert(status, newdata, NULL);
    }

    _size=0;
    double blen_j, mlen_i, patchAngle;
    

    for (j=0;j<bsize;j++)
    {
//        blen_j = norm(B_in.x+3*j);//  B_len[j];

        blen_j = B_mag_in(j);

        lw = queryLower(status, blen_j-mag_th, NULL);
        up = queryUpper(status, blen_j+mag_th, NULL);

//        Vector3 patchCentre(B_in(0,j)/blen_j, B_in(1,j)/blen_j, B_in(2,j)/blen_j);
        Vector3 patchCentre(B_in(0,j), B_in(1,j), B_in(2,j));

        
        for (i=lw; i<=up; i++)
        {
            mlen_i = M_len(i);
            rtree_i = rtreeArrayTmp[i];

            // check d1 and d2 thres here

            if (abs(CD1_in(j) - SD1_in(i)) > 2 * dist_th){
                continue;
            }

            if (abs(CD2_in(j) - SD2_in(i)) > 2 * dist_th){
                continue;
            }


            if( rtree_i == NULL )
            {
                rtree_i = rtreeArrayTmp[i] = new RTree();
                _size++;
            }

            //Any rotation on M[i] will intersect B[j]
            if(blen_j<=EPSILON || mlen_i<=DUMMY_PRECISION)
            {
                rtree_i->setFull();
            }
            else
            {
                // patch angle
                patchAngle = dist_th;
                if(patchAngle<0) //
                {
                    rtree_i->setFull();
                }
                else
                {
//                    mxAssert(fabs(patchCentre.norm()-1.0)<=DUMMY_PRECISION,"");
                    rtree_i->addPatch(patchCentre, patchAngle);
                }
            }
        }
    }

    free_Binarytree(status);

    if(_size>0)
    {

        M.setSize(_size); //M = Eigen::Matrix3Xd::Zero(3, _size);
        tree_array = new RTree*[_size];
        j=0;
        for (i=0; i<msize; i++)
        {
            if (rtreeArrayTmp[i]!=NULL)
            {
                tree_array[j] = rtreeArrayTmp[i];
                M(0,j)=M_tmp(0,i);
                M(1,j)=M_tmp(1,i);
                M(2,j)=M_tmp(2,i);
//                M(0,j)=M_tmp(0,i)/M_len(i);
//                M(1,j)=M_tmp(1,i)/M_len(i);
//                M(2,j)=M_tmp(2,i)/M_len(i);
                j++;
            }
        }
    }
    
    delete []rtreeArrayTmp;
    return _size;
}


template <class SSR>
int MCircIndexation<SSR>::sweep(
        int *matchList, int matchListSize, std::vector<bool> &matches)
{
//    mxAssert(matches.size()>=matchListSize,"");
//    mxAssert(matchListSize>0 && matchListSize<=M_in.n,"");

    int i, j, k;
    int lw, up;

    //-----------------------------
    //        Sort
    //-----------------------------

    // sorting temp. vars.
    double *stmp;

    RTree **rtreeArrayTmp;

    const int bsize = B_in.n ;//B_in.cols();
//    mxAssert(bsize>0,"size of B must be >0");

    rtreeArrayTmp = new RTree*[matchListSize]();

    double *MX_tmp = (double *) malloc(matchListSize*sizeof(double));
    double *MY_tmp = (double *) malloc(matchListSize*sizeof(double));
    double *MZ_tmp = (double *) malloc(matchListSize*sizeof(double));

    //Vector M_len = Vector(matchListSize);
    double *M_len = (double *) malloc(matchListSize*sizeof(double));


    j=0;
    for(i=0; i<matchListSize; i++)
    {
        //Eigen::Vector3d m = M_in.col(matchList[i]).head(3);
        double *m = M_in.x + 3*matchList[i]; //.col(matchList[i])

        MX_tmp[j] = m[0]; //m(0);
        MY_tmp[j] = m[1]; //m(1);
        MZ_tmp[j] = m[2]; //m(2);
        M_len[j] = norm(m); //m.norm();
        j++;
    }
//    mxAssert(j==matchListSize,"Wrong match list size");

    for(i=0; i<matchListSize; i++)
    {
        rtreeArrayTmp[i] = NULL;
    }

    unsigned int *midx = reg::util::sort_index(M_len, matchListSize);
    stmp = (double *)malloc(matchListSize*sizeof(double));

    reg::util::sorted_by_index2(MX_tmp, midx, matchListSize, stmp);
    reg::util::sorted_by_index2(MY_tmp, midx, matchListSize, stmp);
    reg::util::sorted_by_index2(MZ_tmp, midx, matchListSize, stmp);
    free(stmp);

    //-----------------------------
    //        Sweep
    //-----------------------------

    RTree *rtree_i;
    binarytree::payload newdata;
    binarytree::treeNode *status;

    status = NULL;
    for (i=0; i<matchListSize; i++)
    {
        // Create plane sweep status.
        newdata.val = M_len[i];
        newdata.order = i;
        status = Insert(status, newdata, NULL);
    }

    _size=0;
    double blen_j, mlen_i, patchAngle;
    for (j=0;j<bsize;j++)
    {
        blen_j = norm(B_in.x+3*j);//B_len[j];

        lw = queryLower(status, blen_j-mag_th, NULL);
        up = queryUpper(status, blen_j+mag_th, NULL);

        Vector3 patchCentre(
                    B_in(0,j)/blen_j,
                    B_in(1,j)/blen_j,
                    B_in(2,j)/blen_j);

        for (i=lw; i<=up; i++)
        {
            mlen_i = M_len[i];
            rtree_i = rtreeArrayTmp[i];
            if( rtree_i == NULL )
            {
                rtree_i = rtreeArrayTmp[i] = new RTree();
                _size++;
            }

            //Any rotation on M[i] will intersect B[j]
            if(blen_j<=DUMMY_PRECISION || mlen_i<=DUMMY_PRECISION)
            {
                rtree_i->setFull();
            }
            else
            {
                // patch angle
//                patchAngle = circleIntersectionAngle(mlen_i, blen_j, dist_th);  /* R, d, r */
                patchAngle = circleIntersectionAngle(mlen_i, mlen_i, dist_th);  /* R, d, r */
                if(patchAngle<0) //
                {
                    rtree_i->setFull();
                }
                else
                {
//                    mxAssert(fabs(patchCentre.norm()-1.0)<=DUMMY_PRECISION,
//                             "patch is not normalised");
                    rtree_i->addPatch(patchCentre, patchAngle);
                }
            }
        } //for sweep
    }

    free_Binarytree(status);

    unsigned int *orig_order_idx;
    orig_order_idx = reg::util::sort_index(midx, matchListSize);
    free(midx);

//    mxAssert(_size>=0 && _size<=matchListSize, "inconsistent size");
    if(_size>0)
    {
        //M = Eigen::Matrix3Xd::Zero(3, _size);
        M.setSize(_size); //M = Eigen::Matrix3Xd::Zero(3, _size);

        tree_array = new RTree*[_size];
        j=0;

        for (i=0; i<matchListSize; i++)
        {
            // Iterate in original order
            k = orig_order_idx[i]; //k is the idx sorted as original input
//            mxAssert(k>=0 && k <matchListSize, "inconsistent index");

            if (rtreeArrayTmp[k]!=NULL)
            {
                tree_array[j] = rtreeArrayTmp[k];
                M(0,j)=MX_tmp[k]/M_len[k];
                M(1,j)=MY_tmp[k]/M_len[k];
                M(2,j)=MZ_tmp[k]/M_len[k];
                j++;

                matches[i]=true;
            }
            else
            {
                matches[i]=false;
            }
        }
//        mxAssert(j==_size, "inconsistent index");
    }
    else //_size ==0
    {
        for (i=0; i<matchListSize; i++)
        {
            matches[i]=0;
        }
    }

    free(M_len);
    free(MX_tmp);
    free(MY_tmp);
    free(MZ_tmp);
    free(orig_order_idx);
    delete []rtreeArrayTmp;

    return _size;
}

template <class SSR>
int MCircIndexation<SSR>::evalUpperBound(SSR ssr, int lwbnd) const
{
//    mxAssert(lwbnd>=0, "inconsistent lower bound");

    size_t bnd=0;
    const double max_angle = ssrMaxAngle(ssr);

     Matrix3 R ;
     fromAxisAngle(R,centre(ssr));
    
    for (size_t i=0; i<_size && (bnd+_size-i > lwbnd); i++)
    {
        bnd += tree_array[i]->matchPatch(multiply(R,M.x+3*i), max_angle);
    }
    
    return bnd;
}



template <class SSR>
int MCircIndexation<SSR>::evalUpperBound(
        SSR ssr, int lwbnd,
        int *matchList, int matchListSize,
        std::vector<bool> &matches) const
{
//    mxAssert(matchListSize>=0 && matchListSize<=_size, "inconsistent matchlist size");

    size_t bnd=0;
    double max_angle = ssrMaxAngle(ssr);
    //const Eigen::Matrix3d R = centre(ssr);

    Matrix3 R;
    fromAxisAngle(R,centre(ssr));

    for (size_t i=0; i<matchListSize && (bnd+matchListSize-i > lwbnd); i++)
    {
        if(tree_array[matchList[i]]->matchPatch(multiply(R,(M.x+3*matchList[i])), max_angle))
        {
            bnd++;
            matches[i]=1;
        }
        else
        {
            matches[i]=0;
        }
    }

    return bnd;
}

template <class SSR>
int MCircIndexation<SSR>::evalLowerBound(SSR ssr) const
{
    size_t qual=0;
    //Matrix3 R = fromAxisAngle(centre(ssr));
    
    Matrix3 R;
    AxisAngle a = centre(ssr);
    fromAxisAngle(R, a);

    double *RM = (double *)malloc(3*_size*sizeof(double));


    matrixMultipy(R.m, M.x, RM, _size);

    double *rm;
    for (size_t i=0; i<_size; i++)
    {
        //qual += tree_array[i]->matchPoint(RM.col(i));
        rm = RM+3*i;
        qual += tree_array[i]->matchPoint(Vector3(rm[0], rm[1], rm[2]));

    }

    free(RM);

    return qual;
}

template <class SSR>
int MCircIndexation<SSR>::evalLowerBound(
        SSR ssr, int *matchList, int matchlistSize) const
{
//    mxAssert(matchlistSize>=0 && matchlistSize<=_size, "inconsistent matchlist size");

    size_t qual=0;
    Matrix3 R;
    fromAxisAngle(R,centre(ssr));

    for (size_t i=0; i<matchlistSize; i++)
    {
        qual += tree_array[matchList[i]]->matchPoint(multiply(R,M.x+3*matchList[i]));
    }

    return qual;
}


template <class SSR>
int MCircIndexation<SSR>::evalLowerBound(
        SSR ssr, int *matchList, int matchlistSize,
        std::vector<bool> &matches) const
{
//    mxAssert(matchlistSize>=0 && matchlistSize<=_size, "inconsistent matchlist size");

    size_t qual=0;
    Matrix3 R;
    fromAxisAngle(R,centre(ssr));

    for (size_t i=0; i<matchlistSize; i++)
    {
        if( tree_array[matchList[i]]->matchPoint(
                    multiply(R,M.x+3*matchList[i]) ) )
        {
            qual++;
            matches[i]=true;
        }
        else
        {
            matches[i]=false;
        }
    }

//    mxAssert(matches.size()>=qual, "inconsistent qual");

    return qual;
}
