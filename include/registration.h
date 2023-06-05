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

#ifndef REG_REG_SEARCH_
#define REG_REG_SEARCH_

#include "reg_common.h"
#include "state.h"

namespace reg {
    namespace search {
        
        //--------------------------------------------------------------------------
        //     Rotation search
        //--------------------------------------------------------------------------
        
        // on points
        int bnb_rsearch_3dof_mcirc(const Matrix3X &X, const Matrix3X &Y,
                                   const Matrix1X &X_mag, const Matrix1X &Y_mag,
                                   const Matrix1X &IN_d1, const Matrix1X &IN_d2,
                                   const Matrix1X &CAT_d1, const Matrix1X &CAT_d2,
                                   const double dist_th, const double mag_th,
                                   int gap,  int lwbnd, int buckets,
                                   AxisAngle &result);
        
        int bnb_rsearch_3dof_mcirc_ml(const Matrix3X &X, const Matrix3X &Y,
                                      const Matrix1X &X_mag, const Matrix1X &Y_mag,
                                      const Matrix1X &IN_d1, const Matrix1X &IN_d2,
                                      const Matrix1X &CAT_d1, const Matrix1X &CAT_d2,
                                      const double dist_th, const double mag_th,
                                      int gap,  int lwbnd, int buckets,
                                      AxisAngle &result);


        
        // on matches
        
        // ang distance
        int rot3_matches(const Matrix3X &X, const Matrix3X &Y, double th,
                         int gap, int lwbnd, AxisAngle &result);
        
        int rot3_matches_ml(const Matrix3X &X, const Matrix3X &Y, double th,
                            int gap, int lwbnd, AxisAngle &result);
        
        
        // l2 distance
        int rot3_matches_l2(const Matrix3X &X, const Matrix3X &Y, double th,
                            int lwbnd, int gap, AxisAngle &result);
        
        int rot3_matches_l2_ml(const Matrix3X &X, const Matrix3X &Y, double th,
                               int lwbnd,int gap, AxisAngle &result);
        
        
        int rot1_matches_l2(const Matrix3X &X, const Matrix3X &Y, double th,
                            int lwbnd, int gap, double &result);
        
        int rot1_matches_l2_ml(const Matrix3X &X, const Matrix3X &Y, double th,
                               int lwbnd, int gap, double &result);
        
        
        
        //--------------------------------------------------------------------------
        //     Registration
        //--------------------------------------------------------------------------
        
        
        // 6 DoF with matches
        int reg6Matches(const Matrix3X &X, const Matrix3X &Y, double th,
                        int lwbnd, int gap, Transform3 &guessAndResult);
        
        int reg4Matches(const Matrix3X &X, const Matrix3X &Y, double th,
                        int lwbnd, int gap, Transform3 &guessAndResult);
        
        
        
        // Nested (6DoF)
        
        int nestedbnb_search_6dof_mcirc(const Matrix3X &X, const Matrix3X &Y, double th,
                                        int gap, int lwbnd,
                                        TranslationSearchSpaceRegion3DOF &trDom,
                                        int(*bnb_rsearch_3dof)(const Matrix3X&, const Matrix3X&,
                                                               double,int, int, int, AxisAngle&),
                                        Transform3 &guessAndResult,
                                        bool use_local_opt);
        
        int nestedbnb_search_6dof_mcirc_ml(const Matrix3X &X, const Matrix3X &Y, double th,
                                           int gap, int lwbnd,
                                           TranslationSearchSpaceRegion3DOF &trDom,
                                           Transform3 &guessAndResult,
                                           bool use_local_opt);
        
        
        // Local (6DoF)
        int local_search_6dof(const Matrix3X &X, const Matrix3X &Y, double th,
                              int lwbnd,
                              TranslationSearchSpaceRegion3DOF trDom,
                              int(*bnb_rsearch_3dof)(const Matrix3X&, const Matrix3X&,
                                                     double,int, int, int, AxisAngle&),
                              Transform3 &guessAndResult);
        
        
    } // End namespace sarch
} // End namespace reg

#endif

