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


#ifndef REG_SEARCH_
#define REG_SEARCH_

#include "reg_common.h"
#include "data_indexation.h"

namespace reg {
    namespace search {
        
        /**
         * @brief Find optimal quality of a transformation
         * @param psi Indexation of point sets.
         * @param knownSolutionQual Known solution quality.
         * @param gap BnB stop gap.
         * @param guessAndResult Guess search region and final region such that
         *        upbnd-lwbnd<=gap
         * @return Quality of the central transform of the optimal region.
         */
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int bnb_search_queue(const DataIndexation<SSR> &psi,
                             int lwbnd, int gap, SSR &guessAndResult);


        template <class SSR, unsigned int BRANCHING_FACTOR>
        int bnb_search_table(const DataIndexation<SSR> &dsi,
                             const Matrix1X &IN_d1, const Matrix1X &IN_d2,
                             const Matrix1X &CAT_d1, const Matrix1X &CAT_d2,
                             const double mag_th,
                             int lwbnd, int gap, int buckets,
                             SSR &guessAndResult);
        
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int searchTableDF(const DataIndexation<SSR> &dsi,
                          int lwbnd, int gap, int buckets,
                          SSR &guessAndResult);
        
        
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int bnb_search_ml_table(const DataIndexation<SSR> &dsi,
                                const Matrix1X &IN_d1, const Matrix1X &IN_d2,
                                const Matrix1X &CAT_d1, const Matrix1X &CAT_d2,
                                const double mag_th,
                                int lwbnd, int gap, int buckets,
                                SSR &guessAndResult);


        
        
    } // End namespace sarch
} // End namespace reg

#include "search.hpp"

#endif

