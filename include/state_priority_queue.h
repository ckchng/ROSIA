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


#ifndef REG_STATE_PRIORITY_QUEUE_H_
#define REG_STATE_PRIORITY_QUEUE_H_

#include "state.h"
#include <cstddef> //NULL

namespace reg {
namespace search {


template <class SSR, typename Scalar=int >
class StatePriorityQueue
{
public:
    enum OptProblem{MINIMISATION, MAXIMISATION};

private:

    class Node
    {
    public:
        SearchState<SSR, Scalar> *state;
        Node *left, *right;

        Node(): state(NULL), left(NULL), right(NULL) {}
        Node(SearchState<SSR, Scalar> *state):state(state), left(NULL), right(NULL){}
        ~Node() {if (state!=NULL) delete state;}
    };

    const OptProblem optProblem;
    Node *head, *tail;
    unsigned int m_size;

public:
    StatePriorityQueue(OptProblem op=MAXIMISATION);
    ~StatePriorityQueue();

    SearchState<SSR, Scalar> *pop();
    void push(SearchState<SSR, Scalar> *state);

    /**
     * @brief Remove and free states with upper bound lower or equal to lwbnd.
     * @param lwbnd Known lower bound.
     */
    void prune(int curbest);

    unsigned int size() const;
};


} // End namespace search
} // End namespace reg

#include "state_priority_queue.hpp"

#endif
