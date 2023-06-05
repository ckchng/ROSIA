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

#include "state_priority_queue.h"
#include "state_priority_hashtable.h"
#include "data_indexation.h"
#include "geometry.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sys/resource.h>

namespace reg
{
namespace search
{

template <class SSR, unsigned int BRANCHING_FACTOR>
int bnb_search_queue(const DataIndexation<SSR> &dsi,
                     int lwbnd, int gap,
                     SSR &guessAndResult)
{
    int i, np;
    int state_lwbnd, upbnd;

    int count_eval_ubbnd;
    int count_eval_lwbnd;

    SearchState<SSR> *state;
    SSR **ssr_array; //array of pointers

    // Initialise
    count_eval_ubbnd = count_eval_lwbnd = 0;

    upbnd = dsi.evalUpperBound(guessAndResult, lwbnd);
    count_eval_ubbnd++;

    // No better solution than the known one in provided SSR
//    mxAssert( (((lwbnd==0) &&(upbnd>=lwbnd)) || (lwbnd>0)), "Bound error");
    if (upbnd-lwbnd <= gap)
    {
        return lwbnd;
    }

    StatePriorityQueue<SSR> queue;

    
    ssr_array = new SSR*[BRANCHING_FACTOR];
    state = new SearchState<SSR>(guessAndResult, upbnd);
    queue.push(state);

    int iter = 0;
    while (queue.size())
    {
        iter++;

        // Find the state with the highest upper bound
        state = queue.pop();

        // Evaluate lower boud
        state_lwbnd = dsi.evalLowerBound( state->ssr );
        count_eval_lwbnd++;

        // Update solution
        if (state_lwbnd > lwbnd)
        {
            lwbnd = state_lwbnd;
            guessAndResult = state->ssr;
            queue.prune(lwbnd);
        }

        // Stopping criterion
//        mxAssert( state->bnd >= lwbnd,  "Bound error");
        if (state->bnd - lwbnd <= gap)
        {
            delete state;
            break;
        }

        // Branch
        np = reg::search::split(state->ssr, ssr_array);
        delete state;

        for(i=0; i<np; i++)
        {
            upbnd = dsi.evalUpperBound(*(ssr_array[i]), lwbnd);
            if ( upbnd > lwbnd )
            {
                state = new SearchState<SSR>(*(ssr_array[i]), upbnd);
                queue.push(state);
            }
            delete ssr_array[i];
        }
        count_eval_ubbnd += np;
    }

    delete []ssr_array;
    return lwbnd;
}



template <class SSR, unsigned int BRANCHING_FACTOR>
int bnb_search_table(const DataIndexation<SSR> &dsi,
                     const Matrix1X &IN_d1, const Matrix1X &IN_d2,
                     const Matrix1X &CAT_d1, const Matrix1X &CAT_d2,
                     const double mag_th,
                     int lwbnd, int gap, int buckets,
                     SSR &guessAndResult)
{
    int i, np;
    int within_fov_lwbnd, missing_star_lwbnd, within_fov_upbnd, missing_star_upbnd, state_lwbnd, upbnd, matched_star_lwbnd, matched_star_upbnd;

    int count_eval_ubbnd;
    int count_eval_lwbnd;
    
    SearchState<SSR, int> *state;
    SSR **ssr_array; //array of pointers
   
    // Initialise
    count_eval_ubbnd = count_eval_lwbnd = 0;


    upbnd = dsi.evalUpperBound(guessAndResult, lwbnd);

    std::cout << "ubnd" <<upbnd << std::endl;

    count_eval_ubbnd++;

    if (upbnd-lwbnd <= gap)
    {
        return lwbnd;
    }
    
    StatePriorityHashtable<SSR, int, SearchState > table(buckets);

    
    ssr_array = new SSR*[BRANCHING_FACTOR];
    state = new SearchState<SSR>(guessAndResult, upbnd);
    table.push(state);

    int iter = 0;
    
    while (table.size())
    {
        iter++;

        // Find the state with the highest upper bound
        state = table.pop();

        state_lwbnd = dsi.evalLowerBound(state->ssr);

        //compute upperbound here
        count_eval_lwbnd++;

        const double max_angle = ssrMaxAngle(state->ssr);

        // Update solution
        if (state_lwbnd > lwbnd)
        {

            lwbnd = state_lwbnd;
            guessAndResult = state->ssr;
            table.prune(lwbnd);
        }


        // Stopping criterion
//        mxAssert( state->bnd >= lwbnd, "Bound error");
        if (state->bnd - lwbnd <= gap)
        {

            delete state;
            break;
        }

        // Branch
        np = reg::search::split(state->ssr, ssr_array);
        delete state;

        for(i=0; i<np; i++)
        {

            upbnd = dsi.evalUpperBound(*(ssr_array[i]), lwbnd);

            if ( upbnd > lwbnd )
            {

                state = new SearchState<SSR>(*(ssr_array[i]), upbnd);
                table.push(state);
            }

            delete ssr_array[i];
        }
        count_eval_ubbnd += np;
    }

    delete []ssr_array;
    return lwbnd;
}

    
    
//TODO: indexation data structure should be in the template...
template <class SSR, unsigned int BRANCHING_FACTOR>
int searchTableDF(const DataIndexation<SSR> &dsi,
                     int lwbnd, int gap, int buckets,
                     SSR &guessAndResult)
{
    int i, np;
    int state_lwbnd, upbnd;

    int count_eval_ubbnd;
    int count_eval_lwbnd;

    SearchState<SSR, int> *state;
    SSR **ssr_array; //array of pointers

    // Initialise
    count_eval_ubbnd = count_eval_lwbnd = 0;


    upbnd = dsi.evalUpperBound(guessAndResult, lwbnd);

    count_eval_ubbnd++;


    if (upbnd-lwbnd <= gap)
    {
        return lwbnd;
    }

    StatePriorityHashtableDF<SSR, int, SearchState > table(buckets);


    ssr_array = new SSR*[BRANCHING_FACTOR];
    state = new SearchState<SSR>(guessAndResult, upbnd);
    table.push(state);

    int iter = 0;

    while (table.size())
    {
        iter++;

        // Find the state with the highest upper bound
        state = table.pop();

        // Evaluate lower boud
        state_lwbnd = dsi.evalLowerBound( state->ssr );
        count_eval_lwbnd++;

        // Update solution
        if (state_lwbnd > lwbnd)
        {
            lwbnd = state_lwbnd;
            guessAndResult = state->ssr;
            table.prune(lwbnd);
        }


        if(iter%10000000==0)
        {
            std::cout<< "bnb: "<< iter <<" lwbnd "<<lwbnd<<" upbnd "<<state->bnd
            <<" tablesize "<<table.size()<<std::endl;

        }

        // Stopping criterion
//        mxAssert( state->bnd >= lwbnd, "Bound error");
        if (state->bnd - lwbnd <= gap)
        {
            delete state;
            break;
        }

        // Branch
        np = reg::search::split(state->ssr, ssr_array);
        delete state;

        for(i=0; i<np; i++)
        {
            upbnd = dsi.evalUpperBound(*(ssr_array[i]), lwbnd);
            if ( upbnd > lwbnd )
            {
                state = new SearchState<SSR>(*(ssr_array[i]), upbnd);
                table.push(state);
            }
            delete ssr_array[i];
        }
        count_eval_ubbnd += np;
    }

    delete []ssr_array;
    return lwbnd;
}

    
    

// Search using Matching Lists
template <class SSR, unsigned int BRANCHING_FACTOR>
int bnb_search_ml_table(const DataIndexation<SSR> &dsi,
                        const Matrix1X &IN_d1, const Matrix1X &IN_d2,
                        const Matrix1X &CAT_d1, const Matrix1X &CAT_d2,
                        const double mag_th,
                        int lwbnd, int gap, int buckets, SSR &guessAndResult)
{
    int i,j,k, np;
    int state_lwbnd, upbnd;
    
    int count_eval_ubbnd;
    int count_eval_lwbnd;
    int *matchList;
    const size_t dsiSize = dsi.size();


    SearchStateML<SSR> *state, *childState;
    SSR **ssr_array; //array of pointers

    // Initialise
    count_eval_ubbnd = count_eval_lwbnd = 0;

    matchList = new int[dsiSize];
    for(i=0; i<dsiSize; i++)
    {
        matchList[i]=i;
    }

    std::vector<bool> matchesMap(dsiSize);
    upbnd = dsi.evalUpperBound(guessAndResult, lwbnd, matchList, dsiSize, matchesMap);

//    mxAssert(upbnd == dsi.evalUpperBound(guessAndResult,0), "");
    count_eval_ubbnd++;

    if (upbnd-lwbnd <= gap)
    {
        delete matchList;
        return lwbnd;
    }

    StatePriorityHashtable<SSR, int, SearchStateML > table(buckets);

    
    //Create new matchlist
    i=0;
    for(j=0; j<dsiSize; j++)
    {
        if(matchesMap[j])
        {
            matchList[i++]=matchList[j];
        }
    }
//    mxAssert(i<=upbnd,"" ); // equal condition is not valid when worked with offsets

    ssr_array = new SSR*[BRANCHING_FACTOR];

    state = new SearchStateML<SSR>(guessAndResult, upbnd, matchList, i);
    table.push(state);


//    struct rusage usage;
//    getrusage(RUSAGE_SELF, &usage);
//    long peak_memory = usage.ru_maxrss;

    int iter = 0;
    while (table.size())
    {

        iter++;

        // Find the state with the highest upper bound
        state = table.pop();

        // Evaluate lower boud
        state_lwbnd = dsi.evalLowerBound(state->ssr, state->matchList, state->matchListSize);
        count_eval_lwbnd++;

        // Update solution
        if (state_lwbnd > lwbnd)
        {
            lwbnd = state_lwbnd;
            guessAndResult = state->ssr;
            table.prune(lwbnd);
        }

//        if(iter%10000==0)
//        {
//            std::cout<< "bnb ml: "<< iter <<" lwbnd "<<lwbnd<<" upbnd "<<state->bnd
//                     <<" table size "<<table.size()<<" state "<<state->ssr<<std::endl;
////	        break;
//        }

        // Stopping criterion
//        mxAssert( state->bnd >= lwbnd, "Bound error");
        if (state->bnd - lwbnd <= gap)
        {
            delete state;
            break;
        }

        // Branch
        np = reg::search::split(state->ssr, ssr_array);
        for(i=0; i<np; i++)
        {
            upbnd = dsi.evalUpperBound(
                        *(ssr_array[i]), lwbnd,
                        state->matchList, state->matchListSize, matchesMap);
//            mxAssert(dsi.evalUpperBound(*(ssr_array[i]), 0, state->matchList, state->matchListSize, matchesMap)==
//                   dsi.evalUpperBound(*(ssr_array[i]), 0),"");

            //std::cout<<"upbnd "<< upbnd <<" curbest" <<curbest<<std::endl;
            if ( upbnd > lwbnd )
            {
                matchList = new int[upbnd];
                k=0;
                for(j=0; j<state->matchListSize; j++)
                {
                    if(matchesMap[j])
                    {
                        matchList[k++]=state->matchList[j];
                    }
                }
//                mxAssert(k<=upbnd,"");

                childState = new SearchStateML<SSR>(*(ssr_array[i]), upbnd, matchList, k);
                table.push(childState);
            }
            delete ssr_array[i];
        }
        delete state;
        count_eval_ubbnd += np;

//        peak_memory = std::max(peak_memory, usage.ru_maxrss);


    }
    // std::cout<<peak_memory<<std::endl;
    
    // write your memory here
//    std::fstream file;
//    // std::string mem_dir = "/home/demo/Desktop/ROSIA_c_arm/output/memory.txt";
//    file.open(mem_dir, std::ios::out | std::ios::app);
//    file<<"peak_memory: "<<peak_memory<<"\n";
//    file.close();

    delete []ssr_array;

    return lwbnd;
}

} // End namespace reg
} // End namespace search

