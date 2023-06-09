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


#ifndef REG_STATE_PRIORITY_HASHTABLE_H_
#define REG_STATE_PRIORITY_HASHTABLE_H_

#include "state.h"
#include <cstddef>

#ifdef __linux__ //asuming gcc
#include <tr1/unordered_map>
using namespace std::tr1;
#else
#include <unordered_map> //asuming clang
#endif

namespace reg {
namespace search {

//---------------------------------------------
//      State Priority Hashtable
//---------------------------------------------

template <typename SSR, typename Scalar, template <typename SSR_, typename Scalar_> class SS >
class StatePriorityHashtable
{

private:
    class Queue;
    struct Hash;

    #ifdef __linux__ //asuming gcc
    typedef std::tr1::unordered_map<Scalar, Queue*, Hash> Map;
    #else
    typedef std::unordered_map<Scalar, Queue*, Hash> Map;
    #endif

    int m_max_upbnd;
    Queue *m_max_upbnd_queue;
    unsigned int m_size; // number of states

    class Queue
    {
    private:
        class Node
        {
        public:
            Node *next;
            SS<SSR, Scalar> *ss;

            Node(SS<SSR, Scalar> *ss):ss(ss){}
            ~Node()
            {
                delete ss;
            }

        };
    public:
        Queue();
        ~Queue();

        Node *head;
        Node *tail;

        int size;
        SS<SSR, Scalar> *pop();
        void push(SS<SSR, Scalar> *ss);
        void stack(SS<SSR, Scalar> *ss); //not the best design...

        void dump(std::ofstream& ofs) const;
    }; // Queue


    struct Hash
    {
        std::size_t operator()(const int& q) const
        {
            unsigned int key=q;
            key = ((key >> 16) ^ key) * 0x45d9f3b;
            key = ((key >> 16) ^ key) * 0x45d9f3b;
            key = ((key >> 16) ^ key);
            return key;
        }
    };

    Map m_map;


public:
    StatePriorityHashtable(int bucketSize);
    ~StatePriorityHashtable();

    SS<SSR, Scalar> *pop();
    void push(SS<SSR, Scalar> *state);
    void prune(int curbest);
    unsigned int size() const;

    //dump table SSR to file
    void dump(std::ofstream& ofs) const;
};

    
    
    
//---------------------------------------------
//      State Priority Hashtable
//---------------------------------------------
    
template <typename SSR, typename Scalar, template <typename SSR_, typename Scalar_> class SS >
class StatePriorityHashtableDF
{
        
    private:
        class Queue;
        struct Hash;
        
#ifdef __linux__ //asuming gcc
        typedef std::tr1::unordered_map<Scalar, Queue*, Hash> Map;
#else
        typedef std::unordered_map<Scalar, Queue*, Hash> Map;
#endif
        
    int m_max_upbnd;
    Queue *m_max_upbnd_queue;
    unsigned int m_size; // number of states
        
    class Queue
    {
    private:
        class Node
        {
        public:
            Node *next;
            SS<SSR, Scalar> *ss;
            
            Node(SS<SSR, Scalar> *ss):ss(ss){}
            ~Node()
            {
                if (ss!=NULL) delete ss;
            }
            
        };
    public:
        Queue();
        ~Queue();
        
        Node *head;
        Node *tail;
        
        int size;
        SS<SSR, Scalar> *pop();
        void push(SS<SSR, Scalar> *ss);
        void stack(SS<SSR, Scalar> *ss); //not the best design...
        
        void dump(std::ofstream& ofs) const;
    }; // Queue
    
    
    struct Hash
    {
        std::size_t operator()(const int& q) const
        {
            unsigned int key=q;
            key = ((key >> 16) ^ key) * 0x45d9f3b;
            key = ((key >> 16) ^ key) * 0x45d9f3b;
            key = ((key >> 16) ^ key);
            return key;
        }
    };
    
    Map m_map;
    
    
public:
    StatePriorityHashtableDF(int bucketSize);
    ~StatePriorityHashtableDF();
    
    SS<SSR, Scalar> *pop();
    void push(SS<SSR, Scalar> *state);
    void prune(int curbest);
    unsigned int size() const;
    
    //dump table SSR to file
    void dump(std::ofstream& ofs) const;
};
    
} // End namespace search
} // End namespace reg


#include "state_priority_hashtable.hpp"
#include "state_priority_hashtable_df.hpp"

#endif
