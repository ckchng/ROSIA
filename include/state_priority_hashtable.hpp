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
using namespace reg::search;


template <typename SSR, typename Scalar, template <typename, typename> class SS >
StatePriorityHashtable<SSR, Scalar, SS>::Queue::Queue(): size(0)
{}

template <typename SSR, typename Scalar, template <typename, typename> class SS >
StatePriorityHashtable<SSR, Scalar, SS>::Queue::~Queue()
{
    Node *node, *tmp;
    node = head;
    for (int i=0; i<size; i++)
    {
       tmp = node;
       node = node->next;
       delete tmp;
    }
}




template <typename SSR, typename Scalar, template <typename, typename> class SS >
SS<SSR, Scalar> *StatePriorityHashtable<SSR, Scalar, SS>::Queue::pop()
{
//    mxAssert(size>0, "pop in an empty queue");

    Node *tmp = head;
    SS<SSR, Scalar> *ss = head->ss;
    head = head->next;
    tmp->ss=NULL;
    delete tmp;
    size--;

//    mxAssert(size>=0, "Invalid size");
    return ss;
}

template <typename SSR, typename Scalar, template <typename, typename> class SS >
void StatePriorityHashtable<SSR, Scalar, SS>::Queue::push(SS<SSR, Scalar> *ss)
{
    Node *node;

    if (size==0)
    {
        head = tail = new Node(ss);
        size = 1;
        return;
    }

    node = tail;
    // append at the end
    tail = new Node(ss);
    node->next = tail;
    size++;
}


template <typename SSR, typename Scalar, template <typename, typename> class SS >
void StatePriorityHashtable<SSR, Scalar, SS>::Queue::stack(SS<SSR, Scalar> *ss)
{
    Node *node;
    
    if (size==0)
    {
        head = tail = new Node(ss);
        size = 1;
        return;
    }
    
    node = head;
    // stack at the begining
    head = new Node(ss);
    //node->next = tail;
    head->next = node;
    size++;
}



//-------------------------------
//     State Priority Hashtable
//-------------------------------

template <typename SSR, typename Scalar, template <typename, typename> class SS >
StatePriorityHashtable<SSR, Scalar, SS>::StatePriorityHashtable(int bucketSize):
    m_max_upbnd(0), m_size(0)
{
#ifdef __linux__ //asuming gcc
    m_map = std::tr1::unordered_map<Scalar, Queue*, Hash>(bucketSize, Hash());
#else
    m_map = std::unordered_map<Scalar, Queue*, Hash>(bucketSize, Hash());
#endif
}

template <typename SSR, typename Scalar, template <typename,typename> class SS >
StatePriorityHashtable<SSR, Scalar, SS>::~StatePriorityHashtable()
{
    typename Map::iterator it;
    for(it = m_map.begin(); it != m_map.end(); ++it)
    {
        delete it->second;
    }
}



template <typename SSR, typename Scalar, template <typename, typename> class SS >
SS<SSR, Scalar> *StatePriorityHashtable<SSR, Scalar, SS>::pop()
{
    SS<SSR, Scalar> *searchState;
    Queue *queue;
    //SSR ssr;

//    mxAssert(m_size>0 , "pop in an Empty table");
//    mxAssert(m_map.count(m_max_upbnd)>0, "No queue in table");
//    mxAssert(m_map[m_max_upbnd]==m_max_upbnd_queue, "m_max_upbnd_queue");

    queue = m_max_upbnd_queue;
    searchState = queue->pop();
    m_size--;

    // Remove queue if it is empty
    if(queue->size==0)
    {
        m_map.erase(m_max_upbnd);
        delete queue;

        // Update max_upbnd
        if (m_size)
        {
            m_max_upbnd--;
            while( m_map.count(m_max_upbnd)==0)
            {
                m_max_upbnd--;
            }
            
            m_max_upbnd_queue = m_map[m_max_upbnd];
        }
        else
        {
            m_max_upbnd=0;
        }
    }

    return searchState;
}

template <typename SSR, typename Scalar, template <typename, typename> class SS >
void StatePriorityHashtable<SSR, Scalar, SS>::push(SS<SSR, Scalar> *state)
{
//    mxAssert(state->bnd>=0, "inserting state with upbnd<=0");
    Queue *q;

    if ( m_map.count(state->bnd)==0)
    {
        q = new Queue();
        m_map[state->bnd] = q;
    }
    else
    {
        q = m_map[state->bnd];
//        mxAssert(q->size>0, "Invalid queue size in hashtable");
//        mxAssert(m_map.count(state->bnd)>0, "Map error");
    }

    q->push(state);
//    mxAssert(q->size>0, "Invalid queue size after push");

    // Update max_upbnd
    if(state->bnd>m_max_upbnd )
    {
        m_max_upbnd = state->bnd;
        m_max_upbnd_queue = q;
    }

    m_size++;
}

template <typename SSR, typename Scalar, template <typename, typename> class SS >
void StatePriorityHashtable<SSR, Scalar, SS>::prune(int lwbnd)
{
    Queue *queue;
    typename Map::iterator it;

    it = m_map.begin();

    while(it != m_map.end())
    {
        if (it->first <= lwbnd)
        {
            queue = it->second;
            // remove and delete queue
            m_size -= queue->size;
            delete queue;

            it = m_map.erase(it);
        }
        else
        {
            it++;
        }
    }

    // Update max_upbnd
    if (m_size==0)
    {
        m_max_upbnd=0;
    }
}


template <typename SSR, typename Scalar, template <typename, typename> class SS >
void StatePriorityHashtable<SSR, Scalar, SS>::dump(std::ofstream& ofs) const
{
    Queue *queue;
    typename Map::const_iterator it;

    it = m_map.begin();

    while(it != m_map.end())
    {
        queue = it->second;
        //queue->dump();
        queue->dump(ofs);
        it++;
    }
}

template <typename SSR, typename Scalar, template <typename, typename> class SS >
void StatePriorityHashtable<SSR, Scalar, SS>::Queue::dump(std::ofstream& ofs) const
{
    Node *node = head;
    SS<SSR, Scalar> *ss;

    for (int i=0; i<size; i++)
    {
        ss= node->ss;
        ofs<<*ss<<"\n";
        node = node->next;
    }
}

template <typename SSR, typename Scalar, template <typename, typename> class SS >
unsigned int StatePriorityHashtable<SSR, Scalar, SS>::size () const
{
    return m_size;
}
