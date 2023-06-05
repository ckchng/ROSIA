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
#ifndef REG_RTREE_
#define REG_RTREE_

#include <iostream>
#include <queue>
#include <list>
#include "reg_common.h"
#include "geometry.h"

#define RTREE_M 25//10
#define RTREE_m 5

//-------------------------------------
//        Operators
//-------------------------------------

using namespace reg::geometry;

namespace reg
{
class RTree;
}

std::ostream& operator<<(std::ostream &, const reg::RTree *);

namespace reg
{

class RTree
{
private:
    class Node
    {
    public:
        Rectangle MBR[RTREE_M];
        Node *child[RTREE_M];
        bool leaf;
        int n;
        Circle *circle;

        Node(bool leaf);
        ~Node();
        friend std::ostream& operator <<(std::ostream &os, const Node *node)
        {
            os<<"[";
            for (int i=0;i<node->n;i++) os<< node->MBR[i] << ", ";
            os<<"]";
            return os;
        }

        void add(Node *node, Rectangle r); // Add child Node
        void add(Circle c, Rectangle r); //Add leaf node
    };

    Node *root;
    int _depth;
    bool _full;

    Node **_node_path;
    int *_idx_path;

    std::list<Circle> negCirclesList;
    std::list<HalfPlane> posHalfPlanesList;
    std::list<HalfPlane> negHalfPlanesList;

    void freeRTree(Node *);

    void updateMBR(Node **path, int *idx_path);
    void splitInter(Node **node_path, int *idx_path);
    void splitLeaf(Node *leaf, Node **node_path, int *idx_path);
    void printPosCircles(Node *node) const;

    bool intersectsPosHalfPlanes(double x, double y) const;
    bool intersectsNegHalfPlanes(double x, double y) const;
    bool intersectsNegHalfPlanes(Circle c) const;

    bool intersectsPosCircles(Node *node, double x, double y) const;
    bool intersectsPosCircles(Node *node, HalfPlane hp, bool hp_sign) const;
    bool intersectsPosCircles(Node *node, Circle c, bool c_sign) const;

    bool intersectsNegCircles(double x, double y) const;
    bool intersectsNegCircles(Circle c) const;


    void insertPos(Node *node, Circle c, Node **node_path, int *idx_path, int lev);
    void insertNeg(Circle c);

    void insertPos(HalfPlane hp);
    void insertNeg(HalfPlane hp);

public:
    RTree();
    ~RTree();
    friend std::ostream& (::operator <<) (std::ostream &, const RTree *);

    void insertPos(Circle c);
    bool intersectsPosCircles(double x, double y) const { return intersectsPosCircles(root, x, y);}
    bool intersectsPosCircles(Circle c, bool c_sign) const {return intersectsPosCircles(root, c, c_sign);}

    void addPatch(const Vector3 &centre, double angle);

    bool matchPoint(const Vector3 &p) const;
    bool matchPatch(const Vector3 &patchCentre, double patchAngle) const;

    int size() const;
    int depth() const {return _depth;}
    bool full() const {return _full;}
    void setFull() {_full=true;}

    void printNegCircles() const;
    void printPosCircles() const;
};

} // namespace reg
#endif
