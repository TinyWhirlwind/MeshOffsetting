#include "mymesh.h"
enum Principle
{
    SPLIT_MIDDLE = 0,
    SPLIT_EQUAL_COUNT,
    SPLIT_SAH
};
class BVH
{
    class Node
    {
        public:
            Node *left, *right;
            int start, end;
            MyMesh::Point min_bound, max_bound;
    
            Node() : left(nullptr), right(nullptr), start(0), end(0) {}
    };
public:
        BVH(CMeshO &mesh);
        ~BVH();

   
};
