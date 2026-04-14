#include "mymesh.h"
#include"AABB.h"
#include <span>
//#include <nth_element>
enum SplitMethod
{
    SAH = 0, 
    HLBVH, 
    Middle, 
    EqualCounts,
};

struct SplitBucket
{
    int count = 0;
    AABBBox bound;
};
class BVHPrimitive
{
public:
    BVHPrimitive(int i, AABBBox box) :primitiveIndex(i), bound(box), centroid(box.Center())
    {
    }
    ~BVHPrimitive() {};

    int primitiveIndex;
    AABBBox bound;
    Point3m centroid;
};

class BVHNode
{
public:
    BVHNode* childNode[2];
    AABBBox bound;
    unsigned int splitAxis;
    unsigned int childCount;
    unsigned int firstPrimOffset;

    void InitLeaf(int first, int count, const AABBBox& b)
    {
        firstPrimOffset = first;
        childCount = count;
        bound = b;
        childNode[0] = nullptr;
        childNode[1] = nullptr;
    }
    void InitInterior(unsigned int axis, BVHNode* cl, BVHNode* cr)
    {
        childNode[0] = cl;
        childNode[1] = cr;
        splitAxis = axis;
        childCount = 0;
        bound = AABBBox::Merge(cr->bound, cl->bound);
    }

};

template <class Primitive>
class BVH
{
public:
    BVH(std::vector<Primitive> p, int maxPrimitiveNode = 1, SplitMethod type = SplitMethod::SAH);
    BVH(const CMeshO& mesh, int maxPrimitiveNode, SplitMethod type);
    ~BVH();

private:

    BVHNode* buildHLBVH();

    // Build the BVH recursively and return the root node.
    // totalNodes stores the total number of BVH nodes.
    // orderedPrimsOffset stores the offset into the ordered primitive array.
    // orderedBoxs stores the reordered primitives used by the BVH.
    BVHNode* buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive> orderedBoxs);

    void flattenBVH(BVHNode* rootNode, int* offset);
    unsigned int calcSplit(BVHNode* node);

private:
    int maxPrimitiveNode;
    std::vector<Primitive> _primitives;
    BVHNode* nodes = nullptr;
    SplitMethod _splitType;
};
