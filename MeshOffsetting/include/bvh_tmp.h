#include "mymesh.h"
#include"AABB.h"
#include <span>
#include <vcg/space/ray3.h>

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
    BVHPrimitive(int i, AABBBox box) :primitiveIndex(i), bound(box)
    {
    }
    ~BVHPrimitive() {};

    int primitiveIndex;
    AABBBox bound;
    Point3m Centroid() const
    {
        return bound.Center();
    }
};
struct alignas(32) LinearBVHNode
{
    AABBBox bounds;
    union
    {
        int leafOffset;   // leaf
        int interiorOffset; // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
};
class BVHNode
{
public:
    BVHNode* childNode[2];
    AABBBox bound;
    unsigned int splitAxis;
    unsigned int childCount;
    unsigned int leafOffset;

    void InitLeaf(int first, int count, const AABBBox& b)
    {
        leafOffset = first;
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

struct QueryResult
{
    double dist;
    int sign;//0/1/-1
    bool intersected;
    MESHLAB_SCALAR s, t; //barycentric coordinates
    Point3m closestPoint;
    int id;

};

template <class Primitive>
class BVH
{
public:
    BVH(std::vector<Primitive> prims, int maxPrimitiveNode = 1, SplitMethod type = SplitMethod::SAH);
    ~BVH();

    void QueryClosestPoint(const Point3m& p, QueryResult& result);
    //static BVH* Create(std::vector<Primitive> primitives);
    //bool IntersectP(const Ray3m& ray, float t);


private:

    BVHNode* buildHLBVH();

    // Build the BVH recursively and return the root node.
    // totalNodes stores the total number of BVH nodes.
    // orderedPrimsOffset stores the offset into the ordered primitive array.
    // orderedBoxs stores the reordered primitives used by the BVH.
    BVHNode* buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive>& orderedPrims);
    int flattenBVH(BVHNode* node, int* offset);
    unsigned int calcSplit(BVHNode* node);

    float CalcDistancePointToBound(const Point3m& p, const AABBBox& box);
    QueryResult CalcDistancePointToPrimitive(const Point3m& p, const Primitive& prim);

private:
    int maxPrimitiveNode;
    std::vector<Primitive> _primitives;
    LinearBVHNode* nodes = nullptr;
    //int _totalNodes = 0;
    SplitMethod _splitType;
};
