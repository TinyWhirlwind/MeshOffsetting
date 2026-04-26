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
template <class BV>
class BVHPrimitive
{
public:
    BVHPrimitive(int i, BV box,const AABBBox& aabbIn) :primitiveIndex(i), bound(box), aabb(aabbIn)
    {
    }
    ~BVHPrimitive() {};

    int primitiveIndex;
    BV bound;
    AABBBox aabb;
    Point3m Centroid() const
    {
        return aabb.Center();
    }
};

template <class BV>
struct alignas(32) LinearBVHNode
{
    BV bounds;
    union
    {
        int leafOffset;   // leaf
        int interiorOffset; // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
};

template <class BV>
class BVHNode
{
public:
    BVHNode<BV>* childNode[2];
    BV bound;
    unsigned int splitAxis;
    unsigned int childCount;
    unsigned int leafOffset;

    void InitLeaf(int first, int count, const BV& b)
    {
        leafOffset = first;
        childCount = count;
        bound = b;
        childNode[0] = nullptr;
        childNode[1] = nullptr;
    }
    void InitInterior(unsigned int axis, BVHNode<BV>* cl, BVHNode<BV>* cr)
    {
        childNode[0] = cl;
        childNode[1] = cr;
        splitAxis = axis;
        childCount = 0;
        bound = BV::Merge(cr->bound, cl->bound);
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

template <class Primitive, class BV>
class BVH
{
public:
    BVH(std::vector<Primitive> prims, int maxPrimitiveNode = 1, SplitMethod type = SplitMethod::SAH);
    ~BVH();

    void QueryClosestPoint(const Point3m& p, QueryResult& result);
    //static BVH* Create(std::vector<Primitive> primitives);
    //bool IntersectP(const Ray3m& ray, float t);


private:

    BVHNode<BV>* buildHLBVH();

    // Build the BVH recursively and return the root node.
    // totalNodes stores the total number of BVH nodes.
    // orderedPrimsOffset stores the offset into the ordered primitive array.
    // orderedBoxs stores the reordered primitives used by the BVH.
    BVHNode<BV>* buildBVH(std::span<BVHPrimitive<BV>> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive>& orderedPrims);
    int flattenBVH(BVHNode<BV>* node, int* offset);
    unsigned int calcSplit(BVHNode<BV>* node);

    float CalcDistancePointToBV(const Point3m& p, const BV& box);
    QueryResult CalcDistancePointToPrimitive(const Point3m& p, const Primitive& prim);

private:
    int maxPrimitiveNode;
    std::vector<Primitive> _primitives;
    LinearBVHNode<BV>* nodes = nullptr;
    //int _totalNodes = 0;
    SplitMethod _splitType;
};
