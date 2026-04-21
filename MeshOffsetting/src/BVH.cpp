#include "BVH.h"
#include <memory_resource>
#include <thread>
#include <algorithm>
#include <utility>

template <class Primitive>
BVH<Primitive>::BVH(std::vector<Primitive> prims, int maxPrimitiveNode, SplitMethod type) :
    _primitives(std::move(prims)), maxPrimitiveNode(std::min(255, maxPrimitiveNode)), _primitives(std::move(prims)), _splitType(type)
{
    std::vector<BVHPrimitive> bvhPrimitives;
    bvhPrimitives.resize(_primitives.size());
    for (int i = 0; i < _primitives.size(); ++i)
    {
        AABBBox ibox = _primitives[i].GetBBox();
        bvhPrimitives[i] = BVHPrimitive(i, ibox);
        
        /*std::pmr::monotonic_buffer_resource resource;
        std::allocator allocator(&resource);
        using Resource = std::pmr::monotonic_buffer_resource;
        using Allocator = std::pmr::polymorphic_allocator<std::byte>
        std::vector<std::unique_ptr<Resource>> threadBufferResources;*/
        //ignore thread part...

    }

    BVHNode* rootNode;
    std::atomic<int> totalNodes{ 0 };
    std::vector<Primitive> orderedPrims;
    if (type == HLBVH)
    {
        
    }
    else
    {
        std::atomic<int> orderedPrimsOffset{ 0 };
        rootNode = buildBVH(std::span<BVHPrimitive>(bvhPrimitives), &totalNodes, &orderedPrimsOffset, orderedPrims);
    }
    _primitives.swap(orderedPrims);

    bvhPrimitives.resize(0);
    bvhPrimitives.shrink_to_fit();
    nodes = new LinearBVHNode[totalNodes];
    int offset = 0;
    flattenBVH(rootNode, &offset);
}

template <class Primitive>
BVH<Primitive>::~BVH()
{
    delete[] nodes;
}

template <class Primitive>
float BVH<Primitive>::CalcDistancePointToBound(const Point3m& p, const AABBBox& box)
{
    Scalarm dx = 0, dy = 0, dz = 0;

    if (p.X() < box.min.X()) dx = box.min.X() - p.X();
    else if (p.X() > box.max.X()) dx = p.X() - box.max.X();

    if (p.Y() < box.min.Y()) dy = box.min.Y() - p.Y();
    else if (p.Y() > box.max.Y()) dy = p.Y() - box.max.Y();

    if (p.Z() < box.min.Z()) dz = box.min.Z() - p.Z();
    else if (p.Z() > box.max.Z()) dz = p.Z() - box.max.Z();

    return dx * dx + dy * dy + dz * dz;
}

template <class Primitive>
bool BVH<Primitive>::CalcDistancePointToPrimitive(const Point3m& p, const Primitive& prim, QueryResult& result)
{
    if constexpr (std::is_same_v<Primitive, CFaceO>)
    {
        //返回重心坐标
        return true;
    }
    return false;
}

template <class Primitive>
void BVH<Primitive>::QueryClosestPoint(const Point3m& p, QueryResult& result)
{
    float minDistance = FLT_MAX;
    Point3m closestPoint;
    Primitive closestPrim;

    std::stack<int> prims;
    int seed = 0;
    prims.push(seed);
    while (!prims.empty())
    {
        int id = prims.pop();
        const LinearBVHNode* node = nodes[id];
        if (CalcDistancePointToBound(node.bounds) > minDistance)
            continue;
        if (node->nPrimitives == 0)
        {
            for (int i = 0; i < node->leafOffset; ++i)
            { 
                QueryResult curResult = CalcDistancePointToPrimitive(_primitives[i]);
                if (curResult.dist < minDistance)
                {
                    result = curResult;
                }
            }
        }
        else
        {
            auto leftNode = nodes[id + 1];
            auto rightNode = nodes[id + 2];
            auto dis0 = CalcDistancePointToBound(leftNode.bounds);
            auto dis1 = CalcDistancePointToBound(rightNode.bounds);
            if (dis0 > dis1)
            {
                if (dis0 < minDistance)
                {
                    prims.push(id + 1);
                }
                if (dis1 < minDistance)
                {
                    prims.push(id + 2);
                }
            }
            else
            {
                if (dis1 < minDistance)
                {
                    prims.push(id + 2);
                }
                if (dis0 < minDistance)
                {
                    prims.push(id + 1);
                }
            }
        }
    }
    
    result.closestPoint = closestPoint;
    result.dist = minDistance;
    result.intersected = ...;
    
}

template <class Primitive>
BVHNode* BVH<Primitive>::buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive> orderedPrims)
{
    ++*totalNodes;
    BVHNode* node;
    AABBBox boxSum;
    for (const auto& it : primitives)
    {
        boxSum.Add(it->bound);
    }

    //float bestCost = primitives.size() * 1.0f;
    float nodeSA = node->bound.SurfaceArea();

    if (nodeSA < epsilon || primitives.size() == 1)
    {
        int leafOffset = orderedPrimsOffset->fetch_add(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i)
        {
            int index = primitives[i].primitiveIndex;
            orderedPrims[leafOffset + i] = primitives[index];
        }
        node->InitLeaf(leafOffset, primitives.size(), boxSum);
        return node;
    }
    else
    {
        // Find the axis with the largest centroid spread.
        AABBBox centerBox;
        for (const auto& it : primitives)
        {
            centerBox.Add(it.Centroid());
        }
        int dim = centerBox.MaxDim();
        // All centroids collapse to a single point on this axis.
        if (centerBox.min[dim] == centerBox.max[dim])
        {
            int leafOffset = orderedPrimsOffset->fetch_add(primitives.size());
            for (size_t i = 0; i < primitives.size(); ++i)
            {
                int index = primitives[i].primitiveIndex;
                orderedPrims[leafOffset + i] = primitives[index];
            }
            node->InitLeaf(leafOffset, primitives.size(), centerBox);
            return node;
        }
        else
        {
            int mid = primitives.size() / 2;
            switch (_splitType)
            {
            case Middle:
            {
                Scalarm middle = (centerBox.min[dim] + centerBox.max[dim]) * 0.5;
                auto midIter = std::partition(primitives.begin(), primitives.end(), [dim, middle](const BVHPrimitive& pi)
                    {
                        return pi.Centroid()[dim] < middle;
                    });
                mid = midIter - primitives.begin();
                if (midIter != primitives.begin() && midIter != primitives.end())
                    break;
            }
            case EqualCounts:
            {
                Scalarm middle = (centerBox.min[dim] + centerBox.max[dim]) * 0.5;
                std::nth_element(primitives.begin(), primitives.begin() + middle,
                    primitives.end(),
                    [dim](const BVHPrimitive& a, const BVHPrimitive& b) {
                        return a.Centroid()[dim] < b.Centroid()[dim];
                    }); // Partially sort primitives around the median.
                break;
            }
            case SAH:
            default:
            {
                if (primitives.size() <= 2)
                {
                    mid = primitives.size() / 2;
                    std::nth_element(primitives.begin(), primitives.begin() + mid,
                        primitives.end(),
                        [dim](const BVHPrimitive& a, const BVHPrimitive& b) {
                            return a.Centroid()[dim] < b.Centroid()[dim];
                        });
                }
                else
                {
                    // Bucket primitives for SAH evaluation.
                    constexpr int nBucket = 12;
                    SplitBucket buckets[nBucket];
                    for (const auto& prim : primitives)
                    {
                        int b = nBucket * centerBox.NormalizePointToBounds(prim.Centroid())[dim];
                        if (b == nBucket)
                        {
                            b = nBucket - 1;
                        }
                        buckets[b].count++;
                        buckets[b].bound.Add(prim.bound);
                    }

                    // Each split cost is the sum of left and right bucket area-weighted counts.
                    constexpr int nSplits = nBucket - 1;
                    float costs[nSplits] = {};

                    AABBBox leftBound;
                    int leftCount = 0;
                    for (int i = 0; i < nSplits; ++i)
                    {
                        leftBound.Add(buckets[i].bound);
                        leftCount += buckets[i].count;
                        costs[i] += leftCount * leftBound.SurfaceArea();
                    }

                    AABBBox rightBound;
                    int rightCount = 0;
                    for (int i = nSplits; i >= 1; ++i)
                    {
                        rightBound.Add(buckets[i].bound);
                        rightCount += buckets[i].count;
                        costs[i - 1] += rightCount * rightBound.SurfaceArea();
                    }

                    int minCostSplitBucket = -1;
                    float minCost = FLT_MAX;
                    for (int i = 0; i < nSplits; ++i)
                    {
                        if (costs[i] < minCost)
                        {
                            minCost = costs[i];
                            minCostSplitBucket = i;
                        }
                    }

                    //sah
                    // Treating this node as a leaf requires testing all primitives in it.
                    float leafCost = primitives.size();
                    minCost = 1.0 / 2.0 + minCost / boxSum.SurfaceArea();
                    // Split when the leaf would be too large or SAH prefers splitting.
                    if (primitives.size() > maxPrimitiveNode || minCost < leafCost)
                    {
                        auto midIter = std::partition(primitives.begin(), primitives.end(),
                            [=](const BVHPrimitive& perBounds)
                            {
                                int index = nBucket * centerBox.NormalizePointToBounds(perBounds.Centroid())[dim];
                                if (index == nBucket) index == nBucket - 1;
                                return index <= minCostSplitBucket;
                            });
                        mid = midIter - primitives.begin();
                    }
                    else
                    {
                        // Reserve a contiguous range in the global ordered primitive array.
                        int leafOffset = orderedPrimsOffset->fetch_add(primitives.size());
                        for (size_t i = 0; i < primitives.size(); ++i)
                        {
                            int index = primitives[i].primitiveIndex;
                            orderedPrims[leafOffset + i] = primitives[index];
                        }
                        // Each leaf stores the first primitive offset, count, and bounding box.
                        node->InitLeaf(leafOffset, primitives.size(), boxSum);
                        return node;
                    }
                }
                break;
            }
            }

            BVHNode* children[2];
            if (primitives.size() > 1024 * 128)
            {
                // ParallelFor
                for (int i = 0; i < 2; ++i)
                {
                    if (i == 0)
                    {
                        children[0] =
                            buildBVH(primitives.subspan(0, mid),
                                totalNodes, orderedPrimsOffset, orderedPrims);
                    }
                    else
                    {
                        children[1] =
                            buildBVH(primitives.subspan(mid),
                                totalNodes, orderedPrimsOffset, orderedPrims);
                    }
                }
            }
            else
            {
                children[0] =
                    buildBVH(primitives.subspan(0, mid),
                        totalNodes, orderedPrimsOffset, orderedPrims);
                children[1] =
                    buildBVH(primitives.subspan(mid),
                        totalNodes, orderedPrimsOffset, orderedPrims);
            }
            node->InitInterior(dim, children[0], children[1]);
        }
    }
    return node;
}

//preorder
template <class Primitive>
int BVH<Primitive>::flattenBVH(BVHNode* node, int* offset)
{
    LinearBVHNode* linearNode = &nodes[*offset];
    linearNode->bound = node->bound;
    int nodeOffset = (*offset)++;
    if (node->childCount > 0)
    {
        linearNode->leafOffset = node->leafOffset;
        linearNode->nPrimitives = node->childCount;
    }
    else
    {
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVH(node->childNode[0], offset);
    }
    return nodeOffset;
}