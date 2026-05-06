#include "BVH.h"
#include <memory_resource>
#include <thread>
#include <algorithm>
#include <utility>
#include <limits>
#include <type_traits>
#include <cmath>

namespace
{
    template <class Primitive>
    decltype(auto) PrimitiveRef(const Primitive& primitive)
    {
        if constexpr (std::is_pointer_v<Primitive>)
        {
            assert(primitive != nullptr);
            return *primitive;
        }
        else
        {
            return (primitive);
        }
    }
}

template <class Primitive>
BVH<Primitive>::BVH(std::vector<Primitive> prims, int maxPrimitiveNode, SplitMethod type) :
    _primitives(std::move(prims)), maxPrimitiveNode(std::min(255, maxPrimitiveNode)), _splitType(type)
{
    std::vector<BVHPrimitive> bvhPrimitives;
    bvhPrimitives.resize(_primitives.size());
    for (int i = 0; i < _primitives.size(); ++i)
    {
        AABBBox ibox;
        PrimitiveRef(_primitives[i]).GetBBox(ibox);
        bvhPrimitives[i] = BVHPrimitive(i, ibox);
        
        /*std::pmr::monotonic_buffer_resource resource;
        std::allocator allocator(&resource);
        using Resource = std::pmr::monotonic_buffer_resource;
        using Allocator = std::pmr::polymorphic_allocator<std::byte>
        std::vector<std::unique_ptr<Resource>> threadBufferResources;*/
        //ignore thread part...

    }

    BVHNode* rootNode = nullptr;
    std::atomic<int> totalNodes{ 0 };
    std::vector<Primitive> orderedPrims(_primitives.size());
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
QueryResult BVH<Primitive>::CalcDistancePointToPrimitive(const Point3m& p, const Primitive& prim)
{
    QueryResult result{};
    result.dist = std::numeric_limits<double>::max();
    result.id = -1;
    result.sign = 0;
    result.intersected = false;
    result.s = 0;
    result.t = 0;
    using PrimitiveValue = std::remove_cv_t<std::remove_pointer_t<Primitive>>;
    if constexpr (std::is_same_v<PrimitiveValue, CFaceO>)
    {
        const auto& face = PrimitiveRef(prim);
        const Point3m a = face.P(0);
        const Point3m b = face.P(1);
        const Point3m c = face.P(2);

        auto setResult = [&](Scalarm s, Scalarm t)
            {
                result.s = s;
                result.t = t;
                result.closestPoint = a + (b - a) * s + (c - a) * t;
                const Point3m diff = p - result.closestPoint;
                result.dist = diff * diff;
                result.intersected = result.dist <= epsilon;
            };

        const Point3m ab = b - a;
        const Point3m ac = c - a;
        const Point3m ap = p - a;

        const Scalarm d1 = ab * ap;
        const Scalarm d2 = ac * ap;
        if (d1 <= 0 && d2 <= 0)
        {
            setResult(0, 0);
            return result;
        }

        const Point3m bp = p - b;
        const Scalarm d3 = ab * bp;
        const Scalarm d4 = ac * bp;
        if (d3 >= 0 && d4 <= d3)
        {
            setResult(1, 0);
            return result;
        }

        const Scalarm vc = d1 * d4 - d3 * d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0)
        {
            const Scalarm denom = d1 - d3;
            const Scalarm v = denom > epsilon ? d1 / denom : 0;
            setResult(v, 0);
            return result;
        }

        const Point3m cp = p - c;
        const Scalarm d5 = ab * cp;
        const Scalarm d6 = ac * cp;
        if (d6 >= 0 && d5 <= d6)
        {
            setResult(0, 1);
            return result;
        }

        const Scalarm vb = d5 * d2 - d1 * d6;
        if (vb <= 0 && d2 >= 0 && d6 <= 0)
        {
            const Scalarm denom = d2 - d6;
            const Scalarm w = denom > epsilon ? d2 / denom : 0;
            setResult(0, w);
            return result;
        }

        const Scalarm va = d3 * d6 - d5 * d4;
        if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
        {
            const Scalarm denom = (d4 - d3) + (d5 - d6);
            const Scalarm w = denom > epsilon ? (d4 - d3) / denom : 0;
            setResult(1 - w, w);
            return result;
        }

        const Scalarm denom = va + vb + vc;
        if (std::abs(denom) > epsilon)
        {
            const Scalarm v = vb / denom;
            const Scalarm w = vc / denom;
            setResult(v, w);
        }
        else
        {
            auto setSegmentResult = [&](const Point3m& p0, const Point3m& p1, Scalarm s0, Scalarm t0, Scalarm s1, Scalarm t1)
            {
                const Point3m edge = p1 - p0;
                const Scalarm len2 = edge * edge;
                Scalarm u = len2 > epsilon ? ((p - p0) * edge) / len2 : 0;
                u = std::clamp(u, Scalarm(0), Scalarm(1));

                QueryResult candidate{};
                candidate.id = -1;
                candidate.sign = 0;
                candidate.intersected = false;
                candidate.s = s0 + (s1 - s0) * u;
                candidate.t = t0 + (t1 - t0) * u;
                candidate.closestPoint = p0 + edge * u;
                const Point3m diff = p - candidate.closestPoint;
                candidate.dist = diff * diff;
                candidate.intersected = candidate.dist <= epsilon;
                if (candidate.dist < result.dist)
                {
                    result = candidate;
                }
            };

            setSegmentResult(a, b, 0, 0, 1, 0);
            setSegmentResult(a, c, 0, 0, 0, 1);
            setSegmentResult(b, c, 1, 0, 0, 1);
        }
    }
    
    return result;
}

template <class Primitive>
void BVH<Primitive>::QueryClosestPoint(const Point3m& p, QueryResult& result)
{
    assert(nodes != nullptr);
    float minDistance = FLT_MAX;
    Point3m closestPoint;
    Primitive closestPrim;

    std::stack<int> prims;
    int seed = 0;
    prims.push(seed);
    while (!prims.empty())
    {
        int id = prims.top();
        prims.pop();

        const LinearBVHNode* node = &nodes[id];
        if (CalcDistancePointToBound(p, node->bounds) >= minDistance)
            continue;
        if (node->nPrimitives != 0)
        {
            for (int i = 0; i < node->nPrimitives; ++i)
            { 
                QueryResult curResult = CalcDistancePointToPrimitive(p,_primitives[node->leafOffset + i]);
                if (curResult.dist < minDistance)
                {
                    minDistance = curResult.dist;
                    result = curResult;
                    result.id = node->leafOffset + i;
                }
            }
        }
        else
        {
            int left = id + 1;
            int right = node->interiorOffset;
            auto dl = CalcDistancePointToBound(p, nodes[left].bounds);
            auto dr = CalcDistancePointToBound(p, nodes[right].bounds);
            if (dl < dr)
            {
                if (dr < minDistance)
                {
                    prims.push(right);
                }
                if (dl < minDistance)
                {
                    prims.push(left);
                }
            }
            else
            {
                if (dl < minDistance)
                {
                    prims.push(left);
                }
                if (dr < minDistance)
                {
                    prims.push(right);
                }
            }
        }
    }
    
    //result.sign = (p - result.closestPoint) * _primitives[result.id].N() > 0 ? 1 : 0;
    return;
}

template <class Primitive>
BVHNode* BVH<Primitive>::buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive>& orderedPrims)
{
    ++*totalNodes;
    BVHNode* node = new BVHNode();
    AABBBox boxSum;
    for (const auto& it : primitives)
    {
        boxSum.Add(it.bound);
    }

    //float bestCost = primitives.size() * 1.0f;
    float nodeSA = boxSum.SurfaceArea();
    if (nodeSA < epsilon || primitives.size() == 1)
    {
        int leafOffset = orderedPrimsOffset->fetch_add(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i)
        {
            int index = primitives[i].primitiveIndex;
            orderedPrims[leafOffset + i] = _primitives[index];
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
                orderedPrims[leafOffset + i] = _primitives[index];
            }
            node->InitLeaf(leafOffset, primitives.size(), boxSum);
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
                                int index = nBucket* centerBox.NormalizePointToBounds(perBounds.Centroid())[dim];
                                if (index == nBucket) index = nBucket - 1;
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
                            orderedPrims[leafOffset + i] = _primitives[index];
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
    linearNode->bounds = node->bound;
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
        linearNode->interiorOffset = flattenBVH(node->childNode[1], offset);
    }
    return nodeOffset;
}

template class BVH<CFaceO*>;
