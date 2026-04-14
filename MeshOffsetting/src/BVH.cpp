#include "BVH.h"

#include <algorithm>
#include <utility>

template <class Primitive>
BVH<Primitive>::BVH(std::vector<Primitive> p, int maxPrimitiveNode, SplitMethod type)
    : maxPrimitiveNode(std::min(255, maxPrimitiveNode)),
    _primitives(std::move(p)),
    _splitType(type)
{

}

template <class Primitive>
BVH<Primitive>::BVH(const CMeshO& mesh, int maxPrimitiveNode, SplitMethod type) :
    _primitives(std::move(mesh.face)), maxPrimitiveNode(std::min(255, maxPrimitiveNode)), _splitType(type)
{
    std::vector<BVHPrimitive> bvhPrimitives;
    bvhPrimitives.resize(_primitives.size());
    for (int i = 0; i < _primitives.size(); ++i)
    {
        AABBBox ibox;
        ibox = _primitives[i].GetBBox();
        bvhPrimitives[i] = BVHPrimitive(i, ibox);
    }

    BVHNode* rootNode;
    std::atomic<int> totalNodes{ 0 };
    std::vector<Primitive> orderedBoxs;
    if (type == HLBVH)
    {
        
    }
    else
    {
        std::atomic<int> orderedPrimsOffset{ 0 };
        rootNode = buildBVH(std::span<BVHPrimitive>(bvhPrimitives), &totalNodes, &orderedPrimsOffset, orderedBoxs);
    }
    _primitives.swap(orderedBoxs);

    bvhPrimitives.resize(0);
    bvhPrimitives.shrink_to_fit();
    nodes = new BVHNode[totalNodes];
    int offset = 0;
    flattenBVH(rootNode, &offset);
}

template <class Primitive>
BVH<Primitive>::~BVH()
{
    delete[] nodes;
}

template <class Primitive>
BVHNode* BVH<Primitive>::buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive> orderedBoxs)
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
        int firstPrimOffset = orderedPrimsOffset->fetch_add(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i)
        {
            int index = primitives[i].primitiveIndex;
            orderedBoxs[firstPrimOffset + i] = primitives[index];
        }
        node->InitLeaf(firstPrimOffset, primitives.size(), boxSum);
        return node;
    }
    else
    {
        // Find the axis with the largest centroid spread.
        AABBBox centerBox;
        for (const auto& it : primitives)
        {
            centerBox.Add(it.centroid);
        }
        int dim = centerBox.MaxDim();
        // All centroids collapse to a single point on this axis.
        if (centerBox.min[dim] == centerBox.max[dim])
        {
            int firstPrimOffset = orderedPrimsOffset->fetch_add(primitives.size());
            for (size_t i = 0; i < primitives.size(); ++i)
            {
                int index = primitives[i].primitiveIndex;
                orderedBoxs[firstPrimOffset + i] = primitives[index];
            }
            node->InitLeaf(firstPrimOffset, primitives.size(), centerBox);
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
                        return pi.centroid[dim] < middle;
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
                        return a.centroid[dim] < b.centroid[dim];
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
                            return a.centroid[dim] < b.centroid[dim];
                        });
                }
                else
                {
                    // Bucket primitives for SAH evaluation.
                    constexpr int nBucket = 12;
                    SplitBucket buckets[nBucket];
                    for (const auto& prim : primitives)
                    {
                        int b = nBucket * centerBox.NormalizePointToBounds(prim.centroid)[dim];
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
                                int index = nBucket * centerBox.NormalizePointToBounds(perBounds.centroid)[dim];
                                if (index == nBucket) index == nBucket - 1;
                                return index <= minCostSplitBucket;
                            });
                        mid = midIter - primitives.begin();
                    }
                    else
                    {
                        // Reserve a contiguous range in the global ordered primitive array.
                        int firstPrimOffset = orderedPrimsOffset->fetch_add(primitives.size());
                        for (size_t i = 0; i < primitives.size(); ++i)
                        {
                            int index = primitives[i].primitiveIndex;
                            orderedBoxs[firstPrimOffset + i] = primitives[index];
                        }
                        // Each leaf stores the first primitive offset, count, and bounding box.
                        node->InitLeaf(firstPrimOffset, primitives.size(), boxSum);
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
                                totalNodes, orderedPrimsOffset, orderedBoxs);
                    }
                    else
                    {
                        children[1] =
                            buildBVH(primitives.subspan(mid),
                                totalNodes, orderedPrimsOffset, orderedBoxs);
                    }
                }
            }
            else
            {
                children[0] =
                    buildBVH(primitives.subspan(0, mid),
                        totalNodes, orderedPrimsOffset, orderedBoxs);
                children[1] =
                    buildBVH(primitives.subspan(mid),
                        totalNodes, orderedPrimsOffset, orderedBoxs);
            }
            node->InitInterior(dim, children[0], children[1]);
        }
    }
    return node;
}
