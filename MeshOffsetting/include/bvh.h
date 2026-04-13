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
template <class Primitive>
class BVH
{
public:
    struct SplitBucket
    {
        int count = 0;
        AABBBox bound;
    };
    class BVHPrimitive
    {
    public:
        BVHPrimitive(int i, AABBBox box):primitiveIndex(i),bound(box), centroid(box.Center())
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
            bound = Merge(cr->bound, cl->bound);
        }

    };

    BVH(const CMeshO& mesh, int maxPrimitiveNode, SplitMethod type) :
        _primitives(std::move(mesh.face)), maxPrimitiveNode(std::min(255, maxPrimitiveNode)), _splitType(type)
    {
        std::vector<BVHPrimitive> bvhPrimitives;
        bvhPrimitives.resize(_primitives.size());
        for (int i = 0;  i< _primitives.size(); ++i)
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
            //动态变化,频繁建树，适合使用HLBVH算法构建BVH树，构建过程中需要频繁的内存分配和释放
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
    ~BVH();


private:

    BVHNode* buildHLBVH()
    {
    }

    //BVH线性构建，递归构建BVH树，返回根节点指针，totalNodes记录BVH树节点总数，orderedPrimsOffset记录有序的原始图元偏移量，orderedBoxs记录有序的原始图元边界框
    BVHNode* buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive> orderedBoxs)
    {
        ++*totalNodes;
        BVHNode node;
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
            node.InitLeaf(firstPrimOffset, primitives.size(), boxSum);
            return node;
        }
        else
        {
            //找重心分布最开的轴
            AABBBox centerBox;
            for (const auto& it:primitives)
            {   
                centerBox.Add(it.centroid);
            }
            int dim = centerBox.MaxDim();
            if (centerBox.min[dim] == centerBox.max[dim])//重心分布在点上
            {
                int firstPrimOffset = orderedPrimsOffset->fetch_add(primitives.size());
                for (size_t i = 0; i < primitives.size(); ++i)
                {
                    int index = primitives[i].primitiveIndex;
                    orderedBoxs[firstPrimOffset + i] = primitives[index];
                }
                node.InitLeaf(firstPrimOffset, primitives.size(), centerBox);
                return node;
            }
            else
            {
                int mid = primitives.size() / 2;
                switch (_splitType)
                {
                case SAH: 
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
                        //分桶
                        constexpr int nBucket = 12;
                        SplitBucket buckets[nBucket];
                        for (const auto& prim : primitives)
                        {
                            int b = nBucket * centerBox.NormalizePointToBounds(prim.centroid);
                            if (b == nBucket)
                            {
                                b = nBucket - 1;
                            }
                            buckets[b].count++;
                            buckets[b].bound.Add(prim.bound);
                        }

                        //calc cost每条split的cost都是左边+右边所有bucket的box*count;
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
                        float leafCost = primitives.size();//不分裂，直接作为叶子节点,需要测试 N 次 primitive
                        minCost = 1.0 / 2.0 + minCost / boxSum.SurfaceArea();
                        //约束叶子节点过大
                        if (primitives.size() > maxPrimitiveNode || minCost < leafCost)
                        {
                            auto midIter = std::partition(primitives.begin(), primitives.end(),
                                [](const BVHPrimitive & perBounds)
                                {
                                    int index = nBucket * centerBox.NormalizePointToBounds(perBounds.bound);
                                    if (index == nBucket) index == nBucket - 1;
                                    return index <= minCostSplitBucket;
                                });
                            mid = midIter - primitives.begin();
                        }
                        else
                        {
                            int firstPrimOffset = orderedPrimsOffset->fetch_add(primitives.size());//给当前叶子节点分配一段“全局连续数组”的位置
                            for (size_t i = 0; i < primitives.size(); ++i)
                            {
                                int index = primitives[i].primitiveIndex;
                                orderedBoxs[firstPrimOffset + i] = primitives[index];
                            }
                            node.InitLeaf(firstPrimOffset, primitives.size(), boxSum);//每个叶子节点存储第一个地址和prim数量以及包围盒
                            return node;
                        }
                    }
                    break;
                    }
                case Middle:
                {
                    Scalarm middle = (centerBox.min[dim] + centerBox.max[dim]) * 0.5;
                    auto midIter = std::partition(primitives.begin(), primitives.end(), [dim, middle](const BVHPrimitive& pi)
                        {
                            return pi.centroid[dim] < middle;
                        });
                    mid = midIter - primitives.begin();
                    if(midIter!= primitives.begin()&& midIter != primitives.end())
                        break;
                }
                case EqualCounts:
                {
                    Scalarm middle = (centerBox.min[dim] + centerBox.max[dim]) * 0.5;
                    std::nth_element(primitives.begin(), primitives.begin() + middle,
                        primitives.end(),
                        [dim](const BVHPrimitive& a, const BVHPrimitive& b) {
                            return a.centroid[dim] < b.centroid[dim];
                        });//部分排序
                    break;
                }
                default:
                    break;
                }

                BVHNode* children[2];
                if (primitives.size() > 1024 * 128)
                {
                    //并行。。。
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
        /*for (int i = 0; i < 3; ++i)
        {
            std::sort(primitives.begin(), primitives.end(), 
                [i](const BVHPrimitive& a, const BVHPrimitive& b) 
                {
                    return a.centroid[i] < b.centroid[i];
                });
        }*/

        return node;
    }

    void flattenBVH(BVHNode* rootNode, int* offset)
    {

    }
    unsigned int calcSplit(BVHNode* node)
    {

    }

private:
    int maxPrimitiveNode;
    std::vector<Primitive> _primitives;
    BVHNode* nodes = nullptr;
    SplitMethod _splitType;
};
