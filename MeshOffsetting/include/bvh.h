#include "mymesh.h"
#include"AABB.h"
#include <span>
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
            child[0] = nullPtr;
            child[1] = nullPtr;
        }
        void InitInterior(unsigned int axis, BVHNode* cl, BVHNode* cr)
        {
            childNode[0] = cl;
            childNode[1] = cr;
            splitAxis = axis;
            childCount = 0;
            bound = cr->bound.Merge(cl->bound)
        }

    };

    BVH(const CMeshO& mesh, int maxPrimitiveCount, SplitMethod type) :
        _primitives(std::move(mesh.face)), maxPrimitiveCount(std::min(255, maxPrimitiveCount)), _splitType(type)
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
        BVHNode rootNode;
        AABBBox box;
        for (auto& it : primitives)
        {
            box.Merge(it->bound);
        }

        float bestCost = primitives.size() * 1.0f;
        float rootSA = rootNode->bound.SurfaceArea();

        if (rootSA < epsilon || primitives.size() == 1)
        {
            int firstPrimOffset = orderedPrimsOffset->fetch_add(primitives.size());
            for (int i = 0; i < primitives.size(); ++i)
            {
                int index = primitives[i].primitiveIndex;
                orderedBoxs[firstPrimOffset + i] = primitives[index];
            }
            rootNode.InitLeaf(firstPrimOffset, primitives.size(), box);
            return rootNode;
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
            if (centerBox.min[dim] == centerBox.max[dim])//重心分布在各个
            {

            }
        }
        for (int i = 0; i < 3; ++i)
        {
            std::sort(primitives.begin(), primitives.end(), 
                [i](const BVHPrimitive& a, const BVHPrimitive& b) 
                {
                    return a.centroid[i] < b.centroid[i];
                });
        }


    }

    void flattenBVH(BVHNode* rootNode, int* offset)
    {

    }
    unsigned int calcSplit(BVHNode* node)
    {

    }

private:
    int maxPrimitiveCount;
    std::vector<Primitive> _primitives;
    BVHNode* nodes = nullptr;

    std::vector<BVHNode> _nodeList;
    SplitMethod _splitType;
};
