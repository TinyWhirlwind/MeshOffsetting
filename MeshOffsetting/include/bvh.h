#include "mymesh.h"
#include"AABB.h"
#include <atomic>
#include <memory>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <list>
#include <memory_resource>
#include <span>
class AABBBox;
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

        void InitLeaf(int first, int n, const AABBBox& b)
        {
            firstPrimOffset = first;
            childCount = n;//?
            bound = b;
            child[0] = nullPtr;
            child[1] = nullPtr;
        }
        void InitInterior(unsigned int axis, BVHNode* l, BVHNode* cr)
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
        std::vector<BVHPrimitive> bvhPrimitives(_primitives.size());
        for (int i = 0;  i< _primitives.size(); ++i)
        {
            _primitives[i].GetBBox(bvhPrimitives[i].bound);
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
    ~BVH();



    void Set(const std::vector<CFaceO>::iterator& _oBegin,
        const std::vector<CFaceO>::iterator& _oEnd,
        int size)
    {
        if (_oBegin > _oEnd)
            return;
        _bvhPrimitives.clear();
        BVHNode& curNode = _nodeList.back();
        for (auto i = _oBegin; i < _oEnd; ++i)
        {
            AABBBox box;
            box.Add(i->P(0));
            box.Add(i->P(1));
            box.Add(i->P(2));
            _bvhPrimitives.push_back(box);

            curNode.bound.Add(i->P(0));
            curNode.bound.Add(i->P(1));
            curNode.bound.Add(i->P(2));
        }
        curNode.childCount = _oEnd - _oBegin + 1;
        //create tree
        /*if (curNode.childCount > 1)
        {
            int _splitIndex = calcSplit(curNode,int)
        }*/
    }


private:
    BVHNode* buildHLBVH()
    {
    }

    BVHNode* buildBVH(std::span<BVHPrimitive> primitives, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimsOffset, std::vector<Primitive> orderedBoxs)
    {
        BVHNode* rootNode = new BVHNode();
        rootNode->bound=
        ++totalNodes;
        for (auto& it : primitives)
        {
            rootNode->bound.Merge(it->bound);
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
