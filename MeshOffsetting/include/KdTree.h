#include<vector>
#include "mymesh.h"
class KdTree
{
public:
    class Node
    {
    public:
        MESHLAB_SCALAR splitValue;
        int leftChildIndex;
        unsigned int dim = 2;
        unsigned int leaf = 1;
        std::vector<CFaceO*> objectList;
        Box3m box;
    };

public:
    KdTree(CMeshO& mesh, unsigned int maxDepth, unsigned int maxObjectNum) : maxDepth(64), maxObjectNum(64)
    {
        epsilon = std::numeric_limits<float>::epsilon();
        auto a = mesh.face.begin();
        Set(mesh.face.begin(), mesh.face.end(), mesh.face.size());
    }
    ~KdTree();

    void Set(const std::vector<CFaceO>::iterator& _oBegin,
        const std::vector<CFaceO>::iterator& _oEnd,
        int size, bool isSelect)
    {
        nodeList.resize(1);
        Node& node = nodeList.back();
        node.leaf = 0;
        node.box.Offset(Point3m(epsilon, epsilon, epsilon));
        Box3m box;
        if (isSelect)
        {
            for (auto i = _oBegin; i != _oEnd; ++i)
            {
                if (i->IsS())
                {
                    node.objectList.push_back(&(*i));
                    box.Add(i->P(0));
                    box.Add(i->P(1));
                    box.Add(i->P(2));
                }
            }
        }
        else
        {
            for (auto i = _oBegin; i != _oEnd; ++i)
            {
                node.objectList.push_back(&(*i));
                box.Add(i->P(0));
                box.Add(i->P(1));
                box.Add(i->P(2));
            }
        }

        node.box = box;
        numLevel = Create(0, 1);
    }

    void Clear()
    {
        for (auto& iter : nodeList)
        {
            iter.objectList.clear();
        }
        nodeList.clear();
    }

    
protected:
    int Create(unsigned int nodeIndex, unsigned int level)
    {
        Node& node = nodeList[nodeIndex];
        Point3m diag = node.box.max - node.box.min;
        unsigned int dim = diag.X() > diag.Y() ? (diag.X() > diag.Z() ? 0 : 2) : (diag.Y() > diag.Z() ? 1 : 2);

        node.splitValue = (node.box.min[dim] + node.box.max[dim]) * 0.5f;
        node.dim = dim;


    }
 
private:
    std::vector<Node> nodeList;
    unsigned int numLevel;
    MESHLAB_SCALAR epsilon;
    unsigned int maxDepth;
    unsigned int maxObjectNum;

};