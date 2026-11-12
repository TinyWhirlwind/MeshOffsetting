# Fast Intersection-free Offset Surface Generation from Freeform Models with Triangular Meshes

利用mesh的包围盒生成Grid

![](MeshOffsetting\images\92ae510665c771c919b1b24dcd59d4b437cb549c.png)

*A. Swept Sphere Volume Hierarchy (SSVH) Filter*

利用扫掠球包围体加速得到最近三角形，“[7区域参数法](https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf)”计算查询点到三角形的最近距离；
*B. Bounding Box Filter*

将Grid每m个grid cell合成一个block，SSVH计算block中心点到网格的无符号距离场剔除不在窄带内的block；
![](MeshOffsetting\images\67c5338de1a2f46ce137171bc1be27e69b209572.png)

*C. Signed Distance Filter*

利用符号距离从 Bounding Box Filter 留下的候选 Block 中，删除位于错误偏移方向的 Block；
![](MeshOffsetting\images\d8f93a79f1946360e5514660dc57ee644d9f7e59.png)

*D. Octree Filter*

细分剩余 Block ，递归删除更小的无效区域

![](MeshOffsetting\images\1691c6ba43a3b69d6a6ca2806a360b9ceddff565.png)
