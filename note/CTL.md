# Customizable Tree Labeling

Preprocess:
* calculate separator decomposition:算sd
Customization:
* CCH preprocess: 得到图的边集合
* tree hierarchy proprocess: 得到side id和具体每个点hub大小
* Truncated Tree Labeling
* buildCustomizedCTL

# Separator Decomposition
每个vertex有一个order，order小的先收缩

# CCH preprocess
* *rank*: rank[i] = vid, 点vid的rank是i，rank为i的点是vid
* *elimination_tree_parent*: elimination_tree_parent[v] = parent，顶点v的父节点是parent。parent是v被contract时，所有邻居里最后一个被contract的点(rank最大的点)
* *upGraph/downGraph*: 每个边存两遍，upGraph做forward search，downGraph做backward search
* *std::vector<int32_t> firstUpInputEdge*;   // The idx of the 1st upward input edge for each edge.
* *std::vector<int32_t> firstDownInputEdge*; // The idx of the 1st downward input edge for each edge.
* *std::vector<int32_t> upInputEdges*;       // The upward input edges.
* *std::vector<int32_t> downInputEdges*;     // The downward input edges.

# Tree Hierarchy Preprocess
* *std::vector<uint64_t> packedSideIds*: sideid
* *BitVector truncateVertex*: no hub label if true
* *std::vector<uint32_t> firstSepSizeSum*: 每个vertex的偏移值(从firstSepSizeSum[v]可以开始存v的sideid)
* *std::vector<uint32_t> sepSizeSumsOnBranch*: 每个vertex的hub数前缀和序列，最后一位为自己的hub数

## forEachSepDecompNodeInDfsOrder(sd, recurse, backtrack):
* visit all elements in the separator decomposition in pre-order
* recurse(parent, child): when visit child for the first time
* backtrack(child, parent): when finish traversing subtree child

## computeSepDecompDepth(sd)
* compute the depth of sd

## getSepDecompNodesToTruncateAt(sd, theta, maxUntruncatedDepth)
* mark the elements as truncated or not, 返回doTruncate这个bitvector
* numVerticesOnBranch.top(): 当前子树大小-自身所在element大小
* numVerticesChild:真正的子树大小
* 如果子树vertex数小于等于theta, truncated
  * doTruncate[child] = true 且 两个child的子树（如果存在）也肯定是true
* maxUntruncatedDepth: 最大没有truncated的element在sd上的高度,最小为1

## computeVertexLocationInSepDecomp(sd, truncatedAtSdNode, depthofVertex)
* truncatedAtSdNode: getSepDecompNodesToTruncateAt的返回值,由element 给element id返回truncated or not
### 作用: 计算depthOfVertex, packedSideIds, truncateVertex
* *packedSideIds*: 对于每个vertex，存一串“从根到其分离器层级的左右走向”的二进制位,左就0,右就1
  * depth: 当前node高度(sideid长度)
* *truncateVertex*: 标记每个vertex是否被truncated，truncatedAtSdNode里的是sd的node，每个node有多个vertex
* depthOfVertex: 计算有效层级数，从根到这个点到有效层深
* 第一次visit某个孩子时
  * 如果parent没有被truncated，用parent的sideid生成孩子的；如果被truncated，直接继承parent的sideid
  * 需interate所以child的vertices，并标记truncateVertex
* 离开child为root的子树时
  * 回复sideid到parent，如果没有被truncated，回复有效层级

## computeSepSizeSumsOnBranches(sd, truncatedAtSdNode, depthofVertex)
### 作用: 生成v的hub数及其在sd上ancestor的前缀
* *sepSizeSumsOnBranch[firstSepSizeSum[v+1]-1]*:v的hub数
* *firstSepSizeSum*: 每个vertex的偏移值(从firstSepSizeSum[v]可以开始存v的sideid)
* *sepSizeSumsOnBranch*: 每个vertex的hub数前缀和序列，最后一位为自己的hub数
### 方法
* 如果根都被truncated了，那大家前缀和都是0
* 否则,先把所有根上vertex的前缀和写入，sepSizeSumsOnCurBranch推入这个element的长度
* 第一次visit某个child时
  * 如果child被truncated了，对于element里每个vertex，其分离器大小前缀和和parent一样
  * 如果没有，iterate每个vertex，其前缀和谁parent前缀和+hub数(hub: parent里element数+同一element在v之前的vertex数)
* 离开时，回复前缀和


# Truncated Tree Labeling

## init(), 每个hub有up label和down label，每个需要空间K*numHubs，初始化空间和labelOffset

# CCH Metric Customization

# CTLMetric
## CCH Customization
* Perfect Customizaiion
  * buildMinimumWeightedCH: 先给原边加上metric，再计算customize metric
  * computeCustomizedMetricInParallel:
  ```cpp
    void computeCustomizedMetricInParallel() noexcept {
    //按elimination tree从底层到高层遍历vertex u
    cch.forEachVertexBottomUp([&](const int u) {
      //对u的每条向上边lower
      FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), u, lower) {
        //边指向的点为v
        const int v = cch.getUpwardGraph().edgeHead(lower);
        //遍历和边(u,v)构成三角点所有第三个点w, u->v->w, u->w
        //lower: u->v, inter: u->w, upper: v->w
        //三角形不等式去掉不合理边
        //up(v->w)< down(u->v)+up(u->w)
        //down(v->w)< down(u->w)+up(u->v)
        cch.forEachUpperTriangle(u, v, lower, [&](int, const int inter, const int upper) {
          atomicFetchMin(upWeights[upper], downWeights[lower] + upWeights[inter]);
          atomicFetchMin(downWeights[upper], downWeights[inter] + upWeights[lower]);
          return true;
        });
      }
    });
  }
  ```
  ## buildCustomizedCTL
  customizeLabelling
  * 对每个点u构建hub label，顺序是从sd高到低:
    * 如果truncated，不构建
  * numHubsU为u的hub数
  * *uUpLabel*： u的up hub label集合, 初始化hub中的distance集合startUUp为无穷大，hub中edge为invalid
  * iterate u在upGraph中的所有邻接点v,这条边(u->v)的权值为upWeight，v的up hub大小为numHubsV [可知v的hub小于u的hub且v和u的LCH为v]
  * v的up hub label集合为*vUpLabel*