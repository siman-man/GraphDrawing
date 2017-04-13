## GraphDrawing

### Problem Statement

In this problem you will be drawing a graph on a plane. Each vertex will be drawn as a point, 

この問題では平面上にグラフを書いてもらいます。 頂点は点として書かれ、

and each edge will be drawn as a line segment between two vertices. Each edge has a desired length 

辺は各頂点をつなげたものです。 各辺は希望の長さを持っています。

- the length it should have on the drawing. Your task is to assign to each vertex integer coordinates from 

( この長さは図面上のものです ) あなたの目的は平面上に点を配置し、各辺の長さを与えられた長さに近づけることです。

a limited area so that the ratio of the actual edge length to the desired length is as similar across all 
edges as possible.

### Implementation

Your code must implement one method plot(int NV, vector <int> edges):

あなたは `plot(int NV, vector <int> edges)` を実装する必要があります。

NV is the number of vertices in the graph.

NVはgraphに含まれる頂点の数です。

edges describes the desired lengths of the edges in the graph: the i-th edge connects vertices edges[3*i] and edges[3*i+1] and has the desired length of edges[3*i+2].

辺はgraph上の長さを含んでいます。i番目の辺情報には [3*i] - [3*i+1] に点情報、[3*i+2] に長さの情報が含まれています。

Each vertex present in the graph is used in at least one edge. No pair of vertices is connected by two different edges.

各グラフの頂点は、かならず1つの辺に利用されています。 二重辺の存在はありません。

The return from this method will describe the assigned coordinates of the graph vertices. 

返り値として頂点の座標のリストを返します。

It must contain exactly 2*NV elements and give the coordinates of the j-th vertex as (ret[2*j], ret[2*j+1]). 

リストは 2*NV のサイズで、それぞれ [2*j], [2*j+1] の形式で座標が含まれています。

Each coordinate of each vertex must be between 0 and 700, inclusive. All vertices must be placed in distinct points.

各頂点の座標は [0 - 700] の間であり、同じ座標に複数の頂点が来ることはありません。



### スコア

For each test case we will calculate your raw score. If your solution produced an invalid return (contained invalid number of elements,

placed two vertices at the same point etc.), raw score for this test case will be 0. Otherwise, raw score will be calculated as follows. 

For each edge of the graph, its actual length on the drawing and the ratio "actual length / desired length" are calculated. 

Then minimal and maximal ratios across all edges are selected, and the raw score for a test case is calculated as MIN_RATIO/MAX_RATIO. 



### テストケースの生成

To generate the desired lengths of the graph edges, we will use the following procedure:

テストケースは以下の手順で生成されます。

* Assign random integer coordinates from [0, 700] x [0, 700] to each vertex of the graph, so that the coordinates of all vertices are distinct.

[0, 700] の範囲でランダムに座標を生成します (各頂点が重複しないように)

* Calculate actual lengths of edges based on these coordinates.

各頂点同士の長さを計算します。

* Pick a random distortion percentage P between 0 and 99, inclusive.

ねじれ係数を [0, 99] の間で選択します。

* For each edge, multiply its length by 1.1^G with probability P (and keep the length unchanged with probability 1-P). G is sampled from normal distribution with mean 0.0 and standard deviation 1.0. The new lengths are capped to be between 1 and 991, inclusive.



### 注意点

* The time limit is 10 seconds per test case (this includes only the time spent in your code). The memory limit is 1024 megabytes.

制限時間は10秒です。メモリの制限は1GBです。

* There is no explicit code size limit. The implicit source code size limit is around 1 MB (it is not advisable to submit codes of size close to that or larger). Once your code is compiled, the binary size should not exceed 1 MB.

コードサイズの制限は100MBとします。

* The compilation time limit is 30 seconds. You can find information about compilers that we use and compilation options here.

コンパイルタイムの制限時間は30秒です。

* There are 10 example test cases and 100 full submission (provisional) test cases.

テストケースは 10、本番ケースは 100です。

* The match is rated.

この試合はレートがつきます。


### 制限

* The number of vertices NV will be between 10 and 1000, inclusive.

頂点の数は [10 - 1000] です。

* The number of edges will be between NV - 1 and NV * min(10, (NV - 1) / 4), inclusive.

辺の数は [NV-1, NV * (min, (NV - 1) / 4)] です。


## 考察

* 辺情報から点の座標を予測する問題。
* 辺の長いものほど情報が多い (該当する部分が限られてくるため)
* 長さが一致していればいいので、座標が90 * nで回転していてもOK