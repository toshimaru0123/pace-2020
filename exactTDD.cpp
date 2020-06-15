#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>

#include <set>
#include <utility>
#include <vector>

#include <chrono>
#include <stack>
#include <unordered_map>

#include <map>

using namespace std;

using Edge = std::pair<int, int>;
using Edges = std::vector<Edge>;
using Vertices = std::vector<int>;

// TreeDepthComponent.hpp
class TreeDepthComponent {
  private:
    int treeDepth;
    int root;
    Vertices vertices;
    std::vector<TreeDepthComponent *> children;

  public:
    TreeDepthComponent() {}
    TreeDepthComponent(int td, int r) : treeDepth(td), root(r) {}

    TreeDepthComponent(int td, int r, std::vector<TreeDepthComponent *> c)
        : treeDepth(td), root(r), children(c) {}

    int getTD() const { return treeDepth; }
    int getRoot() const { return root; }
    Vertices getVertices() const { return vertices; }
    std::vector<TreeDepthComponent *> getChildren() const { return children; }
};

// Graph.hpp
class Graph {
  protected:
    int const n;
    int const e;
    Vertices const vertices;
    Edges const edges;
    std::vector<std::vector<int>> adjList;
    std::vector<std::vector<bool>> adjMatrix;

  public:
    Graph(int n, int e, Vertices vertices, Edges edges)
        : n(n), e(e), vertices(vertices), edges(edges) {
        // construct adjacent list from input
        adjList.resize(n + 1);
        adjMatrix.resize(n + 1, std::vector<bool>(n + 1));
        for (auto const edge : edges) {
            int u = edge.first, v = edge.second;
            adjList[u].push_back(v);
            adjList[v].push_back(u);
            adjMatrix[u][v] = true;
            adjMatrix[v][u] = true;
        }
    }

    int numNodes() const { return n; }
    int numEdges() const { return e; }
    Vertices getVertices() const { return vertices; }
    Edge getEdge(int i) const { return edges[i]; }
    Edges getEdges() const { return edges; }
    std::vector<std::vector<int>> getAdjList() const { return adjList; }
    std::vector<int> getNeighbors(int v) const { return adjList[v]; }
    bool isEdge(int u, int v) const { return adjMatrix[u][v]; }
};

class CompressedGraph : public Graph {
  private:
    std::vector<std::vector<TreeDepthComponent>> shrankVertex;
    std::vector<int> shrankList;
    std::vector<int> td;
    std::vector<bool> isCliquePart;

  public:
    CompressedGraph(int n, int e, Vertices vertices, Edges edges,
                    std::vector<std::vector<TreeDepthComponent>> sv,
                    std::vector<int> sl, std::vector<int> td,
                    std::vector<bool> isCliquePart)
        : Graph(n, e, vertices, edges), shrankList(sl), td(td),
          isCliquePart(isCliquePart) {
        shrankVertex = sv;
    }

    std::vector<std::vector<TreeDepthComponent>> getShrankVertices() const {
        return shrankVertex;
    }
    std::vector<int> getShrankList() const { return shrankList; }
    bool isShrank(int i) const {
        return (shrankVertex[i].size() != 0) ? true : false;
    }
    std::vector<TreeDepthComponent> getShrank(int i) const {
        return shrankVertex[i];
    }
    int getTD(int i) const { return td[i]; }
};

// TreeDepthDecomposition.hpp
TreeDepthComponent computeStart(const CompressedGraph &G, Vertices V);

// Preprocess.hpp
bool isLength3(const Graph &G, int vertexNum, std::vector<int> neighbers);
CompressedGraph compressGraph(const Graph &G);

// TDClique.hpp
bool isClique(const CompressedGraph &G, std::vector<int> V);
TreeDepthComponent *computeTDonClique(const CompressedGraph &G,
                                      std::vector<int> V);

// Preprocess.cpp
bool isLength2(const Graph &G, int vertexNum, std::vector<int> neighbers) {
    // vertexNumの隣接点
    int neighber = G.getNeighbors(vertexNum)[0];

    // vertexNum（次数1）の隣接点の次数が2じゃないならfalse
    if (neighbers.size() != 2) {
        return false;
    }

    for (auto const &v : neighbers) {
        if (v == vertexNum)
            continue;

        auto const connect = G.getNeighbors(v);
        // neighberのvertexNumではない隣接点が2以下ならfalse
        if (connect.size() <= 2)
            return false;

        for (auto const &v2 : connect) {
            if (v2 == neighber)
                continue;

            auto const &v2_neighber = G.getNeighbors(v2);
            // neighber以外のvの隣接点の次数がどれか一つでも2以上ならスターではない
            if (v2_neighber.size() >= 2) {
                return true;
            }
        }
    }
    return false;
}

bool isLength3(const Graph &G, int vertexNum, std::vector<int> neighbers) {
    if (neighbers.size() == 2) {
        for (auto const &v : neighbers) {
            if (v == vertexNum)
                continue;

            auto const connect2 = G.getNeighbors(v);
            if (connect2.size() == 2) {
                return true;
            } else {
                return false;
            }
        }
    }
    return false;
}

CompressedGraph compressGraph(const Graph &G) {
    auto vertices = G.getVertices();
    auto adjList = G.getAdjList();
    Edges edges;
    std::vector<std::vector<TreeDepthComponent>> shrankVertices(G.numNodes() +
                                                                1);
    std::vector<int> shrankList;
    std::vector<bool> isCliquePart;
    std::vector<int> td(G.numNodes() + 1, 1);
    std::vector<bool> eraseList(G.numNodes() + 1);
    for (unsigned vertexNum = 1; vertexNum < adjList.size(); ++vertexNum) {
        if (eraseList[vertexNum])
            continue;

        // 次数1の頂点の接続先をチェック
        if (adjList[vertexNum].size() == 1) {
            std::vector<int> deletedVertices;
            int shrankVertex = 0;
            // vertexNumの隣接点の隣接点
            auto const connect = G.getNeighbors(adjList[vertexNum][0]);

            // 長さ2の時
            if (isLength2(G, vertexNum, connect)) {
                for (auto const &v : connect) {
                    if (v == vertexNum) {
                        continue;
                    }
                    shrankVertex = v;
                    Vertices::iterator itr_v =
                        std::find(vertices.begin(), vertices.end(), vertexNum);
                    if (itr_v != vertices.end()) {
                        vertices.erase(itr_v);
                    }

                    itr_v = std::find(vertices.begin(), vertices.end(),
                                      adjList[vertexNum][0]);
                    if (itr_v != vertices.end()) {
                        vertices.erase(itr_v);
                    }

                    eraseList[vertexNum] = true;
                    eraseList[adjList[vertexNum][0]] = true;
                    // shrankVertexが圧縮された頂点として登録されてないならListに登録する
                    if (shrankVertices[shrankVertex].size() == 0) {
                        shrankList.push_back(shrankVertex);
                    }
                    auto child = new TreeDepthComponent(1, vertexNum);
                    std::vector<TreeDepthComponent *> children(1, child);
                    auto parent = new TreeDepthComponent(
                        2, adjList[vertexNum][0], children);
                    shrankVertices[shrankVertex].push_back(*parent);
                    if (td[shrankVertex] < 3)
                        td[shrankVertex] = 3;
                }
                continue;
            }

            // 長さ3の時
            if (isLength3(G, vertexNum, connect)) {
                for (auto const &v : connect) {
                    if (v == vertexNum)
                        continue;
                    auto const connect2 = G.getNeighbors(v);
                    for (auto const &v2 : connect2) {
                        if (v2 == adjList[vertexNum][0])
                            continue;
                        shrankVertex = v2;
                        Vertices::iterator itr_v = std::find(
                            vertices.begin(), vertices.end(), vertexNum);
                        if (itr_v != vertices.end()) {
                            vertices.erase(itr_v);
                        }
                        // vertices.erase(vertexNum);

                        itr_v = std::find(vertices.begin(), vertices.end(),
                                          adjList[vertexNum][0]);
                        if (itr_v != vertices.end()) {
                            vertices.erase(itr_v);
                        }
                        // vertices.erase(adjList[vertexNum][0]);

                        itr_v = std::find(vertices.begin(), vertices.end(), v);
                        if (itr_v != vertices.end()) {
                            vertices.erase(itr_v);
                        }
                        // vertices.erase(v);

                        eraseList[vertexNum] = true;
                        eraseList[adjList[vertexNum][0]] = true;
                        eraseList[v] = true;

                        // shrankVertexが圧縮された頂点として登録されてないならListに登録する
                        if (shrankVertices[shrankVertex].size() == 0) {
                            shrankList.push_back(shrankVertex);
                        }
                        auto child1 = new TreeDepthComponent(1, vertexNum);
                        auto child2 = new TreeDepthComponent(1, v);
                        std::vector<TreeDepthComponent *> children;
                        children.push_back(child1);
                        children.push_back(child2);
                        auto tdc = TreeDepthComponent(2, adjList[vertexNum][0],
                                                      children);
                        shrankVertices[shrankVertex].push_back(tdc);
                        if (td[shrankVertex] < 3)
                            td[shrankVertex] = 3;
                    }
                }
                continue;
            }

            if (shrankVertices[adjList[vertexNum][0]].size() != 0) {
                Vertices::iterator itr_v =
                    std::find(vertices.begin(), vertices.end(), vertexNum);
                if (itr_v != vertices.end()) {
                    vertices.erase(itr_v);
                }
                // vertices.erase(vertexNum);
                eraseList[vertexNum] = true;
                shrankVertices[adjList[vertexNum][0]].push_back(
                    TreeDepthComponent(1, vertexNum));
                continue;
            }

            // pathの長さが1の時
            Vertices::iterator itr_v =
                std::find(vertices.begin(), vertices.end(), vertexNum);
            if (itr_v != vertices.end()) {
                vertices.erase(itr_v);
            }
            // vertices.erase(vertexNum);
            eraseList[vertexNum] = true;
            shrankVertex = adjList[vertexNum][0];
            shrankVertices[shrankVertex].push_back(
                TreeDepthComponent(1, vertexNum));
            shrankList.push_back(shrankVertex);
            td[shrankVertex] = 2;
        }
    }
    // TODO: Simplicial vertexの削除

    for (auto const e : G.getEdges()) {
        if (eraseList[e.first] || eraseList[e.second])
            continue;
        edges.push_back(e);
    }
    // std::sort(edges.begin(), edges.end());

    CompressedGraph G_(G.numNodes(), G.numEdges(), vertices, edges,
                       shrankVertices, shrankList, td, isCliquePart);

    return G_;
}

// TDClique.cpp
bool isClique(const CompressedGraph &G, const std::vector<int> V) {
    for (size_t i = 0; i < V.size(); ++i) {
        for (size_t j = i + 1; j < V.size(); ++j) {
            if (!G.isEdge(V[i], V[j])) {
                return false;
            }
        }
    }
    return true;
}

bool pairSecondCompare(const std::pair<int, int> &lhs,
                       const std::pair<int, int> &rhs) {
    // 昇順で比較
    return lhs.second < rhs.second;
}

// Assumption: G[V] is clique.
TreeDepthComponent *computeTDonClique(const CompressedGraph &G,
                                      const std::vector<int> V) {
    // Vの各頂点ごとのtdを格納するpair型の配列 (index, td)
    std::vector<std::pair<int, int>> tmp_tds(V.size());

    // 各頂点vのindexとtdのpairを格納
    for (size_t i = 0; i < V.size(); ++i) {
        // std::cout << V[i] << ": td= " << G.getTD(V[i]) << std::endl;

        tmp_tds[i] = std::make_pair(i, G.getTD(V[i]));
    }

    // 各頂点vのtdを基準に昇順にソート
    std::sort(tmp_tds.begin(), tmp_tds.end(), pairSecondCompare);

    // 一つ前のTDCを保存用
    TreeDepthComponent *pre_tdc = nullptr;

    for (size_t i = 0; i < V.size(); ++i) {
        int index = tmp_tds[i].first;

        if (i == 0) {
            pre_tdc = new TreeDepthComponent(1, V[index]);
            continue;
        }

        // vector<TreeDepthComponet *>に格納
        std::vector<TreeDepthComponent *> children{pre_tdc};
        pre_tdc = new TreeDepthComponent(i + 1, V[index], children);
    }

    // forを抜けた時点でpre_tdcにはV-cliqueのTDCが保存されている．
    return pre_tdc;
}

// TDTree.cpp
class Node {
  public:
    int vertex;
    int parent;                // a pointer to parent node
    std::vector<int> children; // pointers to children

    Node(){};
    Node(int v, int parent) : vertex(v), parent(parent){};
    void setParent(const int v) { parent = v; };
    void addChild(const int v) { children.push_back(v); };
};

void printRootedTree(std::vector<Node> T) {
    std::cout << "Print rooted tree" << std::endl;

    for (size_t i = 0; i < T.size(); i++) {
        std::cout << i << "(" << T[i].vertex << ") ";
        if (!T[i].children.empty()) {
            for (const auto &child : T[i].children) {
                std::cout << child << " ";
            }
        }
        std::cout << std::endl;
    }
}

bool isConnected(const CompressedGraph &G, std::vector<int> V) {
    // DEBUG: ここG.numNodes() + 1 じゃないのん?
    std::vector<bool> visited(G.numNodes() + 1, true);
    for (auto &v : V) {
        visited[v] = false;
    }

    std::stack<int> s;
    s.push(V[0]);
    visited[V[0]] = true;
    while (!s.empty()) { // G[V] で深さ優先探索 （処理が変かも）
        int v = s.top();
        s.pop();
        const auto &Nv = G.getNeighbors(v);
        for (const auto &w : Nv) {
            if (!visited[w]) {
                visited[w] = true;
                s.push(w);
            }
        }
    }

    for (auto &v : V) {
        if (!visited[v])
            return false; // 未到達の点あり
    }
    return true;
}

/*
// counting the number of edges in G[V].
int numEdges(const CompressedGraph &G, std::vector<int> V) {
    int num = 0;
    for (const auto &v : V) {
        for (const auto &w : V) {
            if (G.isEdge(v, w))
                num++;
        }
    }
    return num / 2;
}
bool isTree(const CompressedGraph &G, std::vector<int> V) {
    if (isConnected(G, V) && numEdges(G, V) == V.size() - 1)
        return true;
    else
        return false;
}
*/

std::vector<Node> computeRootedTree(const CompressedGraph &G,
                                    std::vector<int> V) {
    std::vector<Node> rt;
    std::map<int, int> vertex2node;

    std::vector<bool> visited(G.numNodes() + 1, true);
    for (auto &v : V) {
        visited[v] = false;
    }

    std::queue<int> q;
    q.push(V[0]);
    visited[V[0]] = true;
    Node root(V[0], 0);
    vertex2node.insert(std::make_pair(V[0], 0));
    rt.push_back(root);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        int parent = vertex2node[v];

        const auto Nv = G.getNeighbors(v);
        for (const auto &w : Nv) {
            if (!visited[w]) {
                visited[w] = true;
                q.push(w);

                rt[parent].addChild(rt.size());
                vertex2node.insert(std::make_pair(w, rt.size()));
                Node node(w, parent);
                rt.push_back(node);
            }
        }
    }
    // printRootedTree(rt);
    return rt;
}

std::pair<std::vector<int>, std::vector<int>> bfs(std::vector<Node> RootedTree,
                                                  int s) {
    // initialize
    std::vector<int> ordering;
    std::vector<int> prev(RootedTree.size());
    std::vector<bool> visited(RootedTree.size(), false);
    std::queue<int> q;
    q.push(s);
    visited[s] = true;
    ordering.push_back(s);
    prev[s] = -1;

    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (const auto &w : RootedTree[v].children) {
            if (!visited[w]) {
                visited[w] = true;
                q.push(w);
                ordering.push_back(w);
                prev[w] = v;
            }
        }
        int parent = RootedTree[v].parent;
        if (!visited[parent]) {
            visited[parent] = true;
            q.push(parent);
            ordering.push_back(parent);
            prev[parent] = v;
        }
    }

    return std::make_pair(ordering, prev);
}

int computeCenterOnTree(const std::vector<Node> &RootedTree) {
    // printRootedTree(RootedTree);
    std::pair<std::vector<int>, std::vector<int>> ordering1 =
        bfs(RootedTree, 0);
    std::pair<std::vector<int>, std::vector<int>> ordering2 =
        bfs(RootedTree, ordering1.first.back());

    // find the longest path from bfs[0] to bfs[n-1]
    std::vector<int> LP;
    int v = ordering2.first.back();
    LP.push_back(v);
    int prev = ordering2.second[v];
    while (prev != -1) {
        LP.push_back(prev);
        prev = ordering2.second[prev];
    }
    int pos_center = LP.size() / 2;

#ifdef DEBUG_TDTree
    cout << "Longest Path: ";
    for (const auto &v : LP) {
        cout << v << " ";
    }
    cout << endl;
    cout << "Center: " << LP[pos_center] << endl;
#endif

    return LP[pos_center];
}

std::vector<std::vector<Node>>
computeComponents(const std::vector<Node> &RootedTree, int center) {
    std::vector<std::vector<Node>> components;

    std::vector<bool> visited(RootedTree.size(), false);
    visited[center] = true;

    for (size_t i = 0; i < RootedTree.size(); i++) {
        if (visited[i])
            continue;

        std::vector<Node> rt;
        std::map<int, int> node2node;

        std::queue<int> q;
        q.push(i);
        visited[i] = true;
        Node node(RootedTree[i].vertex, 0);
        rt.push_back(node);
        node2node.insert(std::make_pair(i, 0));
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            int v_new_index = node2node[v];
            for (const auto &w : RootedTree[v].children) {
                if (!visited[w]) {
                    visited[w] = true;
                    q.push(w);

                    rt[v_new_index].addChild(rt.size());
                    node2node.insert(std::make_pair(w, rt.size()));
                    Node node(RootedTree[w].vertex, v_new_index);
                    rt.push_back(node);
                }
            }
            // for the parent of v
            const int parent = RootedTree[v].parent;
            if (!visited[parent]) {
                visited[parent] = true;
                q.push(parent);

                rt[v_new_index].addChild(rt.size());
                node2node.insert(std::make_pair(parent, rt.size()));
                Node node(RootedTree[parent].vertex, v_new_index);
                rt.push_back(node);
            }
        }
        components.push_back(rt);
    }

    return components;
}

TreeDepthComponent *computeTDonTree(const std::vector<Node> &RootedTree) {
    // std::cerr << "Computing TDonTree" << std::endl;

    // size 1 case // サイズ 1 の木を返して終了
    if (RootedTree.size() == 1) {
        return new TreeDepthComponent(1, RootedTree[0].vertex);
    }

    std::vector<TreeDepthComponent *> children;
    // compute center v (v is a node index of RootedTree)
    int center = computeCenterOnTree(RootedTree);

    // compute connected components by removing v
    std::vector<std::vector<Node>> Components =
        computeComponents(RootedTree, center);

    // for each component C, compute treedepth decomposition on tree,
    // recursively.
    int depth = 0;
    for (const auto &C : Components) {
#ifdef DEBUG_TDTree
        printRootedTree(C);
#endif
        TreeDepthComponent *tdc = computeTDonTree(C);
        children.push_back(tdc);
        if (depth < tdc->getTD())
            depth = tdc->getTD();
    }

    return new TreeDepthComponent(depth + 1, RootedTree[center].vertex,
                                  children);
}

// Assumption: G[V] is tree.
TreeDepthComponent *computeTDonTree(const CompressedGraph &G,
                                    std::vector<int> V) {
    std::vector<Node> RootedTree = computeRootedTree(G, V);

    // 根は RootedTree[0] とする．
    return computeTDonTree(RootedTree);
}

void printTDC(const TreeDepthComponent root) {
    std::cout << "Root id: " << root.getRoot() << std::endl;
    std::cout << "children: ";
    for (auto const &c : root.getChildren()) {
        std::cout << c->getRoot() << " ";
    }
    std::cout << std::endl;

    for (auto const &c : root.getChildren()) {
        printTDC(*c);
    }
}

// TreeDepthDecomposition.cpp
#define ACTIVE_TDTree
#define ACTIVE_TDClique
#define ACTIVE_PRUNE
#define ACTIVE_D_v

int currentOpt = 0;

class Hash { // ハッシュ関数オブジェクト
  public:
    size_t operator()(const Vertices &V) const {
        const int C = 997; // 適当な素数
        size_t t = 0;
        for (auto const e : V) {
            t = t * C + e;
        }
        return t;
    }
};

std::unordered_map<Vertices, TreeDepthComponent, Hash> table;
std::unordered_map<Vertices, TreeDepthComponent, Hash> untable;

void printTDCs_(const TreeDepthComponent root) {
    std::cout << "Root id: " << root.getRoot() << std::endl;
    std::cout << "children: ";
    for (auto const &c : root.getChildren()) {
        std::cout << c->getRoot() << " ";
    }
    std::cout << std::endl;

    for (auto const &c : root.getChildren()) {
        printTDCs_(*c);
    }
}

std::vector<Vertices> computeConnectedComponents(const Graph &G,
                                                 const Vertices &v) {
    std::vector<int> vertices(v.begin(), v.end());

    std::stack<int> S;
    std::vector<int> component(vertices.size(), 0);
    std::vector<bool> checked(vertices.size(), false);
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (checked[i])
            continue;

        S.push(i);
        while (!S.empty()) {
            int index = S.top();
            S.pop();
            checked[index] = true;
            component[index] = i;
            auto adjList = G.getAdjList();
            for (auto const neighbor : adjList[vertices[index]]) {
                // std::cout << "vertex: " << vertices[index]
                //<< ", Neighbor: " << neighbor << std::endl;
                auto itr =
                    std::find(vertices.begin(), vertices.end(), neighbor);
                if (itr == vertices.end())
                    continue;
                int indexNeighbor = std::distance(vertices.begin(), itr);
                if (!checked[indexNeighbor])
                    S.push(indexNeighbor);
            }
        }
    }

    std::vector<Vertices> tmp(vertices.size());
    for (size_t i = 0; i < component.size(); ++i) {
        tmp[component[i]].push_back(vertices[i]);
    }

    std::vector<Vertices> D;
    for (auto const &d : tmp) {
        if (!d.empty())
            D.push_back(d);
    }

    return D;
}

std::pair<unsigned, bool> numEdgesAndCheckDoTDtree(const CompressedGraph &G,
                                                   const std::vector<int> V) {
    int num = 0;
    bool doTDTree = true;

    for (size_t i = 0; i < V.size(); ++i) {
        if (G.isShrank(V[i]))
            doTDTree = false;
        for (size_t j = i + 1; j < V.size(); ++j) {
            if (G.isEdge(V[i], V[j]))
                num++;
        }
    }

    return std::make_pair(num, doTDTree);
}

// 連結成分の数で降順にソートするための比較関数
bool CompareDsize(const std::pair<int, std::vector<Vertices>> &lhs,
                  const std::pair<int, std::vector<Vertices>> &rhs) {
    // 降順で比較
    return lhs.second.size() > rhs.second.size();
}

bool pairSecondCompareDown(const std::pair<int, int> &lhs,
                           const std::pair<int, int> &rhs) {
    // 降順で比較
    return lhs.second > rhs.second;
}

std::vector<std::pair<int, std::vector<Vertices>>>
computeD_v(const CompressedGraph &G, Vertices &V) {
    std::vector<std::pair<int, std::vector<Vertices>>> D_v(V.size());
    // 連結成分の計算
    int tmp_index = 0;
    for (auto const &v : V) {
        Vertices tmp(V);
        Vertices::iterator itr_v = std::find(tmp.begin(), tmp.end(), v);
        if (itr_v != tmp.end()) {
            tmp.erase(itr_v);
        }
        std::vector<Vertices> D = computeConnectedComponents(G, tmp);
        D_v[tmp_index++] = std::make_pair(v, D);
    }
    // D_v[]はv以外の頂点の連結成分の数を基準に降順でソートされている．
    std::sort(D_v.begin(), D_v.end(), CompareDsize);

    return D_v;
}

TreeDepthComponent *computeOptimalTD(const CompressedGraph &G, Vertices V,
                                     int preV, int level) {
    if (table.find(V) != table.end()) {
        return &table[V];
    }

    // Vの辺の本数を事前に計算（TDTreeとTDCliqueで使う）
    auto numEdge_doTDTree = numEdgesAndCheckDoTDtree(G, V);
    unsigned numEdge = numEdge_doTDTree.first;
    bool doTDTree = numEdge_doTDTree.second;

#ifdef ACTIVE_TDTree
    // Check if G(V) is Tree.
    if (numEdge == V.size() - 1 && doTDTree) {
        // TreeのTreedepth Decompositionが返ってくる．
        auto treeTDC = computeTDonTree(G, V);
        table.insert({V, *treeTDC});
        return treeTDC;
    }
#endif

#ifdef ACTIVE_TDClique
    // Check if G(V) is Clique.
    if (numEdge == (V.size() * (V.size() - 1)) / 2 && V.size() > 2) {
        // CliqueのTreedepth Decompositionが返ってくる．
        auto cliqueTDC = computeTDonClique(G, V);
        table.insert({V, *cliqueTDC});
        return cliqueTDC;
    }
#endif

    int opt_td = G.numNodes() + 1, root = 0;
    std::vector<TreeDepthComponent *> children;

    // D_v[]はv以外の頂点の連結成分の数を基準に降順でソートされている．

    // 連結成分の計算
    std::vector<std::pair<int, std::vector<Vertices>>> D_v(computeD_v(G, V));

    // for (auto const &v : V) {
    for (auto const &D_v_i : D_v) {
        int v = D_v_i.first;
        // Prune unnecessary computing
        if (preV != 0 && preV > v) {
            continue;
        }

        // Compute connected component of G[V \ v]
        std::vector<Vertices> D(D_v_i.second);
        std::vector<TreeDepthComponent *> candidateChildren;
        int maxTD = 0;
        bool pruned = false;
        if (D.size() == 1) {
            auto tdc = computeOptimalTD(G, D[0], v, level + 1);
// Prune rooted v
#ifdef ACTIVE_PRUNE
            if (tdc->getTD() + level > currentOpt) {
                pruned = true;
                continue;
            }
#endif
            candidateChildren.push_back(tdc);
            maxTD = std::max(maxTD, tdc->getTD());
        } else {
            for (auto const C : D) {
                auto tdc = computeOptimalTD(G, C, 0, level + 1);
// Prune rooted v
#ifdef ACTIVE_PRUNE
                if (tdc->getTD() + level > currentOpt) {
                    pruned = true;
                    break;
                }
#endif
                candidateChildren.push_back(tdc);
                maxTD = std::max(maxTD, tdc->getTD());
            }
        }
#ifdef ACTIVE_PRUNE
        if (pruned)
            continue;
#endif

        // If v is shrank, To make a comparison between TD(v) and maxTD + 1
        if (G.isShrank(v)) {
            if (G.getTD(v) > maxTD + 1) {
                maxTD = G.getTD(v) - 1;
            }
        }

        // Update treedepth and root
        if (maxTD + 1 < opt_td) {
            root = v;
            opt_td = maxTD + 1;
            children.clear();
            for (auto const child : candidateChildren) {
                children.push_back(child);
            }
        }

        // Save temporary solution
        if (level == 1) {
            currentOpt = opt_td;
        }
    }

    auto opt = new TreeDepthComponent(opt_td, root, children);

    if (opt->getTD() < G.numNodes() + 1) {
        table.insert({V, *opt});
    }

    return opt;
}

TreeDepthComponent computeStart(const CompressedGraph &G, Vertices V) {
    for (int i = 1; i <= G.numNodes(); ++i) {
        Vertices v(1, i);
        // v.insert(i);

        TreeDepthComponent TDC(G.getTD(i), i);
        table.insert({v, TDC});
    }

    // 暫定解を適当な大きい値に
    currentOpt = G.numNodes() + 1;

    auto tmp = computeOptimalTD(G, V, 0, 1);
    return *tmp;
}

// main.cpp
void printGraph(Graph G) {
    std::cout << "|V| = " << G.numNodes() << ", |E| = " << G.numEdges()
              << std::endl;
    for (auto const e : G.getVertices()) {
        std::cout << e << std::endl;
    }

    for (auto const edge : G.getEdges()) {
        std::cout << edge.first << " " << edge.second << std::endl;
    }
}

void printOptTDC(const TreeDepthComponent root) {
    std::cout << "Root id: " << root.getRoot() << std::endl;
    std::cout << "children: ";
    for (auto const &c : root.getChildren()) {
        std::cout << c->getRoot() << " ";
    }
    std::cout << std::endl;

    for (auto const &c : root.getChildren()) {
        printOptTDC(*c);
    }
}

void calTDC(const CompressedGraph G, const TreeDepthComponent tdc, int parent,
            std::vector<int> &ans) {
    if (G.isShrank(tdc.getRoot())) {
        auto const tmps = G.getShrank(tdc.getRoot());
        for (auto const tmp : tmps) {
            if (tmp.getTD() == 1) {
                ans[tmp.getRoot()] = tdc.getRoot();
                continue;
            }
            if (tmp.getTD() == 2) {
                for (auto const c : tmp.getChildren()) {
                    ans[c->getRoot()] = tmp.getRoot();
                }
                ans[tmp.getRoot()] = tdc.getRoot();
                continue;
            }
        }
    }
    for (auto const &c : tdc.getChildren()) {
        ans[c->getRoot()] = parent;
    }

    for (auto const &c : tdc.getChildren()) {
        calTDC(G, *c, c->getRoot(), ans);
    }
}

void outputSolutionOnPaceFormat(const CompressedGraph G,
                                const TreeDepthComponent tdc) {
    std::cout << tdc.getTD() << std::endl;
    std::vector<int> ans(G.numNodes() + 1, -1);
    calTDC(G, tdc, tdc.getRoot(), ans);
    ans[tdc.getRoot()] = 0;
    for (int i = 1; i <= G.numNodes(); ++i) {
        std::cout << ans[i] << std::endl;
    }
}

/*
void outputSolutionOnPaceFormat(const CompressedGraph G,
                                const TreeDepthComponent tdc, char *file) {
    std::string fileName = file;
    fileName.erase(fileName.size() - 3);
    fileName = fileName + ".tree";
    std::ofstream ofs(fileName);
    std::cout << "Treedepth: " << tdc.getTD() << std::endl;
    ofs << tdc.getTD() << '\n';
    std::vector<int> ans(G.numNodes() + 1, -1);
    calTDC(G, tdc, tdc.getRoot(), ans);
    ans[tdc.getRoot()] = 0;
    for (int i = 1; i <= G.numNodes(); ++i) {
        std::cout << ans[i] << std::endl;
        ofs << ans[i] << '\n';
    }
}
*/

Graph inputGraph() {
    int n = 0, e = 0;
    std::set<int> vertices_;
    Edges edges;
    std::string line;

    while (std::getline(std::cin, line)) {
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ' ')) {
            if (item == "c") {
                break;
            }
            if (item == "p") {
                std::getline(ss, item, ' ');
                std::string numNode;
                std::string numEdge;
                std::getline(ss, numNode, ' ');
                std::getline(ss, numEdge, ' ');
                n = std::stoi(numNode);
                e = std::stoi(numEdge);
                break;
            }
            int u = std::stoi(item);
            std::getline(ss, item, ' ');
            int v = std::stoi(item);
            if (v < u)
                std::swap(u, v); // (u < v)
            Edge edge = std::make_pair(u, v);
            vertices_.insert(u);
            vertices_.insert(v);
            edges.push_back(edge);
            break;
        }
    }

    std::sort(edges.begin(), edges.end());
    std::vector<int> vertices(vertices_.begin(), vertices_.end());
    Graph G(n, e, vertices, edges);

    return G;
}

/*
Graph inputGraph(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Error: Input a graph file(.gr)." << std::endl;
        exit(-1);
    }
    std::ifstream ifs(argv[1]);
    if (!ifs) {
        std::cout << "Error: Failed to open file." << std::endl;
        exit(-1);
    }
    int n = 0, e = 0;
    std::set<int> vertices_;
    Edges edges;
    std::string line;
    while (std::getline(ifs, line)) {
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ' ')) {
            if (item == "c") {
                break;
            }
            if (item == "p") {
                std::getline(ss, item, ' ');
                std::string numNode;
                std::string numEdge;
                std::getline(ss, numNode, ' ');
                std::getline(ss, numEdge, ' ');
                n = std::stoi(numNode);
                e = std::stoi(numEdge);
                break;
            }
            int u = std::stoi(item);
            std::getline(ss, item, ' ');
            int v = std::stoi(item);
            if (v < u)
                std::swap(u, v); // (u < v)
            Edge edge = std::make_pair(u, v);
            vertices_.insert(u);
            vertices_.insert(v);
            edges.push_back(edge);
            break;
        }
    }
    // edge = (u,v)をuを基準に昇順にソート
    std::sort(edges.begin(), edges.end());
    std::vector<int> vertices(vertices_.begin(), vertices_.end());
    Graph G(n, e, vertices, edges);
    return G;
}
*/

int main(int argc, char **argv) {
    // Graph G = inputGraph(argc, argv);
    Graph G = inputGraph();
    // 前処理
    auto G_ = compressGraph(G);

    /*
    std::cout << "Input Graph:" << std::endl;
    printGraph(G);
    std::cout << "Compressed Graph:" << std::endl;
    printGraph(G_);
    */

    // auto start = std::chrono::system_clock::now();
    auto tdc = computeStart(G_, G_.getVertices());
    outputSolutionOnPaceFormat(G_, tdc);
    // outputSolutionOnPaceFormat(G_, tdc, argv[1]);
    // printOptTDC(tdc);
    // auto end = std::chrono::system_clock::now();
    // auto elapsed =
    // std::chrono::duration_cast<std::chrono::milliseconds>(end -
    // start).count();
    // std::cout << "treedepth: " << tdc.getTD() << std::endl;
    // std::cout << "Computational time(Recursion): ";
    // std::cout << elapsed << "ms" << std::endl;
}