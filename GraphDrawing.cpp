// C++11
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <queue>
#include <string.h>
#include <iostream>
#include <vector>

using namespace std;
typedef long long ll;

const int MAX_V = 1000;

const ll CYCLE_PER_SEC = 2400000000;
ll g_startCycle;
double TIME_LIMIT = 10.0;
const int FIELD_SIZE = 701;

unsigned long long xor128() {
    static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
    unsigned long long rt = (rx ^ (rx << 11));
    rx = ry;
    ry = rz;
    rz = rw;
    return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
}

unsigned long long int getCycle() {
    unsigned int low, high;
    __asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
    return ((unsigned long long int) low) | ((unsigned long long int) high << 32);
}

double getTime(unsigned long long int begin_cycle) {
    return (double) (getCycle() - begin_cycle) / CYCLE_PER_SEC;
}

inline double calcDist(int y1, int x1, int y2, int x2) {
    int dy = y1 - y2;
    int dx = x1 - x2;

    assert(sqrt(dy * dy + dx * dx) >= 1.0);
    return sqrt(dy * dy + dx * dx);
}

typedef struct Coord {
    int y;
    int x;

    Coord(int y = -1, int x = -1) {
        this->y = y;
        this->x = x;
    }
} COORD;

typedef struct Vertex {
    int y;
    int x;
    bool fixed;
    vector<int> neighbors;

    Vertex(int y = -1, int x = -1) {
        this->y = y;
        this->x = x;
        this->fixed = false;
    }
} VERTEX;

typedef struct Node {
    int from;
    int to;
    double diff;

    Node(int from = -1, int to = -1, double diff = -1.0) {
        this->from = from;
        this->to = to;
        this->diff = diff;
    }

    bool operator>(const Node &n) const {
        return diff < n.diff;
    }
} NODE;

vector <VERTEX> g_vertexList;
bool g_field[FIELD_SIZE + 1][FIELD_SIZE + 1];
int g_distField[MAX_V + 1][MAX_V + 1];
int g_desireDist[MAX_V][MAX_V];
int g_edgeCount;
int g_vertexCount;

typedef struct Edge {
    int from;
    int to;
    int desire;
    int dist;

    Edge(int from = -1, int to = -1, int desire = -1) {
        this->from = from;
        this->to = to;
        this->desire = desire;
    }

    double diff() {
        return dist / (double) desire;
    }

} EDGE;

vector <Edge> g_edgeList;

class GraphDrawing {
public:
    void init(int N, vector<int> edges) {
        g_vertexCount = N;
        g_edgeCount = edges.size() / 3;
        g_startCycle = getCycle();

        memset(g_desireDist, -1, sizeof(g_desireDist));

        // initialize vertex
        for (int i = 0; i < N; ++i) {
            int y = xor128() % FIELD_SIZE;
            int x = xor128() % FIELD_SIZE;
            g_vertexList.push_back(VERTEX(y, x));
        }

        for (int i = 0; i < g_edgeCount; i++) {
            int from = edges[i * 3];
            int to = edges[i * 3 + 1];
            int dist = edges[i * 3 + 2];

            g_desireDist[from][to] = dist;
            g_desireDist[to][from] = dist;

            g_vertexList[from].neighbors.push_back(to);
            g_vertexList[to].neighbors.push_back(from);
            g_edgeList.push_back(Edge(from, to, dist));
        }

        fprintf(stderr, "VertexCount = %d, EdgeCount = %d\n", g_vertexCount, g_edgeCount);
    }

    void initFieldData() {
        memset(g_field, false, sizeof(g_field));
        memset(g_distField, -1, sizeof(g_distField));

        for (int i = 0; i < g_vertexCount; i++) {
            VERTEX *v = getVertex(i);
            g_field[v->y][v->x] = true;

            int nsize = v->neighbors.size();
            for (int j = 0; j < nsize; j++) {
                int neighborId = v->neighbors[j];
                VERTEX *n = getVertex(neighborId);
                int dist = calcDist(v->y, v->x, n->y, n->x);

                g_distField[i][neighborId] = dist;
                g_distField[neighborId][i] = dist;
            }
        }
    }

    vector<int> plot(int N, vector<int> edges) {
        init(N, edges);

        vector <Vertex> bvl = getBestVertexList();

        vector<int> ret = createAnswer(bvl);
        double score = calcScore();

        return ret;
    }

    vector <Vertex> getBestVertexList() {
        fprintf(stderr, "getBestVertexList =>\n");
        initFieldData();

        ll tryCount = 0;
        //queue <Node> pque;
        priority_queue <Node, vector<Node>, greater<Node>> pque;

        for (int i = 0; i < g_edgeCount; i++) {
            EDGE *edge = getEdge(i);
            int dist = g_distField[edge->from][edge->to];
            double diff = fabs(1.0 - dist / (double) g_desireDist[edge->from][edge->to]);
            pque.push(Node(edge->from, edge->to, diff));
        }

        while (!pque.empty()) {
            //NODE node = pque.front();
            NODE node = pque.top();
            pque.pop();
            int dist = g_distField[node.from][node.to];

            VERTEX *v = getVertex(node.from);
            VERTEX *n = getVertex(node.to);

            if (v->fixed) continue;
            fprintf(stderr, "%d -> %d : (%d/%d) = %f\n", node.from, node.to, dist, g_desireDist[node.from][node.to], node.diff);
            double bestDiff = INT_MAX;
            COORD bestCoord;

            for (int i = 0; i < 1000; i++) {
                COORD coord = createRandomCoord();
                int dist = calcDist(coord.y, coord.x, n->y, n->x);
                double diff = fabs(1.0 - dist / (double) g_desireDist[node.from][node.to]);

                if (bestDiff > diff) {
                    bestDiff = diff;
                    bestCoord = coord;
                }
            }

            if (bestDiff <= 0.1) {
                v->fixed = true;
                n->fixed = true;
            }
            tryCount++;
            moveVertex(node.from, bestCoord.y, bestCoord.x);
            fprintf(stderr, "%d -> %d : (%d/%d) = %f\n", node.from, node.to, g_distField[node.from][node.to], g_desireDist[node.from][node.to], bestDiff);

            int nsize = v->neighbors.size();
            for (int i = 0; i < nsize; i++) {
                int neighborId = v->neighbors[i];
                int dist = g_distField[node.from][neighborId];
                double diff = fabs(1.0 - dist / (double) g_desireDist[node.from][neighborId]);

                pque.push(Node(node.from, neighborId, diff));
            }

            if (tryCount % 1 == 0) {
                double minDiff = DBL_MAX;
                double maxDiff = 0.0;
                priority_queue <Node, vector<Node>, greater<Node>> tque;
                for (int i = 0; i < g_edgeCount; i++) {
                    EDGE *edge = getEdge(i);
                    int dist = g_distField[edge->from][edge->to];
                    double diff = fabs(1.0 - dist / (double) g_desireDist[edge->from][edge->to]);
                    double sdiff = dist / (double) g_desireDist[edge->from][edge->to];
                    tque.push(Node(edge->from, edge->to, sdiff));
                    minDiff = min(minDiff, sdiff);
                    maxDiff = max(maxDiff, sdiff);
                }
                fprintf(stderr,"(%f/%f) = score = %f\n", minDiff, maxDiff, minDiff / maxDiff);

                pque = tque;
            }

            if (tryCount > 1000) {
                break;
            }
        }

        fprintf(stderr, "tryCount = %lld\n", tryCount);
        return g_vertexList;
    }

    void moveVertex(int id, int y, int x) {
        VERTEX *v = getVertex(id);
        g_field[v->y][v->x] = false;
        v->y = y;
        v->x = x;
        g_field[y][x] = true;

        int nsize = v->neighbors.size();
        for (int i = 0; i < nsize; i++) {
            int neighborId = v->neighbors[i];
            VERTEX *n = getVertex(neighborId);
            int dist = calcDist(v->y, v->x, n->y, n->x);

            g_distField[id][neighborId] = dist;
            g_distField[neighborId][id] = dist;
        }
    }

    COORD createRandomCoord() {
        int y, x;

        do {
            y = xor128() % FIELD_SIZE;
            x = xor128() % FIELD_SIZE;
        } while (g_field[y][x]);

        return COORD(y, x);
    }

    vector<int> createAnswer(vector <Vertex> &vertexList) {
        vector<int> ret;

        for (int i = 0; i < g_vertexCount; i++) {
            Vertex v = vertexList[i];
            ret.push_back(v.y);
            ret.push_back(v.x);
        }

        return ret;
    }

    double calcScore() {
        double minDiff = DBL_MAX;
        double maxDiff = 0.0;

        for (int i = 0; i < g_edgeCount; i++) {
            EDGE *edge = getEdge(i);
            int dist = g_distField[edge->from][edge->to];
            double diff = dist / (double) g_desireDist[edge->from][edge->to];

            minDiff = min(minDiff, diff);
            maxDiff = max(maxDiff, diff);
        }

        return minDiff / maxDiff;
    }

    EDGE *getEdge(int id) {
        return &g_edgeList[id];
    }

    VERTEX *getVertex(int id) {
        return &g_vertexList[id];
    }

};


// -------8<------- end of solution submitted to the website -------8<-------
template<class T>
void getVector(vector <T> &v) { for (int i = 0; i < v.size(); ++i) cin >> v[i]; }

int main() {
    TIME_LIMIT = 1.0;
    GraphDrawing gd;
    int N;
    cin >> N;
    int E;
    cin >> E;
    vector<int> edges(E);
    getVector(edges);
    vector<int> ret = gd.plot(N, edges);
    cout << ret.size() << endl;
    for (int i = 0; i < (int) ret.size(); ++i) { cout << ret[i] << endl; }
    cout.flush();
}
