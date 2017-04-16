// C++11
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cmath>
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

typedef struct Coord {
    int y;
    int x;

    Coord(int y, int x) {
        this->y = y;
        this->x = x;
    }
} COORD;

typedef struct Edge {
    int from;
    int to;
    int dist;

    Edge(int from = -1, int to = -1, int dist = -1) {
        this->from = from;
        this->to = to;
        this->dist = dist;
    }
} EDGE;

typedef struct Vertex {
    int y;
    int x;
    double minDiff;
    double maxDiff;
    vector<int> neighbors;

    Vertex(int y = -1, int x = -1, double minDiff = -1.0, double maxDiff = -1.0) {
        this->y = y;
        this->x = x;
        this->minDiff = minDiff;
        this->maxDiff = maxDiff;
    }

    double range() {
        assert(maxDiff - minDiff >= 0);
        return maxDiff - minDiff;
    }

} VERTEX;

struct Diff {
    double min;
    double max;
    double range;

    Diff(double minDiff = -1.0, double maxDiff = -1.0, double range = -1.0) {
        this->min = minDiff;
        this->max = maxDiff;
        this->range = range;
    }
};

struct Node {
    double totalRange;
    double maxDiff;
    double minDiff;

    Node(double totalRange = -1.0, double maxDiff = -1.0, double minDiff = -1.0) {
        this->totalRange = totalRange;
        this->maxDiff = maxDiff;
        this->minDiff = minDiff;
    }

    double value() {
        return 0.0;
    }

    double score() {
        return minDiff / maxDiff;
    }
};

vector <VERTEX> g_vertexList;
vector <Edge> g_edgeList;
bool g_field[FIELD_SIZE + 1][FIELD_SIZE + 1];
int g_desireDist[MAX_V][MAX_V];
int g_edgeCount;
int g_vertexCount;


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

    void initFieldData(vector<VERTEX> &vl) {
        memset(g_field, false, sizeof(g_field));

        for (int i = 0; i < vl.size(); i++) {
            VERTEX *v = &vl[i];
            Diff diff = calcVertexDiff(vl, i);
            v->minDiff = diff.min;
            v->maxDiff = diff.max;
            g_field[v->y][v->x] = true;
        }
    }

    vector<int> plot(int N, vector<int> edges) {
        init(N, edges);

        vector <Vertex> bvl = getBestVertexList(g_vertexList);

        vector<int> ret = createAnswer(bvl);
        Node node = calcScore(bvl);
        fprintf(stderr, "CurrentScore = %f\n", node.score());

        return ret;
    }

    vector <Vertex> getBestVertexList(vector <Vertex> vertexList) {
        vector <Vertex> bestVertexList = vertexList;
        initFieldData(bestVertexList);

        double currentTime = getTime(g_startCycle);
        Node node = calcScore(bestVertexList);
        double bestValue = node.totalRange;

        ll tryCount = 0;

        while (currentTime < TIME_LIMIT) {
            int vertexId = xor128() % g_vertexCount;
            Vertex ov = bestVertexList[vertexId];
            COORD c = createRandomCoord();

            assert(g_field[ov.y][ov.x]);
            g_field[ov.y][ov.x] = false;
            assert(!g_field[c.y][c.x]);
            g_field[c.y][c.x] = true;

            bestVertexList[vertexId].y = c.y;
            bestVertexList[vertexId].x = c.x;
            Diff d = calcVertexDiff(bestVertexList, vertexId);
            double diff = bestValue + (d.range - ov.range());

            if (bestValue > diff) {
                bestValue = diff;
                bestVertexList[vertexId].minDiff = d.min;
                bestVertexList[vertexId].maxDiff = d.max;
            } else {
                assert(!g_field[ov.y][ov.x]);
                g_field[ov.y][ov.x] = true;
                assert(g_field[c.y][c.x]);
                g_field[c.y][c.x] = false;
                bestVertexList[vertexId].y = ov.y;
                bestVertexList[vertexId].x = ov.x;
            }

            currentTime = getTime(g_startCycle);
            tryCount++;
        }

        fprintf(stderr, "tryCount = %lld\n", tryCount);
        return bestVertexList;
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

    Node calcScore(vector <Vertex> &vertexList) {
        double minDiff = DBL_MAX;
        double maxDiff = 0.0;
        double totalRange = 0.0;

        for (int i = 0; i < g_vertexCount; i++) {
            Vertex *v = &vertexList[i];
            minDiff = min(minDiff, v->minDiff);
            maxDiff = max(maxDiff, v->maxDiff);
            totalRange += v->range();
        }

        return Node(totalRange, maxDiff, minDiff);
    }

    Diff calcVertexDiff(vector <Vertex> &vertexList, int fromId) {
        double minDiff = DBL_MAX;
        double maxDiff = 0.0;
        Vertex from = vertexList[fromId];
        int nsize = from.neighbors.size();

        for (int j = 0; j < nsize; j++) {
            int toId = from.neighbors[j];
            Vertex to = vertexList[toId];
            double dist = calcDist(from, to);
            assert(g_desireDist[fromId][toId] > 0);
            double diff = dist / (double) g_desireDist[fromId][toId];
            minDiff = min(minDiff, diff);
            maxDiff = max(maxDiff, diff);
        }

        return Diff(minDiff, maxDiff, maxDiff - minDiff);
    }

    inline double calcDist(Vertex &v1, Vertex &v2) {
        int dy = v1.y - v2.y;
        int dx = v1.x - v2.x;

        return sqrt(dy * dy + dx * dx);
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
