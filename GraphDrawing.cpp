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

struct Edge {
    int from;
    int to;
    int dist;

    Edge(int from = -1, int to = -1, int dist = -1) {
        this->from = from;
        this->to = to;
        this->dist = dist;
    }
};

struct Vertex {
    int y;
    int x;

    Vertex(int y = -1, int x = -1) {
        this->y = y;
        this->x = x;
    }
};

vector <Vertex> g_vertexList;
vector<Edge> g_edgeList;
bool g_field[FIELD_SIZE+1][FIELD_SIZE+1];
int g_desireDist[MAX_V][MAX_V];
int g_edgeCount;
int g_vertexCount;


class GraphDrawing {
public:
    void init(int N, vector<int> edges) {
        g_vertexCount = N;
        g_edgeCount = edges.size() / 3;
        g_startCycle = getCycle();

        for (int i = 0; i < g_edgeCount; i++) {
            int from = edges[i * 3];
            int to = edges[i * 3 + 1];
            int dist = edges[i * 3 + 2];
            g_edgeList.push_back(Edge(from, to, dist));
        }

        memset(g_desireDist, -1, sizeof(g_desireDist));

        fprintf(stderr, "VertexCount = %d, EdgeCount = %d\n", g_vertexCount, g_edgeCount);
    }

    void initField(vector<Vertex> vl) {
        memset(g_field, false, sizeof(g_field));

        for (int i = 0; i < vl.size(); i++) {
            Vertex v = vl[i];
            g_field[v.y][v.x] = true;
        }
    }

    vector<int> plot(int N, vector<int> edges) {
        init(N, edges);

        for (int i = 0; i < N; ++i) {
            int y = xor128() % FIELD_SIZE;
            int x = xor128() % FIELD_SIZE;
            g_vertexList.push_back(Vertex(y, x));
        }

        vector <Vertex> bvl = getBestVertexList(g_vertexList);

        vector<int> ret = createAnswer(bvl);
        double score = calcScore(bvl);
        fprintf(stderr, "CurrentScore = %f\n", score);

        return ret;
    }

    vector <Vertex> getBestVertexList(vector <Vertex> vertexList) {
        vector <Vertex> bestVertexList = vertexList;
        vector <Vertex> goodVertexList = vertexList;
        double currentTime = getTime(g_startCycle);
        double bestScore = calcScore(bestVertexList);
        double goodScore = bestScore;
        double diffScore;

        initField(bestVertexList);

        int R = 1000000;
        double k = 0.1;
        ll tryCount = 0;

        while (currentTime < TIME_LIMIT) {
            double remainTime = TIME_LIMIT - currentTime;
            int i = xor128() % g_vertexCount;
            Vertex ov = goodVertexList[i];
            Vertex v = createRandomVertex();
            assert(g_field[ov.y][ov.x]);
            g_field[ov.y][ov.x] = false;
            assert(!g_field[v.y][v.x]);
            g_field[v.y][v.x] = true;
            goodVertexList[i] = v;
            double score = calcScore(goodVertexList);
            diffScore = bestScore - score;

            if (bestScore < score) {
                bestScore = score;
                bestVertexList = goodVertexList;
            }

            if (goodScore < score || (xor128() % R < R * exp(-diffScore / (k * remainTime)))) {
                goodScore = score;
            } else {
                assert(!g_field[ov.y][ov.x]);
                g_field[ov.y][ov.x] = true;
                assert(g_field[v.y][v.x]);
                g_field[v.y][v.x] = false;
                goodVertexList[i] = ov;
            }

            currentTime = getTime(g_startCycle);
            tryCount++;
        }

        fprintf(stderr, "tryCount = %lld\n", tryCount);
        return bestVertexList;
    }

    Vertex createRandomVertex() {
        int y, x;

        do {
            y = xor128() % FIELD_SIZE;
            x = xor128() % FIELD_SIZE;
        } while (g_field[y][x]);

        return Vertex(y, x);
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

    double calcScore(vector <Vertex> &vertexList) {
        double minDiff = DBL_MAX;
        double maxDiff = 0.0;

        for (int i = 0; i < g_edgeCount; i++) {
            Edge edge = g_edgeList[i];
            Vertex v1 = vertexList[edge.from];
            Vertex v2 = vertexList[edge.to];
            double dist = calcDist(v1, v2);

            double diff = dist / (double) edge.dist;
            minDiff = min(minDiff, diff);
            maxDiff = max(maxDiff, diff);
        }

        return minDiff / maxDiff;
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
    int N; cin >> N;
    int E; cin >> E;
    vector<int> edges(E);
    getVector(edges);
    vector<int> ret = gd.plot(N, edges);
    cout << ret.size() << endl;
    for (int i = 0; i < (int) ret.size(); ++i) { cout << ret[i] << endl; }
    cout.flush();
}
