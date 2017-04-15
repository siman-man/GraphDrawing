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
double TIME_LIMIT = 10.0;
int FIELD_SIZE = 701;

unsigned long long xor128() {
    static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
    unsigned long long rt = (rx ^ (rx << 11));
    rx = ry;
    ry = rz;
    rz = rw;
    return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
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

Vertex g_vertexList[MAX_V];
Edge g_edgeList[MAX_V];
int g_desireDist[MAX_V][MAX_V];
int g_edgeCount;
int g_vertexCount;


class GraphDrawing {
public:
    void init(int N, vector<int> edges) {
        g_vertexCount = N;
        g_edgeCount = edges.size() / 3;

        for (int i = 0; i < g_edgeCount; i++) {
            int from = edges[i * 3];
            int to = edges[i * 3 + 1];
            int dist = edges[i * 3 + 2];
            g_edgeList[i] = Edge(from, to, dist);

            fprintf(stderr, "Edge %d: from = %d, to = %d, d = %d\n", i, from, to, dist);
        }

        memset(g_desireDist, -1, sizeof(g_desireDist));

        fprintf(stderr, "VertexCount = %d, EdgeCount = %d\n", g_vertexCount, g_edgeCount);
    }

    vector<int> plot(int N, vector<int> edges) {
        init(N, edges);

        for (int i = 0; i < N; ++i) {
            int y = xor128() % FIELD_SIZE;
            int x = xor128() % FIELD_SIZE;
            g_vertexList[i] = Vertex(y, x);
        }

        vector<int> ret = createAnswer();
        double score = calcScore();

        fprintf(stderr,"CurrentScore = %f\n", score);

        return ret;
    }

    vector<int> createAnswer() {
        vector<int> ret;

        for (int i = 0; i < g_vertexCount; i++) {
            Vertex v = g_vertexList[i];
            ret.push_back(v.y);
            ret.push_back(v.x);
        }

        return ret;
    }

    double calcScore() {
        double minDiff = DBL_MAX;
        double maxDiff = 0.0;

        for (int i = 0; i < g_edgeCount; i++) {
            Edge edge = g_edgeList[i];
            Vertex v1 = g_vertexList[edge.from];
            Vertex v2 = g_vertexList[edge.to];
            double dist = calcDist(v1, v2);

            double diff = dist / (double) edge.dist;
            minDiff = min(minDiff, diff);
            maxDiff = max(maxDiff, diff);
        }

        return minDiff / maxDiff;
    }

    inline doublet s calcDist(Vertex &v1, Vertex &v2) {
        int dy = v1.y - v2.y;
        int dx = v1.x - v2.x;

        return sqrt(dy * dy + dx * dx);
    }
};
// -------8<------- end of solution submitted to the website -------8<-------

template<class T>
void getVector(vector <T> &v) { for (int i = 0; i < v.size(); ++i) cin >> v[i]; }

int main() {
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
