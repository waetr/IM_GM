//
// Created by asuka on 22-7-14.
//

#ifndef UNTITLED_GRAPH_H
#define UNTITLED_GRAPH_H

#include <random>

#include "preliminaries.h"

/*!
 * @brief triple<int, double, double> to represent an edge
 * first argument : index of outgoing int64
 * second argument : p_{u,v} in IC/WC model
 */
struct Edge {
    int64 v;
    double p;

    Edge() = default;

    Edge(int64 v, double p) : v(v), p(p) {}
};

const std::vector<Edge> edgeList_instance;

class Graph {
public:
    /*!
     * @param n : maximum index of node
     * @param m : number of edges
     * @param g : adjacency list
     */
    int64 n{};
    int64 m{};
    std::vector<std::vector<Edge> > g, gT;
    std::vector<int64> deg_in, deg_out;
    model_type diff_model{};

    /*!
     * @brief Init the graph for a default size.
     */
    Graph() = default;

    /*!
     * @brief A destructor for graph.
    */
    ~Graph() = default;

    /*!
     * @brief add an weighted edge into graph.
     * @param source : source node
     * @param target : destination node
     * @param weight : weight of edge, 1.0 as default
     */
    void add_edge(int64 source, int64 target) {
        n = std::max(n, std::max(source, target) + 1);
        while (g.size() < n) {
            g.emplace_back(edgeList_instance);
            gT.emplace_back(edgeList_instance);
            deg_in.emplace_back(0);
            deg_out.emplace_back(0);
        }
        m++;
        deg_in[target]++;
        deg_out[source]++;
        g[source].emplace_back(Edge(target, 0));
        gT[target].emplace_back(Edge(source, 0));
    }

/*!
 * @brief load a graph through file
 * @param filename : the path of the file
 * @param type : detrmine the graph type is directed or undirected
 */
    void open(const std::string &filename) {
        std::ifstream inFile(filename, std::ios::in);
        if (!inFile.is_open()) {
            std::cerr << "(get error) graph file not found: " << filename << std::endl;
            std::exit(-1);
        }
        int64 x = fast_read(inFile);
        int64 y = fast_read(inFile);
        while (x != -1 || y != -1) {
            add_edge(x, y);
            x = fast_read(inFile);
            y = fast_read(inFile);
        }
        inFile.close();
    }

    Graph(const std::string &filename) : Graph() {
        this->open(filename);
    }

    /*!
     * @brief Set the diffusion model to IC/LT. (Now only IC is supported.)
     * @param new_type : the name of the diffusion model.
     * @param p0 : if p0 > 0, then the activation probability will be set as p0 uniformly.
     */
    void set_diffusion_model(model_type new_type, double p0 = 0) {
        diff_model = new_type;
        if (p0 > 0) {
            printf("p0:%.3f\n", p0);
            for (int64 i = 0; i < n; i++) {
                for (int64 j = 0; j < g[i].size(); j++) {
                    g[i][j].p = p0;
                }
            }
            for (int64 i = 0; i < n; i++) {
                for (int64 j = 0; j < gT[i].size(); j++) {
                    gT[i][j].p = p0;
                }
            }
        } else {
            for (int64 i = 0; i < n; i++) {
                for (int64 j = 0; j < g[i].size(); j++) {
                    g[i][j].p = 1.0 / deg_in[g[i][j].v];
                }
            }
            for (int64 i = 0; i < n; i++) {
                for (int64 j = 0; j < gT[i].size(); j++) {
                    gT[i][j].p = 1.0 / deg_in[i];
                }
            }
        }
    }
};

/*!
 * A data structure for the advanced option that supports access to each active participant's neighboring nodes,
 * and specifying the active participant to which each node belongs during the algorithm.
 */
class CandidateNeigh {
private:
    int64 k{};
    int64 *num_neighbours{};
    std::vector<int64> *f{};//f[v] means in-coming active participant of v
    std::vector<int64> Ap;//active participant
public:
    std::vector<int64> N; //candidate neighbour set

    CandidateNeigh() = default;

    /*!
     * @brief create an instance.
     * @param graph : graph
     * @param A : set of active participant
     * @param k0 : budget
     */
    CandidateNeigh(Graph &graph, std::vector<int64> &A, int64 k0) {
        Ap = A;
        k = k0;
        num_neighbours = new int64[graph.n]();
        f = new std::vector<int64>[graph.n]();
        std::set<int64> N_unique, A_unique;
        for (auto u : A) A_unique.insert(u);

        //push all candidate neighbour to S, and update f
        for (auto u : A) {
            for (auto &edge : graph.g[u]) {
                if (A_unique.find(edge.v) == A_unique.end()) {
                    N_unique.insert(edge.v);
                    f[edge.v].emplace_back(u);
                }
            }
        }
        N.assign(N_unique.begin(), N_unique.end());
        //Disrupt the f-vector to achieve random assignment
        for (auto u : N)
            std::shuffle(f[u].begin(), f[u].end(), std::mt19937(std::random_device()()));
    }

    ~CandidateNeigh() {
        delete[] num_neighbours;
        delete[] f;
    }

    /*!
     * @brief assign function.
     * @param Num : node number of graph
     * @param x : the source class to assign
     */
    void assign(int64 Num, CandidateNeigh &x) {
        if (this == &x) return;
        delete[] num_neighbours;
        delete[] f;
        k = x.k;
        Ap = x.Ap;
        N = x.N;
        num_neighbours = new int64[Num]();
        f = new std::vector<int64>[Num]();
        for (int i = 0; i < Num; i++) {
            num_neighbours[i] = x.num_neighbours[i];
            f[i] = x.f[i];
        }
    }

    /*!
     * @brief assign a source participant of candidate node v
     * @param v : node
     * @return if no AP can be chosen, return -1; otherwise return the index of the Ap
     */
    int64 source_participant(int64 v) {
        ///assignment strategy: random
        while (!f[v].empty() && num_neighbours[f[v][f[v].size() - 1]] >= k) f[v].pop_back();
        return f[v].empty() ? -1 : f[v][f[v].size() - 1];
    }
//    int64 source_participant(int64 v) {
//        ///assignment strategy: Priority is given to the ap with the least number of candidates.
//        if(f[v].empty()) return -1;
//        int u = 0;
//        for(int i = 1; i < f[v].size(); i++) if(num_neighbours[f[v][i]] < num_neighbours[f[v][u]]) u = i;
//        if(num_neighbours[f[v][u]] >= k) return -1;
//        return f[v][u];
//    }

    /*!
     * @brief choose the chosen ap in source_participant(), modify the structure
     * @param u : the ap chosen by source_participant()
     */
    void choose(int64 u) {
        assert(u != -1);
        num_neighbours[u]++;
    }
};


class RRContainer {
private:
    ///@brief temporary array for RI_Gen
    int64 *dist;
    bool *DijkstraVis;
    size_t _sizeOfRRsets = 0;
    bool RIFlag = false; //false : FI sets true: RI sets
    //set which is excluded from the spread propagation (bitwise representation)
    std::vector<bool> excludedNodes;


public:
    ///collection of RI sets
    std::vector<std::vector<int64>> R;

    ///covered[u] marks which RI sets the node u is covered by
    std::vector<int64> *covered;
    ///coveredNum[u] marks how many RI sets the node u is covered by
    int64 *coveredNum;

    RRContainer() {
        dist = nullptr;
        DijkstraVis = nullptr;
        covered = nullptr;
        coveredNum = nullptr;
    }

    explicit RRContainer(Graph &G, std::vector<int64> &A, bool isRRSet) {
        excludedNodes.resize(G.n, false);
        for (auto u : A) excludedNodes[u] = true;
        RIFlag = isRRSet;
        dist = new int64[G.n]();
        DijkstraVis = new bool[G.n]();
        covered = new std::vector<int64>[G.n]();
        coveredNum = new int64[G.n]();
        memset(dist, -1, sizeof(int64) * G.n);
    }

    ~RRContainer() {
        delete[] dist;
        delete[] DijkstraVis;
        delete[] covered;
        delete[] coveredNum;
    }

    /*!
     *
     * @return number of RR sets in this RR's container.
     */
    size_t numOfRRsets() const {
        return R.size();
    }

    /*!
     * @return sum of the size of each RR set in this RR's container.
     */
    size_t sizeOfRRsets() const {
        return _sizeOfRRsets;
    }

    /*!
 * @brief Algorithm for CTIC to generate Reserve-Influence or Forward-Influence set of IMM.
 * @param graph : the graph
 * @param uStart : the starting nodes of this RI/FI set
 * @param RR : returns the RI/FI set as an passed parameter
 */
    void RI_Gen(Graph &graph, std::vector<int64> &uStart, std::vector<int64> &RR) {
        assert(RR.empty());
        auto *edge_list = RIFlag ? &graph.gT : &graph.g;
        if (graph.diff_model == IC) {
            std::deque<int64> Q;
            for (int64 u : uStart)
                if (!excludedNodes[u]) {
                    DijkstraVis[u] = true;
                    Q.push_back(u);
                }
            while (!Q.empty()) {
                int64 u = Q.front();
                Q.pop_front();
                RR.emplace_back(u);
                for (auto &edgeT : (*edge_list)[u]) {
                    if (excludedNodes[edgeT.v] || DijkstraVis[edgeT.v]) continue;
                    bool activate_success = (random_real() < edgeT.p);
                    if (activate_success) {
                        DijkstraVis[edgeT.v] = true;
                        Q.push_back(edgeT.v);
                    }
                }
            }
        }

        for (int64 u : RR) {
            DijkstraVis[u] = false;
        }
    }


    /*!
     * @brief Insert a random RR set into the container.
     * @param G
     */
    void insertOneRandomRRset(Graph &G, std::uniform_int_distribution<int64> &uniformIntDistribution) {
        std::vector<int64> RR;
        int64 v = uniformIntDistribution(mt19937engine);
        std::vector<int64> vStart = {v};
        RI_Gen(G, vStart, RR);
        R.emplace_back(RR);
        _sizeOfRRsets += RR.size();
        for (int64 u : RR) {
            covered[u].emplace_back(R.size() - 1);
            coveredNum[u]++;
        }
    }


    /*!
     * @brief While the required size is larger than the current size, add random RR sets
     * @param G
     * @param size
     */
    void resize(Graph &G, size_t size) {
        assert(R.size() <= size);
        std::uniform_int_distribution<int64> uniformIntDistribution(0, G.n - 1);
        while (R.size() < size) insertOneRandomRRset(G, uniformIntDistribution);
    }

    /*!
     * @brief calculate the coverage of vertex set S on these RR sets
     * @param G : the graph
     * @param vecSeed : vertex set S
     * @return : the number of RR sets that are covered by S
     */
    int64 self_inf_cal(Graph &G, std::vector<int64> &vecSeed) const {
        std::vector<bool> vecBoolVst = std::vector<bool>(R.size());
        for (auto seed : vecSeed) {
            for (auto node : covered[seed]) {
                vecBoolVst[node] = true;
            }
        }
        return std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
    }

    /*!
 * @brief calculate the coverage of vertex set S on these RR sets
 * @param G : the graph
 * @param vecSeed : vertex set S
 * @return : the number of RR sets that are covered by S
 */
    int64 self_inf_cal(Graph &G, std::vector<bi_node> &vecSeed) const {
        std::vector<bool> vecBoolVst = std::vector<bool>(R.size());
        for (auto seed : vecSeed) {
            for (auto node : covered[seed.first]) {
                vecBoolVst[node] = true;
            }
        }
        return std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
    }
};


#endif //UNTITLED_GRAPH_H
