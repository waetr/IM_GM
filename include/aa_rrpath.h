#ifndef AA_RRPATH_H
#define AA_RRPATH_H

#include <random>

#include "preliminaries.h"
#include "graphs.h"

class VRRPath {

public:
    ///@brief temporary array for RI_Gen
    bool *DijkstraVis;
    size_t _sizeOfRRsets = 0;
    std::vector<bool> excludedNodes;
    size_t G_n;
    size_t R_size = 0, all_R_size = 0;

    std::vector<std::vector<int64>> R;
    std::vector<std::vector<int64>> R_edge;

    ///covered[u] marks which RI sets the node u is covered by
    std::vector<int64> *covered, **edge_covered;
    ///coveredNum[u] marks how many RI sets the node u is covered by
    int64 *coveredNum, **edge_coveredNum;

    VRRPath() {
        DijkstraVis = nullptr;
        covered = nullptr;
        coveredNum = nullptr;
        edge_covered = nullptr;
        edge_coveredNum = nullptr;
    }

    VRRPath(Graph &G, std::vector<int64> &S) {
        G_n = G.n;
        excludedNodes.resize(G.n, false);
        for (auto u: S) excludedNodes[u] = true;
        DijkstraVis = new bool[G.n]();
        covered = new std::vector<int64>[G.n]();
        coveredNum = new int64[G.n]();
        edge_covered = new std::vector<int64> *[G.n]();
        edge_coveredNum = new int64 *[G.n]();
        for (int i = 0; i < G.n; ++i) {
            edge_covered[i] = new std::vector<int64>[G.deg_in[i]]();
            edge_coveredNum[i] = new int64[G.deg_in[i]]();
        }
    }

    ~VRRPath() {
        delete[] DijkstraVis;
        delete[] covered;
        delete[] coveredNum;
        for (int i = 0; i < G_n; ++i) {
            delete[] edge_covered[i];
            delete[] edge_coveredNum[i];
        }
        delete[] edge_covered;
        delete[] edge_coveredNum;
    }

    /*!
     *
     * @return number of RR sets in this RR's container.
     */
    size_t numOfRRsets() const {
        return R_size;
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
    bool RI_Gen(Graph &graph, int64 u, std::vector<int64> &RR, std::vector<int64> &RR_edge) {
        bool flag = false;
        auto *edge_list = &graph.gT;
        while (graph.deg_in[u] > 0) {
            RR.emplace_back(u);
            DijkstraVis[u] = true;
            std::vector<int64> in_neighbours;
            std::vector<double> weights;
            for (auto &edgeT: (*edge_list)[u]) {
                in_neighbours.emplace_back(edgeT.v);
                weights.emplace_back(edgeT.p);
            }
            std::discrete_distribution<> dist(weights.begin(), weights.end());
            int64 sampled_num = dist(mt19937engine);
            int64 sampled_v = in_neighbours[sampled_num];
            RR_edge.emplace_back(sampled_num);
            if (excludedNodes[sampled_v] || DijkstraVis[sampled_v]) {
                if (excludedNodes[sampled_v]) flag = true;
                break;
            }
            u = sampled_v;
        }

        for (int64 v: RR) {
            DijkstraVis[v] = false;
        }
        return flag;
    }


    /*!
     * @brief Insert a random RR set into the container.
     * @param G
     */
    void insertOneRandomRRset(Graph &G, std::uniform_int_distribution<int64> &uniformIntDistribution) {
        std::vector<int64> RR;
        std::vector<int64> RR_edge;
        bool success = false;
        int tot = 0;
        while (!success) {
            tot++;
            RR.clear();
            RR_edge.clear();
            int64 v = uniformIntDistribution(mt19937engine);
            while (excludedNodes[v]) v = uniformIntDistribution(mt19937engine);
            success = RI_Gen(G, v, RR, RR_edge);
        }
        R_size += 1;
        all_R_size += tot;
        _sizeOfRRsets += RR.size() + RR_edge.size();
        for (int64 u: RR) {
            covered[u].emplace_back(R_size - 1);
            coveredNum[u]++;
        }
        for (int i = 0; i < RR_edge.size(); i++) {
            edge_covered[RR[i]][RR_edge[i]].emplace_back(R_size - 1);
            edge_coveredNum[RR[i]][RR_edge[i]]++;
        }
        R.emplace_back(RR);
        R_edge.emplace_back(RR_edge);
    }


    /*!
     * @brief While the required size is larger than the current size, add random RR sets
     * @param G
     * @param size
     */
    void resize(Graph &G, size_t size) {
        std::uniform_int_distribution<int64> uniformIntDistribution(0, G.n - 1);
        while (R_size < size) insertOneRandomRRset(G, uniformIntDistribution);
    }

    void resize1(Graph &G, size_t size) {
        std::uniform_int_distribution<int64> uniformIntDistribution(0, G.n - 1);
        while (all_R_size < size) insertOneRandomRRset(G, uniformIntDistribution);
    }


    // seed.first = node;  seed.second = -1 ? NULL : (index of edge)
    int64 self_inf_cal(std::vector<std::pair<int64, int64>> &vecSeed) const {
        std::vector<bool> vecBoolVst = std::vector<bool>(R_size);
        for (auto seed: vecSeed) {
            if (seed.second == -1) {
                for (auto node: covered[seed.first]) {
                    vecBoolVst[node] = true;
                }
            } else {
                for (auto node: edge_covered[seed.first][seed.second]) {
                    vecBoolVst[node] = true;
                }
            }
        }
        return std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
    }
};


#endif //AA_RRPATH_h
