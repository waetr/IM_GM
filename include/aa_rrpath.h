#ifndef AA_RRPATH_H
#define AA_RRPATH_H

#include <random>

#include "preliminaries.h"
#include "graphs.h"


class VRRPath {
private:
    ///@brief temporary array for RI_Gen
    bool *DijkstraVis;
    size_t _sizeOfRRsets = 0;
    std::vector<bool> excludedNodes;
    size_t G_n;


public:
    ///collection of RI sets
    std::vector<std::vector<int64>> R;

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
            edge_covered[i] = new std::vector<int64>[G.deg_out[i]]();
            edge_coveredNum[i] = new int64[G.deg_out[i]]();
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
        auto *edge_list = &graph.gT;
        auto u = uStart[0];
        while (true) {
            RR.emplace_back(u);
            std::vector<int64> in_neighbours;
            std::vector<double> weights;
            for (auto &edgeT: (*edge_list)[u]) {
                in_neighbours.emplace_back(edgeT.v);
                weights.emplace_back(edgeT.p);
            }
            std::discrete_distribution<> dist(weights.begin(), weights.end());
            int64 sampled_v = in_neighbours[dist(mt19937engine)];
            if (excludedNodes[sampled_v] || DijkstraVis[sampled_v]) break;
            u = sampled_v;
            DijkstraVis[sampled_v] = true;
        }

        for (int64 v: RR) {
            DijkstraVis[v] = false;
        }
    }


    /*!
     * @brief Insert a random RR set into the container.
     * @param G
     */
    void insertOneRandomRRset(Graph &G, std::uniform_int_distribution<int64> &uniformIntDistribution) {
        std::vector<int64> RR;
        while (RR.empty() || !excludedNodes[RR[RR.size() - 1]]) {
            RR.clear();
            int64 v = uniformIntDistribution(mt19937engine);
            while (excludedNodes[v]) v = uniformIntDistribution(mt19937engine);
            std::vector<int64> vStart = {v};
            RI_Gen(G, vStart, RR);
        }
        R.emplace_back(RR);
        _sizeOfRRsets += RR.size();
        for (int64 u: RR) {
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
        for (auto seed: vecSeed) {
            for (auto node: covered[seed]) {
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
        for (auto seed: vecSeed) {
            for (auto node: covered[seed.first]) {
                vecBoolVst[node] = true;
            }
        }
        return std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
    }
};


#endif //AA_RRPATH_h
