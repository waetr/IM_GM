#ifndef MRIM_RRSET_H
#define MRIM_RRSET_H

#include <random>

#include "preliminaries.h"


class MultiRRContainer {
private:
    ///@brief temporary array for RI_Gen
    bool *DijkstraVis;
    size_t _sizeOfRRsets = 0;
    //set which is excluded from the spread propagation (bitwise representation)
    std::vector<bool> excludedNodes;
    int64 G_n = 0;
    int64 round = 0;


public:
    ///collection of multi-RI sets
    std::vector<std::vector<std::vector<int64>>> multi_R;

    ///covered[u] marks which RI sets the node u is covered by
    std::vector<int64> *covered;
    ///coveredNum[u] marks how many RI sets the node u is covered by
    int64 *coveredNum;

    MultiRRContainer() {
        DijkstraVis = nullptr;
        covered = nullptr;
        coveredNum = nullptr;
    }

    explicit MultiRRContainer(Graph &G, int64 mrim_round) {
        round = mrim_round;
        excludedNodes.resize(G.n, false);
        DijkstraVis = new bool[G.n]();
        covered = new std::vector<int64>[G.n * round]();
        coveredNum = new int64[G.n * round]();
        G_n = G.n;
    }

    ~MultiRRContainer() {
        delete[] DijkstraVis;
        delete[] covered;
        delete[] coveredNum;
    }

    /*!
     *
     * @return number of RR sets in this RR's container.
     */
    size_t numOfRRsets() const {
        return multi_R.size();
    }

    /*!
     * @return sum of the size of each RR set in this RR's container.
     */
    size_t sizeOfRRsets() const {
        return _sizeOfRRsets;
    }

    size_t idx(size_t node_idx, size_t round_idx) const {
        return round_idx * G_n + node_idx;
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
        if (graph.diff_model == IC) {
            std::deque<int64> Q;
            for (int64 u: uStart)
                if (!excludedNodes[u]) {
                    DijkstraVis[u] = true;
                    Q.push_back(u);
                }
            while (!Q.empty()) {
                int64 u = Q.front();
                Q.pop_front();
                RR.emplace_back(u);
                for (auto &edgeT: (*edge_list)[u]) {
                    if (excludedNodes[edgeT.v] || DijkstraVis[edgeT.v]) continue;
                    bool activate_success = (random_real() < edgeT.p);
                    if (activate_success) {
                        DijkstraVis[edgeT.v] = true;
                        Q.push_back(edgeT.v);
                    }
                }
            }
        }

        for (int64 u: RR) {
            DijkstraVis[u] = false;
        }
    }


    /*!
     * @brief Insert a random RR set into the container.
     * @param G
     */
    void insertOneRandomRRset(Graph &G, std::uniform_int_distribution<int64> &uniformIntDistribution) {
        std::vector<int64> RR;
        std::vector<std::vector<int64>> multi_RR;
        int64 v = uniformIntDistribution(mt19937engine);
        std::vector<int64> vStart = {v};
        for (int i = 0; i < round; ++i) {
            RR.clear();
            RI_Gen(G, vStart, RR);
            multi_RR.emplace_back(RR);
            _sizeOfRRsets += RR.size();
            for (int64 u: RR) {
                covered[idx(u, i)].emplace_back(multi_R.size());
                coveredNum[idx(u, i)]++;
            }
        }
        multi_R.emplace_back(multi_RR);
    }


    /*!
     * @brief While the required size is larger than the current size, add random RR sets
     * @param G
     * @param size
     */
    void resize(Graph &G, size_t size) {
        assert(multi_R.size() <= size);
        std::uniform_int_distribution<int64> uniformIntDistribution(0, G.n - 1);
        while (multi_R.size() < size) insertOneRandomRRset(G, uniformIntDistribution);
    }

    /*!
* @brief calculate the coverage of vertex set S on these RR sets
* @param G : the graph
* @param vecSeed : vertex set S
* @return : the number of RR sets that are covered by S
*/
    int64 self_inf_cal_multi(std::vector<bi_node> &vecSeed) const {
        assert(!multi_R.empty());
        std::vector<bool> vecBoolVst = std::vector<bool>(multi_R.size());
        for (auto seed: vecSeed) {
            for (auto node: covered[idx(seed.first, seed.second)]) {
                vecBoolVst[node] = true;
            }
        }
        return std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
    }
};


#endif //MRIM_RRSET_H
