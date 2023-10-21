#ifndef EXP_SIMULATION_H
#define EXP_SIMULATION_H

#include "graphs.h"


/*!
 * @brief generate a random node set of graph.
 * @param graph : the graph
 * @param A : stores the node set. Suppose it is initialized as empty.
 * @param size : the size of the node set
 */
void generate_ap(Graph &graph, std::vector<int64> &A, int64 size = 1) {
    A.clear();
    std::set<int64> S;
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    for (int64 i = 0; i < size; i++) {
        int64 v = uniformIntDistribution(mt19937engine);
        while (1) {
            bool flag = false;
            if(S.find(v) != S.end()) flag = true;
            else {
                for (auto &edge: graph.g[v]) {
                    if(S.find(edge.v) != S.end()) {
                        flag = true;
                        break;
                    }
                }
            }
            if(!flag) break;
            v = uniformIntDistribution(mt19937engine);
        }
        A.emplace_back(v);
        S.insert(v);
        for (auto &edge: graph.g[v]) {
            S.insert(edge.v);
        }
    }
}



/*!
 * @brief Calculate the degree of neighbor overlap at active participant.
 * @param graph : the graph
 * @param seeds : the active participant set
 * @param iteration_rounds : The number of selection
 * @return the mean overlap ratio
 */
double estimate_neighbor_overlap(Graph &graph, std::vector<int64> &seeds) {
    auto *num = new int64[graph.n]();
    int64 tot = 0, overlap = 0;
    for (int64 u : seeds)
        for (auto &e : graph.g[u])
            num[e.v]++;
    for (int64 u : seeds) {
        for (auto &e : graph.g[u]) {
            if (num[e.v] > 0) tot++;
            if (num[e.v] > 1) overlap++;
            num[e.v] = 0;
        }
    }
    delete[] num;
    return (double) overlap / tot;
}

/*!
 * @brief New calculation of influence spread, in which the effect of ap is excluded
 * @param graph : the graph that define propagation models(IC-M)
 * @param S : the seed set
 * @param A : the active participant
 * @return the estimated value of influence spread
 */
double MC_simulation(Graph &graph, std::vector<int64> &S, std::vector<int64> &A, int64 it_rounds = -1) {
    RRContainer RRI(graph, A, false);
    double res = 0;
    int64 it_ = (it_rounds == -1) ? MC_iteration_rounds : it_rounds;
    std::vector<int64> RR;
    for (int i = 1; i <= it_; i++) {
        RRI.RI_Gen(graph, S, RR);
        res += RR.size();
        RR.clear();
    }
    return res / it_;
}

/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
double effic_inf(Graph &graph, std::vector<bi_node> &S, std::vector<int64> &A) {
    const double delta = 1e-3, eps = 0.01, c = 2.0 * (exp(1.0) - 2.0);
    const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
    size_t numHyperEdge = 0, numCoverd = 0;
    std::vector<bool> vecBoolSeed(graph.n), exclusive(graph.n);
    for (auto seed : S) vecBoolSeed[seed.first] = true;
    for (auto ap : A) exclusive[ap] = true;
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    std::vector<int64> nodes;
    std::queue<int64> Q;
    while (numCoverd < LambdaL) {
        numHyperEdge++;
        const auto uStart = uniformIntDistribution(mt19937engine);
        if (exclusive[uStart]) {
            continue;
        }
        if (vecBoolSeed[uStart]) {
            // Stop, this sample is covered
            numCoverd++;
            continue;
        }
        exclusive[uStart] = true;
        Q.push(uStart);
        nodes.emplace_back(uStart);
        bool break_flag = false;
        while (!Q.empty()) {
            int64 u = Q.front();
            Q.pop();
            for (auto &edgeT : graph.gT[u]) {
                if (exclusive[edgeT.v]) continue;
                if (random_real() < edgeT.p) {
                    if (vecBoolSeed[edgeT.v]) {
                        numCoverd++;
                        break_flag = true;
                        break;
                    }
                    exclusive[edgeT.v] = true;
                    Q.push(edgeT.v);
                    nodes.emplace_back(edgeT.v);
                }
            }
            if (break_flag) break;
        }
        for (auto e : nodes) exclusive[e] = false;
        nodes.clear();
        Q = std::queue<int64>(); // clear the queue
    }
    return 1.0 * numCoverd * graph.n / numHyperEdge;
}

#endif //EXP_SIMULATION_H