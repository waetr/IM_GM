
#ifndef SUBSET_IM_IA_H
#define SUBSET_IM_IA_H

#include "IMs.h"
#include "OPIM_new.h"

void IA_simulation(Graph &graph, std::vector<int64> &S, std::vector<double> &IA) {
    assert(IA.size() == S.size());
    int64 it_rounds = 1000;
    std::set<int64> source[graph.n];
    std::vector<bool> vis(graph.n);

    double source_active_time = 0;

    //An size-n index hash map
    std::vector<int64> seed_index(graph.n);
    for (int i = 0; i < S.size(); i++) seed_index[S[i]] = i;

    for (int i = 1; i <= it_rounds; i++) {
        std::cout << i << "\n";
        for (int64 e: S) vis[e] = true;
        for (int64 e: S) source[e].insert(e);
        std::vector<int64> A = S, prev_new = S;
        while (!prev_new.empty()) {
            std::vector<int64> New;
            for (int64 u: prev_new) {
                for (auto &e: graph.g[u]) {
                    int64 v = e.v;
                    if (vis[v]) continue;
                    if (random_real() < e.p) {
                        if (source[v].empty()) New.push_back(v);
                        source[v].insert(source[u].begin(), source[u].end());
                    }
                }
            }
            for (int64 e: New) vis[e] = true;
            A.insert(A.end(), New.begin(), New.end());
            prev_new.assign(New.begin(), New.end());
        }
        double one_source_active_time = 0;
        for (int64 u: A) {
            one_source_active_time += source[u].size();
            for (int64 s: source[u]) {
                IA[seed_index[s]] += 1.0 / source[u].size();
            }
            vis[u] = false;
            source[u].clear();
        }
        source_active_time += one_source_active_time / A.size();
    }
    for (auto &e: IA) e /= it_rounds;
    source_active_time /= it_rounds;
    std::cout << "active time: " << source_active_time << "\n";
}

#endif //SUBSET_IM_IA_H
