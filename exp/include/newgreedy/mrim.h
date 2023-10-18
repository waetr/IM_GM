#ifndef EXP_MRIM_H
#define EXP_MRIM_H

#include "IMs.h"
#include "OPIM_new.h"
#include "newgreedy\Heap.h"


void Cross_Round_Node_Selection(Graph &G, RRContainer &RRI, int64 T, int64 k, std::vector<bi_node> &seeds) {
    coveredNum_tmp = new int64[G.n * T];
    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * T * sizeof(int64));
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    std::vector<int> cardinalities(T, 0);
    std::vector<bool> vis(G.n * T, false);
    for (int i = 0; i < T * k; ++i) {
        int64 covered_value = -1, u0, t0;
        for (int t = 0; t < T; ++t) {
            if (cardinalities[t] >= k) continue;
            for (int u = 0; u < G.n; ++u) {
                if (vis[t * G.n + u]) continue;
                if (coveredNum_tmp[t * G.n + u] > covered_value) {
                    covered_value = coveredNum_tmp[t * G.n + u];
                    u0 = u;
                    t0 = t;
                }
            }
        }
        vis[t0 * G.n + u0] = true;
        seeds.emplace_back(u0, t0);
        cardinalities[t0] += 1;
        for (auto RRIndex: RRI.covered[t0 * G.n + u0]) {
            if (RRSetCovered[RRIndex]) continue;
            for (int t = 0; t < T; ++t) {
                for (auto u: RRI.multi_R[RRIndex][t]) {
                    coveredNum_tmp[t * G.n + u]--;
                }
            }
            RRSetCovered[RRIndex] = true;
        }
    }
    delete[] coveredNum_tmp;
}

double M_CGreedy_Partition(Graph &G, RRContainer &RRI, int64 T, int64 k, int64 t_max, std::vector<bi_node> &bi_seeds) {
    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0;
    std::vector<int64> vis_time(G.n, T * t_max + 1);
    std::vector<double> q_R(RRI.numOfRRsets(), 1);
    //priority queue for lazy sampling
    auto Q = new HeapArray<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>[T];
    for (size_t i = 0; i < T; i++) {
        Q[i].val = new std::pair<double, int64>[G.n];
        Q[i].ids = new std::pair<unsigned, unsigned>[G.n];
        Q[i].k = 0;
    }

    for (int t = 1; t <= t_max; t++) {
        for (int i = 0; i < T; i++) {
            Q[i].k = 0;
            for (int j = 0; j < G.n; j++) {
                double value_v = 0;
                for (auto rr: RRI.covered[i * G.n + j]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[i * G.n + j]);
                Q[i].k++;
                heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>(Q[i].k, Q[i].val, Q[i].ids,
                                                                                         std::make_pair(value_v,
                                                                                                        vis_time[j]),
                                                                                         std::make_pair(0, j));
            }
        }
        std::vector<int64> cardinality(T, 0);
        for (int i = 0; i < T * k; i++) {
            int round = i % T;
            if (cardinality[round] >= k) continue;
            while (Q[round].k != 0) {
                int64 u = Q[round].ids[0].second;
                int64 it_round = Q[round].ids[0].first;
                heap_pop<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>(Q[round].k, Q[round].val,
                                                                                        Q[round].ids);
                Q[round].k -= 1;
                if (it_round == i) {
                    //choose
                    cardinality[round] += 1;
                    for (long rr: RRI.covered[round * G.n + u]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_x[round * G.n + u] - 1) /
                                   (double) (t_max - frac_x[round * G.n + u]);
                        Fx += q_R_old - q_R[rr];
                    }
                    frac_x[round * G.n + u] += 1;
                    vis_time[u] -= 1;
                    break;
                } else {
                    double value_v = 0;
                    for (long rr: RRI.covered[round * G.n + u]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_x[round * G.n + u]);
                    Q[round].k++;
                    heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>(Q[round].k, Q[round].val,
                                                                                             Q[round].ids,
                                                                                             std::make_pair(value_v,
                                                                                                            vis_time[u]),
                                                                                             std::make_pair(i, u));
                }
            }
        }
    }

    double tight_bound = Fx;
    std::vector<std::vector<double>> value(T);
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < G.n; ++j) {
            double value_v = 0;
            for (long rr: RRI.covered[i * G.n + j]) {
                value_v += q_R[rr];
            }
            value[i].push_back(value_v);
        }
        std::nth_element(value[i].begin(), value[i].begin() + k - 1, value[i].end(), std::greater<>());
        for (int j = 0; j < k; ++j) {
            tight_bound += value[i][j];
        }
    }

    //advanced-rounding
    for (int i = 0; i < T; i++) {
        int x = 0;
        while (x < G.n && (frac_x[G.n * i + x] == 0 || frac_x[G.n * i + x] == t_max)) x++;
        if (x == G.n) continue;
        int y = x + 1;
        while (y < G.n) {
            while (y < G.n && (frac_x[G.n * i + y] == 0 || frac_x[G.n * i + y] == t_max)) y++;
            double dx = 0, dy = 0;
            for (long rr: RRI.covered[G.n * i + x]) {
                dx += q_R[rr];
            }
            dx *= (double) t_max / (double) (t_max - frac_x[G.n * i + x]);
            for (long rr: RRI.covered[G.n * i + y]) {
                dy += q_R[rr];
            }
            dy *= (double) t_max / (double) (t_max - frac_x[G.n * i + y]);
            if (dx < dy) std::swap(x, y);
            if (t_max - frac_x[G.n * i + x] > frac_x[G.n * i + y]) {
                for (long rr: RRI.covered[G.n * i + x]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] *= (double) (t_max - frac_x[G.n * i + x] - frac_x[G.n * i + y]) /
                               (double) (t_max - frac_x[G.n * i + x]);
                    Fx += q_R_old - q_R[rr];
                }
                frac_x[G.n * i + x] += frac_x[G.n * i + y];
                for (long rr: RRI.covered[G.n * i + y]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] *= (double) t_max / (double) (t_max - frac_x[G.n * i + y]);
                    Fx += q_R_old - q_R[rr];
                }
                frac_x[G.n * i + y] = 0;
                y = std::max(x, y) + 1;
            } else {
                for (long rr: RRI.covered[G.n * i + y]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] *= (double) (t_max - (frac_x[G.n * i + x] + frac_x[G.n * i + y] - t_max)) /
                               (double) (t_max - frac_x[G.n * i + y]);
                    Fx += q_R_old - q_R[rr];
                }
                frac_x[G.n * i + y] = frac_x[G.n * i + x] + frac_x[G.n * i + y] - t_max;
                for (long rr: RRI.covered[G.n * i + x]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] = 0;
                    Fx += q_R_old;
                }
                frac_x[G.n * i + x] = t_max;
                int t = x;
                x = y;
                y = std::max(t, y) + 1;
            }
            if (frac_x[G.n * i + x] == 0 || frac_x[G.n * i + x] == t_max) {
                x = y;
                while (x < G.n && (frac_x[G.n * i + x] == 0 || frac_x[G.n * i + x] == t_max)) x++;
                if (x >= G.n) break;
                y = x + 1;
            }
        }
    }

    for (int j = 0; j < T; j++) {
        for (int l = 0; l < G.n; l++) {
            if (frac_x[G.n * j + l] == t_max) bi_seeds.emplace_back(l, j);
        }
    }

    for (size_t i = 0; i < T; i++) {
        delete[] Q[i].val;
        delete[] Q[i].ids;
    }
    delete[] Q;
    return tight_bound;
}


#endif //EXP_MRIM_H
