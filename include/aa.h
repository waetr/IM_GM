//
// Created by admin on 10/11/2023.
//

#ifndef EXP_AA_H
#define EXP_AA_H

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"
#include "aa_rrpath.h"

double calc_bound_AA(Graph &G, int64 k_N, int64 k_E, VRRPath &RRI, std::vector<double> &q_R) {
    double sum = 0;
    std::vector<double> value_N, value_E;
    for (int i = 0; i < G.n; ++i) {
        double value_v = 0;
        for (auto rr: RRI.covered[i]) {
            value_v += q_R[rr];
        }
        value_N.emplace_back(value_v);
    }
    for (int i = 0; i < G.n; ++i) {
        for (int j = 0; j < G.deg_in[i]; ++j) {
            double value_v = 0;
            for (auto rr: RRI.edge_covered[i][j]) {
                value_v += q_R[rr];
            }
            value_E.emplace_back(value_v);
        }
    }
    std::nth_element(value_N.begin(), value_N.begin() + k_N - 1, value_N.end(), std::greater<>());
    std::nth_element(value_E.begin(), value_E.begin() + k_E - 1, value_E.end(), std::greater<>());
    for (int j = 0; j < k_N; ++j) {
        sum += value_N[j];
    }
    for (int j = 0; j < k_E; ++j) {
        sum += value_E[j];
    }
    return sum;
}

void rounding_AA(Graph &G, std::vector<std::vector<bi_node>> &bases, std::vector<int64> &frac_N,
                 std::vector<std::vector<int64>> &frac_E,
                 std::vector<double> &q_R, VRRPath &RRI, int64 t_max, std::vector<bi_node> &seeds) {
    std::vector<bool> pos_N(G.n, false), pos1_N(G.n, false);
    std::vector<std::vector<bool>> pos_E(G.n), pos1_E(G.n);
    for (int i = 0; i < G.n; ++i) {
        pos_E[i].resize(G.deg_in[i], false);
        pos1_E[i].resize(G.deg_in[i], false);
    }
    for (auto &i: bases[0]) {
        //first: node second: edge-idx
        if (i.second == -1) pos_N[i.first] = true;
        else pos_E[i.first][i.second] = true;
    }
    for (int i = 0; i < bases.size() - 1; ++i) {
        for (auto &j: bases[0]) {
            //first: node second: edge-idx
            if (j.second == -1) pos_N[j.first] = true;
            else pos_E[j.first][j.second] = true;
        }
        int x = 0, y = 0;
        while (x < pos_N.size() && y < pos1_N.size()) {
            while (x < pos_N.size() && !(pos_N[x] && !pos1_N[x])) x++;
            while (y < pos_N.size() && !(!pos_N[y] && pos1_N[y])) y++;
            if (x == pos_N.size() || y == pos1_N.size()) break;
            //x, y: node

            // round x and y
            double dx = 0, dy = 0;
            for (auto rr: RRI.covered[x]) {
                dx += q_R[rr];
            }
            dx *= (double) t_max / (double) (t_max - frac_N[x]);
            for (auto rr: RRI.covered[y]) {
                dy += q_R[rr];
            }
            dy *= (double) t_max / (double) (t_max - frac_N[y]);
            int x_ = x, y_ = y, swap_flag = 1;
            if (dx < dy) {
                std::swap(x_, y_);
                pos_N[x] = false, pos_N[y] = true;
                swap_flag = i + 1;
            } else {
                pos1_N[y] = false, pos1_N[x] = true;
            }
            for (auto rr: RRI.covered[y_]) {
                q_R[rr] *= (double) (t_max - frac_N[y_] + swap_flag) /
                           (double) (t_max - frac_N[y_]);
            }
            frac_N[y_] -= swap_flag;
            for (auto rr: RRI.covered[x_]) {
                q_R[rr] *= (double) (t_max - frac_N[x_] - swap_flag) /
                           (double) (t_max - frac_N[x_]);
            }
            frac_N[x_] += swap_flag;
            //end of rounding x and y
        }

        int xi = 0, xj = 0, yi = 0, yj = 0;
        while (true) {
            for (; xi < G.n; xi++) {
                bool flag = false;
                for (; xj < G.deg_in[xi]; ++xj) {
                    if (pos_E[xi][xj] && !pos1_E[xi][xj]) {
                        flag = true;
                        break;
                    }
                }
                if (flag) break;
                xj = 0;
            }
            for (; yi < G.n; yi++) {
                bool flag = false;
                for (; yj < G.deg_in[yi]; ++yj) {
                    if (!pos_E[yi][yj] && pos1_E[yi][yj]) {
                        flag = true;
                        break;
                    }
                }
                if (flag) break;
                yj = 0;
            }
            if (xi >= G.n || yi >= G.n) break;

            //start of rounding
            double dx = 0, dy = 0;
            for (auto rr: RRI.edge_covered[xi][xj]) {
                dx += q_R[rr];
            }
            dx *= (double) t_max / (double) (t_max - frac_E[xi][xj]);
            for (auto rr: RRI.edge_covered[yi][yj]) {
                dy += q_R[rr];
            }
            dy *= (double) t_max / (double) (t_max - frac_E[yi][yj]);
            int xi_ = xi, xj_ = xj, yi_ = yi, yj_ = yj, swap_flag = 1;
            if (dx < dy) {
                std::swap(xi_, yi_);
                std::swap(xj_, yj_);
                pos_E[xi][xj] = false, pos_E[yi][yj] = true;
                swap_flag = i + 1;
            } else {
                pos1_E[yi][yj] = false, pos1_E[xi][xj] = true;
            }
            for (auto rr: RRI.edge_covered[yi_][yj_]) {
                q_R[rr] *= (double) (t_max - frac_E[yi_][yj_] + swap_flag) /
                           (double) (t_max - frac_E[yi_][yj_]);
            }
            frac_E[yi_][yj_] -= swap_flag;
            for (auto rr: RRI.edge_covered[xi_][xj_]) {
                q_R[rr] *= (double) (t_max - frac_E[xi_][xj_] - swap_flag) /
                           (double) (t_max - frac_E[xi_][xj_]);
            }
            frac_E[xi_][xj_] += swap_flag;
        }
    }
    for (int i = 0; i < pos_N.size(); ++i) {
        if (pos_N[i]) seeds.emplace_back(i, -1);
    }
    for (int i = 0; i < pos_E.size(); ++i) {
        for (int j = 0; j < pos_E[i].size(); ++j) {
            if (pos_E[i][j]) seeds.emplace_back(i, j);
        }
    }
}

double CGreedy_AA(Graph &G, VRRPath &RRI, int64 k_N, int64 k_E, int64 t_max, std::vector<bi_node> &bi_seeds) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_N(G.n, 0);
    std::vector<std::vector<int64>> frac_E(G.n);
    for (int i = 0; i < G.n; ++i) frac_E[i].resize(G.deg_in[i]);
    //temporary varible
    double Fx = 0;
    double tight_bound = RRI.numOfRRsets();
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    //priority queue for lazy sampling
    HeapArray<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>> Q{};
    Q.val = new double[G.n + G.m];
    Q.ids = new std::pair<std::pair<unsigned, int>, unsigned>[G.n + G.m];
    Q.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);
    double cur = clock();
    for (int t = 1; t <= t_max; t++) {
        Q.k = 0;
        tight_bound = std::min(tight_bound, Fx + calc_bound_AA(G, k_N, k_E, RRI, q_R));
        for (int i = 0; i < G.n; i++) {
            if (RRI.excludedNodes[i]) continue; //is seed
            double value_v = 0;
            for (auto rr: RRI.covered[i]) {
                value_v += q_R[rr];
            }
            value_v *= (double) t_max / (double) (t_max - frac_N[i]);
            Q.k++;
            heap_push<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                   value_v, std::make_pair(
                            std::make_pair(i, -1), 0));
            for (int j = 0; j < G.deg_in[i]; j++) {
                value_v = 0;
                for (auto rr: RRI.edge_covered[i][j]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_E[i][j]);
                Q.k++;
                heap_push<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                       value_v, std::make_pair(
                                std::make_pair(i, j), 0));
            }
        }
        int64 c_N = 0, c_E = 0;
        for (int i = 0; i < k_N + k_E; i++) {
            while (Q.k != 0) {
                //j: node l: edge-index
                int64 j = Q.ids[0].first.first;
                int64 l = Q.ids[0].first.second;
                int64 it_round = Q.ids[0].second;
                heap_pop<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q.k, Q.val, Q.ids);
                Q.k -= 1;
                if ((c_N == k_N && l == -1) || (c_E == k_E && l != -1)) continue;
                if (it_round == i) {
                    //choose
                    if (l == -1) {
                        for (auto rr: RRI.covered[j]) {
                            double q_R_old = q_R[rr];
                            q_R[rr] *= (double) (t_max - frac_N[j] - 1) /
                                       (double) (t_max - frac_N[j]);
                            Fx += q_R_old - q_R[rr];
                        }
                        frac_N[j] += 1;
                        c_N++;
                    } else {
                        for (auto rr: RRI.edge_covered[j][l]) {
                            double q_R_old = q_R[rr];
                            q_R[rr] *= (double) (t_max - frac_E[j][l] - 1) /
                                       (double) (t_max - frac_E[j][l]);
                            Fx += q_R_old - q_R[rr];
                        }
                        frac_E[j][l] += 1;
                        c_E++;
                    }
                    bases[t - 1].emplace_back(j, l);
                    break;
                } else {
                    double value_v = 0;
                    if (l == -1) {
                        for (auto rr: RRI.covered[j]) {
                            value_v += q_R[rr];
                        }
                        value_v *= (double) t_max / (double) (t_max - frac_N[j]);
                    } else {
                        for (auto rr: RRI.edge_covered[j][l]) {
                            value_v += q_R[rr];
                        }
                        value_v *= (double) t_max / (double) (t_max - frac_E[j][l]);
                    }
                    Q.k++;
                    heap_push<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                           value_v, std::make_pair(
                                    std::make_pair(j, l), i));
                }
            }
            if(i == (k_E+k_N)/2)tight_bound = std::min(tight_bound, Fx + calc_bound_AA(G, k_N, k_E, RRI, q_R));
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    tight_bound = std::min(tight_bound, Fx + calc_bound_AA(G, k_N, k_E, RRI, q_R));

    rounding_AA(G, bases, frac_N, frac_E, q_R, RRI, t_max, bi_seeds);
    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}

#endif //EXP_AA_H
