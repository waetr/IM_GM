//
// Created by admin on 10/11/2023.
//

#ifndef EXP_AA_H
#define EXP_AA_H

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"
#include "aa_rrpath.h"

void TRGreedy_AA(Graph &G, VRRPath &RRI, int64 k_N, int64 k_E, double eps, std::vector<bi_node> &seeds) {
    coveredNum_tmp = new int64[G.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * sizeof(int64));
    auto **edgeCoveredNum_tmp = new int64 *[G.n]();
    for (int i = 0; i < G.n; ++i) {
        edgeCoveredNum_tmp[i] = new int64[G.deg_in[i]]();
        memcpy(edgeCoveredNum_tmp[i], RRI.edge_coveredNum[i], G.deg_in[i] * sizeof(int64));
    }
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    int64 c_N = 0, c_E = 0;
    std::vector<bool> vis(G.n, false);
    std::vector<std::vector<bool>> vis_edge(G.n);
    for (int i = 0; i < G.n; ++i) {
        vis_edge[i].resize(G.deg_in[i], false);
    }
    double w0 = 0;
    for (int i = 0; i < G.n; ++i) {
        w0 = std::max(w0, 1.0 * coveredNum_tmp[i]);
        for (int j = 0; j < G.deg_in[i]; j++) w0 = std::max(w0, 1.0 * edgeCoveredNum_tmp[i][j]);
    }
    for (double w = w0; w > eps * w0 / (k_N + k_E); w *= 1.0 - eps) {
        for (int i = 0; i < G.n; ++i) {
            if (c_N < k_N && !vis[i] && coveredNum_tmp[i] >= w) {
                vis[i] = true;
                seeds.emplace_back(i, -1);
                c_N++;
                for (auto RRIndex: RRI.covered[i]) {
                    if (RRSetCovered[RRIndex]) continue;
                    for (int l = 0; l < RRI.R[RRIndex].size(); ++l) {
                        auto u0 = RRI.R[RRIndex][l];
                        auto e0 = RRI.R_edge[RRIndex][l];
                        coveredNum_tmp[u0]--;
                        edgeCoveredNum_tmp[u0][e0]--;
                    }
                    RRSetCovered[RRIndex] = true;
                }
            }
            for (int j = 0; j < G.deg_in[i]; ++j) {
                if (c_E < k_E && !vis_edge[i][j] && edgeCoveredNum_tmp[i][j] >= w) {
                    vis_edge[i][j] = true;
                    seeds.emplace_back(i, j);
                    c_E++;
                    for (auto RRIndex: RRI.edge_covered[i][j]) {
                        if (RRSetCovered[RRIndex]) continue;
                        for (int l = 0; l < RRI.R[RRIndex].size(); ++l) {
                            auto u0 = RRI.R[RRIndex][l];
                            auto e0 = RRI.R_edge[RRIndex][l];
                            coveredNum_tmp[u0]--;
                            edgeCoveredNum_tmp[u0][e0]--;
                        }
                        RRSetCovered[RRIndex] = true;
                    }
                }
            }
            if (c_N == k_N && c_E == k_E) break;
        }
        if (c_N == k_N && c_E == k_E) break;
    }
    delete[] coveredNum_tmp;
    for (int i = 0; i < G.n; ++i) {
        delete[] edgeCoveredNum_tmp[i];
    }
    delete[] edgeCoveredNum_tmp;
}


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

double CGreedy_AA(Graph &G, VRRPath &RRI, int64 k_N, int64 k_E, int64 t_max, std::vector<bi_node> &bi_seeds,
                  bool bound_flag = false) {

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
        if (bound_flag) tight_bound = std::min(tight_bound, Fx + calc_bound_AA(G, k_N, k_E, RRI, q_R));
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
        }
    }

    if (bound_flag) tight_bound = std::min(tight_bound, Fx + calc_bound_AA(G, k_N, k_E, RRI, q_R));

    rounding_AA(G, bases, frac_N, frac_E, q_R, RRI, t_max, bi_seeds);
    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}

double CGreedy_AA_PM(Graph &G, VRRPath &RRI, int64 k_N, int64 k_E, int64 t_max, std::vector<bi_node> &bi_seeds,
                     bool bound_flag = false) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_N(G.n, 0);
    std::vector<std::vector<int64>> frac_E(G.n);
    for (int i = 0; i < G.n; ++i) frac_E[i].resize(G.deg_in[i]);
    //temporary varible
    double Fx = 0;
    double tight_bound = RRI.numOfRRsets();
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    //priority queue for lazy sampling
    HeapArray<CMax<double, std::pair<unsigned, unsigned>>> Q_N{};
    Q_N.val = new double[G.n];
    Q_N.ids = new std::pair<unsigned, unsigned>[G.n];
    Q_N.k = 0;

    HeapArray<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>> Q_E{};
    Q_E.val = new double[G.m];
    Q_E.ids = new std::pair<std::pair<unsigned, int>, unsigned>[G.m];
    Q_E.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);
    double cur = clock();
    for (int t = 1; t <= t_max; t++) {
        Q_N.k = 0;
        Q_E.k = 0;
        std::vector<double> value_N, value_E;
        for (int i = 0; i < G.n; i++) {
            if (RRI.excludedNodes[i]) continue; //is seed
            double value_v = 0;
            for (auto rr: RRI.covered[i]) {
                value_v += q_R[rr];
            }
            if (bound_flag) value_N.emplace_back(value_v);
            value_v *= (double) t_max / (double) (t_max - frac_N[i]);
            Q_N.k++;
            heap_push<CMax<double, std::pair<unsigned, unsigned>>>(Q_N.k, Q_N.val, Q_N.ids,
                                                                   value_v, std::make_pair(i, 0));
            for (int j = 0; j < G.deg_in[i]; j++) {
                value_v = 0;
                for (auto rr: RRI.edge_covered[i][j]) {
                    value_v += q_R[rr];
                }
                if (bound_flag) value_E.emplace_back(value_v);
                value_v *= (double) t_max / (double) (t_max - frac_E[i][j]);
                Q_E.k++;
                heap_push<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q_E.k, Q_E.val, Q_E.ids,
                                                                                       value_v, std::make_pair(
                                std::make_pair(i, j), 0));
            }
        }
        if (bound_flag) {
            int64 sum = 0;
            std::nth_element(value_N.begin(), value_N.begin() + k_N - 1, value_N.end(), std::greater<>());
            std::nth_element(value_E.begin(), value_E.begin() + k_E - 1, value_E.end(), std::greater<>());
            for (int j = 0; j < k_N; ++j) {
                sum += value_N[j];
            }
            for (int j = 0; j < k_E; ++j) {
                sum += value_E[j];
            }
            tight_bound = std::min(tight_bound, Fx + sum);
        }
        for (int i = 0; i < k_N; i++) {
            while (Q_N.k != 0) {
                //j: node l: edge-index
                int64 j = Q_N.ids[0].first;
                int64 it_round = Q_N.ids[0].second;
                heap_pop<CMax<double, std::pair<unsigned, unsigned>>>(Q_N.k, Q_N.val, Q_N.ids);
                Q_N.k -= 1;
                if (it_round == i) {
                    //choose
                    for (auto rr: RRI.covered[j]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_N[j] - 1) /
                                   (double) (t_max - frac_N[j]);
                        Fx += q_R_old - q_R[rr];
                    }
                    frac_N[j] += 1;
                    bases[t - 1].emplace_back(j, -1);
                    break;
                } else {
                    double value_v = 0;
                    for (auto rr: RRI.covered[j]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_N[j]);
                    Q_N.k++;
                    heap_push<CMax<double, std::pair<unsigned, unsigned>>>(Q_N.k, Q_N.val, Q_N.ids,
                                                                           value_v, std::make_pair(j, i));
                }
            }
        }
        for (int i = 0; i < k_E; i++) {
            while (Q_E.k != 0) {
                //j: node l: edge-index
                int64 j = Q_E.ids[0].first.first;
                int64 l = Q_E.ids[0].first.second;
                int64 it_round = Q_E.ids[0].second;
                heap_pop<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q_E.k, Q_E.val, Q_E.ids);
                Q_E.k -= 1;
                if (it_round == i) {
                    //choose
                    for (auto rr: RRI.edge_covered[j][l]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_E[j][l] - 1) /
                                   (double) (t_max - frac_E[j][l]);
                        Fx += q_R_old - q_R[rr];
                    }
                    frac_E[j][l] += 1;
                    bases[t - 1].emplace_back(j, l);
                    break;
                } else {
                    double value_v = 0;
                    for (auto rr: RRI.edge_covered[j][l]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_E[j][l]);
                    Q_E.k++;
                    heap_push<CMax<double, std::pair<std::pair<unsigned, int>, unsigned>>>(Q_E.k, Q_E.val, Q_E.ids,
                                                                                           value_v, std::make_pair(
                                    std::make_pair(j, l), i));
                }
            }
        }
    }

    if (bound_flag) tight_bound = std::min(tight_bound, Fx + calc_bound_AA(G, k_N, k_E, RRI, q_R));

    rounding_AA(G, bases, frac_N, frac_E, q_R, RRI, t_max, bi_seeds);
    delete[] Q_E.val;
    delete[] Q_E.ids;
    delete[] Q_N.val;
    delete[] Q_N.ids;
    return tight_bound;
}

double get_lower_bound(Graph &G, std::vector<int64> &A, int64 k_N, int64 k_E) {
    const double d0 = log(12.0 * G.n);
    VRRPath R(G, A);
    R.resize1(G, (size_t) 512);
    std::vector<bi_node> bi_seeds;
    CGreedy_AA_PM(G, R, k_N, k_E, 1, bi_seeds);
    auto lowerC = (double) R.self_inf_cal(bi_seeds);
    double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
    //std:: cout << "lower = " << lower << " : " << lowerC << " value = " << lower * (G.n - A.size()) / R.all_R_size << "\n";
    return lower * (G.n - A.size()) / R.all_R_size;
}


double OPIM_AA(Graph &G, std::vector<int64> &A, int64 k_N, int64 k_E, std::vector<bi_node> &bi_seeds, double eps) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / G.n;
    const double approx = 1.0 - 1.0 / exp(1);
    const double approx1 = approx - eps / 2;
    double opt_lower_bound = k_N + k_E;
    int64 slope = G.n - A.size();
    VRRPath R1(G, A), R2(G, A);

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = logcnk(slope, k_N) + logcnk(G.m, k_E);
    double C_max = 8.0 * slope * sqr(
            approx1 * sqrt(log(6.0 / delta)) + sqrt(approx1 * (sum_log + log(6.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = 8.0 * sqr(
            approx1 * sqrt(log(6.0 / delta)) + sqrt(approx1 * (sum_log + log(6.0 / delta)))) / opt_lower_bound;
    cur = clock();
    R1.resize1(G, (size_t) C_0);
    R2.resize1(G, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    int64 a = 2;
    while (1.0 / pow(1.0 + 1.0 / a, a) > 1.0 / exp(1) + eps / 2.0) a++;

    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
//        double upperC_1 = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
//        bi_seeds.clear();
        double upperC = CGreedy_AA_PM(G, R1, k_N, k_E, a, bi_seeds, true);
        double upperC1 = (double) R1.self_inf_cal(bi_seeds) / approx1;
        upperC = std::min(upperC, upperC1);
        auto lowerC = (double) R2.self_inf_cal(bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
       // printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 >= approx - eps || R1.all_R_size >= C_max) break;
        cur = clock();
        int up_rate = a0 < 0.01 ? 32 : ((a0 < (approx - eps) / 2) ? 8 : 2);
        R1.resize1(G, R1.all_R_size * up_rate);
        R2.resize1(G, R2.all_R_size * up_rate);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

double MGGreedy_AA(Graph &G, VRRPath &RRI, int64 k_N, int64 k_E, std::vector<bi_node> &seeds) {
    coveredNum_tmp = new int64[G.n];
    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * sizeof(int64));
    auto **edgeCoveredNum_tmp = new int64 *[G.n]();
    for (int i = 0; i < G.n; ++i) {
        edgeCoveredNum_tmp[i] = new int64[G.deg_in[i]]();
        memcpy(edgeCoveredNum_tmp[i], RRI.edge_coveredNum[i], G.deg_in[i] * sizeof(int64));
    }
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    int64 c_N = 0, c_E = 0;
    std::vector<bool> vis(G.n, false);
    std::vector<std::vector<bool>> vis_edge(G.n);
    for (int i = 0; i < G.n; ++i) {
        vis_edge[i].resize(G.deg_in[i], false);
    }
    int64 influence = 0;
    for (int t = 0; t < k_N + k_E; t++) {
        int i0 = 0, j0 = -1, value_tmp = -1;
        if (c_N < k_N) {
            for (int i = 0; i < G.n; ++i) {
                if (!vis[i] && coveredNum_tmp[i] > value_tmp) {
                    i0 = i;
                    j0 = -1;
                    value_tmp = coveredNum_tmp[i];
                }
            }
        }
        if (c_E < k_E) {
            for (int i = 0; i < G.n; ++i) {
                for (int j = 0; j < G.deg_in[i]; ++j) {
                    if (c_E < k_E && !vis_edge[i][j] && edgeCoveredNum_tmp[i][j] > value_tmp) {
                        i0 = i;
                        j0 = j;
                        value_tmp = edgeCoveredNum_tmp[i][j];
                    }
                }
            }
        }
        assert(c_E + c_N < k_E + k_N);
        influence += value_tmp;
        if (j0 == -1) {
            vis[i0] = true;
            seeds.emplace_back(i0, -1);
            c_N++;
            for (auto RRIndex: RRI.covered[i0]) {
                if (RRSetCovered[RRIndex]) continue;
                for (int l = 0; l < RRI.R[RRIndex].size(); ++l) {
                    auto u0 = RRI.R[RRIndex][l];
                    auto e0 = RRI.R_edge[RRIndex][l];
                    coveredNum_tmp[u0]--;
                    edgeCoveredNum_tmp[u0][e0]--;
                }
                RRSetCovered[RRIndex] = true;
            }
        } else {
            vis_edge[i0][j0] = true;
            seeds.emplace_back(i0, j0);
            c_E++;
            for (auto RRIndex: RRI.edge_covered[i0][j0]) {
                if (RRSetCovered[RRIndex]) continue;
                for (int l = 0; l < RRI.R[RRIndex].size(); ++l) {
                    auto u0 = RRI.R[RRIndex][l];
                    auto e0 = RRI.R_edge[RRIndex][l];
                    coveredNum_tmp[u0]--;
                    edgeCoveredNum_tmp[u0][e0]--;
                }
                RRSetCovered[RRIndex] = true;
            }
        }
    }
    delete[] coveredNum_tmp;
    for (int i = 0; i < G.n; ++i) {
        delete[] edgeCoveredNum_tmp[i];
    }
    delete[] edgeCoveredNum_tmp;
    return (double) influence / RRI.numOfRRsets();
}


double IMM_AA(Graph &G, std::vector<int64> &A, int64 k_N, int64 k_E, std::vector<bi_node> &bi_seeds, double eps) {
    double epsilon1 = eps * sqrt(2);
    double iota = 1.0 + log(2) / log(G.n);
    double LB = 1;
    int64 n_ = G.n - A.size();
    double sum_log = logcnk(n_, k_N) + logcnk(G.m, k_E);
    auto End = (int) (log2(G.n) + 1e-9 - 1);
    VRRPath RRI(G, A);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 1; i <= End; i++) {
        auto ci = (int64) ((2.0 + 2.0 * epsilon1 / 3) * (sum_log + iota * log(n_) + log(log2(n_))) / sqr(epsilon1) *
                           pow(2.0, i));
        //std::cout << "ci:" << ci;
        RRI.resize(G, ci);
        bi_seeds.clear();
        CGreedy_AA(G, RRI, k_N, k_E, 1, bi_seeds);
        double ept = RRI.self_inf_cal(bi_seeds) / RRI.numOfRRsets();
        //std::cout << " aa:" << ept << " ee:" << (1.0 + epsilon1) / pow(2.0, i) << "\n";
        if (ept > (1.0 + epsilon1) / pow(2.0, i)) {
            LB = ept * n_ / (1.0 + epsilon1);
            break;
        }
    }
    double e = exp(1);
    double alpha = sqrt(iota * log(G.n) + log(2));
    double beta = sqrt(0.5 * (sum_log + iota * log(G.n) + log(2)));
    auto C = (int64) (2.0 * G.n * sqr(0.5 * alpha + beta) / LB / sqr(eps));
    std::cout << "C:" << C << "\n";
    RRI.resize(G, C);
    bi_seeds.clear();
    CGreedy_AA(G, RRI, k_N, k_E, 1, bi_seeds);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}


#endif //EXP_AA_H
