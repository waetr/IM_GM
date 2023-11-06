#ifndef EXP_GREEDY_H
#define EXP_GREEDY_H

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"


double calc_bound(std::vector<std::vector<int64>> &candidates, std::vector<double> &q_R, int64 k, RRContainer &RRI) {
    double sum = 0;
    std::vector<std::vector<double>> value(candidates.size());
    for (int i = 0; i < candidates.size(); i++) {
        for (int j = 0; j < candidates[i].size(); ++j) {
            auto v = candidates[i][j];
            double value_v = 0;
            for (long rr: RRI.covered[v]) {
                value_v += q_R[rr];
            }
            value[i].push_back(value_v);
        }
        int64 k_max = std::min((int64) value[i].size(), k);
        std::nth_element(value[i].begin(), value[i].begin() + k_max - 1, value[i].end(), std::greater<>());
        for (int j = 0; j < k_max; ++j) {
            sum += value[i][j];
        }
    }
    return sum;
}

void swap_rounding(std::vector<std::vector<int64>> &candidates, std::vector<std::vector<bi_node>> &bases,
                   std::vector<bi_node> &seeds) {
    std::vector<std::vector<bool>> pos(candidates.size());
    std::vector<std::vector<bool>> pos1(candidates.size());
    for (int i = 0; i < candidates.size(); i++) {
        pos[i].resize(candidates[i].size(), false);
        pos1[i].resize(candidates[i].size(), false);
    }
    for (auto &i: bases[0]) {
        pos[i.first][i.second] = true;
    }
    for (int i = 0; i < bases.size() - 1; ++i) {
        for (int j = 0; j < pos1.size(); j++) {
            for (int k = 0; k < pos1[j].size(); ++k) {
                pos1[j][k] = false;
            }
        }

        for (int j = 0; j < bases[i + 1].size(); ++j) {
            pos1[bases[i + 1][j].first][bases[i + 1][j].second] = true;
        }
        for (int j = 0; j < pos.size(); ++j) {
            int x = 0, y = 0;
            while (x < pos[j].size() && y < pos[j].size()) {
                while (x < pos[j].size() && !(pos[j][x] == true && pos1[j][x] == false)) x++;
                while (y < pos[j].size() && !(pos[j][y] == false && pos1[j][y] == true)) y++;
                if (x == pos[j].size() || y == pos[j].size()) break;
                double xx = random_real();
                if (xx < (double) (i + 1) / (i + 2)) pos1[j][y] = false, pos1[j][x] = true;
                else pos[j][x] = false, pos[j][y] = true;
            }
        }
    }
    for (int i = 0; i < pos.size(); ++i) {
        for (int j = 0; j < pos[i].size(); ++j) {
            if (pos[i][j]) seeds.emplace_back(candidates[i][j], i);
        }
    }
}

void advanced_rounding(std::vector<std::vector<int64>> &candidates, std::vector<std::vector<int64>> &frac_x,
                       std::vector<double> &q_R, RRContainer &RRI,int64 t_max,std::vector<bi_node> &seeds) {
    //the number of partition
    const int64 d = candidates.size();
    //advanced-rounding
    for (int i = 0; i < d; i++) {
        int x = 0;
        while (x < frac_x[i].size() && (frac_x[i][x] == 0 || frac_x[i][x] == t_max)) x++;
        if (x == frac_x[i].size()) continue;
        int y = x + 1;
        while (y < frac_x[i].size()) {
            while (y < frac_x[i].size() && (frac_x[i][y] == 0 || frac_x[i][y] == t_max)) y++;
            double dx = 0, dy = 0;
            for (long rr: RRI.covered[candidates[i][x]]) {
                dx += q_R[rr];
            }
            dx *= (double) t_max / (double) (t_max - frac_x[i][x]);
            for (long rr: RRI.covered[candidates[i][y]]) {
                dy += q_R[rr];
            }
            dy *= (double) t_max / (double) (t_max - frac_x[i][y]);
            if (dx < dy) std::swap(x, y);
            if (t_max - frac_x[i][x] > frac_x[i][y]) {
                for (long rr: RRI.covered[candidates[i][x]]) {
                    q_R[rr] *= (double) (t_max - frac_x[i][x] - frac_x[i][y]) /
                               (double) (t_max - frac_x[i][x]);
                }
                frac_x[i][x] += frac_x[i][y];
                for (long rr: RRI.covered[candidates[i][y]]) {
                    q_R[rr] *= (double) t_max / (double) (t_max - frac_x[i][y]);
                }
                frac_x[i][y] = 0;
                y = std::max(x, y) + 1;
            } else {
                for (long rr: RRI.covered[candidates[i][y]]) {
                    q_R[rr] *= (double) (t_max - (frac_x[i][x] + frac_x[i][y] - t_max)) /
                               (double) (t_max - frac_x[i][y]);
                }
                frac_x[i][y] = frac_x[i][x] + frac_x[i][y] - t_max;
                for (long rr: RRI.covered[candidates[i][x]]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] = 0;
                }
                frac_x[i][x] = t_max;
                int t = x;
                x = y;
                y = std::max(t, y) + 1;
            }
            if (frac_x[i][x] == 0 || frac_x[i][x] == t_max) {
                x = y;
                while (x < frac_x[i].size() && (frac_x[i][x] == 0 || frac_x[i][x] == t_max)) x++;
                if (x >= frac_x[i].size()) break;
                y = x + 1;
            }
        }
    }

    for (int j = 0; j < d; j++) {
        for (int l = 0; l < candidates[j].size(); l++) {
            if (frac_x[j][l] == t_max) seeds.emplace_back(candidates[j][l], j);
        }
    }
}

double
Combined_Greedy(Graph &G, std::vector<std::vector<int64>> &candidates, int64 k, RRContainer &RRI, int64 t_max,
                std::vector<bi_node> &bi_seeds, double delta) {
    //the number of partition
    const int64 d = candidates.size();

    //the fractional solution (use integers to avoid float error)
    std::vector<std::vector<int64>> frac_x(d);
    for (int i = 0; i < d; i++) {
        frac_x[i].resize(candidates[i].size(), 0);
    }
    //temporary varible
    double Fx = 0;
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    int total_size = 0;
    for (int i = 0; i < d; ++i) {
        total_size += candidates[i].size();
    }

    HeapArray<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>> Q{};
    Q.val = new std::pair<double, int64>[total_size];
    Q.ids = new std::pair<unsigned, std::pair<unsigned, unsigned >>[total_size];
    Q.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);

    for (int t = 1; t <= t_max; t++) {
        Q.k = 0;
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < candidates[i].size(); j++) {
                double value_v = 0;
                for (long rr: RRI.covered[candidates[i][j]]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[i][j]);
                Q.k++;
                heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>>(Q.k,
                                                                                                               Q.val,
                                                                                                               Q.ids,
                                                                                                               std::make_pair(
                                                                                                                       value_v,
                                                                                                                       candidates[i][j]),
                                                                                                               std::make_pair(
                                                                                                                       0,
                                                                                                                       std::make_pair(
                                                                                                                               i,
                                                                                                                               j)));
            }
        }
        std::vector<int64> cardinality(d, 0);
        for (int i = 0; i < d * k; i++) {
            while (Q.k != 0) {
                int64 j = Q.ids[0].second.first;
                int64 l = Q.ids[0].second.second;
                int64 it_round = Q.ids[0].first;
                heap_pop<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>>(Q.k,
                                                                                                              Q.val,
                                                                                                              Q.ids);
                Q.k -= 1;
                if (cardinality[j] >= std::min((size_t) k, candidates[j].size())) continue;
                if (it_round == i) {
                    //choose
                    cardinality[j] += 1;
                    for (long rr: RRI.covered[candidates[j][l]]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_x[j][l] - 1) /
                                   (double) (t_max - frac_x[j][l]);
                        Fx += q_R_old - q_R[rr];
                    }

                    frac_x[j][l] += 1;
                    bases[t - 1].emplace_back(j, l);
                    break;
                } else {
                    double value_v = 0;
                    for (long rr: RRI.covered[candidates[j][l]]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_x[j][l]);
                    Q.k++;
                    heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>>(Q.k,
                                                                                                                   Q.val,
                                                                                                                   Q.ids,
                                                                                                                   std::make_pair(
                                                                                                                           value_v,
                                                                                                                           candidates[j][l]),
                                                                                                                   std::make_pair(
                                                                                                                           i,
                                                                                                                           std::make_pair(
                                                                                                                                   j,
                                                                                                                                   l)));
                }
            }
        }
    }

    double tight_bound = Fx + calc_bound(candidates, q_R, k, RRI);

    std::vector<bi_node> temp_seeds;
    int64 coverage = 0;
    int64 D = ceil(32.0 * log(1.0 / delta) * t_max * t_max / (G.n * Fx / RRI.numOfRRsets()));

    for (int i = 0; i < D; ++i) {
        temp_seeds.clear();
        swap_rounding(candidates, bases, temp_seeds);
        int64 temp_coverage = RRI.self_inf_cal(G, temp_seeds);
        if (temp_coverage > coverage) {
            bi_seeds.assign(temp_seeds.begin(), temp_seeds.end());
            coverage = temp_coverage;
        }
    }

    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}
double
Combined_Greedy_Partition(Graph &G, std::vector<std::vector<int64>> &candidates, int64 k, RRContainer &RRI,
                           int64 t_max,
                           std::vector<bi_node> &bi_seeds) {
    //the number of partition
    const int64 d = candidates.size();

    //the fractional solution (use integers to avoid float error)
    std::vector<std::vector<int64>> frac_x(d);
    for (int i = 0; i < d; i++) {
        frac_x[i].resize(candidates[i].size(), 0);
    }
    //temporary varible
    double Fx = 0;
    std::vector<double> q_R(RRI.numOfRRsets(), 1);
    //priority queue for lazy sampling
    auto Q = new HeapArray<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>[d];
    for (size_t i = 0; i < d; i++) {
        Q[i].val = new std::pair<double, int64>[frac_x[i].size()];
        Q[i].ids = new std::pair<unsigned, unsigned>[frac_x[i].size()];
        Q[i].k = 0;
    }

    std::vector<int> xxx(G.n);
    for (int i = 0; i < G.n; ++i) {
        xxx[i] = i;
    }

    for (int t = 1; t <= t_max; t++) {
        for (int i = 0; i < d; i++) {
            Q[i].k = 0;
            for (int j = 0; j < candidates[i].size(); j++) {
                double value_v = 0;
                for (auto rr: RRI.covered[candidates[i][j]]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[i][j]);
                Q[i].k++;
                heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>(Q[i].k, Q[i].val, Q[i].ids,
                                                                                         std::make_pair(value_v,
                                                                                                        xxx[candidates[i][j]]),
                                                                                         std::make_pair(0, j));
            }
        }
        std::vector<int64> cardinality(d, 0);
        for (int i = 0; i < d * k; i++) {
            int j = i % d;
            if (cardinality[j] >= std::min((size_t) k, candidates[j].size()) || Q[j].k == 0) continue;
            while (Q[j].k != 0) {
                int64 l = Q[j].ids[0].second;
                int64 it_round = Q[j].ids[0].first;
                heap_pop<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>(Q[j].k, Q[j].val, Q[j].ids);
                Q[j].k -= 1;
                if (it_round == i) {
                    //choose
                    cardinality[j] += 1;
                    for (long rr: RRI.covered[candidates[j][l]]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_x[j][l] - 1) /
                                   (double) (t_max - frac_x[j][l]);
                        Fx += q_R_old - q_R[rr];
                    }
                    frac_x[j][l] += 1;
                    break;
                } else {
                    double value_v = 0;
                    for (long rr: RRI.covered[candidates[j][l]]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_x[j][l]);
                    Q[j].k++;
                    heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, unsigned>>>(Q[j].k, Q[j].val, Q[j].ids,
                                                                                             std::make_pair(value_v,
                                                                                                            xxx[candidates[j][l]]),
                                                                                             std::make_pair(i, l));
                }
            }
        }
    }

    double tight_bound = Fx + calc_bound(candidates, q_R, k, RRI);

    advanced_rounding(candidates, frac_x, q_R, RRI, t_max, bi_seeds);

    for (size_t i = 0; i < d; i++) {
        delete[] Q[i].val;
        delete[] Q[i].ids;
    }
    delete[] Q;
    return tight_bound;
}

double OPIM_Matroid(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &bi_seeds, double eps) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / graph.n;
    const double approx = 1.0 - 1.0 / exp(1) - eps / 2.0;
    std::vector<bi_node> pre_seeds;
    method_random(graph, k, A, pre_seeds);
    int64 opt_lower_bound = pre_seeds.size();
    RRContainer R1(graph, A, true), R2(graph, A, true);

    std::vector<std::vector<int64>> candidates(A.size());

    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = 0;
    for (auto &candidate: candidates) sum_log += logcnk(candidate.size(), k);
    double C_max = 8.0 * graph.n * sqr(
            approx * sqrt(log(8.0 / delta)) + sqrt(approx * (sum_log + log(8.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = 256 * C_max * eps * eps / graph.n;
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(4.0 * i_max / delta);

    int64 a = ceil(1.0 / eps);
    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
        double d1 = std::min(delta / (4.0 * i_max), delta * eps * eps / graph.n);
//        double upperC_1 = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
//        bi_seeds.clear();
        double upperC = Combined_Greedy(graph, candidates, k, R1, a, bi_seeds, d1);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 >= approx - eps / 2.0 || i == i_max) break;
        cur = clock();
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

double OPIM_Partition(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &bi_seeds, double eps) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / graph.n;
    const double approx = 1.0 - 1.0 / exp(1) - 3.0 * eps / 4.0;
    std::vector<bi_node> pre_seeds;
    method_random(graph, k, A, pre_seeds);
    int64 opt_lower_bound = pre_seeds.size();
    RRContainer R1(graph, A, true), R2(graph, A, true);

    std::vector<std::vector<int64>> candidates(A.size());

    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = 0;
    for (auto &candidate: candidates) sum_log += logcnk(candidate.size(), k);
    double C_max = 32.0 * graph.n * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = C_max * eps * eps * 64 / graph.n;
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    int64 a = 1;
    while (1.0 / pow(1.0 + 1.0 / a, a) > 1.0 / exp(1) + 3.0 * eps / 4.0) a++;
    std::cout << "a=" << a << "\n";
    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
//        double upperC_1 = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
//        bi_seeds.clear();
        double upperC = Combined_Greedy_Partition(graph, candidates, k, R1, a, bi_seeds);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 >= approx - eps / 4.0 || i == i_max) break;
        cur = clock();
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

double OPIM_Partition_1(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &bi_seeds, double eps) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / graph.n;
    const double approx = 1.0 - 1.0 / exp(1) - eps / 2.0;
    std::vector<bi_node> pre_seeds;
    method_random(graph, k, A, pre_seeds);
    int64 opt_lower_bound = pre_seeds.size();
    RRContainer R1(graph, A, true), R2(graph, A, true);

    std::vector<std::vector<int64>> candidates(A.size());

    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = 0;
    for (auto &candidate: candidates) sum_log += logcnk(candidate.size(), k);
    double C_max = 32.0 * graph.n * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = C_max * eps * eps * 64 / graph.n;
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    int64 a = 1;
    while (1.0 / pow(1.0 + 1.0 / a, a) > 1.0 / exp(1) + eps / 2.0) a++;
    std::cout << "a=" << a << "\n";
    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
        double upperC = Combined_Greedy_Partition(graph, candidates, k, R1, a, bi_seeds);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 >= approx - eps / 2.0 || i == i_max) break;
        cur = clock();
        R1.resize(graph, R1.R.size() * 2ll);
        R2.resize(graph, R2.R.size() * 2ll);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

double
naive_RR(Graph &G, std::vector<std::vector<int64>> &candidates, int64 k, RRContainer &RRI,
         std::vector<bi_node> &bi_seeds) {
    //the number of partition
    const int64 d = candidates.size();

    //the fractional solution (use integers to avoid float error)
    std::vector<std::vector<bool>> vis(d);
    for (int i = 0; i < d; i++) {
        vis[i].resize(candidates[i].size(), false);
    }
    //temporary varible
    auto coveredNum_tmp0 = new int64[G.n];
    memcpy(coveredNum_tmp0, RRI.coveredNum, G.n * sizeof(int64));
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    std::vector<int64> cardinalities(d, 0);
    int64 current_influence = 0, tight_bound = INT64_MAX;

    for (int i = 0; i < d * k; ++i) {
        int64 covered_value = -1, u0, t0 = i % d;
        if (t0 == 0) {
            int64 tight_mg = 0;
            std::vector<int64> tight_list;
            for (int j = 0; j < d; j++) {
                tight_list.clear();
                for (auto e: candidates[j]) tight_list.emplace_back(coveredNum_tmp0[e]);
                int64 k_max = std::min((int64) tight_list.size(), k);
                std::nth_element(tight_list.begin(), tight_list.begin() + k_max - 1, tight_list.end(),
                                 std::greater<>());
                for (int l = 0; l < k_max; l++) {
                    tight_mg += tight_list[l];
                }
            }
            tight_bound = std::min(tight_bound, current_influence + tight_mg);
        }
        if (cardinalities[t0] == k || cardinalities[t0] == candidates[t0].size()) continue;
        for (int u = 0; u < candidates[t0].size(); ++u) {
            if (vis[t0][u]) continue;
            if (coveredNum_tmp0[candidates[t0][u]] > covered_value) {
                covered_value = coveredNum_tmp0[candidates[t0][u]];
                u0 = u;
            }
        }
        vis[t0][u0] = true;
        current_influence += coveredNum_tmp0[candidates[t0][u0]];
        bi_seeds.emplace_back(candidates[t0][u0], t0);
        cardinalities[t0] += 1;
        for (auto RRIndex: RRI.covered[candidates[t0][u0]]) {
            if (RRSetCovered[RRIndex]) continue;
            for (auto u: RRI.R[RRIndex]) {
                coveredNum_tmp0[u]--;
            }
            RRSetCovered[RRIndex] = true;
        }
    }
    delete[] coveredNum_tmp0;
    return tight_bound;
}

double OPIM_RR(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &bi_seeds, double eps) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / graph.n;
    const double approx = 0.5;
    std::vector<bi_node> pre_seeds;
    method_random(graph, k, A, pre_seeds);
    int64 opt_lower_bound = pre_seeds.size();
    RRContainer R1(graph, A, true), R2(graph, A, true);

    std::vector<std::vector<int64>> candidates(A.size());

    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = 0;
    for (auto &candidate: candidates) sum_log += logcnk(candidate.size(), k);
    double C_max = 2.0 * graph.n * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = 2.0 * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / opt_lower_bound;
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
//        double upperC_1 = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
//        bi_seeds.clear();
        double upperC = Combined_Greedy_Partition(graph, candidates, k, R1, 1, bi_seeds);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 >= approx - eps || i == i_max) break;
        cur = clock();
        C_0 *= 2;
        R1.resize(graph, (size_t) C_0);
        R2.resize(graph, (size_t) C_0);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

double OPIM_MG(Graph &graph, int64 k, std::vector<int64> &A, std::vector<bi_node> &bi_seeds, double eps) {
    assert(bi_seeds.empty());
    const double delta = 1.0 / graph.n;
    const double approx = 0.5;
    std::vector<bi_node> pre_seeds;
    method_random(graph, k, A, pre_seeds);
    int64 opt_lower_bound = pre_seeds.size();
    RRContainer R1(graph, A, true), R2(graph, A, true);

    std::vector<std::vector<int64>> candidates(A.size());

    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: graph.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = 0;
    for (auto &candidate: candidates) sum_log += logcnk(candidate.size(), k);
    double C_max = 2.0 * graph.n * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / eps / eps /
                   opt_lower_bound;
    double C_0 = 2.0 * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (sum_log + log(6.0 / delta)))) / opt_lower_bound;
    cur = clock();
    R1.resize(graph, (size_t) C_0);
    R2.resize(graph, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    for (int64 i = 1; i <= i_max; i++) {
        bi_seeds.clear();
        cur = clock();
//        double upperC_1 = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
//        bi_seeds.clear();
        double upperC = Combined_Greedy(graph, candidates, k, R1, 1, bi_seeds, 1.0 / graph.n);
        auto lowerC = (double) R2.self_inf_cal(graph, bi_seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.R.size(), lowerC, upperC);
        if (a0 >= approx - eps || i == i_max) break;
        cur = clock();
        C_0 *= 2;
        R1.resize(graph, (size_t) C_0);
        R2.resize(graph, (size_t) C_0);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

#endif //EXP_GREEDY_H
