#ifndef EXP_RM_H
#define EXP_RM_H

#include <random>

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"
#include "rm_rrset.h"


void TGreedy_RM(Graph &G, RMRRContainer &RRI, int64 T, double eps, std::vector<bi_node> &seeds) {
    coveredNum_tmp = new int64[G.n * T];
    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * T * sizeof(int64));
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    std::vector<bool> vis(G.n, false);
    double w0 = 0;
    for (int i = 0; i < G.n * T; ++i) w0 = std::max(w0, 1.0 * coveredNum_tmp[i]);
    for (double w = w0; w > eps * w0 / G.n; w *= 1.0 - eps) {
        for (int u = 0; u < G.n; ++u) {
            if (vis[u]) continue;
            for (int t = 0; t < T; ++t) {
                if (coveredNum_tmp[t * G.n + u] >= w) {
                    vis[u] = true;
                    seeds.emplace_back(u, t);
                    for (auto RRIndex: RRI.covered[t * G.n + u]) {
                        if (RRSetCovered[RRIndex]) continue;
                        auto t0 = RRI.multi_R[RRIndex].first;
                        for (auto u0: RRI.multi_R[RRIndex].second) {
                            coveredNum_tmp[t0 * G.n + u0]--;
                        }
                        RRSetCovered[RRIndex] = true;
                    }
                    break;
                }
            }
        }
    }
    delete[] coveredNum_tmp;
}

double calc_bound_RM(Graph &G, int64 T, RMRRContainer &RRI, std::vector<double> &q_R) {
    double sum = 0;
    for (int i = 0; i < G.n; i++) {
        double value = 0;
        for (int j = 0; j < T; ++j) {
            double value_v = 0;
            for (auto rr: RRI.covered[j * G.n + i]) {
                value_v += q_R[rr];
            }
            value = std::max(value, value_v);
        }
        sum += value;
    }
    return sum;
}

void rounding_RM(Graph &G, std::vector<std::vector<bi_node>> &bases, std::vector<int64> &frac_x,
                 std::vector<double> &q_R, RMRRContainer &RRI, int64 t_max, std::vector<bi_node> &seeds) {
    std::vector<int64> pos(G.n);
    std::vector<int64> pos1(G.n);
    for (auto &i: bases[0]) {
        //first: round second:node
        pos[i.second] = i.first;
    }
    for (int i = 0; i < bases.size() - 1; ++i) {
        for (auto &j: bases[i + 1]) {
            //first: round second:node
            pos1[j.second] = j.first;
        }
        for (int j = 0; j < pos.size(); ++j) {
            //j: node
            int x = pos[j], y = pos1[j];
            //x, y: round
            if (x != y) {
                // round (x,j) and (y,j)
                double dx = 0, dy = 0;
                for (auto rr: RRI.covered[G.n * x + j]) {
                    dx += q_R[rr];
                }
                dx *= (double) t_max / (double) (t_max - frac_x[G.n * x + j]);
                for (auto rr: RRI.covered[G.n * y + j]) {
                    dy += q_R[rr];
                }
                dy *= (double) t_max / (double) (t_max - frac_x[G.n * y + j]);
                int x_ = x, y_ = y, swap_flag = 1;
                if (dx < dy) { //|| (dx == dy && random_real() > (double) (i + 1) / (i + 2))
                    std::swap(x_, y_);
                    pos[j] = y;
                    swap_flag = i + 1;
                } else {
                    pos1[j] = x;
                }
                for (auto rr: RRI.covered[G.n * y_ + j]) {
                    q_R[rr] *= (double) (t_max - frac_x[G.n * y_ + j] + swap_flag) /
                               (double) (t_max - frac_x[G.n * y_ + j]);
                }
                frac_x[G.n * y_ + j] -= swap_flag;
                for (auto rr: RRI.covered[G.n * x_ + j]) {
                    q_R[rr] *= (double) (t_max - frac_x[G.n * x_ + j] - swap_flag) /
                               (double) (t_max - frac_x[G.n * x_ + j]);
                }
                frac_x[G.n * x_ + j] += swap_flag;
                //end of rounding x and y
            }
        }
    }
    for (int i = 0; i < pos.size(); ++i) {
        // first : node second : round
        seeds.emplace_back(i, pos[i]);
    }
}

double CGreedy_RM(Graph &G, RMRRContainer &RRI, int64 T, int64 t_max, std::vector<bi_node> &bi_seeds, bool bound_flag = false) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0, tight_bound = RRI.numOfRRsets();
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    //priority queue for lazy sampling
    HeapArray<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>> Q{};
    Q.val = new double[G.n * T];
    Q.ids = new std::pair<std::pair<unsigned, unsigned>, unsigned>[G.n * T];
    Q.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);

    double cur = clock();
    for (int t = 1; t <= t_max; t++) {
        Q.k = 0;
        double sum = 0;
        for (int j = 0; j < G.n; j++) {
            double value_max = 0;
            for (int i = 0; i < T; i++) {
                double value_v = 0;
                for (auto rr: RRI.covered[G.n * i + j]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[G.n * i + j]);
                value_max = std::max(value_max, value_v);
                Q.k++;
                heap_push<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                            value_v, std::make_pair(
                                std::make_pair(i, j), 0));
            }
            sum += value_max;
        }
        if (bound_flag) tight_bound = std::min(tight_bound, Fx + sum);
        std::vector<bool> node_selected(G.n, false);
        for (int i = 0; i < G.n; i++) {
            while (Q.k != 0) {
                //j: round l:node
                int64 j = Q.ids[0].first.first;
                int64 l = Q.ids[0].first.second;
                int64 it_round = Q.ids[0].second;
                heap_pop<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k, Q.val, Q.ids);
                Q.k -= 1;
                if (node_selected[l]) continue;
                if (it_round == i) {
                    //choose
                    node_selected[l] = true;
                    for (auto rr: RRI.covered[G.n * j + l]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_x[G.n * j + l] - 1) /
                                   (double) (t_max - frac_x[G.n * j + l]);
                        Fx += q_R_old - q_R[rr];
                    }
                    frac_x[G.n * j + l] += 1;
                    bases[t - 1].emplace_back(j, l);
                    break;
                } else {
                    double value_v = 0;
                    for (auto rr: RRI.covered[G.n * j + l]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_x[G.n * j + l]);
                    Q.k++;
                    heap_push<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k, Q.val, Q.ids,
                                                                                                value_v, std::make_pair(
                                    std::make_pair(j, l), i));
                }
            }
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    if (bound_flag) tight_bound = std::min(tight_bound, Fx + calc_bound_RM(G, T, RRI, q_R));

    rounding_RM(G, bases, frac_x, q_R, RRI, t_max, bi_seeds);
    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}

double CGreedy_RM_PM(Graph &G, RMRRContainer &RRI, int64 T, int64 t_max, std::vector<bi_node> &bi_seeds, bool bound_flag = false) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0, tight_bound = RRI.numOfRRsets();
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    std::vector<std::vector<bi_node>> bases(t_max);

    double cur = clock();
    for (int t = 1; t <= t_max; t++) {
        if (bound_flag) tight_bound = std::min(tight_bound, Fx + calc_bound_RM(G, T, RRI, q_R));
        for (int i = 0; i < G.n; i++) {
            double value = -1;
            int selected_round = -1;
            for (int j = 0; j < T; j++) {
                double value_v = 0;
                for (auto rr: RRI.covered[G.n * j + i]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[G.n * j + i]);
                if (value_v > value) selected_round = j, value = value_v;
            }
            for (auto rr: RRI.covered[G.n * selected_round + i]) {
                double q_R_old = q_R[rr];
                q_R[rr] *= (double) (t_max - frac_x[G.n * selected_round + i] - 1) /
                           (double) (t_max - frac_x[G.n * selected_round + i]);
                Fx += q_R_old - q_R[rr];
            }
            frac_x[G.n * selected_round + i] += 1;
            bases[t - 1].emplace_back(selected_round, i);
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    if (bound_flag) tight_bound = std::min(tight_bound, Fx + calc_bound_RM(G, T, RRI, q_R));

    rounding_RM(G, bases, frac_x, q_R, RRI, t_max, bi_seeds);
    return tight_bound;
}

double get_lower_bound(Graph &G, int64 T) {
    const double d0 = log(12.0 * G.n);
    RMRRContainer R(G, T);
    R.resize(G, (size_t) 256);
    std::vector<bi_node> bi_seeds;
    CGreedy_RM_PM(G, R, T, 1, bi_seeds, true);
    auto lowerC = (double) R.self_inf_cal_multi(bi_seeds);
    double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
    std:: cout << "lower = " << lower << " : " << lowerC << " value = " << lower * (G.n * T) / R.numOfRRsets() << "\n";
    return lower * (G.n * T) / R.numOfRRsets();
}

double OPIM_RM(Graph &G, int64 T, double eps, std::vector<bi_node> &seeds) {
    const double delta = 1.0 / G.n;
    const double approx = 1.0 - 1.0 / exp(1);
    const double approx1 = approx - eps / 2;
    double opt_lower_bound = get_lower_bound(G, T);
    int64 slope = T*G.n;
    RMRRContainer R1(G, T), R2(G, T);

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double sum_log = G.n * log(T);
    double C_max = 8.0 * slope * sqr(
            approx1 * sqrt(log(12.0 / delta)) + sqrt(approx1 * (sum_log + log(12.0 / delta)))) / eps / eps / opt_lower_bound;
    double C_0 = 8.0 * sqr(
            approx1 * sqrt(log(12.0 / delta)) + sqrt(approx1 * (sum_log + log(12.0 / delta)))) / opt_lower_bound;
    cur = clock();
    R1.resize(G, (size_t) C_0);
    R2.resize(G, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    int64 a = 2;
    while (1.0 / pow(1.0 + 1.0 / a, a) > 1.0 / exp(1) + eps / 2.0) a++;
    for (int64 i = 1; i <= i_max; i++) {
        seeds.clear();
        cur = clock();
        double upperC = CGreedy_RM_PM(G, R1, T, a, seeds, true);
        double upperC1 = (double) R1.self_inf_cal_multi(seeds) / approx1;
        upperC = std::min(upperC, upperC1);
        auto lowerC = (double) R2.self_inf_cal_multi(seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        printf(" a0:%.3f theta0:%zu upperOPT: %.3f lowerCur: %.3f\n", a0, R1.numOfRRsets(), upperC,
               lowerC);
        if (a0 >= approx - eps || R1.numOfRRsets() >= C_max) break;
        cur = clock();
        int up_rate = a0 < 0.01 ? 32 : ((a0 < (approx - eps) / 2) ? 8 : 2);
        R1.resize(G, R1.numOfRRsets() * up_rate);
        R2.resize(G, R2.numOfRRsets() * up_rate);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

void MGGreedy_RM(Graph &G, RMRRContainer &RRI, int64 T, std::vector<bi_node> &seeds) {
    coveredNum_tmp = new int64[G.n * T];
    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * T * sizeof(int64));
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    std::vector<bool> vis(G.n, false);
    for (int i = 0; i < G.n; ++i) {
        int t_, u_, value_tmp = -1;
        for (int u = 0; u < G.n; ++u) {
            if (vis[u]) continue;
            for (int t = 0; t < T; ++t) {
                if (coveredNum_tmp[t * G.n + u] > value_tmp) {
                    t_ = t;
                    u_ = u;
                    value_tmp = coveredNum_tmp[t * G.n + u];
                }
            }
        }
        vis[u_] = true;
        seeds.emplace_back(u_, t_);
        for (auto RRIndex: RRI.covered[t_ * G.n + u_]) {
            if (RRSetCovered[RRIndex]) continue;
            auto t0 = RRI.multi_R[RRIndex].first;
            for (auto u0: RRI.multi_R[RRIndex].second) {
                coveredNum_tmp[t0 * G.n + u0]--;
            }
            RRSetCovered[RRIndex] = true;
        }
    }
    delete[] coveredNum_tmp;
}

void RM_without_oracle(Graph &G, int64 T, double eps, std::vector<bi_node> &seeds) {
    const double delta = 1.0 / G.n;
    const double approx = 0.5;
    RMRRContainer R1(G, T), R2(G, T);

    double C_max = 8.0 * G.n * sqr(
            approx * sqrt(log(16.0 / delta)) + sqrt(approx * (T * G.n + log(16.0 / delta)))) / eps / eps ;
    double C_0 = 1;
    R1.resize(G, (size_t) C_0);
    R2.resize(G, (size_t) C_0);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log((2.0+ T)  / delta / i_max);

    for (int64 i = 1; i <= i_max; i++) {
        seeds.clear();
        MGGreedy_RM(G, R1, T, seeds);
        double upperC = R1.self_inf_cal_multi(seeds) / approx;
        auto lowerC = (double) R2.self_inf_cal_multi(seeds);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        printf(" a0:%.3f theta0:%zu upperOPT: %.3f lowerCur: %.3f\n", a0, R1.numOfRRsets(), upperC,
               lowerC);
        if (a0 >= approx - eps || R1.numOfRRsets() >= C_max) break;
        int up_rate = a0 < 0.01 ? 32 : ((a0 < (approx - eps) / 2) ? 8 : 2);
        R1.resize(G, R1.numOfRRsets() * up_rate);
        R2.resize(G, R2.numOfRRsets() * up_rate);
    }
}


#endif //EXP_RM_H
