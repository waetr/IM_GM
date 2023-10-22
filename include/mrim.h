#ifndef EXP_MRIM_H
#define EXP_MRIM_H

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"
#include "mrim_rrset.h"


void Cross_Round_Node_Selection(Graph &G, MultiRRContainer &RRI, int64 T, int64 k, std::vector<bi_node> &seeds) {
    coveredNum_tmp = new int64[G.n * T];
    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * T * sizeof(int64));
    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
    std::vector<int> cardinalities(T, 0);
    std::vector<bool> vis(G.n * T, false);
    std::vector<int64> vis_time(G.n, T + 1);
    for (int i = 0; i < T * k; ++i) {
        int64 covered_value = -1, u0, t0;
        for (int t = 0; t < T; ++t) {
            if (cardinalities[t] >= k) continue;
            for (int u = 0; u < G.n; ++u) {
                if (vis[t * G.n + u]) continue;
                if (coveredNum_tmp[t * G.n + u] > covered_value ||
                    (coveredNum_tmp[t * G.n + u] == covered_value && vis_time[u] > vis_time[u0])) {
                    covered_value = coveredNum_tmp[t * G.n + u];
                    u0 = u;
                    t0 = t;
                }
            }
        }
        vis[t0 * G.n + u0] = true;
        vis_time[u0] -= 1;
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


double M_calc_bound(Graph &G, int64 T, int64 k, MultiRRContainer &RRI, std::vector<double> &q_R) {
    double sum = 0;
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
            sum += value[i][j];
        }
    }
    return sum;
}

void M_swap_rounding(Graph &G, int64 T, std::vector<std::vector<bi_node>> &bases, std::vector<bi_node> &seeds) {
    std::vector<std::vector<bool>> pos(T);
    std::vector<std::vector<bool>> pos1(T);
    for (int i = 0; i < T; i++) {
        pos[i].resize(G.n, false);
        pos1[i].resize(G.n, false);
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
            if (pos[i][j]) seeds.emplace_back(j, i);
        }
    }
}

double
M_CGreedy(Graph &G, MultiRRContainer &RRI, int64 T, int64 k, int64 t_max, std::vector<bi_node> &bi_seeds, double delta) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0;
    std::vector<int64> vis_time(G.n, T * t_max + 1);
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    HeapArray<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>> Q{};
    Q.val = new std::pair<double, int64>[G.n * T];
    Q.ids = new std::pair<unsigned, std::pair<unsigned, unsigned >>[G.n * T];
    Q.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);

    for (int t = 1; t <= t_max; t++) {
        Q.k = 0;
        for (int i = 0; i < T; i++) {
            for (int j = 0; j < G.n; j++) {
                double value_v = 0;
                for (long rr: RRI.covered[G.n * i + j]) {
                    value_v += q_R[rr];
                }
                value_v *= (double) t_max / (double) (t_max - frac_x[G.n * i + j]);
                Q.k++;
                heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>>(Q.k,
                                                                                                               Q.val,
                                                                                                               Q.ids,
                                                                                                               std::make_pair(
                                                                                                                       value_v,
                                                                                                                       vis_time[j]),
                                                                                                               std::make_pair(
                                                                                                                       0,
                                                                                                                       std::make_pair(
                                                                                                                               i,
                                                                                                                               j)));
            }
        }
        std::vector<int64> cardinality(T, 0);
        for (int i = 0; i < T * k; i++) {
            while (Q.k != 0) {
                //j: round l:node
                int64 j = Q.ids[0].second.first;
                int64 l = Q.ids[0].second.second;
                int64 it_round = Q.ids[0].first;
                heap_pop<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>>(Q.k,
                                                                                                              Q.val,
                                                                                                              Q.ids);
                Q.k -= 1;
                if (cardinality[j] >= k) continue;
                if (it_round == i) {
                    //choose
                    cardinality[j] += 1;
                    for (long rr: RRI.covered[G.n * j + l]) {
                        double q_R_old = q_R[rr];
                        q_R[rr] *= (double) (t_max - frac_x[G.n * j + l] - 1) /
                                   (double) (t_max - frac_x[G.n * j + l]);
                        Fx += q_R_old - q_R[rr];
                    }
                    vis_time[l] -= 1;
                    frac_x[G.n * j + l] += 1;
                    bases[t - 1].emplace_back(j, l);
                    break;
                } else {
                    double value_v = 0;
                    for (long rr: RRI.covered[G.n * j + l]) {
                        value_v += q_R[rr];
                    }
                    value_v *= (double) t_max / (double) (t_max - frac_x[G.n * j + l]);
                    Q.k++;
                    heap_push<CMax<std::pair<double, int64>, std::pair<unsigned, std::pair<unsigned, unsigned >>>>(Q.k,
                                                                                                                   Q.val,
                                                                                                                   Q.ids,
                                                                                                                   std::make_pair(
                                                                                                                           value_v,
                                                                                                                           vis_time[l]),
                                                                                                                   std::make_pair(
                                                                                                                           i,
                                                                                                                           std::make_pair(
                                                                                                                                   j,
                                                                                                                                   l)));
                }
            }
        }
    }

    std::cout << "[before rounding:" << Fx << "]";
    double tight_bound = Fx + M_calc_bound(G, T, k, RRI, q_R);

    std::vector<bi_node> temp_seeds;
    int64 coverage = 0;
    int64 D = ceil(32.0 * log(1.0 / delta) * t_max * t_max / (G.n * Fx / RRI.numOfRRsets()));
    std::cout << "D:" << D << " ";

    for (int i = 0; i < D; ++i) {
        temp_seeds.clear();
        M_swap_rounding(G, T, bases, temp_seeds);
        int64 temp_coverage = RRI.self_inf_cal_multi(temp_seeds);
        if (temp_coverage > coverage) {
            bi_seeds.assign(temp_seeds.begin(), temp_seeds.end());
            coverage = temp_coverage;
        }
    }

    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}


double M_CGreedy_Partition(Graph &G, MultiRRContainer &RRI, int64 T, int64 k, int64 t_max, std::vector<bi_node> &bi_seeds) {
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

    double tight_bound = Fx + M_calc_bound(G, T, k, RRI, q_R);

    //advanced-rounding
    for (int i = 0; i < T; i++) {
        std::vector<std::pair<double, int>> gradient;
        for (int j = 0; j < G.n; ++j) {
            if (frac_x[G.n * i + j] > 0 && frac_x[G.n * i + j] < t_max) {
                gradient.emplace_back(0, j);
            }
        }
        while (!gradient.empty()) {
            for (int j = 0; j < gradient.size(); ++j) {
                double dx = 0;
                for (long rr: RRI.covered[G.n * i + gradient[j].second]) {
                    dx += q_R[rr];
                }
                dx *= (double) t_max / (double) (t_max - frac_x[G.n * i + gradient[j].second]);
                gradient[j].first = dx;
            }
            std::sort(gradient.begin(), gradient.end());
            int x = gradient[gradient.size() - 1].second, y = gradient[0].second;
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
            }
            if (frac_x[G.n * i + x] == t_max) {
                gradient.pop_back();
            }
            if (frac_x[G.n * i + y] == 0) {
                std::swap(gradient[0], gradient[gradient.size() - 1]);
                gradient.pop_back();
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

double CR_NAIMM(Graph &G, int64 T, int64 k, double eps, std::vector<bi_node> &seeds) {
    double iota = 1 + log(2) / log(G.n), eps0 = eps * sqrt(2);
    double alpha = sqrt(iota * log(G.n) + log(2)), beta = sqrt(0.5 * T * logcnk(G.n, k) + sqr(alpha));
    double theta = 2.0 * (2.0 + eps0 * 2.0 / 3.0) * (T * logcnk(G.n, k) + iota * log(G.n) + log(log2(G.n))) / eps0 /
                   eps0;
    double lambda1 = 2.0 * G.n * T * sqr((1.0 - 1.0 / exp(1)) * alpha + beta) / eps / eps;
    double LB = 1;
    int i_max = (int) log2(G.n - 1);

    auto start_time = std::chrono::high_resolution_clock::now();

    MultiRRContainer R(G, T);
//    for (int i = 1; i <= i_max; ++i) {
//        seeds.clear();
//        R.resize(G, (int64) theta);
//        Cross_Round_Node_Selection(G, R, T, k, seeds);
//        std::cout << "theta:" << theta << " inf_call:" << 1.0 * R.self_inf_cal_multi(seeds) / R.numOfRRsets() << " pow:"
//                  << (1.0 + eps0) / pow(2, i) << std::endl;
//        if (1.0 * R.self_inf_cal_multi(seeds) / R.numOfRRsets() >= (1.0 + eps0) / pow(2, i)) {
//            LB = R.self_inf_cal_multi(seeds) * G.n / R.numOfRRsets() / (1.0 + eps0);
//            break;
//        }
//        theta *= 2;
//    }
//    std::cout << "final C=" << (int64) (lambda1 / LB) << "\n";
//    R.resize(G, (int64) (lambda1 / LB));
    R.resize(G, 1984);
    seeds.clear();
    Cross_Round_Node_Selection(G, R, T, k, seeds);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

double CR_OPIM_Partition(Graph &G, int64 T, int64 k, double eps, std::vector<bi_node> &seeds) {
    const double delta = 1.0 / G.n;
    const double approx = 1.0 - 1.0 / exp(1) - 3.0 * eps / 4.0;
    MultiRRContainer R1(G, T), R2(G, T);

    auto start_time = std::chrono::high_resolution_clock::now();
    double time1 = 0, time2 = 0, cur;
    double C_max = 32.0 * G.n * sqr(
            approx * sqrt(log(6.0 / delta)) + sqrt(approx * (T * logcnk(G.n, k) + log(6.0 / delta)))) / eps / eps /
                   (T * k);
    double C_0 = C_max * 8 * eps * eps / G.n;
    cur = clock();
    R1.resize(G, (size_t) C_0);
    R2.resize(G, (size_t) C_0);
    time1 += time_by(cur);
    auto i_max = (int64) (log2(C_max / C_0) + 1);
    double d0 = log(3.0 * i_max / delta);

    int64 a = 2;
    while (1.0 / pow(1.0 + 1.0 / a, a) > 1.0 / exp(1) + 3.0 * eps / 4.0) a++;
    std::cout << "a=" << a << "\n";
    for (int64 i = 1; i <= i_max; i++) {
        std::cout << "i=" << i << "\n";
        seeds.clear();
        cur = clock();
//        double upperC_1 = RR_OPIM_Selection(graph, A, k, bi_seeds, R1, true);
//        bi_seeds.clear();
        double upperC = M_CGreedy_Partition(G, R1, T, k, a, seeds);
//                M_CGreedy_Partition(G, R1, T, k, a, seeds);
//        double upperC = R1.self_inf_cal_multi(seeds) / (1.0 - 1.0 / exp(1) - 3.0 * eps / 4.0);
        auto lowerC = (double) R2.self_inf_cal_multi(seeds);
        time2 += time_by(cur);
        double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
        double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
        double a0 = lower / upper;
        //printf("a0:%.3f theta0:%zu lowerC: %.3f upperC: %.3f\n", a0, R1.numOfRRsets(), lowerC, upperC);
        if (a0 >= approx - eps / 4.0 || i == i_max) break;
        cur = clock();
        R1.resize(G, R1.numOfRRsets() * 2ll);
        R2.resize(G, R2.numOfRRsets() * 2ll);
        time1 += time_by(cur);
    }
    printf("time1: %.3f time2: %.3f size: %zu\n", time1, time2, R1.numOfRRsets());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
double effic_inf_multi(Graph &graph, std::vector<bi_node> &S, int64 T) {
    const double delta = 1e-3, eps = 0.01, c = 2.0 * (exp(1.0) - 2.0);
    const double LambdaL = 20000;//1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
    size_t numHyperEdge = 0, numCoverd = 0;
    std::vector<bool> exclusive(graph.n);
    std::vector<std::vector<bool>> vecBoolSeed(T);
    for (int i = 0; i < T; i++) vecBoolSeed[i].resize(graph.n, false);
    for (auto seed: S) vecBoolSeed[seed.second][seed.first] = true;
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    std::vector<int64> nodes;
    std::queue<int64> Q;
    while (numCoverd < LambdaL) {
        bool flag = false;
        numHyperEdge++;
        const auto uStart = uniformIntDistribution(mt19937engine);
        for (int i = 0; i < T; ++i) {
            if (vecBoolSeed[i][uStart]) {
                flag = true;
                break;
            }
        }
        if (flag) {
            // Stop, this sample is covered
            numCoverd++;
            continue;
        }
        for (int i = 0; i < T; ++i) {
            exclusive[uStart] = true;
            Q.push(uStart);
            nodes.emplace_back(uStart);
            while (!Q.empty()) {
                int64 u = Q.front();
                Q.pop();
                for (auto &edgeT: graph.gT[u]) {
                    if (exclusive[edgeT.v]) continue;
                    if (random_real() < edgeT.p) {
                        if (vecBoolSeed[i][edgeT.v]) {
                            numCoverd++;
                            flag = true;
                            break;
                        }
                        exclusive[edgeT.v] = true;
                        Q.push(edgeT.v);
                        nodes.emplace_back(edgeT.v);
                    }
                }
                if (flag) break;
            }
            for (auto e: nodes) exclusive[e] = false;
            nodes.clear();
            Q = std::queue<int64>(); // clear the queue
            if (flag) break;
        }
    }
    return 1.0 * numCoverd * graph.n / numHyperEdge;
}


#endif //EXP_MRIM_H
