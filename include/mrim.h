#ifndef EXP_MRIM_H
#define EXP_MRIM_H

#include <random>

#include "IMs.h"
#include "OPIM_new.h"
#include "Heap.h"
#include "mrim_rrset.h"


//void Cross_Round_Node_Selection(Graph &G, MultiRRContainer &RRI, int64 T, int64 k, std::vector<bi_node> &seeds) {
//    coveredNum_tmp = new int64[G.n * T];
//    memcpy(coveredNum_tmp, RRI.coveredNum, G.n * T * sizeof(int64));
//    std::vector<bool> RRSetCovered(RRI.numOfRRsets(), false);
//    std::vector<int> cardinalities(T, 0);
//    std::vector<bool> vis(G.n * T, false);
//    for (int i = 0; i < T * k; ++i) {
//        int64 covered_value = -1, u0, t0;
//        for (int t = 0; t < T; ++t) {
//            if (cardinalities[t] >= k) continue;
//            for (int u = 0; u < G.n; ++u) {
//                if (vis[t * G.n + u]) continue;
//                if (coveredNum_tmp[t * G.n + u] > covered_value ||
//                    (coveredNum_tmp[t * G.n + u] == covered_value && t * G.n + u > t0 * G.n + u0)) {
//                    covered_value = coveredNum_tmp[t * G.n + u];
//                    u0 = u;
//                    t0 = t;
//                }
//            }
//        }
//        vis[t0 * G.n + u0] = true;
//        seeds.emplace_back(u0, t0);
//        cardinalities[t0] += 1;
//        for (auto RRIndex: RRI.covered[t0 * G.n + u0]) {
//            if (RRSetCovered[RRIndex]) continue;
//            for (int t = 0; t < T; ++t) {
//                for (auto u: RRI.multi_R[RRIndex][t]) {
//                    coveredNum_tmp[t * G.n + u]--;
//                }
//            }
//            RRSetCovered[RRIndex] = true;
//        }
//    }
//    delete[] coveredNum_tmp;
//}

double calc_bound_MRIM(Graph &G, int64 T, int64 k, MultiRRContainer &RRI, std::vector<double> &q_R) {
    double sum = 0;
    std::vector<std::vector<double>> value(T);
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < G.n; ++j) {
            double value_v = 0;
            for (auto rr: RRI.covered[i * G.n + j]) {
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

void rounding_MRIM(Graph &G, int64 T, std::vector<std::vector<bi_node>> &bases, std::vector<int64> &frac_x,
                      std::vector<double> &q_R, MultiRRContainer &RRI, int64 t_max, std::vector<bi_node> &seeds) {
    std::vector<std::vector<bool>> pos(T);
    std::vector<std::vector<bool>> pos1(T);
    for (int i = 0; i < T; i++) {
        pos[i].resize(G.n, false);
        pos1[i].resize(G.n, false);
    }
    for (auto &i: bases[0]) {
        //first: round second:node
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
            //j: round
            int x = 0, y = 0;
            while (x < pos[j].size() && y < pos[j].size()) {
                while (x < pos[j].size() && !(pos[j][x] == true && pos1[j][x] == false)) x++;
                while (y < pos[j].size() && !(pos[j][y] == false && pos1[j][y] == true)) y++;
                if (x == pos[j].size() || y == pos[j].size()) break;
                //x, y: node

                // round x and y
                double dx = 0, dy = 0;
                for (auto rr: RRI.covered[G.n * j + x]) {
                    dx += q_R[rr];
                }
                dx *= (double) t_max / (double) (t_max - frac_x[G.n * j + x]);
                for (auto rr: RRI.covered[G.n * j + y]) {
                    dy += q_R[rr];
                }
                dy *= (double) t_max / (double) (t_max - frac_x[G.n * j + y]);
                int x_ = x, y_ = y, swap_flag = 1;
                if (dx < dy) { //|| (dx == dy && random_real() > (double) (i + 1) / (i + 2))
                    std::swap(x_, y_);
                    pos[j][x] = false, pos[j][y] = true;
                    swap_flag = i + 1;
                } else {
                    pos1[j][y] = false, pos1[j][x] = true;
                }
                for (auto rr: RRI.covered[G.n * j + y_]) {
                    q_R[rr] *= (double) (t_max - frac_x[G.n * j + y_] + swap_flag) /
                            (double) (t_max - frac_x[G.n * j + y_]);
                }
                frac_x[G.n * j + y_] -= swap_flag;
                for (auto rr: RRI.covered[G.n * j + x_]) {
                    q_R[rr] *= (double) (t_max - frac_x[G.n * j + x_] - swap_flag) /
                            (double) (t_max - frac_x[G.n * j + x_]);
                }
                frac_x[G.n * j + x_] += swap_flag;
                //end of rounding x and y
            }
        }
    }
    for (int i = 0; i < pos.size(); ++i) {
        for (int j = 0; j < pos[i].size(); ++j) {
            if (pos[i][j]) seeds.emplace_back(j, i);
        }
    }
}

double CGreedy_MRIM(Graph &G, MultiRRContainer &RRI, int64 T, int64 k, int64 t_max, std::vector<bi_node> &bi_seeds) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0;
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    //priority queue for lazy sampling
    HeapArray<CMax<double, std::pair<std::pair<unsigned, unsigned >, unsigned>>> Q{};
    Q.val = new double[G.n * T];
    Q.ids = new std::pair<std::pair<unsigned, unsigned>, unsigned>[G.n * T];
    Q.k = 0;

    std::vector<std::vector<bi_node>> bases(t_max);

    double cur = clock();
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
                heap_push<CMax<double, std::pair<std::pair<unsigned, unsigned >, unsigned>>>(Q.k,Q.val,Q.ids,
                                                                                             value_v, std::make_pair(std::make_pair(i,j),0));
            }
        }
        std::vector<int64> cardinality(T, 0);
        for (int i = 0; i < T * k; i++) {
            while (Q.k != 0) {
                //j: round l:node
                int64 j = Q.ids[0].first.first;
                int64 l = Q.ids[0].first.second;
                int64 it_round = Q.ids[0].second;
                heap_pop<CMax<double, std::pair<std::pair<unsigned, unsigned >, unsigned>>>(Q.k,Q.val,Q.ids);
                Q.k -= 1;
                if (cardinality[j] >= k) continue;
                if (it_round == i) {
                    //choose
                    cardinality[j] += 1;
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
                    heap_push<CMax<double, std::pair<std::pair<unsigned, unsigned>, unsigned>>>(Q.k,Q.val,Q.ids,value_v,std::make_pair(std::make_pair(j,l),i));
                }
            }
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    double tight_bound = Fx + calc_bound_MRIM(G, T, k, RRI, q_R);

    rounding_MRIM(G, T, bases, frac_x, q_R, RRI, t_max, bi_seeds);
    delete[] Q.val;
    delete[] Q.ids;
    return tight_bound;
}

double CGreedy_PM_MRIM(Graph &G, MultiRRContainer &RRI, int64 T, int64 k, int64 t_max, std::vector<bi_node> &bi_seeds) {

    //the fractional solution (use integers to avoid float error)
    std::vector<int64> frac_x(G.n * T, 0);
    //temporary varible
    double Fx = 0;
    std::vector<double> q_R(RRI.numOfRRsets(), 1);

    //priority queue for lazy sampling
    auto Q = new HeapArray<CMax<double, std::pair<unsigned, unsigned >>>[T];
    for (size_t i = 0; i < T; i++) {
        Q[i].val = new double[G.n];
        Q[i].ids = new std::pair<unsigned, unsigned>[G.n];
        Q[i].k = 0;
    }

    std::vector<std::vector<bi_node>> bases(t_max);

    double cur = clock();
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
                heap_push<CMax<double, std::pair<unsigned, unsigned >>>(Q[i].k, Q[i].val, Q[i].ids,
                                                                                         value_v,
                                                                                         std::make_pair(j, 0));
            }
        }
        std::vector<int64> cardinality(T, 0);
        for (int i = 0; i < T * k; i++) {
            int64 j = i % T;
            if (cardinality[j] >= k) continue;
            while (Q[j].k != 0) {
                //j: round l:node
                int64 l = Q[j].ids[0].first;
                int64 it_round = Q[j].ids[0].second;
                heap_pop<CMax<double, std::pair<unsigned, unsigned >>>(Q[j].k,Q[j].val,Q[j].ids);
                Q[j].k -= 1;
                if (it_round == i) {
                    //choose
                    cardinality[j] += 1;
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
                    Q[j].k++;
                    heap_push<CMax<double, std::pair<unsigned, unsigned >>>(Q[j].k,Q[j].val,Q[j].ids,value_v,std::make_pair(l,i));
                }
            }
        }
    }
    std::cout << "alg time = " << (clock() - cur) / CLOCKS_PER_SEC;

    double tight_bound = Fx + calc_bound_MRIM(G, T, k, RRI, q_R);

    rounding_MRIM(G, T, bases, frac_x, q_R, RRI, t_max, bi_seeds);
    for (size_t i = 0; i < T; i++) {
        delete[] Q[i].val;
        delete[] Q[i].ids;
    }
    delete[] Q;
    return tight_bound;
}


//double CR_NAIMM(Graph &G, int64 T, int64 k, double eps, std::vector<bi_node> &seeds) {
//    double iota = 1 + log(2) / log(G.n), eps0 = eps * sqrt(2);
//    double alpha = sqrt(iota * log(G.n) + log(2)), beta = sqrt(0.5 * T * logcnk(G.n, k) + sqr(alpha));
//    double theta = 2.0 * (2.0 + eps0 * 2.0 / 3.0) * (T * logcnk(G.n, k) + iota * log(G.n) + log(log2(G.n))) / eps0 /
//                   eps0;
//    double lambda1 = 2.0 * G.n * T * sqr((1.0 - 1.0 / exp(1)) * alpha + beta) / eps / eps;
//    double LB = 1;
//    int i_max = (int) log2(G.n - 1);
//
//    auto start_time = std::chrono::high_resolution_clock::now();
//
//    MultiRRContainer R(G, T);
//    for (int i = 1; i <= i_max; ++i) {
//        seeds.clear();
//        R.resize(G, (int64) theta);
//        Cross_Round_Node_Selection(G, R, T, k, seeds);
////        std::cout << "theta:" << theta << " inf_call:" << 1.0 * R.self_inf_cal_multi(seeds) / R.numOfRRsets() << " pow:"
////                  << (1.0 + eps0) / pow(2, i) << std::endl;
//        if (1.0 * R.self_inf_cal_multi(seeds) / R.numOfRRsets() >= (1.0 + eps0) / pow(2, i)) {
//            LB = R.self_inf_cal_multi(seeds) * G.n / R.numOfRRsets() / (1.0 + eps0);
//            break;
//        }
//        theta *= 2;
//    }
////    std::cout << "final C=" << (int64) (lambda1 / LB) << "\n";
//    R.resize(G, (int64) (lambda1 / LB));
//    seeds.clear();
//    Cross_Round_Node_Selection(G, R, T, k, seeds);
//    auto end_time = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = end_time - start_time;
//    return elapsed.count();
//}

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
        double upperC = CGreedy_MRIM(G, R1, T, k, a, seeds);
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
