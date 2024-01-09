#include "IMs.h"
#include "OPIM_new.h"
#include "aa.h"
#include <map>

using namespace std;

double MC_sim(Graph &G, vector<int64> &seeds, int MC_rounds, vector<bi_node> &excludes) {
    vector<bool> exclude_node(G.n, false);
    vector<vector<bool>> exclude_edge(G.n);
    for (int i = 0; i < G.n; ++i) exclude_edge[i].resize(G.deg_out[i], false);
    for (auto &exclude: excludes) {
        if (exclude.second == -1) exclude_node[exclude.first] = true;
        else {
            int j;
            int u = G.gT[exclude.first][exclude.second].v;
            double p = G.gT[exclude.first][exclude.second].p;
            for (j = 0; j < G.g[u].size(); ++j) {
                if (G.g[u][j].v == exclude.first && G.g[u][j].p == p) {
                    exclude_edge[u][j] = true;
                    break;
                }
            }
        }
    }
    double total = 0;
    for (int times = 0; times < MC_rounds; times++) {
        vector<bool> vis(G.n, false);
        vector<double> threshold(G.n);
        for (int i = 0; i < G.n; ++i) threshold[i] = random_real();
        vector<double> threshold_sum(G.n, 0);
        vector<int64> A = seeds;
        vector<int64> PrevNew = A;
        for (auto u: A) {
            vis[u] = true;
        }
        while (!PrevNew.empty()) {
            vector<int64> New;
            for (auto u: PrevNew) {
                for (int i = 0; i < G.g[u].size(); ++i) {
                    if (exclude_edge[u][i]) continue;
                    int64 v = G.g[u][i].v;
                    if (vis[v] || exclude_node[v]) continue;
                    threshold_sum[v] += G.g[u][i].p;
                    assert(G.g[u][i].p == 1.0 / G.deg_in[v]);
                }
            }
            for (auto u: PrevNew) {
                for (int i = 0; i < G.g[u].size(); ++i) {
                    if (exclude_edge[u][i]) continue;
                    int64 v = G.g[u][i].v;
                    if (vis[v] || exclude_node[v]) continue;
                    if (threshold_sum[v] > 1.001) {
                        cerr << "threshold surpassed:" << threshold_sum[v] << "\n";
                        exit(0);
                    }
                    if (threshold_sum[v] > threshold[v]) {
                        vis[v] = true;
                        New.push_back(v);
                    }
                }
            }
            for (auto u: New) A.push_back(u);
            PrevNew = New;
        }
        total += A.size();
    }
    return total / MC_rounds;
}

int main(int argc, char const *argv[]) {
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    vector<bi_node> seeds, emptyset;
    auto k_N = atoi(argv[2]);
    auto k_T = atoi(argv[3]);
    auto partial_or_all = atoi(argv[4]);

    std::minstd_rand fixed_engine(2024);

    //set S as the top-1 largest degree nodes
    set<int64> S_;
    vector<pair<int, int>> degree_order;
    for (int i = 0; i < G.n; ++i) degree_order.emplace_back(G.deg_out[i], i);
    std::sort(degree_order.begin(), degree_order.end(), std::greater<>());
    for (int i = 0; i < 1; ++i) S_.insert(degree_order[i].second);

    for (auto u: S_) {
        for (auto e : G.g[u]) {
            if (distrib(fixed_engine) <= e.p) S_.insert(e.v);
        }
    }
    printf("size(A) = %zu\n", S_.size());
    vector<int64> A;
    A.assign(S_.begin(), S_.end());

    VRRPath R_judge(G, A);
    R_judge.resize(G, 100000);
    printf("Judge set generated.\n");

    if (partial_or_all == 0) {
        //pre-definition
        auto start_time = std::chrono::high_resolution_clock::now();
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        vector<double> time_[9][7], spread[9][7],rrset_time[9];
        for (int i = 0; i < 5; ++i) {
            double rrset_time_tmp = 0;
            printf("%d-th round:\n", i);
            VRRPath R(G, A);
            R.resize(G, 2);
            for (int j = 0; j < 9; j++) {
                rrset_time[j].emplace_back(rrset_time_tmp);
                printf("\t# of RR sets = %zu | time = %.3f\n", R.numOfRRsets(), average(rrset_time[j]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_AA(G, R, k_N, k_T, 2, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][0].emplace_back(elapsed.count()), spread[j][0].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tCG-MG-2 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][0]), SD(time_[j][0]),average(spread[j][0]), SD(spread[j][0]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_AA_PM(G, R, k_N, k_T, 2, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][1].emplace_back(elapsed.count()), spread[j][1].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tCG-MGPM-2 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][1]), SD(time_[j][1]),average(spread[j][1]), SD(spread[j][1]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_AA(G, R, k_N, k_T, 8, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][2].emplace_back(elapsed.count()), spread[j][2].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tCG-MG-8 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][2]), SD(time_[j][2]),average(spread[j][2]), SD(spread[j][2]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_AA_PM(G, R, k_N, k_T, 8, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][3].emplace_back(elapsed.count()), spread[j][3].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tCG-MGPM-8 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][3]), SD(time_[j][3]),average(spread[j][3]), SD(spread[j][3]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_AA(G, R, k_N, k_T, 1, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][4].emplace_back(elapsed.count()), spread[j][4].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tMG time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][4]), SD(time_[j][4]),average(spread[j][4]), SD(spread[j][4]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_AA_PM(G, R, k_N, k_T, 1, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][5].emplace_back(elapsed.count()), spread[j][5].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tMG-PM time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][5]), SD(time_[j][5]),average(spread[j][5]), SD(spread[j][5]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                TRGreedy_AA(G, R, k_N, k_T, 0.05, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][6].emplace_back(elapsed.count()), spread[j][6].emplace_back(1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size);
                printf("\tTR time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][6]), SD(time_[j][6]),average(spread[j][6]), SD(spread[j][6]));
                if(j<8) {
                    start_time = std::chrono::high_resolution_clock::now();
                    R.resize(G, R.numOfRRsets() * 4);
                    end_time = std::chrono::high_resolution_clock::now();
                    elapsed = end_time - start_time;
                    rrset_time_tmp += elapsed.count();
                }
            }
        }
    } else {
        vector<double> eps_batch = {0.5, 0.4, 0.3, 0.2, 0.1};
        vector<double> time_[5][2], spread[5][2];
        for (int i = 0; i < 5; ++i) {
            printf("%d-th round:\n", i);
            for (int j = 0; j < eps_batch.size(); ++j) {
                auto eps = eps_batch[j];
                printf("\teps=%.1f\n", eps);
                seeds.clear();
                auto x = OPIM_AA(G, A, k_N, k_T, seeds, eps), y = 1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size;
                time_[j][0].push_back(x), spread[j][0].push_back(y);
                printf("\tOURS time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][0]), SD(time_[j][0]),average(spread[j][0]), SD(spread[j][0]));
//                seeds.clear();
//                x = IMM_AA(G, A, k_N, k_T, seeds, eps), y = 1.0 * (G.n - A.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size;
//                time_[j][1].push_back(x), spread[j][1].push_back(y);
//                printf("\tIMM time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][1]), SD(time_[j][1]),average(spread[j][1]), SD(spread[j][1]));
            }
        }
    }

    return 0;
}