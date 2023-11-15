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
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds, emptyset;
    auto k_N = atoi(argv[2]);
    auto k_T = atoi(argv[3]);
    auto k_seed = atoi(argv[4]);

    //set S as the top-100 largest degree nodes
    vector<int64> S;
    vector<pair<int, int>> degree_order;
    for (int i = 0; i < G.n; ++i) degree_order.emplace_back(G.deg_out[i], i);
    std::sort(degree_order.begin(), degree_order.end(), std::greater<>());
    for (int i = 0; i < k_seed; ++i) S.emplace_back(degree_order[i].second);

    VRRPath R_judge(G, S);
    R_judge.resize(G, 100000);

    VRRPath R1(G, S), R2(G, S);
    R1.resize(G, 500);
    R2.resize(G, 500);

    cout << "RR set generated!\n";

    double d0 = log(3.0*log(G.n*100)*G.n);

    for (int i = 1; i <= 11; ++i) {
        R1.resize(G, R1.numOfRRsets() * 2);
        R2.resize(G, R2.numOfRRsets() * 2);
        cout << "# of RR sets = " << R1.numOfRRsets() << endl;
        for(int t = 1; t <= 4; t *= 4) {
            seeds.clear();
            cur = clock();
            auto upperC = CGreedy_AA(G, R1, k_N, k_T, t, seeds);
            double lowerC = 1.0 * R2.self_inf_cal(seeds);
            double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
            double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
            cout <<  " CG t = " << t << " time = " << (clock() - cur) / CLOCKS_PER_SEC;
            cout << " spread = " << R_judge.self_inf_cal(seeds) << " upperC:" << upperC << " cov:" << R1.self_inf_cal(seeds) << " a:" << lower / upper << endl;
        }
    }

//    for (int i = 1; i <= 1; ++i) {
//        R.resize(G, R.numOfRRsets() * 2);
//        cout << "# of RR sets = " << R.numOfRRsets() << " real = " << R.all_R_size << endl;
//        for (int t = 1; t <= 4; t *= 2) {
//            seeds.clear();
//            cur = clock();
//            auto x = CGreedy_AA(G, R, k_N, k_T, t, seeds);
//            cout << " CG t = " << t << " time = " << (clock() - cur) / CLOCKS_PER_SEC;
//            cout << "debug value = " << R.self_inf_cal(seeds) << "|" << x << endl;
//            //cout << " value = " << 1.0 * (G.n - S.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size << endl;
//        }
//    }

//    CGreedy_AA(G, R, k_N, k_T, 4, seeds);
//
//    double A1 = MC_sim(G, S, 500, emptyset) - MC_sim(G, S, 500, seeds);
//    double A2 = 1.0 * (G.n - S.size()) * R_judge.self_inf_cal(seeds) / R_judge.all_R_size;
//    cout << "\n" << A1 << " " << A2 << endl;1

    return 0;
}