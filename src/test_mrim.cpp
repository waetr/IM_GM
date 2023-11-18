#include "IMs.h"
#include "OPIM_new.h"
#include "mrim.h"


using namespace std;

double MC_sim(Graph &G, vector<bi_node> &seeds, int T, int MC_rounds) {
    double total = 0;
    vector<vector<int64>> seeds_ordered(T);
    for (auto e: seeds) {
        seeds_ordered[e.second].push_back(e.first);
    }
    for (int times = 0; times < MC_rounds; times++) {
        set<int64> A_total;
        for (int t = 0; t < T; ++t) {
            vector<bool> vis(G.n, false);
            vector<int64> A = seeds_ordered[t];
            vector<int64> PrevNew = A;
            for (auto u: A) {
                vis[u] = true;
            }
            while (!PrevNew.empty()) {
                vector<int64> New;
                for (auto u: PrevNew) {
                    for (auto edge: G.g[u]) {
                        int64 v = edge.v;
                        if (vis[v]) continue;
                        if (random_real() < edge.p) {
                            vis[v] = true;
                            New.push_back(v);
                        }
                    }
                }
                for (auto u: New) A.push_back(u);
                PrevNew = New;
            }
            for (auto u: A) A_total.insert(u);
        }
        total += A_total.size();
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
    vector<bi_node> seeds;
    auto k = atoi(argv[2]);
    auto T = atoi(argv[3]);
    //vector<int> size_batch = {atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])};

    MultiRRContainer R_judge(G, T), R(G, T);
    R_judge.resize(G, 5);
    cout << "Judge set generated: " << R_judge.numOfRRsets() << endl;

    vector<double> eps_batch = {0.5, 0.4, 0.3, 0.2, 0.1};
    for (int i = 0; i < 5; i++) {
        cout << "eps = " << eps_batch[i] << endl;
        seeds.clear();
        OPIM_MRIM(G,T,k,eps_batch[i], seeds);
        cout << " spread_MG = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        R.resize(G, size_batch[i]);
//        cout << "eps = " << eps_batch[i] << endl;
//        seeds.clear();
//        CGreedy_MRIM(G, R, T, k, 1, seeds);
//        cout << " spread_MG = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_PM_MRIM(G, R, T, k, 1, seeds);
//        cout << " spread_Local = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        TGreedy_MRIM(G, R, T, k, 0.05, seeds);
//        cout << " spread_Threshold = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_PM_MRIM(G, R, T, k, 4, seeds);
//        cout << " spread_Local1 = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
    }

    return 0;
}