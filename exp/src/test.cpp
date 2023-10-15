#include "IMs.h"
#include "OPIM_new.h"


using namespace std;

double
Combined_Greedy(std::vector<std::vector<int64>> &candidates,
                std::vector<std::vector<int64>> &frac_x, int64 k, RRContainer &RRI, int64 t_max) {
    ///temporary varible
    double Fx = 0;
    vector<double> q_R(RRI.numOfRRsets(), 1);
    const int d = candidates.size();

    for (int t = 1; t <= t_max; t++) {
        vector<int64> cardinality(d, 0);
        std::vector<std::vector<bool>> vis(d);
        for (int i = 0; i < d; i++) vis[i].resize(candidates[i].size(), false);
        for (int i = 0; i < d * k; i++) {
            double value = -1;
            int64 candidate_num = -1, candidate_idx = -1;
            //for (int j = 0; j < d; j++)
            {
                int j = i % d;
                if (cardinality[j] >= min((size_t)k, candidates[j].size())) continue;
                for (int l = 0; l < candidates[j].size(); l++) {
                    if (vis[j][l]) continue;
                    double value_v = 0;
                    for (long rr : RRI.covered[candidates[j][l]]) {
                        value_v += q_R[rr];
                    }
                    value_v /= (double) (t_max - frac_x[j][l]) / (double) t_max;
                    if (value_v > value) {
                        candidate_num = j;
                        candidate_idx = l;
                        value = value_v;
                    }
                }
            }
            if (candidate_num == -1) {
                continue;
            }
            cardinality[candidate_num] += 1;
            vis[candidate_num][candidate_idx] = true;
            for (long rr : RRI.covered[candidates[candidate_num][candidate_idx]]) {
                double q_R_old = q_R[rr];
                q_R[rr] *= (double) (t_max - frac_x[candidate_num][candidate_idx] - 1) /
                           (double) (t_max - frac_x[candidate_num][candidate_idx]);
                Fx += q_R_old - q_R[rr];
            }
            frac_x[candidate_num][candidate_idx] += 1;
        }
    }

    for (int i = 0; i < d; i++) {
        int x = 0;
        while (x < frac_x[i].size() && (frac_x[i][x] == 0 || frac_x[i][x] == t_max)) x++;
        if (x == frac_x[i].size()) continue;
        int y = x + 1;
        while (y < frac_x[i].size()) {
            while (y < frac_x[i].size() && (frac_x[i][y] == 0 || frac_x[i][y] == t_max)) y++;
            double dx = 0, dy = 0;
            for (long rr : RRI.covered[candidates[i][x]]) {
                dx += q_R[rr];
            }
            dx /= (double) (t_max - frac_x[i][x]) / t_max;
            for (long rr : RRI.covered[candidates[i][y]]) {
                dy += q_R[rr];
            }
            dy /= (double) (t_max - frac_x[i][y]) / t_max;
            if (dx < dy) std::swap(x, y);
            if (t_max - frac_x[i][x] > frac_x[i][y]) {
                for (long rr : RRI.covered[candidates[i][x]]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] *= (double) (t_max - frac_x[i][x] - frac_x[i][y]) /
                               (double) (t_max - frac_x[i][x]);
                    Fx += q_R_old - q_R[rr];
                }
                frac_x[i][x] += frac_x[i][y];
                for (long rr : RRI.covered[candidates[i][y]]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] *= (double) t_max / (double) (t_max - frac_x[i][y]);
                    Fx += q_R_old - q_R[rr];
                }
                frac_x[i][y] = 0;
                y = max(x, y) + 1;
            } else {
                for (long rr : RRI.covered[candidates[i][y]]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] *= (double) (t_max - (frac_x[i][x] + frac_x[i][y] - t_max)) /
                               (double) (t_max - frac_x[i][y]);
                    Fx += q_R_old - q_R[rr];
                }
                frac_x[i][y] = frac_x[i][x] + frac_x[i][y] - t_max;
                for (long rr : RRI.covered[candidates[i][x]]) {
                    double q_R_old = q_R[rr];
                    q_R[rr] = 0;
                    Fx += q_R_old;
                }
                frac_x[i][x] = t_max;
                int t = x;
                x = y;
                y = max(t, y) + 1;
            }
            if (frac_x[i][x] == 0 || frac_x[i][x] == t_max) {
                x = y;
                while (x < frac_x[i].size() && (frac_x[i][x] == 0 || frac_x[i][x] == t_max)) x++;
                if (x >= frac_x[i].size()) break;
                y = x + 1;
            }
        }
    }
    return Fx;
}

int main(int argc, char const *argv[]) {
    //blog-catalog 5 0.1 0
    if (argc != 5) {
        printf("Usage: [dataset_name] k eps greedy_option");
        return 0;
    }
    string dataset_name = argv[1];
    string graphFilePath = "../data/" + dataset_name + ".txt";
    auto args_k = (int64) atoi(argv[2]);
    auto args_eps = (double) atof(argv[3]);
    auto args_greedy = atoi(argv[4]);
    printf("Start--\n");
    double cur = clock();
    Graph G(graphFilePath);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;
    auto k = args_k;

    vector<int64> A;
    generate_ap(G, A, 1000);
//    vector<int64> A;
//    for (auto e: AB)
//        if (G.deg_out[e] > k) A.push_back(e);
    cout << "overlap:" << estimate_neighbor_overlap(G, A) << endl;
    RRContainer R(G, A, true);
    R.resize(G, 64000);
    cout << "RR sets generated!\n";
    coveredNum_tmp = new int64[G.n];
    nodeRemain = new bool[G.n];


    MG_OPIM_Selection(G, A, k, seeds, R, false);
    cout << R.self_inf_cal(G, seeds) << endl;
    seeds.clear();

    RR_OPIM_Selection(G, A, k, seeds, R, false);
    cout << R.self_inf_cal(G, seeds) << endl;

    std::vector<std::vector<int64>> candidates(A.size());
    std::vector<std::vector<int64>> frac_x(A.size());
    for (int i = 0; i < A.size(); i++) {
        for (auto edge: G.g[A[i]]) {
            if (std::find(A.begin(), A.end(), edge.v) == A.end()) candidates[i].push_back(edge.v);
        }
        assert(candidates.size() >= k);
        frac_x[i].resize(candidates[i].size(), 0);
    }

    for (int t = 1; t <= 5; t++) {
        double fx = Combined_Greedy(candidates, frac_x, k, R, t);
        cout << "t=" << t << " fx=" << fx << " ";

        vector<bi_node> sseed;
        for (int i = 0; i < A.size(); i++) {
            int kk = 0;
            for (int j = 0; j < frac_x[i].size(); j++) {
                if (frac_x[i][j] == t) {
                    kk++;
                    sseed.emplace_back(candidates[i][j], A[i]);
                }
            }
            if (kk != min((size_t)k, candidates[i].size())) {
                cout << "i:" << i << " kk:" << kk << " " << k << endl;
            }
        }
        cout << " fx_=" << R.self_inf_cal(G, sseed) << endl;

        for (int i = 0; i < A.size(); i++) {
            for (int j = 0; j < frac_x[i].size(); j++)
                frac_x[i][j] = 0;
        }
    }

    cout << endl;

    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return 0;
}