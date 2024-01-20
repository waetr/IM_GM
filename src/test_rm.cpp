#include "rm.h"


using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    vector<bi_node> seeds;
    auto T = atoi(argv[2]);
    auto partial_or_all = atoi(argv[3]);

    RMRRContainer R_judge(G, T);
    R_judge.resize(G, 100000);

    if (partial_or_all == 0) {
        //pre-definition
        auto start_time = std::chrono::high_resolution_clock::now();
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        vector<double> time_[9][7], spread[9][7], rrset_time[9];
        for (int i = 0; i < 5; ++i) {
            double rrset_time_tmp = 0;
            printf("%d-th experiment round:\n", i);
            RMRRContainer R(G, T);
            R.resize(G, 32);
            for (int j = 0; j < 9; j++) {
                rrset_time[j].emplace_back(rrset_time_tmp);
                printf("\t# of RR sets = %zu | time = %.3f\n", R.numOfRRsets(), average(rrset_time[j]));
                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_RM(G, R, T, 2, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][0].emplace_back(elapsed.count()), spread[j][0].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
                printf("\tCG-MG-2 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][0]), SD(time_[j][0]),average(spread[j][0]), SD(spread[j][0]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_RM_PM(G, R, T, 2, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][1].emplace_back(elapsed.count()), spread[j][1].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
                printf("\tCG-MGPM-2 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][1]), SD(time_[j][1]),average(spread[j][1]), SD(spread[j][1]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_RM(G, R, T, 8, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][2].emplace_back(elapsed.count()), spread[j][2].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
                printf("\tCG-MG-8 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][2]), SD(time_[j][2]),average(spread[j][2]), SD(spread[j][2]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_RM_PM(G, R, T, 8, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][3].emplace_back(elapsed.count()), spread[j][3].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
                printf("\tCG-MGPM-8 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][3]), SD(time_[j][3]),average(spread[j][3]), SD(spread[j][3]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_RM(G, R, T, 1, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][4].emplace_back(elapsed.count()), spread[j][4].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
                printf("\tMG time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][4]), SD(time_[j][4]),average(spread[j][4]), SD(spread[j][4]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                CGreedy_RM_PM(G, R, T, 1, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][5].emplace_back(elapsed.count()), spread[j][5].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
                printf("\tMG-PM time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][5]), SD(time_[j][5]),average(spread[j][5]), SD(spread[j][5]));

                seeds.clear();
                start_time = std::chrono::high_resolution_clock::now();
                TGreedy_RM(G, R, T, 0.05, seeds);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                time_[j][6].emplace_back(elapsed.count()), spread[j][6].emplace_back(1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets());
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
                auto x = OPIM_RM(G, T, eps,seeds), y = 1.0 * R_judge.self_inf_cal_multi(seeds) * G.n / R_judge.numOfRRsets();
                time_[j][0].push_back(x), spread[j][0].push_back(y);
                printf("\tOURS time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][0]), SD(time_[j][0]),average(spread[j][0]), SD(spread[j][0]));
                seeds.clear();
                x = RM_without_oracle(G, T, eps,seeds), y = 1.0 * R_judge.self_inf_cal_multi(seeds) * G.n / R_judge.numOfRRsets();
                time_[j][1].push_back(x), spread[j][1].push_back(y);
                printf("\tOPIM0 time = %.3f(%.3f) spread = %.3f(%.3f)\n", average(time_[j][1]), SD(time_[j][1]),average(spread[j][1]), SD(spread[j][1]));
            }
        }
    }

    return 0;
}