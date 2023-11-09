#include "IMs.h"
#include "OPIM_new.h"
#include "rm.h"


using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;
    auto T = atoi(argv[2]);

    RMRRContainer R_judge(G, T);
    R_judge.resize(G, 1000000);

    RMRRContainer R(G, T);
    R.resize(G, 500);

    cout << "RR set generated!\n";

    for (int i = 1; i <= 10; ++i) {
        R.resize(G, R.numOfRRsets() * 2);
        cout << "# of RR sets = " << R.numOfRRsets() << endl;
        for(int t = 1; t <= 16; t *= 2) {
            seeds.clear();
            cur = clock();
            CGreedy_RM(G, R, T, t, seeds);
            cout <<  " CG t = " << t << " time = " << (clock() - cur) / CLOCKS_PER_SEC;
            cout << " value = " << R.self_inf_cal_multi(seeds) << "|" << R_judge.self_inf_cal_multi(seeds) << endl;
            vector<int> num(G.n, 0);
            for (int j = 0; j < seeds.size(); ++j) {
                num[seeds[j].first]++;
            }
            for (int j = 0; j < G.n; ++j) {
                assert(num[i] <= 1);
            }
        }
        for(int t = 1; t <= 16; t *= 2) {
            seeds.clear();
            cur = clock();
            CGreedy_RM_PM(G, R, T, t, seeds);
            cout <<  " CG-PM t = " << t << " time = " << (clock() - cur) / CLOCKS_PER_SEC;
            cout << " value = " << R.self_inf_cal_multi(seeds) << "|" << R_judge.self_inf_cal_multi(seeds) << endl;
            vector<int> num(G.n, 0);
            for (int j = 0; j < seeds.size(); ++j) {
                num[seeds[j].first]++;
            }
            for (int j = 0; j < G.n; ++j) {
                assert(num[i] <= 1);
            }
        }
    }

    return 0;
}