#include "IMs.h"
#include "OPIM_new.h"
#include "aa.h"


using namespace std;

int main(int argc, char const *argv[]) {
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;
    auto k_N = atoi(argv[2]);
    auto k_T = atoi(argv[3]);

    //set S as the top-100 largest degree nodes
    vector<int64> S;
    vector<pair<int, int>> degree_order;
    for (int i = 0; i < G.n; ++i) degree_order.emplace_back(G.deg_out[i], i);
    std::sort(degree_order.begin(), degree_order.end(), std::greater<>());
    for (int i = 0; i < 50; ++i) S.emplace_back(degree_order[i].second);

    VRRPath R_judge(G, S);
    R_judge.resize(G, 100000);

    VRRPath R(G, S);
    R.resize(G, 50);

    cout << "RR set generated!\n";

    for (int i = 1; i <= 10; ++i) {
        R.resize(G, R.numOfRRsets() * 2);
        cout << "# of RR sets = " << R.numOfRRsets() << endl;
        for(int t = 1; t <= 4; t *= 2) {
            seeds.clear();
            cur = clock();
            CGreedy_AA(G, R, k_N, k_T, t, seeds);
            cout <<  " CG t = " << t << " time = " << (clock() - cur) / CLOCKS_PER_SEC;
            cout << " value = " << R_judge.self_inf_cal(seeds) << endl;
        }
    }

    return 0;
}