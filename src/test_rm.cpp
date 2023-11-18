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
//    vector<int> size_batch = {atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7])};

    RMRRContainer R_judge(G, T), R(G, T);
    R_judge.resize(G, 1000000);
    cout << "Judge set generated: " << R_judge.numOfRRsets() << endl;

    vector<double> eps_batch = {0.5, 0.4, 0.3, 0.2, 0.1};
    for (int i = 0; i < 5; i++) {
        cout << "eps = " << eps_batch[i] << endl;
        seeds.clear();
        OPIM_RM(G,T,eps_batch[i], seeds);
        cout << " spread_MG = " << 1.0 * T * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        R.resize(G, size_batch[i]);
//        cout << "eps = " << eps_batch[i] << endl;
//        seeds.clear();
//        CGreedy_RM(G, R, T, 1, seeds);
//        cout << " spread_MG = " << 1.0 * T * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_RM_PM(G, R, T, 1, seeds);
//        cout << " spread_Local = " << 1.0 * T * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        TGreedy_RM(G, R, T, 0.05, seeds);
//        cout << " spread_Threshold = " << 1.0 * T * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_RM_PM(G, R, T,4, seeds);
//        cout << " spread_Local1 = " << 1.0 * T * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
    }

    return 0;
}