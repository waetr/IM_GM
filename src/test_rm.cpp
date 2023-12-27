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
    auto partial_or_all = atoi(argv[3]);

    if (partial_or_all == 0) {

    } else {
        vector<double> eps_batch = {0.5, 0.4, 0.3, 0.2, 0.1};
        OPIM_RM(G,T,0.3,seeds);
        RM_without_oracle(G,T,0.3,seeds);
    };



//    RMRRContainer R_judge(G, T), R(G, T);
//    R_judge.resize(G, 1000000);
//    R.resize(G, 1);
//    cout << "Judge set generated: " << R_judge.numOfRRsets() << endl;
//
//    for (int i = 0; i < 23; i++) {
//        R.resize(G, R.numOfRRsets() * 2);
//        cout << "# of RR sets = " << R.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_RM(G, R, T, 4, seeds);
//        cout << " spread_CG-MG = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_RM_PM(G, R, T, 4, seeds);
//        cout << " spread_CG-PM = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_RM(G, R, T, 1, seeds);
//        cout << " spread_MG = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        CGreedy_RM_PM(G, R, T, 1, seeds);
//        cout << " spread_Local = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//        seeds.clear();
//        TGreedy_RM(G, R, T, 0.05, seeds);
//        cout << " spread_Threshold = " << 1.0 * G.n * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
//
//    }

    return 0;
}