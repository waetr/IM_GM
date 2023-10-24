#include "IMs.h"
#include "OPIM_new.h"
#include "mrim.h"


using namespace std;

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

    MultiRRContainer R_judge(G, T);
    R_judge.resize(G, 100000);


    MultiRRContainer R(G, T);
    R.resize(G, 50);
    for (int i = 0; i < 10; ++i) {
        R.resize(G, R.numOfRRsets() * 2);
        cout << "RR set num=" << R.numOfRRsets() << endl;
        cur = clock();
        Cross_Round_Node_Selection(G, R, T, k, seeds);
        cout << "RR:" << R.self_inf_cal_multi(seeds) << "|" << R_judge.self_inf_cal_multi(seeds) << " time=" << clock() - cur
             << endl;
        seeds.clear();
        for (int j = 1; j <= 8; j *= 2) {
            cur = clock();
            M_CGreedy(G, R, T, k, j, seeds, 1.0/G.n);
            cout << "M-CG(t=" << j << "):" << R.self_inf_cal_multi(seeds) << "|" << R_judge.self_inf_cal_multi(seeds) << " time=" << clock() - cur
                 << endl;
            seeds.clear();
        }
        for (int j = 1; j <= 8; j *= 2) {
            cur = clock();
            M_CGreedy_Partition(G, R, T, k, j, seeds);
            cout << "CG(t=" << j << "):" << R.self_inf_cal_multi(seeds) << "|" << R_judge.self_inf_cal_multi(seeds) << " time=" << clock() - cur
                 << endl;
            seeds.clear();
        }
        for (int j = 1; j <= 8; j *= 2) {
            cur = clock();
            M_CGreedy_Partition1(G, R, T, k, j, seeds);
            cout << "CG1(t=" << j << "):" << R.self_inf_cal_multi(seeds) << "|" << R_judge.self_inf_cal_multi(seeds) << " time=" << clock() - cur
                 << endl;
            seeds.clear();
        }
        cout << endl;
    }

    return 0;
}