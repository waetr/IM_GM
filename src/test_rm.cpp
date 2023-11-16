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
    cout << "Judge set generated: " << R_judge.numOfRRsets() << endl;

    vector<double> eps_batch = {0.5, 0.4, 0.3, 0.2, 0.1};
    for (double eps : eps_batch) {
        cout << "eps = " << eps << endl;
        OPIM_RM(G, T, eps, seeds);
        cout << "final spread = " << 1.0 * G.n * T * R_judge.self_inf_cal_multi(seeds) / R_judge.numOfRRsets() << endl;
    }

    return 0;
}