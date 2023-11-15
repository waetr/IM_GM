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

    RMRRContainer R1(G, T), R2(G, T);
    R1.resize(G, 500);
    R2.resize(G, 500);

    cout << "RR set generated!\n";

    double d0 = log(3.0*log(G.n*100)*G.n);

    for (int i = 1; i <= 14; ++i) {
        R1.resize(G, R1.numOfRRsets() * 2);
        R2.resize(G, R2.numOfRRsets() * 2);
        cout << "# of RR sets = " << R1.numOfRRsets() << endl;
        for(int t = 1; t <= 4; t *= 4) {
            seeds.clear();
            cur = clock();
            auto upperC = CGreedy_RM(G, R1, T, t, seeds);
            double lowerC = 1.0 * R2.self_inf_cal_multi(seeds);
            double lower = sqr(sqrt(lowerC + 2.0 * d0 / 9.0) - sqrt(d0 / 2.0)) - d0 / 18.0;
            double upper = sqr(sqrt(upperC + d0 / 2.0) + sqrt(d0 / 2.0));
            double upper1 = sqr(sqrt(1.0 * R1.self_inf_cal_multi(seeds) + d0 / 2.0) + sqrt(d0 / 2.0));
            cout <<  " CG t = " << t << " time = " << (clock() - cur) / CLOCKS_PER_SEC;
            cout << " spread = " << R_judge.self_inf_cal_multi(seeds) << " cov:" << R1.self_inf_cal_multi(seeds) << " a:" << lower / upper << " a1:" << lower / upper1 << endl;
        }
    }

    return 0;
}