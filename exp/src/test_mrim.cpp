#include "IMs.h"
#include "OPIM_new.h"
#include "newgreedy/mrim.h"


using namespace std;

int main(int argc, char const *argv[]) {
    //blog-catalog 5 0.1 0
    if (argc != 5) {
        printf("Usage: [dataset_name] k eps greedy_option");
        return 0;
    }
    double cur = clock();
    string dataset_name = argv[1];
    string graphFilePath = "../data/" + dataset_name + ".txt";
    auto args_k = (int64) atoi(argv[2]);
    auto args_eps = (double) atof(argv[3]);
    auto args_greedy = atoi(argv[4]);
    printf("Start--\n");
    Graph G(graphFilePath);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;
    auto k = args_k;
    auto eps = args_eps;

    int64 T = 10;

    RRContainer R(G, T);
    R.resize(G, 5000);
    RRContainer R1(G, T);
    R1.resize(G, 10000);

    cout << "RR set generated!\n";

    cur = clock();
    Cross_Round_Node_Selection(G, R, T, k, seeds);

    cout << "naive greedy:" << R1.self_inf_cal_multi(seeds);
    cout << " time:" << clock() - cur << endl;

    for (int t = 1; t <= 5; ++t) {
        seeds.clear();
        cout << "cgreedy:" << M_CGreedy_Partition(G, R, T, k, t, seeds);
        cout << " time:" << clock() - cur;
        cout << " value:" << R.self_inf_cal_multi(seeds) << "|" << R1.self_inf_cal_multi(seeds) << " t:" << t << endl;
    }

    return 0;
}