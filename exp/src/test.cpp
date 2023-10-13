#include "IMs.h"
#include "OPIM_new.h"

using namespace std;

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

    vector<int64> A;
    generate_ap(G, A, 5000);
    RRContainer R(G, A, true);
    R.resize(G, 500000);
    cout << "RR sets generated!\n";
    coveredNum_tmp = new int64[G.n];
    nodeRemain = new bool[G.n];

    auto k = args_k;
    MG_OPIM_Selection(G, A, k, seeds, R, false);
    cout << R.self_inf_cal(G, seeds) << endl;
    seeds.clear();

    RR_OPIM_Selection(G, A, k, seeds, R, false);
    cout << R.self_inf_cal(G, seeds) << endl;

    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return 0;
}