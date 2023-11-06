#include "IMs.h"
#include "OPIM_new.h"
#include "greedy.h"


using namespace std;


int main(int argc, char const *argv[]) {
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds, seeds0;
    auto k = atoi(argv[2]);
    auto A_size = atoi(argv[3]);
    double eps_ = atof(argv[4]);

    vector<int64> A;
    generate_ap(G, A, A_size);

    RRContainer R_judge(G, A, true);
    R_judge.resize(G, 100000);

    std::vector<std::vector<int64>> candidates(A.size());
    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: G.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }


    RRContainer R(G, A, true);
    R.resize(G, 500);
    for (int i = 0; i < 10; ++i) {
        cout << "R_size=" << R.numOfRRsets() << endl;
        R.resize(G, R.numOfRRsets() * 2);
        for (int j = 1; j <= 32; j *= 2) {
            Combined_Greedy(G, candidates, k, R, j, seeds, seeds0, 1.0 / G.n);
            cout << "CG-M(t=" << j << "):" << R_judge.self_inf_cal(G, seeds) << "||" << R_judge.self_inf_cal(G, seeds0) << endl;
            seeds.clear();
            seeds0.clear();
        }
        cout << endl;
    }

    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return 0;
}