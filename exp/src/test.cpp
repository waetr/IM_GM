#include "IMs.h"
#include "OPIM_new.h"
#include "newgreedy/greedy.h"


using namespace std;

double effic_inf2(Graph &graph, std::vector<bi_node> &S, std::vector<int64> &A) {
    RRContainer R0(graph, A, true);
    R0.resize(graph, 10000);
    return R0.self_inf_cal(graph, S);
}

double effic_inf2(Graph &graph, std::vector<int64> &S, std::vector<int64> &A) {
    RRContainer R0(graph, A, true);
    R0.resize(graph, 10000);
    return R0.self_inf_cal(graph, S);
}

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

    vector<int64> A;
    generate_ap(G, A, 50);
    cout << "overlap:" << estimate_neighbor_overlap(G, A) << endl;
//    RRContainer R(G, A, true);
//    R.resize(G, 100000);
//    cout << "RR sets generated!\n";
//    coveredNum_tmp = new int64[G.n];
//    nodeRemain = new bool[G.n];
//
////    cur = clock();
////    cout << "mg coverage:" << MG_OPIM_Selection(G, A, k, seeds, R, false) << " time=" << clock() - cur;
////    cout << " spread:" << effic_inf2(G, seeds, A) << " size:" << seeds.size() << endl;
////    seeds.clear();
//
//    cur = clock();
//    cout << "RR bound: " << RR_OPIM_Selection(G, A, k, seeds, R, true) << " time=" << clock() - cur;
//    cout << " spread:" << effic_inf2(G, seeds, A) << " size:" << seeds.size() << endl;
//
//    std::vector<std::vector<int64>> candidates(A.size());
//
//    std::set<int64> A_reorder(A.begin(), A.end());
//    for (int i = 0; i < A.size(); i++) {
//        for (auto e: G.g[A[i]])
//            if (A_reorder.find(e.v) == A_reorder.end()) {
//                candidates[i].push_back(e.v);
//            }
//    }
//
//    for (int t = 1; t <= 32; t *= 2) {
//        seeds.clear();
//        cur = clock();
//        cout << "Matroid bound=" << Combined_Greedy(G, candidates, k, R, t, seeds, 1.0 / G.n);
//
//        cout << " CG t=" << t << " coverage=" << R.self_inf_cal(G, seeds) << " time=" << clock() - cur;
//        cout << " spread=" << effic_inf2(G, seeds, A) << " size:" << seeds.size() << endl;
//    }
//
//    cout << endl;
//
//    for (int t = 1; t <= 32; t *= 2) {
//        seeds.clear();
//        cur = clock();
//        cout << "Partition bound=" << Combined_Greedy_Partition(G, candidates, k, R, t, seeds);
//
//        cout << " CG t=" << t << " coverage=" << R.self_inf_cal(G, seeds) << " time=" << clock() - cur;
//        cout << " spread=" << effic_inf2(G, seeds, A) << " size:" << seeds.size() << endl;
//    }
//
//    cout << endl;
//
//    delete[] coveredNum_tmp;
//    delete[] nodeRemain;

    vector<double> eps_batch = {0.4,0.2,0.1,0.05,0.02,0.01};

    for (auto eps_ : eps_batch) {
        cout << "eps=" << eps_ << endl;
        printf("RR-OPIM+:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, eps_, "RR+"));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

        printf("OPIM-Matroid:\n\ttime = %.3f\n", OPIM_Matroid(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

        printf("OPIM-Partition:\n\ttime = %.3f\n", OPIM_Partition(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
    }
    return 0;
}