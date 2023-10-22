#include "IMs.h"
#include "OPIM_new.h"
#include "greedy.h"


using namespace std;

double effic_inf2(Graph &graph, std::vector<bi_node> &S, std::vector<int64> &A) {
    RRContainer R0(graph, A, true);
    R0.resize(graph, 50000);
    return R0.self_inf_cal(graph, S);
}

double effic_inf2(Graph &graph, std::vector<int64> &S, std::vector<int64> &A) {
    RRContainer R0(graph, A, true);
    R0.resize(graph, 10000);
    return R0.self_inf_cal(graph, S);
}

void max_degree(Graph &G, std::vector<int64> &S, int64 size) {
    vector<pair<int, int>> degrees;
    for (int i = 0; i < G.n; i++) degrees.emplace_back(G.deg_out[i], i);
    std::sort(degrees.begin(), degrees.end());
    reverse(degrees.begin(), degrees.end());
    for (int i = 0; i < size; i++) S.emplace_back(degrees[i].second);
}

void generate_ap1(Graph &graph, std::vector<int64> &A, int64 size = 1) {
    A.clear();
    std::set<int64> S;
    std::uniform_int_distribution<int64> uniformIntDistribution(0, graph.n - 1);
    for (int64 i = 0; i < size; i++) {
        int64 v = uniformIntDistribution(mt19937engine);
        while (std::find(A.begin(), A.end(), v) != A.end()) v = uniformIntDistribution(mt19937engine);
        A.emplace_back(v);
    }
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
    //max_degree(G, A, 1000);
    generate_ap1(G, A, 2000);
    cout << "overlap:" << estimate_neighbor_overlap(G, A) << endl;
//    RRContainer R1(G, A, true);
//    R1.resize(G, 10000);
//    RRContainer R(G, A, true);
//    R.resize(G, 1000);
//    cout << "RR sets generated!\n";
//    coveredNum_tmp = new int64[G.n];
//    nodeRemain = new bool[G.n];
//
//    cur = clock();
//    cout << "mg coverage:" << MG_OPIM_Selection(G, A, k, seeds, R, false) << " time=" << clock() - cur;
//    cout << " spread:" << effic_inf2(G, seeds, A) << " size:" << seeds.size() << endl;
//    seeds.clear();
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
//    for (int t = 1; t <= 5; t ++) {
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
//    for (int t = 1; t <= 5; t += 1) {
//        seeds.clear();
//        cur = clock();
//        cout << "1Partition bound=" << Combined_Greedy_Partition1(G, candidates, k, R, t, seeds);
//
//        cout << " CG t=" << t << " coverage=" << R.self_inf_cal(G, seeds) << " time=" << clock() - cur;
//        cout << " spread=" << effic_inf2(G, seeds, A) << " size:" << seeds.size() << endl;
//    }
//
//    cout << endl;
//
//    delete[] coveredNum_tmp;
//    delete[] nodeRemain;

    vector<double> eps_batch = {0.1};

    for (auto eps_: eps_batch) {
        cout << "eps=" << eps_ << endl;
        printf("RR-OPIM+:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, eps_, "RR+"));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

//        printf("OPIM-Matroid:\n\ttime = %.3f\n", OPIM_Matroid(G, k, A, seeds, eps_));
//        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
//        seeds.clear();

        printf("OPIM-Partition:\n\ttime = %.3f\n", OPIM_Partition(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

        printf("OPIM-Partition1:\n\ttime = %.3f\n", OPIM_Partition1(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
    }
    return 0;
}