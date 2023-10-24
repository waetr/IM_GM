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
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);
    printf("Evaluating influence in [0.99,1.01]*EPT with prob. 99.9%%.\n");
    vector<bi_node> seeds;
    auto k = atoi(argv[2]);
    auto A_size = atoi(argv[3]);
    double eps_ = atof(argv[4]);

    vector<int64> A;
    generate_ap(G, A, A_size);

    RRContainer R_judge(G, A, true);
    R_judge.resize(G, 1000000);

    std::vector<std::vector<int64>> candidates(A.size());
    std::set<int64> A_reorder(A.begin(), A.end());
    for (int i = 0; i < A.size(); i++) {
        for (auto e: G.g[A[i]])
            if (A_reorder.find(e.v) == A_reorder.end()) {
                candidates[i].push_back(e.v);
            }
    }
    coveredNum_tmp = new int64[G.n];
    nodeRemain = new bool[G.n];


    RRContainer R(G, A, true);
    R.resize(G, 500);
    for (int i = 0; i < 10; ++i) {
        R.resize(G, R.numOfRRsets() * 2);
        cout << "RR set num=" << R.numOfRRsets() << endl;
        cur = clock();
        RR_OPIM_Selection(G, A, k, seeds, R, false);
        cout << "RR:" << R.self_inf_cal(G, seeds) << "|" << R_judge.self_inf_cal(G, seeds) << " time=" << clock() - cur
             << endl;
        seeds.clear();
        cur = clock();
        naive_RR(G, candidates, k, R, seeds);
        cout << "RR-:" << R.self_inf_cal(G, seeds) << "|" << R_judge.self_inf_cal(G, seeds) << " time=" << clock() - cur
             << endl;
        seeds.clear();
        for (int j = 1; j <= 8; j *= 2) {
            cur = clock();
            Combined_Greedy(G, candidates, k, R, j, seeds, 1.0 / G.n);
            cout << "CG-M(t=" << j << "):" << R.self_inf_cal(G, seeds) << "|" << R_judge.self_inf_cal(G, seeds)
                 << " time=" << clock() - cur << endl;
            seeds.clear();
        }
        for (int j = 1; j <= 8; j *= 2) {
            cur = clock();
            Combined_Greedy_Partition1(G, candidates, k, R, j, seeds);
            cout << "CG(t=" << j << "):" << R.self_inf_cal(G, seeds) << "|" << R_judge.self_inf_cal(G, seeds)
                 << " time=" << clock() - cur << endl;
            seeds.clear();
        }
        for (int j = 1; j <= 8; j *= 2) {
            cur = clock();
            Combined_Greedy_Partition(G, candidates, k, R, j, seeds);
            cout << "CG+(t=" << j << "):" << R.self_inf_cal(G, seeds) << "|" << R_judge.self_inf_cal(G, seeds)
                 << " time=" << clock() - cur << endl;
            seeds.clear();
        }
        cout << endl;
    }

    delete[] coveredNum_tmp;
    delete[] nodeRemain;
    return 0;
}