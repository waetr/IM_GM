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

    cout << "overlap:" << estimate_neighbor_overlap(G, A) << endl;

    {
        cout << "eps=" << eps_ << endl;
        printf("RR-OPIM+:\n\ttime = %.3f\n", method_FOPIM(G, k, A, seeds, eps_, "RR+"));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

        printf("OPIM-Partition:\n\ttime = %.3f\n", OPIM_Partition(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

        printf("OPIM-Partition1:\n\ttime = %.3f\n", OPIM_Partition1(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();

        printf("RR-OPIM-:\n\ttime = %.3f\n", OPIM_RR(G, k, A, seeds, eps_));
        printf("\tsize = %zu\n\tspread = %.3f\n", seeds.size(), effic_inf(G, seeds, A));
        seeds.clear();
    }
    return 0;
}