#include <random>

#include "IA.h"

using namespace std;

int main(int argc, char const *argv[]) {
    //blog-catalog
    if (argc != 2) {
        printf("Usage: [dataset_name]");
        return 0;
    }
    string dataset_name = argv[1];
    string graphFilePath = "../data/" + dataset_name + ".txt";
    string ApFilePath = "../data/" + dataset_name + ".ap";
    printf("Start--\n");
    Graph G(graphFilePath);
    G.set_diffusion_model(IC); // Only support IC

    vector<pair<int, int>> ss;
    for (int i = 0; i < G.n; i++) {
        ss.emplace_back(G.deg_out[i], i);
    }
    sort(ss.begin(), ss.end());
    reverse(ss.begin(), ss.end());
    vector<int64> seed, empty(0);
    for (int i = 0; i < 50; i++) seed.emplace_back(ss[i].second);
    shuffle(seed.begin(), seed.end(), std::mt19937(std::random_device()()));
    seed.resize(10);
    for(auto u : seed) cout << u << " ";
    cout << endl;
    vector<double> IA(seed.size(), 0);

    IA_simulation(G, seed, IA);

//    //overall spread
//    cout << "Overall: " << MC_simulation(G, seed, empty) << endl;
//    //single spread
//    double sum = 0;
//    vector<int64> one_seed(1);
//    for(auto u : seed) {
//        one_seed[0] = u;
//        sum += MC_simulation(G, one_seed, empty);
//    }
//    cout << "Single spread: " << sum << endl;
//    //degree
//    sum = 0;
//    for(auto u : seed) sum += G.deg_out[u];
//    cout << "Degree: " << sum << endl;
    //pagerank
//    double sum = 0;
//    std::vector<double> pi(G.n, 0);
//    power_iteration(G, pi, 0.2);
//    for(int i = 0; i < G.n; i++) {
//        sum += pi[i];
//    }
//    cout << "Page Rank: " << sum << endl;

    return 0;
}