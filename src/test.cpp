#include "IMs.h"
#include "OPIM_new.h"
#include "greedy.h"


using namespace std;

std::minstd_rand mt19937engine1(123);

int main(int argc, char const *argv[]) {
    double cur = clock();
    printf("Start--\n");
    Graph G(argv[1]);
    G.set_diffusion_model(IC); // Only support IC
    printf("read time = %.3f n=%ld m=%ld\n", time_by(cur), G.n, G.m);

    std::uniform_int_distribution<int64> r1(0, G.n-1);
    std::uniform_real_distribution<double> r2(0.0, 1.0);

    ofstream fout("twitter_weight.txt");
    for (int64 i = 0; i < 30000000; ++i) {
            fout << r1(mt19937engine1) % 100000 + 1000000000 << " " << r1(mt19937engine1) + 1000000000 << " " << r2(mt19937engine1) << endl;
    }
    fout.close();
}