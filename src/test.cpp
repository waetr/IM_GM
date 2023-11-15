#include "IMs.h"
#include "OPIM_new.h"
#include "greedy.h"


using namespace std;

int main(int argc, char const *argv[]) {
    int64 a = 1;
    while (1.0 / pow(1.0 + 1.0 / a, a) > 1.0 / exp(1) + 0.1) a++;
    cout << a;
    return 0;
}