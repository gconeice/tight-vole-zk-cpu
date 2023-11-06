#include <iostream>
#include <random>
#include "zkcpu.h"

using namespace std;

int main() {
    std::random_device rd; // obtain a random number from hardware
    auto xxx = rd();
    
    ZKCPU zkcpu(4,4,xxx);
    zkcpu.rand_cpu();
    //zkcpu.print();
    //zkcpu.test_eval_f5();
    zkcpu.test_proof_f5(4);
    return 0;
}