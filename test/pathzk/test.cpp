#include <random>
#include "zkcpu.h"
#include "emp-zk/emp-zk.h"
#include <iostream>
#include "emp-tool/emp-tool.h"
#if defined(__linux__)
	#include <sys/time.h>
	#include <sys/resource.h>
#elif defined(__APPLE__)
	#include <unistd.h>
	#include <sys/resource.h>
	#include <mach/mach.h>
#endif

using namespace emp;
using namespace std;

int main() {
    std::random_device rd; // obtain a random number from hardware
    auto xxx = rd();
    
    ZKCPU zkcpu(10,15,xxx);
    zkcpu.rand_cpu();
    //zkcpu.print();
    //zkcpu.test_eval_f5();
    zkcpu.test_proof_f61(100);
    return 0;
}