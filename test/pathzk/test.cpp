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

int port, party;
const int threads = 1;

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, size_t branch_size, size_t reg_size, size_t step) {

    // initial zk exec
    auto init_start = clock_start();
    setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);
    auto init_time = time_from(init_start);


    // set up randomized CPU (with circuits as instructions)
    std::__1::random_device::result_type cir_seed;
    if (party == ALICE) {
        std::random_device rd; // obtain a random number from hardware
        cir_seed = rd();
        ZKFpExec::zk_exec->send_data(&cir_seed, sizeof(std::__1::random_device::result_type));
    } else {
        ZKFpExec::zk_exec->recv_data(&cir_seed, sizeof(std::__1::random_device::result_type));
    }

    ZKCPU zkcpu(branch_size, reg_size, cir_seed);
    zkcpu.rand_cpu();
    //zkcpu.print();

    std::cout << "CPU has been randomized." << std::endl;

    // prover selects random witness
    if (party == ALICE) {
        

    } else { // V:Bob

    }

    // proof start sync
    if (party == ALICE) {
        int ack = 0;
        ZKFpExec::zk_exec->send_data(&ack, sizeof(int));
    } else {
        int ack;
        ZKFpExec::zk_exec->recv_data(&ack, sizeof(int));
    }
	
	auto start = clock_start();

	IntFp final_res = IntFp(0, PUBLIC);
	batch_reveal_check_zero(&final_res, 1);

	finalize_zk_arith<BoolIO<NetIO>>();
	auto timeuse = time_from(start);	
	cout << init_time+timeuse << " us\t" << party << " " << endl;
	std::cout << std::endl;


#if defined(__linux__)
	struct rusage rusage;
	if (!getrusage(RUSAGE_SELF, &rusage))
		std::cout << "[Linux]Peak resident set size: " << (size_t)rusage.ru_maxrss << std::endl;
	else std::cout << "[Linux]Query RSS failed" << std::endl;
#elif defined(__APPLE__)
	struct mach_task_basic_info info;
	mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
	if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS)
		std::cout << "[Mac]Peak resident set size: " << (size_t)info.resident_size_max << std::endl;
	else std::cout << "[Mac]Query RSS failed" << std::endl;
#endif
}

int main(int argc, char** argv) {
    // std::random_device rd; // obtain a random number from hardware
    // auto xxx = rd();
    
    // ZKCPU zkcpu(10,15,xxx);
    // zkcpu.rand_cpu();
    // //zkcpu.print();
    // //zkcpu.test_eval_f5();
    // zkcpu.test_proof_f61(100);
    // return 0;

	// if(argc < 6) {
	// 	std::cout << "usage: a.out PARTY(1/2) PORT IP DIMENSION #BRANCH" << std::endl;
	// 	return -1;
	// }
	
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO>* ios[threads];
	for(int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE?nullptr:argv[3],port+i), party==ALICE);

	std::cout << std::endl << "------------ circuit zero-knowledge proof test ------------" << std::endl << std::endl;;

	// int num = 0;
	// int branch = 0;
	// if(argc < 3) {
	// 	std::cout << "usage: bin/arith/matrix_mul_arith PARTY PORT DIMENSION" << std::endl;
	// 	return -1;
	// } else if (argc == 3) {
	// 	num = 10;
	// 	branch = 10;
	// } else {
	// 	num = atoi(argv[4]);
	// 	branch = atoi(argv[5]);
	// }
	

	test_circuit_zk(ios, party, 4, 4, 10);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;

}