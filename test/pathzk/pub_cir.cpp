#include <random>
#include "zkcpu.h"
#include "emp-zk/emp-zk.h"
#include <iostream>
#include "emp-tool/emp-tool.h"
#include "zk-ram/zk-ram.h"
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

uint64_t comm(BoolIO<NetIO> *ios[threads]) {
	uint64_t c = 0;
	for(int i = 0; i < threads; ++i)
		c += ios[i]->counter;
	return c;
}
uint64_t comm2(BoolIO<NetIO> *ios[threads]) {
	uint64_t c = 0;
	for(int i = 0; i < threads; ++i)
		c += ios[i]->io->counter;
	return c;
}

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, string cpu_cfg_file, string input_file, string final_reg_file, string cid_file) {

    // testing communication
    uint64_t com1 = comm(ios);
	uint64_t com11 = comm2(ios);

    // initial zk exec
    auto init_start = clock_start();
    setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);
    auto init_time = time_from(init_start);

    ZKCPU zkcpu(0, 0, 0);
    zkcpu.load_cpu_from_file(cpu_cfg_file);
    size_t reg_size = zkcpu.m;

    size_t size;
    ifstream fin;
    fin.open(input_file);
    std::vector<uint64_t> in_val;
    fin >> size; in_val.resize(size);
    for (int i = 0; i < size; i++) fin >> in_val[i];
    fin.close();
    fin.open(final_reg_file);
    std::vector<uint64_t> final_reg;
    fin >> size; final_reg.resize(size);
    for (int i = 0; i < size; i++) fin >> final_reg[i];
    fin.close();
    fin.open(cid_file);
    std::vector<size_t> cids;
    fin >> size; cids.resize(size);
    for (int i = 0; i < size; i++) fin >> cids[i];
    fin.close();

    cout << in_val.size() << endl;

    size_t cur_in = 0;
    IntFp last_input(0, PUBLIC);
	
	auto start = clock_start();
    std::vector<IntFp> wire; wire.clear();
    for (int i = 0; i < reg_size; i++) wire.push_back(IntFp(0, PUBLIC));
    for (int i = 0; i < cids.size(); i++) {
        size_t cid = cids[i];
        for (int j = reg_size; j < zkcpu.br[cid].bank.size(); j++) {
            if (zkcpu.br[cid].bank[j].op == OPTYPE::INPUT) {
                last_input = IntFp(in_val[cur_in++], ALICE);
                wire.push_back(last_input);
            } else if (zkcpu.br[cid].bank[j].op == OPTYPE::ADD) { 
                wire.push_back(wire[zkcpu.br[cid].bank[j].l] + wire[zkcpu.br[cid].bank[j].r]);
            } else {
                wire.push_back(wire[zkcpu.br[cid].bank[j].l] * wire[zkcpu.br[cid].bank[j].r]);
            }
        }
        IntFp correct_fetch = wire[reg_size - 1] + (PR - cid);
        batch_reveal_check_zero(&correct_fetch, 1);
        std::vector<IntFp> new_wire; new_wire.clear();
        for (int j = wire.size() - reg_size + 1; j < wire.size(); j++) new_wire.push_back(wire[j]);
        new_wire.push_back(last_input);
        wire = new_wire;
    }

    for (int i = 0; i < wire.size(); i++) {
        IntFp check_zero = wire[i] + (PR - final_reg[i]);
        batch_reveal_check_zero(&check_zero, 1);
    }
    
	finalize_zk_arith<BoolIO<NetIO>>();
	auto timeuse = time_from(start);	
	cout << init_time+timeuse << " us\t" << party << " " << endl;
	std::cout << std::endl;

	uint64_t com2 = comm(ios) - com1;
	uint64_t com22 = comm2(ios) - com11;
	std::cout << "communication (B): " << com2 << std::endl;
	std::cout << "communication (B): " << com22 << std::endl;

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
	

	test_circuit_zk(ios, party, string(argv[4]), string(argv[5]), string(argv[6]), string(argv[7]));

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;

}