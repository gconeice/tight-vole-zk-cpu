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
    std::vector<f61> l_val; std::vector<f61> r_val;
    std::vector<f61> in_val;
    std::vector<size_t> cids;
    std::vector<f61> final_reg;
    std::size_t path_length; // this is the only public knowledge
    if (party == ALICE) {
        zkcpu.rand_f61_witness_by_step(step, cids, in_val, l_val, r_val, final_reg);
        path_length = in_val.size();
        ZKFpExec::zk_exec->send_data(&path_length, sizeof(size_t));
        for (int i = 0; i < final_reg.size(); i++) ZKFpExec::zk_exec->send_data(&final_reg[i], sizeof(f61));
    } else { // V:Bob
        ZKFpExec::zk_exec->recv_data(&path_length, sizeof(size_t));
        final_reg.resize(reg_size);
        for (int i = 0; i < final_reg.size(); i++) ZKFpExec::zk_exec->recv_data(&final_reg[i], sizeof(f61));        
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

    // P:Alice commits l,r,o,in
    std::vector<IntFp> in(path_length);
    std::vector<IntFp> l(path_length);
    std::vector<IntFp> r(path_length);
    std::vector<IntFp> o(path_length);

    for (int i = 0; i < path_length; i++) in[i] = IntFp(party == ALICE ? in_val[i].val : 0, ALICE);
    for (int i = 0; i < path_length; i++) l[i] = IntFp(party == ALICE ? l_val[i].val : 0, ALICE);
    for (int i = 0; i < path_length; i++) r[i] = IntFp(party == ALICE ? r_val[i].val : 0, ALICE);
    for (int i = 0; i < path_length; i++) o[i] = l[i] * r[i];

    // P:Alice commits 0-1 program p and index indicator id
    std::vector<IntFp> p(path_length);
    std::vector<IntFp> id(path_length);
    if (party == ALICE) {
        size_t scan_i = 0;
        for (int i = 0; i < cids.size(); i++) {
            uint64_t cid = cids[i];
            size_t sub_path_length = zkcpu.br[cid].m + zkcpu.br[cid].n + 1;
            while (sub_path_length--) {
                p[scan_i] = IntFp(0, ALICE);
                id[scan_i++] = IntFp(cid, ALICE);
            }
            p[scan_i] = IntFp(1, ALICE);
            id[scan_i++] = IntFp(cid, ALICE);
        }

    } else {
        for (int i = 0; i < path_length; i++) {
            p[i] = IntFp(0, ALICE);
            id[i] = IntFp(0, ALICE);
        }
    }

    // V: issue uniform \chi to generate topology vectors
    uint64_t chi;
    if (party == ALICE) {
		ZKFpExec::zk_exec->recv_data(&chi, sizeof(uint64_t));
        chi = chi % PR; // to prevent cheating V
    } else {
		PRG().random_data(&chi, sizeof(uint64_t));
		chi = chi % PR;
		ZKFpExec::zk_exec->send_data(&chi, sizeof(uint64_t));	        
    }
    f61 f61_chi(chi);
    std::vector<f61> topo_vec[branch_size];
    for (int i = 0; i < branch_size; i++) zkcpu.br[i].comp_topo_vec_pselect<f61>(f61_chi, i, topo_vec[i]);    

    // P: commits the (entire) topology vector    

    // TODO: prove o \otimes p = \vec{0}

	// batch_reveal_check_zero(&tmp[0], tmp.size());

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