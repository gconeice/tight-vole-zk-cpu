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

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, size_t branch_size, size_t reg_size, size_t step, size_t cir_size) {

    // testing communication
    uint64_t com1 = comm(ios);
	uint64_t com11 = comm2(ios);

    // initial zk exec
    auto init_start = clock_start();
    setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);
    auto init_time = time_from(init_start);


    // set up randomized CPU (with circuits as instructions)
    std::random_device::result_type cir_seed;
    if (party == ALICE) {
        std::random_device rd; // obtain a random number from hardware
        cir_seed = rd();
        ZKFpExec::zk_exec->send_data(&cir_seed, sizeof(std::random_device::result_type));
    } else {
        ZKFpExec::zk_exec->recv_data(&cir_seed, sizeof(std::random_device::result_type));
    }

    ZKCPU zkcpu(branch_size, reg_size, cir_seed);
    zkcpu.rand_uniform_cpu(cir_size);
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

	// Unit for constant offsets
	// IntFp one = IntFp(1, PUBLIC);
	uint64_t delta = ZKFpExec::zk_exec->get_delta();    
    f61 f61_delta(delta);
	
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
    f61 f61_chi2 = f61_chi * f61_chi;
    std::vector<f61> topo_vec[branch_size];
    for (int i = 0; i < branch_size; i++) zkcpu.br[i].comp_topo_vec_pselect<f61>(f61_chi, i, topo_vec[i]);    

    // P: commits the (entire) topology vector cv
    std::vector<IntFp> cv(path_length * 2);
    if (party == ALICE) {
        size_t scan_i = 0;
        for (int i = 0; i < cids.size(); i++) {
            uint64_t cid = cids[i];
            for (int j = 0; j < topo_vec[cid].size(); j++) cv[scan_i++] = IntFp(topo_vec[cid][j].val, ALICE);
        }
    } else {
        for (int i = 0; i < path_length * 2; i++) cv[i] = IntFp(0, ALICE);
    }

    // prove p \times (1-p) = \vec{0} && last{p} = 1; I.e., p is a 0-1 program
    IntFp notlastp = p[path_length-1].negate() + 1;
    batch_reveal_check_zero(&notlastp, 1);
    if (party == ALICE) {
		uint64_t chal; // random challenge
		ZKFpExec::zk_exec->recv_data(&chal, sizeof(uint64_t));
        chal = chal % PR;

        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 C0 = f61::zero(), C1 = f61::zero();
        
        for (int i = 0; i < path_length; i++) {
            IntFp notp = p[i].negate() + 1;
            f61 f61_p = f61(LOW64(p[i].value));
            f61 f61_notp = f61(LOW64(notp.value));
            C0 += coeff * f61_p * f61_notp;
            if (HIGH64(p[i].value) == 0) {
                C1 += coeff * f61_p;
            } else { // == 1
                C1 += coeff * f61_notp;
            }
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        C1 += f61(HIGH64(random_mask));
        C1 = f61::minor(C1.val);
        C0 += f61(LOW64(random_mask));

        ZKFpExec::zk_exec->send_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->send_data(&C1, sizeof(f61));        

    } else {
		uint64_t chal; // random challenge
		PRG().random_data(&chal, sizeof(uint64_t));
		chal = chal % PR;
		ZKFpExec::zk_exec->send_data(&chal, sizeof(uint64_t));	
        
        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 acc = f61::zero();

        for (int i = 0; i < path_length; i++) {
            IntFp notp = p[i].negate() + 1;
            acc += coeff * f61(LOW64(p[i].value)) * f61(LOW64(notp.value));
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        acc += f61(LOW64(random_mask));
        
        f61 C0, C1;
        ZKFpExec::zk_exec->recv_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->recv_data(&C1, sizeof(f61));  

        // std::cout << (f61_delta*C1 + C0).val << std::endl;
        // std::cout << acc.val << std::endl;
        if ((f61_delta*C1 + C0).val == acc.val) std::cout << "[Check]: p is a 0-1 program" << std::endl;
        else {
            std::cout << "[Cheat]: p is not a 0-1 program" << std::endl;
            exit(-1);
        }
    }

    // ZKUROM (tail-heavy)

    // prove, for all odd positions, (cv-1) \otimes p = \vec{0}; I.e., tail heavy
    if (party == ALICE) {
		uint64_t chal; // random challenge
		ZKFpExec::zk_exec->recv_data(&chal, sizeof(uint64_t));
        chal = chal % PR;

        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 C0 = f61::zero(), C1 = f61::zero();
        
        for (int i = 0; i < path_length; i++) {
            IntFp cvm1 = cv[2*i+1] + (PR - 1);
            f61 f61_p = f61(LOW64(p[i].value));
            f61 f61_cvm1 = f61(LOW64(cvm1.value));
            C0 += coeff * f61_p * f61_cvm1;
            C1 += coeff * f61_p * f61(HIGH64(cvm1.value)) + coeff * f61_cvm1 * f61(HIGH64(p[i].value));
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        C1 += f61(HIGH64(random_mask));
        C1 = f61::minor(C1.val);
        C0 += f61(LOW64(random_mask));

        ZKFpExec::zk_exec->send_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->send_data(&C1, sizeof(f61));        
    } else {
		uint64_t chal; // random challenge
		PRG().random_data(&chal, sizeof(uint64_t));
		chal = chal % PR;
		ZKFpExec::zk_exec->send_data(&chal, sizeof(uint64_t));	
        
        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 acc = f61::zero();

        for (int i = 0; i < path_length; i++) {
            IntFp cvm1 = cv[2*i+1] + (PR - 1);
            acc += coeff * f61(LOW64(p[i].value)) * f61(LOW64(cvm1.value));
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        acc += f61(LOW64(random_mask));
        
        f61 C0, C1;
        ZKFpExec::zk_exec->recv_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->recv_data(&C1, sizeof(f61));  

        // std::cout << (f61_delta*C1 + C0).val << std::endl;
        // std::cout << acc.val << std::endl;
        if ((f61_delta*C1 + C0).val == acc.val) std::cout << "[Check]: each fragment ends with 0" << std::endl;
        else {
            std::cout << "[Cheat]: some fragment ends with non-zero" << std::endl;
            exit(-1);
        }
    }        

    // V issues uniform \gamma to compress vectors to tokens (macs)   
    uint64_t gamma;
    if (party == ALICE) {
		ZKFpExec::zk_exec->recv_data(&gamma, sizeof(uint64_t));
        gamma = gamma % PR; // to prevent cheating V
    } else {
		PRG().random_data(&gamma, sizeof(uint64_t));
		gamma = gamma % PR;
		ZKFpExec::zk_exec->send_data(&gamma, sizeof(uint64_t));	        
    }
    f61 f61_gamma(gamma);
    f61 f61_gamma2 = f61_gamma * f61_gamma;
    std::vector<uint64_t> init_mac(branch_size);
    // compute token
    for (int i = 0; i < branch_size; i++) {
        f61 coeff = f61::unit();
        f61 acc = f61::zero();
        for (int j = 0; j < topo_vec[i].size(); j++) {
            acc += coeff * topo_vec[i][j];
            coeff *= f61_gamma;            
        }
        init_mac[i] = acc.val;
        // cout << "mac " << i << ":" << init_mac[i] << endl;
    }
    ZKROM macrom(branch_size);
    macrom.Public_Setup(init_mac);
    std::vector<IntFp> mac(path_length);
    for (int i = 0; i < path_length; i++) mac[i] = macrom.Access(id[i]);
    macrom.Teardown_Batch_Public(party, 16);
    std::cout << "[Check]: zkrom is performed correctly" << std::endl;

    // prove zkurom is formed correctly
    // linear scan to generate (1,gamma,...,gamma^?,1,...,gamma^?,1...)
    std::vector<IntFp> s_cont(path_length);
    IntFp cur = IntFp(1, PUBLIC);
    for (int i = 0; i < path_length; i++) {
        // s_jump[2*i] = s_jump[2*i+1] = cur;
        s_cont[i] = cur;
        // update the next cur
        IntFp next_cur = p[i] + (p[i].negate() + 1) * cur * f61_gamma2.val;
        cur = next_cur;
    }

    // history: brute force
    // IntFp acc_load_mac = IntFp(0, PUBLIC);
    // IntFp acc_calc_mac = IntFp(0, PUBLIC);    
    // for (int i = 0; i < path_length; i++) {
    //     acc_load_mac = acc_load_mac + mac[i] * p[i];
    //     acc_calc_mac = acc_calc_mac + cv[2*i] * s_cont[i] + cv[2*i+1] * s_cont[i] * gamma;
    //     IntFp diff = p[i] * (acc_load_mac + acc_calc_mac.negate());
    //     batch_reveal_check_zero(&diff, 1);
    // }    

    // switch to more effient version
    if (party == ALICE) {
		uint64_t chal; // random challenge
		ZKFpExec::zk_exec->recv_data(&chal, sizeof(uint64_t));
        chal = chal % PR;

        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 C0 = f61::zero(), C1 = f61::zero(), C2 = f61::zero();
        f61 A0 = f61::zero(), A1 = f61::zero(), A2 = f61::zero();
        
        for (int i = 0; i < path_length; i++) {
            // obtain high low
            f61 h64_mac_i = f61(HIGH64(mac[i].value)), l64_mac_i = f61(LOW64(mac[i].value));
            f61 h64_p_i = f61(HIGH64(p[i].value)), l64_p_i = f61(LOW64(p[i].value));
            f61 h64_scont_i = f61(HIGH64(s_cont[i].value)), l64_scont_i = f61(LOW64(s_cont[i].value));
            f61 h64_cv_first = f61(HIGH64(cv[2*i].value)), l64_cv_first = f61(LOW64(cv[2*i].value));
            f61 h64_cv_second = f61(HIGH64(cv[2*i+1].value)), l64_cv_second = f61(LOW64(cv[2*i+1].value));
            // subtract mac[i] * p[i]
            A0 += f61::minor( (l64_mac_i * l64_p_i).val );
            A1 += f61::minor( (l64_mac_i * h64_p_i + l64_p_i * h64_mac_i).val );
            A2 += f61::minor( (h64_mac_i * h64_p_i).val );
            // add cv[2*i] * s_cont[i]
            A0 += l64_cv_first * l64_scont_i;
            A1 += l64_cv_first * h64_scont_i + l64_scont_i * h64_cv_first;
            A2 += h64_cv_first * h64_scont_i;
            // add cv[2*i+1] * s_cont[i] * gamma
            A0 += l64_cv_second * l64_scont_i * f61_gamma;
            A1 += l64_cv_second * h64_scont_i * f61_gamma + l64_scont_i * h64_cv_second * f61_gamma;
            A2 += h64_cv_second * h64_scont_i * f61_gamma;

            // case on p[i] = 0 or p[i] = 1
            if (h64_p_i.val == 0) { // case 0
                f61 shift = coeff * l64_p_i;
                C2 += shift * A2;
                C1 += shift * f61::minor( A1.val );
                C0 += shift * A0;
            } else { // case 1
                f61 mcoeff = f61::minor( coeff.val );
                C2 += coeff * A1;
                C1 += mcoeff * (A1 * l64_p_i + A0);
                C0 += coeff * A0 * l64_p_i;                
            }
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        size_t block_size = 3;
        uint64_t random_pad[block_size]; memset(random_pad, 0, sizeof(random_pad));
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        // -P.HIGH64 \Delta + P.LOW64 = V.LOW64
        random_pad[0] = LOW64(random_mask);
        random_pad[1] = PR - HIGH64(random_mask);
        uint64_t tmp[block_size+1]; 
        for (int i = 1; i < block_size-1; i++) {
            random_mask = ZKFpExec::zk_exec->get_one_role();                
            tmp[0] = 0;
            uint64_t x = PR - HIGH64(random_mask), M = LOW64(random_mask);
            for (int j = 0; j <= i; j++) tmp[j+1] = mult_mod(x, random_pad[j]);
            for (int j = 0; j <= i; j++) random_pad[j] = add_mod(tmp[j], mult_mod(M, random_pad[j])); 
            random_pad[i+1] = tmp[i+1];  
            random_mask = ZKFpExec::zk_exec->get_one_role();
            random_pad[0] = add_mod(random_pad[0], LOW64(random_mask));
            random_pad[1] = add_mod(random_pad[1], PR - HIGH64(random_mask));
        }
        C0 += f61(random_pad[0]);
        C1 += f61(random_pad[1]);
        C2 += f61(random_pad[2]);

        ZKFpExec::zk_exec->send_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->send_data(&C1, sizeof(f61));        
        ZKFpExec::zk_exec->send_data(&C2, sizeof(f61));    
    } else {
		uint64_t chal; // random challenge
		PRG().random_data(&chal, sizeof(uint64_t));
		chal = chal % PR;
		ZKFpExec::zk_exec->send_data(&chal, sizeof(uint64_t));	
        
        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 acc = f61::zero();
        f61 exp_proof = f61::zero();

        for (int i = 0; i < path_length; i++) {
            // obtain high low
            f61 l64_mac_i = f61(LOW64(mac[i].value));
            f61 l64_p_i = f61(LOW64(p[i].value));
            f61 l64_scont_i = f61(LOW64(s_cont[i].value));
            f61 l64_cv_first = f61(LOW64(cv[2*i].value));
            f61 l64_cv_second = f61(LOW64(cv[2*i+1].value));  
            // subtract mac[i] * p[i]          
            acc += f61::minor( (l64_mac_i * l64_p_i).val );
            // add cv[2*i] * s_cont[i]
            acc += l64_cv_first * l64_scont_i;
            // add cv[2*i+1] * s_cont[i] * gamma
            acc += l64_cv_second * l64_scont_i * gamma;
            // update the expected proof
            exp_proof += coeff * acc * l64_p_i;
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        size_t block_size = 3;
        uint64_t random_pad;
         __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        // V.LOW64
        random_pad = LOW64(random_mask); 
        for (int i = 1; i < block_size-1; i++) {
            random_mask = ZKFpExec::zk_exec->get_one_role();
            random_pad = mult_mod(random_pad, LOW64(random_mask));
            random_mask = ZKFpExec::zk_exec->get_one_role();
            random_pad = add_mod(random_pad, LOW64(random_mask));
        }
        exp_proof += f61(random_pad);

        f61 C0, C1, C2;
        ZKFpExec::zk_exec->recv_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->recv_data(&C1, sizeof(f61));  
        ZKFpExec::zk_exec->recv_data(&C2, sizeof(f61));  

        // std::cout << (f61_delta*f61_delta*C2 + f61_delta*C1 + C0).val << std::endl;
        // std::cout << exp_proof.val << std::endl;
        if ((f61_delta*f61_delta*C2 + f61_delta*C1 + C0).val == exp_proof.val) std::cout << "[Check]: each fragment is from zkurom" << std::endl;
        else {
            std::cout << "[Cheat]: some fragment is wrong" << std::endl;
            exit(-1);
        }
    }    



    // prove o \otimes p = \vec{0}, i.e., every instruction's checking output is 0 
    if (party == ALICE) {
		uint64_t chal; // random challenge
		ZKFpExec::zk_exec->recv_data(&chal, sizeof(uint64_t));
        chal = chal % PR;

        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 C0 = f61::zero(), C1 = f61::zero();
        
        for (int i = 0; i < path_length; i++) {
            f61 f61_p = f61(LOW64(p[i].value));
            f61 f61_o = f61(LOW64(o[i].value));
            C0 += coeff * f61_p * f61_o;
            C1 += coeff * f61_p * f61(HIGH64(o[i].value)) + coeff * f61_o * f61(HIGH64(p[i].value));
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        C1 += f61(HIGH64(random_mask));
        C1 = f61::minor(C1.val);
        C0 += f61(LOW64(random_mask));

        ZKFpExec::zk_exec->send_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->send_data(&C1, sizeof(f61));        
    } else {
		uint64_t chal; // random challenge
		PRG().random_data(&chal, sizeof(uint64_t));
		chal = chal % PR;
		ZKFpExec::zk_exec->send_data(&chal, sizeof(uint64_t));	
        
        f61 f61_chal(chal);
        f61 coeff = f61::unit();
        f61 acc = f61::zero();

        for (int i = 0; i < path_length; i++) {
            acc += coeff * f61(LOW64(p[i].value)) * f61(LOW64(o[i].value));
            coeff *= f61_chal;
        }

        // mask the proofs with random_mask
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        acc += f61(LOW64(random_mask));
        
        f61 C0, C1;
        ZKFpExec::zk_exec->recv_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->recv_data(&C1, sizeof(f61));  

        // std::cout << (f61_delta*C1 + C0).val << std::endl;
        // std::cout << acc.val << std::endl;
        if ((f61_delta*C1 + C0).val == acc.val) std::cout << "[Check]: each instruction checking's output is 0" << std::endl;
        else {
            std::cout << "[Cheat]: some instruction checking's output is non-zero" << std::endl;
            exit(-1);
        }
    }    

    // prove (1,chi,\ldots) \times M \times (in,o,in,o,\ldots) = (1,chi,\ldots) \times (l,r,l,r,\ldots)
    // This is performed by showing () \times (1,0,0,...) + (chi^? \odot (1,1,1,...,chi^?,...) \otimes cv) \times (in,o,in,o,...) = (1,chi,chi^2,...) \times (l,r,l,r,...,1,1,1,reg[0],1,reg[1],...)    
    
    // compute add_l and shift_l
    f61 add_l = f61::unit() + f61_chi;
    f61 shift_l = f61_chi2;
    for (int i = 0; i < reg_size; i++) {
        add_l += shift_l;
        shift_l *= f61_chi2;
    }

    // generate shift_l \odot (1,1,1,...,chi^?,...) \otimes cv
    // linear scan to generate (1,1,1,...,chi^?,...)
    std::vector<IntFp> s_jump(path_length);
    cur = IntFp(1, PUBLIC);
    f61 tmp_chi = f61_chi2;
    for (int i = 0; i < path_length; i++) {
        // s_jump[2*i] = s_jump[2*i+1] = cur;
        s_jump[i] = cur * shift_l.val;
        // update the next cur
        IntFp next_cur = p[i] * tmp_chi.val + (p[i].negate() + 1) * cur;
        cur = next_cur;
        tmp_chi *= f61_chi2;
    }

    // compute sum_r (inner-product on the right hand side)
    IntFp sum_r(0, PUBLIC);
    tmp_chi = f61::unit();
    for (int i = 0; i < path_length; i++) {
        sum_r = sum_r + l[i] * tmp_chi.val; tmp_chi *= f61_chi;
        sum_r = sum_r + r[i] * tmp_chi.val; tmp_chi *= f61_chi;
    }
    // add by add_r
    sum_r = sum_r + tmp_chi.val; tmp_chi *= f61_chi;
    sum_r = sum_r + tmp_chi.val; tmp_chi *= f61_chi;
    for (int i = 0; i < reg_size; i++) {
        sum_r = sum_r + tmp_chi.val; tmp_chi *= f61_chi;
        sum_r = sum_r + (tmp_chi * final_reg[i]).val; tmp_chi *= f61_chi;
    }
    sum_r = sum_r + (PR - add_l.val);
    sum_r = sum_r.negate();

    // history   
    // IntFp sum_l(0, PUBLIC);
    // for (int i = 0; i < path_length; i++) {
    //     sum_l = sum_l + s_jump[i] * cv[2*i] * in[i];
    //     sum_l = sum_l + s_jump[i] * cv[2*i+1] * o[i];
    // }            
    // IntFp iszero = sum_l + sum_r;
    // batch_reveal_check_zero(&iszero, 1);
    // (generalized) inner product optimization 
    if (party == ALICE) {
        f61 C0 = f61::zero(), C1 = f61::zero(), C2 = f61::zero();        
        for (int i = 0; i < path_length; i++) {
            // obtain high low
            f61 h64_sjump_i = f61(HIGH64(s_jump[i].value)), l64_sjump_i = f61(LOW64(s_jump[i].value));
            f61 h64_cv_first = f61(HIGH64(cv[2*i].value)), l64_cv_first = f61(LOW64(cv[2*i].value));
            f61 h64_cv_second = f61(HIGH64(cv[2*i+1].value)), l64_cv_second = f61(LOW64(cv[2*i+1].value));
            f61 h64_in_i = f61(HIGH64(in[i].value)), l64_in_i = f61(LOW64(in[i].value));
            f61 h64_o_i = f61(HIGH64(o[i].value)), l64_o_i = f61(LOW64(o[i].value));
            // add s_jump[i] * cv_first * in[i]
            C0 += l64_sjump_i * l64_cv_first * l64_in_i;
            C1 += h64_sjump_i * l64_cv_first * l64_in_i + h64_cv_first * l64_sjump_i * l64_in_i + h64_in_i * l64_sjump_i * l64_cv_first;
            C2 += l64_sjump_i * h64_cv_first * h64_in_i + l64_cv_first * h64_sjump_i * h64_in_i + l64_in_i * h64_sjump_i * h64_cv_first;            
            // add s_jump[i] * cv_second * o[i]
            C0 += l64_sjump_i * l64_cv_second * l64_o_i;
            C1 += h64_sjump_i * l64_cv_second * l64_o_i + h64_cv_second * l64_sjump_i * l64_o_i + h64_o_i * l64_sjump_i * l64_cv_second;
            C2 += l64_sjump_i * h64_cv_second * h64_o_i + l64_cv_second * h64_sjump_i * h64_o_i + l64_o_i * h64_sjump_i * h64_cv_second;
        }
        C1 = f61::minor( C1.val );
        C2 += f61(LOW64(sum_r.value));

        // mask the proofs with random_mask
        size_t block_size = 3;
        uint64_t random_pad[block_size]; memset(random_pad, 0, sizeof(random_pad));
        __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        // -P.HIGH64 \Delta + P.LOW64 = V.LOW64
        random_pad[0] = LOW64(random_mask);
        random_pad[1] = PR - HIGH64(random_mask);
        uint64_t tmp[block_size+1]; 
        for (int i = 1; i < block_size-1; i++) {
            random_mask = ZKFpExec::zk_exec->get_one_role();                
            tmp[0] = 0;
            uint64_t x = PR - HIGH64(random_mask), M = LOW64(random_mask);
            for (int j = 0; j <= i; j++) tmp[j+1] = mult_mod(x, random_pad[j]);
            for (int j = 0; j <= i; j++) random_pad[j] = add_mod(tmp[j], mult_mod(M, random_pad[j])); 
            random_pad[i+1] = tmp[i+1];  
            random_mask = ZKFpExec::zk_exec->get_one_role();
            random_pad[0] = add_mod(random_pad[0], LOW64(random_mask));
            random_pad[1] = add_mod(random_pad[1], PR - HIGH64(random_mask));
        }
        C0 += f61(random_pad[0]);
        C1 += f61(random_pad[1]);
        C2 += f61(random_pad[2]);

        ZKFpExec::zk_exec->send_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->send_data(&C1, sizeof(f61));        
        ZKFpExec::zk_exec->send_data(&C2, sizeof(f61));    
    } else {        
        f61 exp_proof = f61::zero();
        for (int i = 0; i < path_length; i++) {
            // obtain high low
            f61 l64_sjump_i = f61(LOW64(s_jump[i].value));
            f61 l64_cv_first = f61(LOW64(cv[2*i].value));
            f61 l64_cv_second = f61(LOW64(cv[2*i+1].value));
            f61 l64_in_i = f61(LOW64(in[i].value));
            f61 l64_o_i = f61(LOW64(o[i].value));
            // add s_jump[i] * cv_first * in[i]
            exp_proof += l64_sjump_i * l64_cv_first * l64_in_i;
            // add s_jump[i] * cv_second * o[i]
            exp_proof += l64_sjump_i * l64_cv_second * l64_o_i;
        }
        exp_proof += f61(LOW64(sum_r.value)) * f61_delta * f61_delta;

        // mask the proofs with random_mask
        size_t block_size = 3;
        uint64_t random_pad;
         __uint128_t random_mask = ZKFpExec::zk_exec->get_one_role();
        // V.LOW64
        random_pad = LOW64(random_mask); 
        for (int i = 1; i < block_size-1; i++) {
            random_mask = ZKFpExec::zk_exec->get_one_role();
            random_pad = mult_mod(random_pad, LOW64(random_mask));
            random_mask = ZKFpExec::zk_exec->get_one_role();
            random_pad = add_mod(random_pad, LOW64(random_mask));
        }
        exp_proof += f61(random_pad);

        f61 C0, C1, C2;
        ZKFpExec::zk_exec->recv_data(&C0, sizeof(f61));
        ZKFpExec::zk_exec->recv_data(&C1, sizeof(f61));  
        ZKFpExec::zk_exec->recv_data(&C2, sizeof(f61));  

        // std::cout << (f61_delta*f61_delta*C2 + f61_delta*C1 + C0).val << std::endl;
        // std::cout << exp_proof.val << std::endl;
        if ((f61_delta*f61_delta*C2 + f61_delta*C1 + C0).val == exp_proof.val) std::cout << "[Check]: the entire circuit is satisfied" << std::endl;
        else {
            std::cout << "[Cheat]: the entire circuit is not satisfied" << std::endl;
            exit(-1);
        }
    }    

    

	// batch_reveal_check_zero(&tmp[0], tmp.size());

	finalize_zk_arith<BoolIO<NetIO>>();
	auto timeuse = time_from(start);	
    cout << "Path Length: " << path_length << std::endl;
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
    size_t branch_size = atoi(argv[4]);
    size_t reg_size = atoi(argv[5]);
    size_t step = atoi(argv[6]);
    size_t cir_size = atoi(argv[7]);
	

	test_circuit_zk(ios, party, branch_size, reg_size, step, cir_size);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;

}