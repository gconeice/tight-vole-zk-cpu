#ifndef IMPROVE_ZK_RAM
#define IMPROVE_ZK_RAM

#include "emp-zk/emp-zk.h"
#include <iostream>
#include "emp-tool/emp-tool.h"
#include "inset_rom.h"
// #include "ram_util.h"
#include <vector>
#if defined(__linux__)
	#include <sys/time.h>
	#include <sys/resource.h>
#elif defined(__APPLE__)
	#include <unistd.h>
	#include <sys/resource.h>
	#include <mach/mach.h>
#endif

using namespace std;
using namespace emp;

struct RAMTuple{
    IntFp idx; // the address or index
    IntFp timestamp; // which time it is
    IntFp val; // the value inside
};

class IZKRAM {
public:
    uint64_t N; // size of the ram
    uint64_t T; // number of access
    std::vector<RAMTuple> read_list;
    std::vector<RAMTuple> write_list;
    // inset ROM for bound_check
    InsetZKROM *bound_check_rom;
    // Alice/P uses latest_pos to find the latest place of an index
    std::vector<uint64_t> latest_pos; 
    IntFp one;

    IZKRAM(uint64_t N, uint64_t total_T) {
        this->N = N;
        this->T = 0;
        bound_check_rom = new InsetZKROM(total_T);        
        one = IntFp(1, PUBLIC);
    }

    void Setup(std::vector<uint64_t>& init_val) {        
        RAMTuple tmp;
        for (int i = 0; i < N; i++) {
            tmp.idx = IntFp(i, PUBLIC);
            tmp.timestamp = IntFp(0, PUBLIC);
            tmp.val = IntFp(init_val[i], ALICE);
            write_list.push_back(tmp);
            latest_pos.push_back(i);
        }
        // Setup bound_check inset_rom
        bound_check_rom->Setup();
    }

    IntFp Access(IntFp& id, IntFp& val, IntFp& rw) {
        uint64_t addr = HIGH64(id.value);
        RAMTuple tmp_w, tmp_r;
        tmp_r.idx = tmp_w.idx = id;
        uint64_t pos = latest_pos[addr];
        // two inputs
        IntFp old_timestamp = IntFp(HIGH64(write_list[pos].timestamp.value), ALICE);
        IntFp old_val = IntFp(HIGH64(write_list[pos].val.value), ALICE);
        tmp_r.timestamp = old_timestamp;
        tmp_r.val = old_val;
        tmp_w.timestamp = IntFp(++T, PUBLIC);
        tmp_w.val = old_val + rw * (val + old_val.negate());
        read_list.push_back(tmp_r);
        write_list.push_back(tmp_w);
        // boundcheck access
        IntFp diff = tmp_w.timestamp + tmp_r.timestamp.negate();
        bound_check_rom->Access(diff);
        // Update the lastest pos and increase the T's counter
        latest_pos[addr] = N + T - 1;
        return old_val;
    }

    void Teardown_Basic(int &party) {
        RAMTuple tmp;
        for (int i = 0; i < N; i++) {
            tmp.idx = IntFp(i, PUBLIC);
            uint64_t pos = latest_pos[i]; // the last pos
            tmp.timestamp = IntFp(HIGH64(write_list[pos].timestamp.value), ALICE);            
            tmp.val = IntFp(HIGH64(write_list[pos].val.value), ALICE);
            read_list.push_back(tmp);
            // boundcheck access
            // bound_check_rom->Access(IntFp(T+1, PUBLIC) - tmp.timestamp);
            // Dicuss this with Dave
        }

        // bound check
        bound_check_rom->Teardown_Basic(party);

        uint64_t A0, A1, A2; // to compress the tuple into one IntFp
        uint64_t X; // evaluate at this value
        if (party == ALICE) {
            ZKFpExec::zk_exec->recv_data(&A0, sizeof(uint64_t));
            ZKFpExec::zk_exec->recv_data(&A1, sizeof(uint64_t));
            ZKFpExec::zk_exec->recv_data(&A2, sizeof(uint64_t));
            ZKFpExec::zk_exec->recv_data(&X, sizeof(uint64_t));
        } else {
            PRG tmpprg;
		    tmpprg.random_data(&A0, sizeof(uint64_t));
		    A0 = A0 % PR;
		    tmpprg.random_data(&A1, sizeof(uint64_t));
		    A1 = A1 % PR;
		    tmpprg.random_data(&A2, sizeof(uint64_t));
		    A2 = A2 % PR;
            tmpprg.random_data(&X, sizeof(uint64_t));
            X = X % PR;
		    ZKFpExec::zk_exec->send_data(&A0, sizeof(uint64_t));  
            ZKFpExec::zk_exec->send_data(&A1, sizeof(uint64_t));  
            ZKFpExec::zk_exec->send_data(&A2, sizeof(uint64_t));  
            ZKFpExec::zk_exec->send_data(&X, sizeof(uint64_t));  
        }
	    
        IntFp prod_read = IntFp(1, PUBLIC);
        IntFp prod_write = IntFp(1, PUBLIC);
        for (int i = 0; i < N+T; i++) {
            prod_read = prod_read * (read_list[i].idx * A0 + read_list[i].timestamp * A1 + read_list[i].val * A2 + X);
            prod_write = prod_write * (write_list[i].idx * A0 + write_list[i].timestamp * A1 + write_list[i].val * A2 + X);
        }
        // check they are equal, namely, the same is zero
        IntFp final_zero = prod_read + prod_write.negate();
	    batch_reveal_check_zero(&final_zero, 1);
        
    }

    void Teardown_Batch(int &party, int &block_size) {
        RAMTuple tmp;
        for (int i = 0; i < N; i++) {
            tmp.idx = IntFp(i, PUBLIC);
            uint64_t pos = latest_pos[i]; // the last pos
            tmp.timestamp = IntFp(HIGH64(write_list[pos].timestamp.value), ALICE);            
            tmp.val = IntFp(HIGH64(write_list[pos].val.value), ALICE);
            read_list.push_back(tmp);
            // boundcheck access
            // bound_check_rom->Access(IntFp(T+1, PUBLIC) - tmp.timestamp);
            // Dicuss this with Dave
        }

        // bound check
        bound_check_rom->Teardown_Batch(party, block_size);

        uint64_t A0, A1, A2; // to compress the tuple into one IntFp
        uint64_t X; // evaluate at this value
        if (party == ALICE) {
            ZKFpExec::zk_exec->recv_data(&A0, sizeof(uint64_t));
            ZKFpExec::zk_exec->recv_data(&A1, sizeof(uint64_t));
            ZKFpExec::zk_exec->recv_data(&A2, sizeof(uint64_t));
            ZKFpExec::zk_exec->recv_data(&X, sizeof(uint64_t));
        } else {
            PRG tmpprg;
		    tmpprg.random_data(&A0, sizeof(uint64_t));
		    A0 = A0 % PR;
		    tmpprg.random_data(&A1, sizeof(uint64_t));
		    A1 = A1 % PR;
		    tmpprg.random_data(&A2, sizeof(uint64_t));
		    A2 = A2 % PR;
            tmpprg.random_data(&X, sizeof(uint64_t));
            X = X % PR;
		    ZKFpExec::zk_exec->send_data(&A0, sizeof(uint64_t));  
            ZKFpExec::zk_exec->send_data(&A1, sizeof(uint64_t));  
            ZKFpExec::zk_exec->send_data(&A2, sizeof(uint64_t));  
            ZKFpExec::zk_exec->send_data(&X, sizeof(uint64_t));  
        }
	    
        IntFp prod_read = IntFp(1, PUBLIC);
        IntFp prod_write = IntFp(1, PUBLIC);
        // vector<IntFp> final_tree_r;
        // vector<IntFp> final_tree_w;

        // for the combine (poly's) term proof
        block seed; 
        if (party == ALICE) {
            ZKFpExec::zk_exec->recv_data(&seed, sizeof(block));
        } else {
            PRG().random_block(&seed, 1);
            ZKFpExec::zk_exec->send_data(&seed, sizeof(block));
        }
        PRG prg(&seed);      
        uint64_t acc_C[block_size];      
        memset(acc_C, 0, sizeof(acc_C));
        uint64_t acc_K = 0;

        uint64_t delta = ZKFpExec::zk_exec->get_delta();
        int now_i = 0;
        uint64_t power_delta = 1;
        for (int i = 1; i < block_size; i++) power_delta = mult_mod(power_delta, delta);
        
        while (now_i < N+T) {
            uint64_t product_r = 1;
            uint64_t product_w = 1;

            uint64_t C_r[block_size+1];
            memset(C_r, 0, sizeof(C_r)); C_r[0] = 1;
            uint64_t K_r = 1;
            uint64_t C_w[block_size+1];
            memset(C_w, 0, sizeof(C_w)); C_w[0] = 1;
            uint64_t K_w = 1;            
            
            uint64_t M, x;      
            uint64_t tmp[block_size+1]; 
            for (int i = 0; i < block_size; i++, now_i++) {
                IntFp tmp_r = read_list[now_i].idx * A0 + read_list[now_i].timestamp * A1 + read_list[now_i].val * A2 + X;
                IntFp tmp_w = write_list[now_i].idx * A0 + write_list[now_i].timestamp * A1 + write_list[now_i].val * A2 + X;
                if (party == ALICE) { // Alice computes the coeffs to prove the commited product is well-formed
                    product_r = mult_mod(product_r, HIGH64(tmp_r.value));
                    product_w = mult_mod(product_w, HIGH64(tmp_w.value)); 
                    // for read
                    M = LOW64(tmp_r.value); x = PR - HIGH64(tmp_r.value);
                    tmp[0] = 0;
                    for (int j = 0; j <= i; j++) tmp[j+1] = mult_mod(x, C_r[j]);
                    for (int j = 0; j <= i; j++) C_r[j] = add_mod(tmp[j], mult_mod(M, C_r[j])); 
                    C_r[i+1] = tmp[i+1];
                    // for write
                    M = LOW64(tmp_w.value); x = PR - HIGH64(tmp_w.value);
                    tmp[0] = 0;
                    for (int j = 0; j <= i; j++) tmp[j+1] = mult_mod(x, C_w[j]);
                    for (int j = 0; j <= i; j++) C_w[j] = add_mod(tmp[j], mult_mod(M, C_w[j])); 
                    C_w[i+1] = tmp[i+1];                    
                } else { // Bob computes the expected proofs
                    K_r = mult_mod(K_r, LOW64(tmp_r.value));
                    K_w = mult_mod(K_w, LOW64(tmp_w.value));
                }
            }
            IntFp combine_r_term = IntFp(product_r, ALICE);
            IntFp combine_w_term = IntFp(product_w, ALICE);
            // final_tree_r.push_back(combine_r_term);
            // final_tree_w.push_back(combine_w_term);
            prod_read = prod_read * combine_r_term;
            prod_write = prod_write * combine_w_term;
            // checking the combine term is correctly formed
            if (party == ALICE) {
                //std::cout << C_r[block_size] << ' ' << HIGH64(combine_r_term.value) << std::endl;
                //std::cout << C_w[block_size] << ' ' << HIGH64(combine_w_term.value) << std::endl;
                C_r[block_size-1] = add_mod(C_r[block_size-1], LOW64(combine_r_term.value));
                C_w[block_size-1] = add_mod(C_w[block_size-1], LOW64(combine_w_term.value));
                uint64_t random_c;
                prg.random_data(&random_c, sizeof(uint64_t));
		        random_c %= PR;
                for (int i = 0; i < block_size; i++) acc_C[i] = add_mod(acc_C[i], mult_mod(C_r[i], random_c));
                prg.random_data(&random_c, sizeof(uint64_t));
		        random_c %= PR;
                for (int i = 0; i < block_size; i++) acc_C[i] = add_mod(acc_C[i], mult_mod(C_w[i], random_c));
            } else {
                K_r = add_mod(K_r, mult_mod(power_delta, LOW64(combine_r_term.value)));
                K_w = add_mod(K_w, mult_mod(power_delta, LOW64(combine_w_term.value)));
                uint64_t random_c;
                prg.random_data(&random_c, sizeof(uint64_t));
		        random_c %= PR;
                acc_K = add_mod(acc_K, mult_mod(random_c, K_r));                
                prg.random_data(&random_c, sizeof(uint64_t));
		        random_c %= PR;
                acc_K = add_mod(acc_K, mult_mod(random_c, K_w));                
            }
        }

        // EpsilonScan(party, final_tree_r, final_tree_r.size(), block_size, prg, acc_C, acc_K, prod_read);
        // EpsilonScan(party, final_tree_w, final_tree_w.size(), block_size, prg, acc_C, acc_K, prod_write);

        // check the polynomial proof
        // add random mask for ZK
        if (party == ALICE) {
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
            for (int i = 0; i < block_size; i++) acc_C[i] = add_mod(acc_C[i], random_pad[i]);
        } else {
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
            acc_K = add_mod(acc_K, random_pad);
        }

        if (party == ALICE) ZKFpExec::zk_exec->send_data(&acc_C, sizeof(acc_C));
        else ZKFpExec::zk_exec->recv_data(&acc_C, sizeof(acc_C));
        // check the proof
        power_delta = 1;
        if (party == BOB) {
            uint64_t proof = 0;
            for (int i = 0; i < block_size; i++) {
                proof = add_mod(proof, mult_mod(power_delta, acc_C[i]));
                power_delta = mult_mod(power_delta, delta);
            }
            if (proof != acc_K) error("Prover cheat!");
        }

        // check they are equal, namely, the same is zero
        IntFp final_zero = prod_read + prod_write.negate();
	    batch_reveal_check_zero(&final_zero, 1);
        
    }

};

#endif