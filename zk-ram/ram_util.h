#ifndef RAM_UTIL
#define RAM_UTIL

#include "emp-zk/emp-zk.h"
#include <iostream>
#include "emp-tool/emp-tool.h"
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

void EpsilonScan(int &party, vector<IntFp> &a, int n, int block_size, PRG &prg, uint64_t *acc_C, uint64_t &acc_K, IntFp &final_zero) {

    while (n >= block_size) {

        if (n%block_size) {
            cout << "NO=" << n%block_size << endl;
        }

        int m = 0;

        uint64_t delta = ZKFpExec::zk_exec->get_delta();
        int now_i = 0;
        uint64_t power_delta = 1;
        for (int i = 1; i < block_size; i++) power_delta = mult_mod(power_delta, delta);
        
        while (now_i < n) {
            uint64_t product_r = 1;

            uint64_t C_r[block_size+1];
            memset(C_r, 0, sizeof(C_r)); C_r[0] = 1;
            uint64_t K_r = 1;
            
            uint64_t M, x;      
            uint64_t tmp[block_size+1]; 
            for (int i = 0; i < block_size; i++, now_i++) {
                IntFp tmp_r = a[now_i];
                if (party == ALICE) { // Alice computes the coeffs to prove the commited product is well-formed
                    product_r = mult_mod(product_r, HIGH64(tmp_r.value));
                    // for read
                    M = LOW64(tmp_r.value); x = PR - HIGH64(tmp_r.value);
                    tmp[0] = 0;
                    for (int j = 0; j <= i; j++) tmp[j+1] = mult_mod(x, C_r[j]);
                    for (int j = 0; j <= i; j++) C_r[j] = add_mod(tmp[j], mult_mod(M, C_r[j])); 
                    C_r[i+1] = tmp[i+1];
                } else { // Bob computes the expected proofs
                    K_r = mult_mod(K_r, LOW64(tmp_r.value));
                }
            }
            IntFp combine_r_term = IntFp(product_r, ALICE);
            // checking the combine term is correctly formed
            if (party == ALICE) {
                //std::cout << C_r[block_size] << ' ' << HIGH64(combine_r_term.value) << std::endl;
                //std::cout << C_w[block_size] << ' ' << HIGH64(combine_w_term.value) << std::endl;
                C_r[block_size-1] = add_mod(C_r[block_size-1], LOW64(combine_r_term.value));
                uint64_t random_c;
                prg.random_data(&random_c, sizeof(uint64_t));
		        random_c %= PR;
                for (int i = 0; i < block_size; i++) acc_C[i] = add_mod(acc_C[i], mult_mod(C_r[i], random_c));
            } else {
                K_r = add_mod(K_r, mult_mod(power_delta, LOW64(combine_r_term.value)));
                uint64_t random_c;
                prg.random_data(&random_c, sizeof(uint64_t));
		        random_c %= PR;
                acc_K = add_mod(acc_K, mult_mod(random_c, K_r));                
            }
            a[m++] = combine_r_term;
        }


        n = m;
    }
    for (int i = 0; i < n; i++) final_zero = final_zero * a[i];
    //batch_reveal_check_zero(&final_zero, 1);
}

#endif