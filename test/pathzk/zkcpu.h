#ifndef ZKCPU_H__
#define ZKCPU_H__

#include <vector>
#include <random>
#include <iostream>

#include "emp-zk/emp-zk.h"
#include "emp-tool/emp-tool.h"
#if defined(__linux__)
	#include <sys/time.h>
	#include <sys/resource.h>
#elif defined(__APPLE__)
	#include <unistd.h>
	#include <sys/resource.h>
	#include <mach/mach.h>
#endif

class f5 {
public:
    int val;

    f5() {val = 0;}
    f5(int val) : val(((val%5)+5)%5) {}

    static f5 minor(int x) {
        x = -x;
        x = ((x % 5) + 5) % 5;
        return x;
    }

    f5 & operator+= (f5 const & rhs) { 
        val += rhs.val; 
        val %= 5;
        return *this; 
    }
    f5 & operator+= (f5 && rhs) { 
        val += rhs.val;
        val %= 5;
        return *this;
    }
    f5 & operator*= (f5 const & rhs) {
        val *= rhs.val;
        val %= 5;
        return *this;
    }
    f5 & operator*= (f5 && rhs) {
        val *= rhs.val;
        val %= 5;
        return *this;
    }

    f5 operator+(f5 rhs) &  { return rhs += *this; }
    f5 operator+(f5 rhs) && { return rhs += std::move(*this); }
    f5 operator*(f5 rhs) &  { return rhs *= *this; }
    f5 operator*(f5 rhs) && { return rhs *= std::move(*this); }

    static f5 unit() {
        return f5(1);
    }

    static f5 zero() {
        return f5(0);
    }
};

class f61 {
public:
    uint64_t val;

    f61(uint64_t val) : val(val) {}
    f61() : val(0) {}

    static f61 unit() {
        return f61(1);
    }
    static f61 zero() {
        return f61(0);
    }
    static f61 minor(uint64_t x) {
        return f61(PR - x);
    }

    f61 & operator+= (f61 const & rhs) { 
        val = add_mod(val, rhs.val);         
        return *this; 
    }
    f61 & operator+= (f61 && rhs) { 
        val = add_mod(val, rhs.val);
        return *this;
    }
    f61 & operator*= (f61 const & rhs) {
        val = mult_mod(val, rhs.val);
        return *this;
    }
    f61 & operator*= (f61 && rhs) {
        val = mult_mod(val, rhs.val);
        return *this;
    }
    f61 operator+(f61 rhs) &  { return rhs += *this; }
    f61 operator+(f61 rhs) && { return rhs += std::move(*this); }
    f61 operator*(f61 rhs) &  { return rhs *= *this; }
    f61 operator*(f61 rhs) && { return rhs *= std::move(*this); }    
};

enum OPTYPE {
    INPUT, ADD, MUL
};

class BaseOp {
public:
    OPTYPE op;
    std::size_t l,r;
    BaseOp(OPTYPE op, std::size_t l, std::size_t r) : op(op), l(l), r(r) {}
};

class Instruction {
public:
    std::size_t m; // number of register
    std::size_t n; // number of ``useful'' multiplication
    std::vector<BaseOp> bank;
    std::vector<std::size_t> lid;
    std::vector<std::size_t> rid;
    std::vector<std::size_t> oid;    

    Instruction(std::size_t m, std::size_t n) : m(m), n(n) {}

    void rand_inst(std::mt19937 &gen) {
        std::uniform_int_distribution<> distr(0, n*n);
        for (int i = 0; i < 2*m+n; i++) bank.push_back(BaseOp(OPTYPE::INPUT, 0, 0));
        int cur = n;
        while (cur) {
            std::size_t last = bank.size();
            std::size_t l = distr(gen) % last;
            std::size_t r = distr(gen) % last;
            if (l > r) std::swap(l, r);
            OPTYPE optype = distr(gen) % 2 ? OPTYPE::ADD : OPTYPE::MUL;
            bank.push_back(BaseOp(optype, l, r));
            if (optype == OPTYPE::MUL) {                
                cur--;
            }
        }
    }

    void rand_inst_pselect(std::mt19937 &gen) {
        std::uniform_int_distribution<> distr(0, n*n);
        for (int i = 0; i < 2*m+n+2; i++) bank.push_back(BaseOp(OPTYPE::INPUT, 0, 0));
        int cur = n;
        while (cur) {
            std::size_t last = bank.size();
            std::size_t l = distr(gen) % last;
            std::size_t r = distr(gen) % last;
            if (l > r) std::swap(l, r);
            OPTYPE optype = distr(gen) % 2 ? OPTYPE::ADD : OPTYPE::MUL;
            bank.push_back(BaseOp(optype, l, r));
            if (optype == OPTYPE::MUL) {
                lid.push_back(l);
                rid.push_back(r);
                oid.push_back(last);                
                cur--;
            }
        }
    }    

    void print() {
        for (int i = 0; i < bank.size(); i++) {
            if (bank[i].op == OPTYPE::INPUT) {
                std::cout << i << ": INPUT " << std::endl;
            } else if (bank[i].op == OPTYPE::ADD) {
                std::cout << i << ": ADD " << bank[i].l << ' ' << bank[i].r << std::endl;
            } else {
                std::cout << i << ": MUL " << bank[i].l << ' ' << bank[i].r << std::endl;
            }
        }
    }

    template <typename T> 
    void Eval(std::vector<T> &wire, std::vector<T> &output) {
        wire.resize(bank.size());
        for (int i = 0; i < bank.size(); i++) {
            if (bank[i].op == OPTYPE::INPUT) continue;
            if (bank[i].op == OPTYPE::ADD) {
                wire[i] = wire[bank[i].l] + wire[bank[i].r];
            } else { // MUL
                wire[i] = wire[bank[i].l] * wire[bank[i].r];
            }
        }
        output.clear();
        for (int i = wire.size() - m; i < wire.size(); i++) output.push_back(wire[i]);
    }

    template <typename T> 
    void EvalWithExtWit(std::vector<T> &wire, std::vector<T> &output, std::vector<T> &lwire, std::vector<T> &rwire) {
        wire.resize(bank.size());
        lwire.clear(); rwire.clear();
        for (int i = 0; i < bank.size(); i++) {
            if (bank[i].op == OPTYPE::INPUT) continue;
            if (bank[i].op == OPTYPE::ADD) {
                wire[i] = wire[bank[i].l] + wire[bank[i].r];
            } else { // MUL
                wire[i] = wire[bank[i].l] * wire[bank[i].r];
                lwire.push_back(wire[bank[i].l]);
                rwire.push_back(wire[bank[i].r]);
            }
        }
        output.clear();
        for (int i = wire.size() - m + 1; i < wire.size(); i++) output.push_back(wire[i]);
    }    

    template <typename T>
    void comp_topo_vec_pselect(T chi, int bid, std::vector<T> &topo_vec) {
        std::vector<T> wire(bank.size());
        T coeff = T::unit();
        T constant = T::zero();

        // set up initial values on l and r
        for (int i = 0; i < n; i++) {
            // setup lid
            wire[lid[i]] += coeff;
            coeff *= chi;            
            // setup rid
            wire[rid[i]] += coeff;
            coeff *= chi;
        }
        // set up initial values on C's output wire: left (1), right (reg[m-1]-bid)
        constant += coeff; coeff *= chi;
        constant += T::minor(bid) * coeff; wire[m-1] += coeff; coeff *= chi;
        // set up inital values on ``special'' 1*1=1 tuple
        constant += coeff; coeff *= chi;
        constant += coeff; coeff *= chi;
        // set up initial values on l/r of register (in the next instruction), for the first m-1 registers
        for (int i = bank.size() - m + 1; i < bank.size(); i++) {
            constant += coeff; coeff *= chi;
            wire[i] += coeff; coeff *= chi;
        }
        // set up initial values on l/r of register (in the next instruction), for the final register
        constant += coeff; coeff *= chi;
        wire[n+2*m+1] += coeff; coeff *= chi;

        // evaluate the circuit backward
        for (int i = bank.size() - 1; i >= 0; i--) 
            if (bank[i].op == OPTYPE::ADD) {
                wire[bank[i].l] += wire[i];
                wire[bank[i].r] += wire[i];
            }

        std::vector<T> tmp_o; tmp_o.clear();
        tmp_o.push_back(constant);
        for (int i = 0; i < m; i++) tmp_o.push_back(wire[i]);
        for (int i = 0; i < n; i++) tmp_o.push_back(wire[oid[i]]);
        tmp_o.push_back(T::unit());

        // merge into topo_vec
        topo_vec.clear();
        for (int i = 0; i < n+m+2; i++) {
            topo_vec.push_back(wire[m+i]);
            topo_vec.push_back(tmp_o[i]);
        }
    }    

};

class ZKCPU {
public:
    std::size_t B; // number of branches
    std::size_t m; // number of regsiter
    std::mt19937 gen;
    std::vector<Instruction> br;
    
    ZKCPU(std::size_t B, std::size_t m, std::__1::random_device::result_type seed) : B(B), m(m) {
        gen = std::mt19937(seed);
    }

    void rand_cpu() {
        for (int i = 0; i < B; i++) {
            br.push_back(Instruction(m, i*100));
            auto last = br.size();
            //br[last - 1].rand_inst(gen);
            br[last - 1].rand_inst_pselect(gen);
        }
    }

    void print() {
        for (int i = 0; i < B; i++) {
            br[i].print();
            std::cout << "=========" << std::endl;
        }
    }

    void test_eval_f5() {
        std::vector<f5> reg[2];
        std::uniform_int_distribution<> distr(0, 4);

        std::size_t now = 0; // rotate        
        for (int i = 0; i < m; i++) reg[now].push_back(f5());        
        
        for (int i = 0; i < 3; i++, now = 1 - now) {
            std::size_t cid = reg[now][m-1].val;
            std::size_t n = br[cid].n;
            for (int j = 0; j < n+m; j++) {
                f5 x;
                x.val = distr(gen);
                reg[now].push_back(x);
            }
            // spec output
            std::cout << i << ":" << std::endl;
            for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;
            // execute
            br[cid].Eval<f5>(reg[now], reg[1-now]);
            // use mod to fetch next instrution
            reg[1-now][m-1].val %= B;            
        }

        // final state
        std::cout << "Final:" << std::endl;
        for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;
    }

    void test_proof_f5(std::size_t step) {
        std::vector<f5> reg[2];
        std::uniform_int_distribution<> distr(0, 4);
        std::uniform_int_distribution<> distrb(0, B-1);

        std::size_t now = 0; // rotate        
        for (int i = 0; i < m; i++) reg[now].push_back(f5());        

        std::vector<f5> in;
        std::vector<f5> l, r, o;
        std::vector<std::size_t> cids;
        
        for (int i = 0; i < step; i++, now = 1 - now) {
            std::size_t cid = reg[now][m-1].val;
            cids.push_back(cid);
            std::size_t n = br[cid].n;
            for (int j = 0; j < n+m+1; j++) {
                f5 x;
                x.val = distr(gen);
                reg[now].push_back(x);
                in.push_back(x);
            }
            f5 next_cid;
            next_cid.val = distrb(gen);
            reg[now].push_back(next_cid);
            in.push_back(next_cid);

            // spec output
            // std::cout << i << ":" << std::endl;
            // for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;
            
            // execute
            std::vector<f5> lwire, rwire;
            br[cid].EvalWithExtWit<f5>(reg[now], reg[1-now], lwire, rwire);

            // setup l,r,o
            // the fisrt 1*1=1
            l.push_back(f5(1)); r.push_back(f5(1)); o.push_back(f5(1));
            // to capture m reg
            for (int j = 0; j < m; j++) {
                l.push_back(f5(1));
                r.push_back(reg[now][j]);
                o.push_back(reg[now][j]);
            }
            // multi tuples in the middle
            for (int j = 0; j < lwire.size(); j++) {
                l.push_back(lwire[j]);
                r.push_back(rwire[j]);
                o.push_back(lwire[j] * rwire[j]);
            }
            // the last 1*0=0
            l.push_back(f5(1)); r.push_back(f5(0)); o.push_back(f5(0));

            // p selects the next instruction
            reg[1-now].push_back(next_cid);
        }

        //std::cout << in.size() << ' ' << l.size() << ' ' << r.size() << ' ' << o.size() << std::endl;

        // final state
        // std::cout << "Final:" << std::endl;
        // for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;        

        f5 chi(distr(gen));
        //if (chi.val == 0) std::cout << "HHHH" << std::endl;
        std::vector<f5> topo_vec[B];
        for (int i = 0; i < B; i++) br[i].comp_topo_vec_pselect<f5>(chi, i, topo_vec[i]);

        f5 sum_l = f5::unit() + chi;        
        f5 coeff = chi * chi;
        for (int i = 0; i < m; i++) {
            sum_l += coeff; coeff *= chi;
            coeff *= chi;
        }
        std::vector<f5> ll; ll.clear();
        for (int i = 0; i < step; i++) {
            std::size_t cid = cids[i];
            f5 tmp_coeff = coeff;
            for (int j = 0; j < topo_vec[cid].size(); j++) {
                ll.push_back(topo_vec[cid][j]*coeff);
                tmp_coeff *= chi;
            }
            coeff = tmp_coeff;
        }
        for (int i = 0; i < in.size(); i++) {
            sum_l += ll[i*2] * in[i];
            sum_l += ll[i*2+1] * o[i];
        }

        f5 sum_r = f5::zero();
        coeff = f5::unit();
        for (int i = 0; i < l.size(); i++) {
            sum_r += coeff * l[i]; coeff *= chi;
            sum_r += coeff * r[i]; coeff *= chi;
        }
        sum_r += coeff; coeff *= chi;
        sum_r += coeff; coeff *= chi;
        for (int i = 0; i < reg[now].size(); i++) {
            sum_r += coeff; coeff *= chi;
            sum_r += coeff * reg[now][i]; coeff *= chi;
        }
        
        std::cout << sum_l.val-sum_r.val << std::endl;
    }

    uint64_t gen_rand_f61(PRG &prg_s) {
        block tmp;
        prg_s.random_block(&tmp, 1);
        return LOW64(tmp);
    }

    void test_proof_f61(std::size_t step) {
        std::vector<f61> reg[2];
        block s_seed; 
        PRG().random_block(&s_seed, 1);
        PRG prg_s(&s_seed); 

        std::size_t now = 0; // rotate        
        for (int i = 0; i < m; i++) reg[now].push_back(f61());        

        std::vector<f61> in;
        std::vector<f61> l, r, o;
        std::vector<std::size_t> cids;
        
        for (int i = 0; i < step; i++, now = 1 - now) {
            std::size_t cid = reg[now][m-1].val;
            cids.push_back(cid);
            std::size_t n = br[cid].n;
            for (int j = 0; j < n+m+1; j++) {
                f61 x(gen_rand_f61(prg_s) % PR);
                reg[now].push_back(x);
                in.push_back(x);
            }
            f61 next_cid(gen_rand_f61(prg_s) % B);
            reg[now].push_back(next_cid);
            in.push_back(next_cid);

            // spec output
            // std::cout << i << ":" << std::endl;
            // for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;
            
            // execute
            std::vector<f61> lwire, rwire;
            br[cid].EvalWithExtWit<f61>(reg[now], reg[1-now], lwire, rwire);

            // setup l,r,o
            // the fisrt 1*1=1
            l.push_back(f61(1)); r.push_back(f61(1)); o.push_back(f61(1));
            // to capture m reg
            for (int j = 0; j < m; j++) {
                l.push_back(f61(1));
                r.push_back(reg[now][j]);
                o.push_back(reg[now][j]);
            }
            // multi tuples in the middle
            for (int j = 0; j < lwire.size(); j++) {
                l.push_back(lwire[j]);
                r.push_back(rwire[j]);
                o.push_back(lwire[j] * rwire[j]);
            }
            // the last 1*0=0
            l.push_back(f61(1)); r.push_back(f61(0)); o.push_back(f61(0));

            // p selects the next instruction
            reg[1-now].push_back(next_cid);
        }

        //std::cout << in.size() << ' ' << l.size() << ' ' << r.size() << ' ' << o.size() << std::endl;

        // final state
        // std::cout << "Final:" << std::endl;
        // for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;        

        f61 chi(gen_rand_f61(prg_s) % PR);
        std::vector<f61> topo_vec[B];
        for (int i = 0; i < B; i++) br[i].comp_topo_vec_pselect<f61>(chi, i, topo_vec[i]);

        f61 sum_l = f61::unit() + chi;        
        f61 coeff = chi * chi;
        for (int i = 0; i < m; i++) {
            sum_l += coeff; coeff *= chi;
            coeff *= chi;
        }
        std::vector<f61> ll; ll.clear();
        for (int i = 0; i < step; i++) {
            std::size_t cid = cids[i];
            f61 tmp_coeff = coeff;
            for (int j = 0; j < topo_vec[cid].size(); j++) {
                ll.push_back(topo_vec[cid][j]*coeff);
                tmp_coeff *= chi;
            }
            coeff = tmp_coeff;
        }
        for (int i = 0; i < in.size(); i++) {
            sum_l += ll[i*2] * in[i];
            sum_l += ll[i*2+1] * o[i];
        }

        f61 sum_r = f61::zero();
        coeff = f61::unit();
        for (int i = 0; i < l.size(); i++) {
            sum_r += coeff * l[i]; coeff *= chi;
            sum_r += coeff * r[i]; coeff *= chi;
        }
        sum_r += coeff; coeff *= chi;
        sum_r += coeff; coeff *= chi;
        for (int i = 0; i < reg[now].size(); i++) {
            sum_r += coeff; coeff *= chi;
            sum_r += coeff * reg[now][i]; coeff *= chi;
        }
        
        std::cout << sum_l.val-sum_r.val << std::endl;
    }


    void rand_f61_witness_by_step(std::size_t step, std::vector<std::size_t> &cids, std::vector<f61> &in, std::vector<f61> &l, std::vector<f61> &r, std::vector<f61> &final_reg) {
        std::vector<f61> reg[2];
        block s_seed; 
        PRG().random_block(&s_seed, 1);
        PRG prg_s(&s_seed); 

        std::size_t now = 0; // rotate        
        for (int i = 0; i < m; i++) reg[now].push_back(f61());        

        cids.clear();
        in.clear();
        l.clear(); r.clear(); // o.clear();
        final_reg.clear();
        
        for (int i = 0; i < step; i++, now = 1 - now) {
            std::size_t cid = reg[now][m-1].val;
            cids.push_back(cid);
            std::size_t n = br[cid].n;
            for (int j = 0; j < n+m+1; j++) {
                f61 x(gen_rand_f61(prg_s) % PR);
                reg[now].push_back(x);
                in.push_back(x);
            }
            f61 next_cid(gen_rand_f61(prg_s) % B);
            reg[now].push_back(next_cid);
            in.push_back(next_cid);

            // execute
            std::vector<f61> lwire, rwire;
            br[cid].EvalWithExtWit<f61>(reg[now], reg[1-now], lwire, rwire);

            // setup l,r,o
            // the fisrt 1*1=1
            l.push_back(f61(1)); r.push_back(f61(1)); // o.push_back(f61(1));
            // to capture m reg
            for (int j = 0; j < m; j++) {
                l.push_back(f61(1));
                r.push_back(reg[now][j]);
                // o.push_back(reg[now][j]);
            }
            // multi tuples in the middle
            for (int j = 0; j < lwire.size(); j++) {
                l.push_back(lwire[j]);
                r.push_back(rwire[j]);
                // o.push_back(lwire[j] * rwire[j]);
            }
            // the last 1*0=0
            l.push_back(f61(1)); r.push_back(f61(0)); // o.push_back(f61(0));

            // p selects the next instruction
            reg[1-now].push_back(next_cid);
        }     
        for (int i = 0; i < reg[now].size(); i++) final_reg.push_back(reg[now][i]);
    }

};

#endif