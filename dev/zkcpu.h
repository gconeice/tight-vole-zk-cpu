#ifndef ZKCPU_H__
#define ZKCPU_H__

#include <vector>
#include <random>
#include <iostream>

class f5 {
public:
    std::size_t val;

    f5() {val = 0;}
    f5(std::size_t val) : val(val) {}

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
            if (optype == OPTYPE::MUL) cur--;
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
            if (optype == OPTYPE::MUL) cur--;
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
            br.push_back(Instruction(m, 10));
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
        std::uniform_int_distribution<> distrb(0, m-1);

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
            std::cout << i << ":" << std::endl;
            for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;
            
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

        // std::cout << in.size() << ' ' << l.size() << ' ' << r.size() << ' ' << o.size() << std::endl;

        // final state
        std::cout << "Final:" << std::endl;
        for (int j = 0; j < reg[now].size(); j++) std::cout << reg[now][j].val << std::endl;        
    }
};

#endif