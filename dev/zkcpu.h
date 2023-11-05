#ifndef ZKCPU_H__
#define ZKCPU_H__

#include <vector>
#include <random>
#include <iostream>

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
            br[last - 1].rand_inst(gen);
        }
    }

    void print() {
        for (int i = 0; i < B; i++) {
            br[i].print();
            std::cout << "=========" << std::endl;
        }
    }
};

#endif