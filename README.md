# VOLE-Based Tight ZK CPU

This is the artifact for the paper: **Tight ZK CPU: Batched ZK Branching with Cost Proportional to Evaluated Instruction** to be presented on ACM CCS 2024.

ePrint link: https://eprint.iacr.org/2024/456

Base
=====
Our tight ZK CPU is based on (1) QuickSilver's (https://eprint.iacr.org/2021/076) repository, which is part of the EMP Toolkit: https://github.com/emp-toolkit/emp-zk; (2) Batchman's (https://github.com/gconeice/stacking-vole-zk) repository; and (3) the state-of-the-art ZK ROM/RAM (https://github.com/gconeice/improved-zk-ram) repository. In particular, we forked the Batchman repo and developed based on it. We also tweak some of EMP's libraries.

MIT license is included as part of each repository.

Hardware and OS
=====
Our experiments were performed over **two** AWS EC2 `m5.2xlarge` instances.
One particular experiment was performed over **two** AWS EC `m5.8xlarge` instances because of our baseline needs more RAM.
We used two machines to emulate ZK prover P and ZK verifier V.
(It can also be executed over one single machine to emulate two parties via localhost.)
We tested our code on a clean installation of `Ubuntu 22.04`, while we believe it also works on OS X and other versions of Ubuntu.

The detailed hardware configuration can be found on: https://aws.amazon.com/ec2/instance-types/m5/.


Setup Environment
=====
`sudo bash setup.sh`

Build and Install
=====
`bash install.sh`

Browsing the Code
=====
`/test/pathzk` contains our core code.

`/emp-zk` contains the EMP Toolkit's ZK library.

`/test/arith` contains the Batchman (baseline) implementations.

`/zk-ram` contains the state-of-the-art ZK (balanced) ROM/RAM implementations.

Expected Executable
=====
After compiling, the executable would show up in `build/bin`, including the following executables:

` test_pathzk_comp_batchman_balance`: This is used to test our tight ZK CPU with balanced instructions and without the rounding optimization, compared to Batchman.

`test_pathzk_comp_batchman_balance_opt`: This is used to test our tight ZK CPU with balanced instructions and the rounding optimization, compared to Batchman.

`test_pathzk_comp_batchman_unbalance`: This is used to test our tight ZK CPU with unbalanced instructions and without the rounding optimization, compared to Batchman.

`test_pathzk_comp_with_pub_cir`: This is used to test our tight ZK CPU with uniformly distributed sizes of instructions and without the rounding optimization, compared to the *insecure execution*. The execution would generate files that can be re-executed by the *insecure execution* with a same CPU configuration.

`test_pathzk_comp_with_pub_cir_opt`: This is used to test our tight ZK CPU with uniformly distributed sizes of instructions and the rounding optimization, compared to the *insecure execution*. The execution would generate files that can be re-executed by the *insecure execution* with a same CPU configuration.

`test_pathzk_fine_grain`: This is used to perform fine-grain analysis for our tight ZK CPU to generate the microbenchmarks.

`test_pathzk_pub_cir`: This is the *insecure execution* baseline, where the entire execution path is revealed to the verifier.

`test_pathzk_test`: This is the clean tight ZK CPU implementation, without any modification for comparison.

Toy Test Example
=====
You can use the following toy example to ensure that the compilation is successful.

Please `cd` to the `build` folder.

On P's machine, execute: `./bin/`

On v's machine, execute: `./bin/`

Here $IP is the P's IP address. If everything goes through, you should see the execution results on P and V. Starting from here, you can reproduce our results. (Details are listed below.)

Benmark Summary
=====
`benchmark_summary.xlsx` recorded our experiments, which were used to plot our tables/figures in the paper.

Reproduce the Result
=====
Please refer to `reproduce.pdf` to see how to get the numbers in `benchmark_summary.xlsx`.

# How to Simulate Network Setting

**We use `tc` command to simulate the network setting.**

     DEV=lo
     
     sudo tc qdisc del dev $DEV root
     
     sudo tc qdisc add dev $DEV root handle 1: tbf rate 1Gbit burst 100000 limit 10000
     
     sudo tc qdisc add dev $DEV parent 1:1 handle 10: netem delay 2msec

Change DEV to the network card you need (e.g., ens5).
Note that `sudo tc qdisc del dev $DEV root` needs to be executed before resetting the network.
It is used to clean the `tc` setting.

Both P and V need to restrict the network, note that for 30ms latency, both parties should be set to `delay 15msec`.

**You can use `iperf` to test the network throughput. Namely:**

P: `iperf -s`

V: `iperf -c [ip addr]`

**You can use `ping` to test the network latency. Namely:**

V: `ping [ip addr]`
