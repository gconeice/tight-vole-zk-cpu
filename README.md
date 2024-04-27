# TIGHT ZK CPU

This repository implements the tight ZK CPU protocol.
See our paper for details.

Basis on EMP
=====
Our tight ZK CPU is based on (1) QuickSilver's (https://eprint.iacr.org/2021/076) repository, which is part of the EMP Toolkit: https://github.com/emp-toolkit/emp-zk; (2) Batchman's (https://github.com/gconeice/stacking-vole-zk) repository; and (3) the state-of-the-art ZK ROM/RAM (https://github.com/gconeice/improved-zk-ram) repository. In particular, we forked their repositories and developed based on them. We also tweak some of EMP's libraries.
In our final open source version, we will further clarify all changes.

MIT license is included as part of each repository.

Setup environment
=====
`sudo bash setup.sh`

Build and install
=====
`bash install.sh`

Browsing the code
=====
`/test/pathzk` contains our core code.

`/emp-zk` contains the EMP Toolkit's ZK library.

`/test/arith` contains the Batchman (baseline) implementations.

`/zk-ram` contains the state-of-the-art ZK (balanced) ROM/RAM implementations.

Expected executable
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

Benmark summary
=====
`benchmark_summary.xlsx` recorded our experiments, which were used to plot our tables/figures in the paper.

Test
=====
We will further explain how one can reproduce all results in the paper when we de-anonymize the repository.