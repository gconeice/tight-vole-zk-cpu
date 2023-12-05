# TIGHT ZK CPU

This repository implements the tight ZK CPU protocol.
See our paper for details.

Basis on EMP
=====
Our tight ZK CPU is based on (1) QuickSilver's (https://eprint.iacr.org/2021/076) repository, which is part of the EMP Toolkit: https://github.com/emp-toolkit/emp-zk; (2) Batchman's (https://github.com/gconeice/stacking-vole-zk) repository; and (3) the state-of-the-art ZK ROM/RAM (https://github.com/gconeice/improved-zk-ram) repository. In particular, we forked their repositories and developed based on them. We also tweak some of EMP's libraries.
In our final open source version, we will further clarify all changes.

MIT license is included as part of each repository.

Setup Environment
=====
`sudo bash setup.sh`

Build and Install
=====
`bash install.sh`

Browsing the code
=====
`/test/pathzk` contains our core code.
`/emp-zk` contains the EMP Toolkit's ZK library.
`/test/arith` contains the Batchman (baseline) implementations.
`/zk-ram` contains the state-of-the-art ZK (balanced) ROM/RAM implementations.

Test
=====
We will further detail how one can reproduce all results in the paper when we e-anonymize the repository.