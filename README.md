Implementation of the [BGW88] Multi-Party Compuation protocol for the semi-honest case.
We also implemented a much efficient way for evaluating multiplication gates by [DIK10].
Secret sharing reconstruction is implemented using the recent exact-repair problem algorithm suggested by [GW17].

The project implements many pure python modules that may be used as-is:
GF - Finite fields for any general prime power.
Polynomials - a generic polynomial and numerous operations on it.
Matrix - a generic simple and powerful matrix module.
Circuit - implementation of an arithmetic circuit.
Transformation - a matrix-based linear transformation.
There are zero dependencies for this project, developed for python 2.7.8

Run main.py for an example run. Different parameters may be confiugred in config.py or directly in main.py.

[BGW88] - Michael Ben-Or, Shafi Goldwasser and Avi Wigderson. Completeness theorems for non-cryptographic fault-tolerant distributed computation.

[DIK10] - Ivan Damgard, Yuval Ishai and Mikkel Kroigaard. Perfectly secure multiparty computation and the computational overhead of cryptography.

[GW17]  - Venkatesan Guruswami and Mary Wootters. Repairing Reed-Solomon Codes.
