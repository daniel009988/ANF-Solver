# DPLL-based SAT-Solver Specifically Designed for solving ANF Instances

This repository is an (extended and enhanced) implementation of the solver presented in the paper "Logical Cryptanalysis with WDSat" by Monika Trimoska, Gilles Dequen & Sorina Ionica, International Conference on Theory and Applications of Satisfiability Testing, SAT 2021: Theory and Applications of Satisfiability Testing – SAT 2021 pp 545–561. Over the last decade, there have been significant efforts in developing efficient XOR-enabled SAT solvers for cryptographic applications. Previously, the authors proposed a solver specialised to cryptographic problems, and more precisely to instances arising from the index calculus attack on the discrete logarithm problem for elliptic curve-based cryptosystems. Its most prominent feature is the module that performs an enhanced version of Gaussian Elimination is concentrated on the theoretical aspects of the new tool, but the running time-per-conflict results suggest that this module uses efficient implementation techniques as well. Thus, the first goal of this paper is to give a comprehensive exposition of the implementation details of the solver. In addition, they show that the approach can be extended to other cryptographic applications, mainly all attacks that involve solving dense Boolean polynomial systems. They give complexity analysis for such systems and we compare different state-of-the-art SAT solvers experimentally, concluding that their solver gives the best results. As a second contribution, they provide an original and economical implementation of a module for handling OR-clauses of any size, as their solver currently handles OR-clauses comprised of up to four literals. They finally provide experimental results showing that this new approach does not impair the performance of the solver.

## Requirements

We have tested the implementation on Ubuntu Linux and MacOS Monterey and expect to work on any modern Linux or Mac based platform. 

## Compilation

To compile the solver, type

```make``` 


After compilation, the file "q_solver" is being produced in the main directory.

## Running the solver

You can run the solver with default parameters:

```./q_solver -i <INSTANCE.ANF>```

This will load the ANF file with the aim to solve it. When a solution is found, the variable assignments are being printed.

## Parameters

The following command shows all available command line options:

```./q_solver -h```

These are the main options:

```-x : to enable Gaussian Elimination (default: enabled)```

```-i file    : input file```

```-v file    : input varmap (ANF_VARNAME->ANF_VARNUM mapping)```

```-w file    : input cnfvarmap (CNF_VARNAME->ANF_VARNUM mapping)```

```-q file    : input cube file (march_cu compatible)```

```-h file    : input hess dump file (hessanf compatible)```

```-r : reverse variable order to by tried```

```-g mvc    : where mvc is a string of comma-separated variables that defines statically the branching order```

```-c cubes  : where cubes is a string of comma-sperated variables that set the cube```


## Source code settings

The file ```config.h``` sets the upper boundaries for the variables. These should be adjusted to the requirements of the ANF formulations which are to be solved. 

## ANF File Format

This input describes XOR-clauses that can contain both unary literals and conjunctions of two or more literals. When a term is unary, only the corresponding literal is written. When a term is a conjunction, first we write .d, where d is the number of literals that comprise the conjunction and then, we list the literals. There can not be negative literals in the formula, but we can add a ⊤ constant instead.

Example:
```x1 ⊕ x2 ⊕ x4 ⊕ x5 ⊕ ⊤```

```x1 ⊕ (x1 ∧ x2) ⊕ x6```

```(x2 ∧ x3 ∧ x6) ⊕ x6 ⊕ ⊤```

is written as

```p cnf 6 3```

```x 1 2 4 5 T 0```

```x 1 .2 1 2 6 0```

```x .3 2 3 6 6 T 0```

<sub>@InProceedings{10.1007/978-3-030-80223-3_37,
author="Trimoska, Monika
and Dequen, Gilles
and Ionica, Sorina",
editor="Li, Chu-Min
and Many{\`a}, Felip",
title="Logical Cryptanalysis with WDSat",
booktitle="Theory and Applications of Satisfiability Testing -- SAT 2021",
year="2021",
publisher="Springer International Publishing",
address="Cham",
pages="545--561",
abstract="Over the last decade, there have been significant efforts in developing efficient XOR-enabled SAT solvers for cryptographic applications. In [22] we proposed a solver specialised to cryptographic problems, and more precisely to instances arising from the index calculus attack on the discrete logarithm problem for elliptic curve-based cryptosystems. Its most prominent feature is the module that performs an enhanced version of Gaussian Elimination. [22] is concentrated on the theoretical aspects of the new tool, but the running time-per-conflict results suggest that this module uses efficient implementation techniques as well. Thus, the first goal of this paper is to give a comprehensive exposition of the implementation details of WDSat. In addition, we show that the WDSat approach can be extended to other cryptographic applications, mainly all attacks that involve solving dense Boolean polynomial systems. We give complexity analysis for such systems and we compare different state-of-the-art SAT solvers experimentally, concluding that WDSat gives the best results. As a second contribution, we provide an original and economical implementation of a module for handling OR-clauses of any size, as WDSat currently handles OR-clauses comprised of up to four literals. We finally provide experimental results showing that this new approach does not impair the performance of the solver.",
isbn="978-3-030-80223-3"
}</sub>