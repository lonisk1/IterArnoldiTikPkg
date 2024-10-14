# IterArnoldiTikPkg
The package contains the basic function for the Iterated Arnoldi-Tikhonov (IAT) method. To reproduce some of the examples from the paper, the software package IR-Tools is needed. A test function in the package requires IR-Tools.

IterArnoldiTikPkg
=============

Software associated with the paper:

A. Buccini, L. Onisk, and L. Reichel, An Arnoldi-based preconditioner for iterated Tikhonov regularization. Numerical Algorithms 92, 223--245 (2023)

---

DESCRIPTION
-----------

This package provides (1) a function cpable of solving linear discrete ill-posed problems terminated according to the discrepancy principle and
(2) a demo based upon the software package IR Tools (by P. C. Hansen, S. Gazzola, and J. G. Nagy).

Release: 1.0, March 2023

Programming language: Matlab 9.9 (R2020b)

License: (see license.md)

---

INSTALLATION
------------

Download and extract the IterArnoldiTikPkg package. Additional files to run the demo may be needed from the software package IR Tools by
P. C. Hansen, S. Gazzola, and J. G. Nagy.

---

PACKAGE USE
------------

As long as the IR Tools package is installed and this pacakge is installed and you are in the correct directory, there should be no issue running the associated demo.

---

PACKAGE STRUCTURE
-----------------

The following is a list with a brief description of the contents of the
IterArnoldiTikPkg package.

* README.md                 : This file.
* ArnoldiIterTik.m          : Modular function for Arnoldi iterative Tikhonov based on the  Arnoldi-Tikhonov method published in 23' in Numerical Algorithms. Method terminates according to the discrepancy principle.
* TestScript_IATmethod.m    : This test script requires that you have the IR-Tools toolbox downloaded and pointed to. You must also have the 'ArnoldiIterTik.m' code mapped to.
