# Discrete-logarithm-solver

A C program that uses the GMP library to compute the Pohlig-Hellman reduction and solve a discrete logarithm instance. Firstly the modulus is factorized, 
then the Pohlig-Hellman reduction is applied and each remaining discrete log problems are solved.

Once they are all solved they are assembled to give the result of the initial problem.


# Installation

First you need to install the GMP library that allows us to handle large numbers in C:

```
sudo apt-get install libgmp3-dev
```

Then the installation is straightforward:

```
git clone https://github.com/timpaquatte/Discrete-logarithm-solver.git
cd Discrete-logarithm-solver
make
```

# Usage

Let's say you want to solve a problem of this form:  g^x = a [N]

Then the command to use is simply:
```
./bin/solve_dlog g a N
```
