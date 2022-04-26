# Quadratic-Sieve
C Implementation of the Quadratic Sieve Algorithm

## How to Compile

```
gcc -O3 main.c quad_sieve.c -o quad_sieve.out -lgmp -lm && ./quad_sieve.out
```

## Requirements
- [GNU Multiple Precision Arithmetic Library](https://gmplib.org)
- A local [Magma Computer Algebra](http://magma.maths.usyd.edu.au/magma/) installation (to be able to factorize larger numbers) 
- [Python3](https://www.python.org/downloads/) for _encoder.py_. If you have a local Magma installed, there is no need for _encoder.py_
