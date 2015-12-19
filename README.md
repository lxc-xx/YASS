# [YASS](http://lxc-xx.github.io/YASS/)
# Yet Another SDP Solver
YASS is a Semidefinite programming (SDP) solver using interior-point method. Currently it has two implementions using Python and C++ respectively. It is aimed to creat a lightweight and fast solver, which can be called from any device and environment.

# What is the semidefinite programming.
Semidefinite programming is a convex optimization problem. Check this [doc](http://sdpa.sourceforge.net/whatissdp.pdf) for introduction.

## Python implementation
Python implementation is based on [Numpy](http://www.numpy.org), a popular Python numerical computing package. 

## C++ implementation
C++ implementation is based on [Eigen] (http://eigen.tuxfamily.org/index.php?title=Main_Page), a C++ template for linear algebra. Eigen is included in the package. So YASS can be compiled does not depend on any non-standard C++ library. It is easy to compile it at any platform where a C++ compiler is available.
# Usage
Here is the example for solving truss1, a problem from SDPLIB
## Python
```bash
$ ./yass.py
Usage: yass.py problem_file init_file output
$ ./yass.py ../example/truss1.dat-s ../example/truss1.dat-ini solution
```
## C++
Compile
```bash
$cd cpp
$make
```
Run
```bash
./yass
Usage: sdp problem_file init_file output_file
$ ./yass ../example/truss1.dat-s ../example/truss1.dat-ini solution
```

# Speed
YASS is keep evolving. At this stage, we are focused on the function development. Here is the current speed test on the sample problems from [SDPLIB](http://euler.nmt.edu/~brian/sdplib/sdplib.html). We compare it with another heavy-weight SDP solver [SDPA](http://sdpa.sourceforge.net).

We tested the YASS C++, YASS Python, SDPA on problems of SDPLIB The experiment is running on a AWS EC2 [g2.8xlarge] (https://aws.amazon.com/ec2/instance-types/) instance.

| Problem  | YASS C++ | YASS Python | SDPA   |
|----------|----------|-------------|--------|
| truss1   | 0.524s   | 0.722s      | 0.117s |
| truss4   | 3.696s   | 3.162s      | 0.066s |
| truss5   | 53.484s  | 25.418s     | 1.802s |
| control1 | N/A      | 1.532s      | 0.482s |
| control2 | N/A      | 24.446s     | 0.665s |
| control3 | N/A      | 4m46.772s   | 1.259s |







# Reference

[1]Freund, Robert M. "Introduction to semidefinite programming (SDP)." Massachusetts Institute of Technology (2004).
