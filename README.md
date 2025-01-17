# triadgraph
Python script for the analysis of energy transfer via triad interactions (triadgraph.py)  
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/smaeyama/triadgraph.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/smaeyama/triadgraph/context:python)

## Overview
Energy transfer via quadratic nonlinearity is described by triad interactions.  

  <img src="https://latex.codecogs.com/gif.latex?\frac{dE_k}{dt}=\sum_p\sum_qS_k^{p,q}" />  

The triad transfer function satisfies the symmetry <img src="https://latex.codecogs.com/gif.latex?S_k^{p,q}=S_k^{q,p}" />, and the detailed balance <img src="https://latex.codecogs.com/gif.latex?S_k^{p,q}+S_p^{q,k}+S_q^{k,p}=0" />.

**triadgraph** provide a tool for symmetrization, directional representation, and network graph visualization of triad interactions. For more details, see Ref. \[1\].

## Usage
**triadgraph** requires external packages: numpy, matplotlib, pygraphviz.

```
from triadgraph import symmetrize_triadtransfer, directional_triadtransfer, \
                       triadgraph_symmetric_all, triadgraph_directional_all
```
- ***symmetrize_triadtransfer*** returns the symmetrized triad transfer function <img src="https://latex.codecogs.com/gif.latex?S_k^{p,q}" /> from the non-symmetric triad transfer function <img src="https://latex.codecogs.com/gif.latex?A_k^{p,q}" /> .
```
S_kpq = symmetrize_triadtransfer(A_kpq, time_axis=0)
```
where S_kpq\[time,k,p,q\] and A_kpq\[time,k,p,q\] are 4D arrays.
- ***directional_triadtransfer*** returns the directional representation <img src="https://latex.codecogs.com/gif.latex?D_{k&space;\leftarrow&space;q}^p" title="D_{k \leftarrow q}^p" /> from the symmetrized triad transfer function.
```
D_kpq = directional_triadtransfer(S_kpq,time_axis=0)
```
where D_kpq\[time,k,p,q\] is a 4D array, represent energy transfer from q to k via coupling with a mediator p.
- ***triadgraph_symmetric_all*** plots a network graph of the symmetrized triad transfer function <img src="https://latex.codecogs.com/gif.latex?S_k^{p,q}" />.
```
t=0
triadgraph_symmetric_all(S_kpq[t,:,:,:])
```
draws a network graph of the symmetrized transfer at time t=0.
- ***triadgraph_directional_all*** plots a network graph of the directional representation <img src="https://latex.codecogs.com/gif.latex?D_{k&space;\leftarrow&space;q}^p" title="D_{k \leftarrow q}^p" />.
```
t=0
triadgraph_directional_all(D_kpq[t,:,:,:])
```
draws a network graph of the directional representation at t=0.

For more details, see help of each functions and example of usage in DEMO_triadgraph.ipynb.

## Reference
[S. Maeyama, M. Sasaki, K. Fujii, T. Kobayashi, R. O. Dendy, Y. Kawachi, H. Arakawa, S. Inagaki,
"On the triad transfer analysis of plasma turbulence: symmetrization, coarse-graining, and directional representation",
New J. Phys. 23, 043049 (2021).](https://doi.org/10.1088/1367-2630/abeffc)
[![doi](https://img.shields.io/badge/doi-10.1088/1367--2630/abeffc-5077AB.svg)](https://doi.org/10.1088/1367-2630/abeffc)
