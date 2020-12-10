# triadgraph
Python script for the analysis of energy transfer via triad interactions

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
- ***symmetrize_triadtransfer*** returns the symmetrized triad transfer function <img src="https://latex.codecogs.com/gif.latex?S_k^{p,q}" /> from the non-symmetric triad transfer function.
- ***directional_triadtransfer*** returns the directional representation <img src="https://latex.codecogs.com/gif.latex?D_{k&space;\leftarrow&space;q}^p" title="D_{k \leftarrow q}^p" /> from the symmetrized triad transfer function.
- ***triadgraph_symmetric_all*** plots a network graph of the symmetrized triad transfer function <img src="https://latex.codecogs.com/gif.latex?S_k^{p,q}" />.
- ***triadgraph_directional_all*** plots a network graph of the directional representation <img src="https://latex.codecogs.com/gif.latex?D_{k&space;\leftarrow&space;q}^p" title="D_{k \leftarrow q}^p" />.

For more details, see help of each functions and example of usage in DEMO_triadgraph.ipynb.

## Reference
\[1\] S. Maeyama, M. Sasaki, K. Fujii, T. Kobayashi, R. O. Dendy, Y. Kawachi, H. Arakawa, S. Inagaki,
"On the triad transfer analysis of plasma turbulence: symmetrization, coarse-graining, and directional representation",
under review.
