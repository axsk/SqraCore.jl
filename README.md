# SqraCore.jl
Approximation of the generator of a Diffusion process by pointwise potential evaluations.

## Theory

Using the Square-Root approximation (SQRA) as in [6]:

The generetor of the Koopman operator of the diffusion in a potential $U$: $dX = -\nabla U dt + \sigma dB_t$ 
on a Voronoi tesselation of the state space is approximated by the matrix $Q$:

$$Q_{ij} = C_{ij} \sqrt\frac{\pi_j}{\pi_i}$$
for $i\neq j$ and $Q_{ii} = -\sum_{j\neq i} Q_{ij}$ with $$C_{ij} = \beta^{-1}\frac{A_{ij}}{h_{ij} V_i}$$ where $A_{ij}$ and $V_{ij}$ are the areas resp. volumes of the cells (common boundary), $h_{ij}$ the distance between the centers of cells $i,j$ and $$\pi_i = \exp (-\beta U(x_i))$$ the stationary distribution at the cell centers.

## References

- [1] Weber (2010). "A subspace approach to Molecular Markov State Models via a New Infinitesimal Generator"
- [2] Lie, Fackeldey, M. Weber (2013)
- [4] Donati, Heida, Keller, Weber (2018). "Estimation of the infinitesimal generator by square-root approximation."
- [3] Donati, Weber, Keller (2021). "Markov models from the square root approximation of the Fokker–Planck equation: calculating the grid-dependent flux."
- [5] Donati, Weber, Keller (2022). "A review of Girsanov reweighting and of square root approximation for building molecular Markov state models."
- [6] Schütte, Klus, Hartmann (2022). "Overcoming the Timescale Barrier in Molecular Dynamics: Transfer Operators, Variational Principles, and Machine Learning."
