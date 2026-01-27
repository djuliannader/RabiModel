# RabiModel
Static and dynamical properties of states in the Asymmetric Quantum Rabi Model with an small modulation in the qubit frequency

clone as:


```bash
git clone https://github.com/djuliannader/RabiModel.git
```



## Data 

This repository provides codes and data for the papers

- Controllable distribution of nonclassicality in qubit-oscillator systems, by D. J. Nader

- Sensitivity enhanced by qubit modulation, by D. J. Nader, P. Stransky, J. Chavez, P. Cejnar and R. Filip.

- Wigner negativity as a precursor of dynamical phase transitions, by D. J. Nader, P. Stransky, J. Novotny, P. Cejnar and R. Filip.


The repository is currently being developed; please await its completion.
The Mathematica notebooks, necessary to recreate the Figures of the articles can be found in the directory Figures/


## Environment Requirements  

To run this repository, please make sure the following environment is available:

- LinearAlgebra (stdlib)  
- DifferentialEquations
- HCubature
- QuantumOptics
- Convex 
- SCS
- PyPlot

## Usage

To run the code, navigate to the /src folder and execute the corresponding file:

- Stationary and Floquet states

Edit the parameters in input.dat and run:

```bash
julia main.jl
```

- Survival amplitude of quench dynamics

Edit the parameters directly in the file main_dqpts.jl and run:


```bash
julia main_dqpts.jl
```

- Survival amplitude of quench dynamics with an initial thermal state

Edit the parameters directly in the file  main_dqpts_thermal.jl and run:

```bash
julia main_dqpts_thermal.jl
```

- Adiabatic ramp for the preparation of Floquet states

Edit the parameters directly in the file Floquet_ramp.jl and run:

```bash
julia Floquet_ramp.jl
```

- Adiabatic ramp for the transference of quantumness

Edit the parameters directly in the file Adiabatic_ramp.jl and run:

```bash
julia Adiabatic_ramp.jl
```
