# RabiModel
Static and dynamical properties of states in the Asymmetric Quantum Rabi Model with an small modulation in the qubit frequency

clone as:


```bash
git clone https://github.com/djuliannader/RabiModel.git
```



## Data 

This repository provides data and codes for the paper Quantum non-Gaussian criticality in the quantum a Rabi model with a small modulation of the
qubit frequency, by D. J. Nader, P. Stransky, J. Chavez, P. Cejnar and R. Filip.
The paper is available on arXiv.
The repository is currently being developed; please await its completion. Some data for the figures can be found in folder /data1.  The files contain Mathematica notebooks, necessary to recreate the Figures of the article.


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

Edit the parameters directly in the file main_adiabaticramp.jl and run:

```bash
julia main_adiabaticramp.jl
```
