# RabiModel
Static and dynamical properties of states in the Quantum Rabi Model and variations

clone as:


```bash
git clone https://github.com/djuliannader/RabiModel.git
```

## Data 

This repository provides data and codes for the paper Quantum non-Gaussian criticality in the quantum a Rabi model with a small modulation of the
qubit frequency, by D. J. Nader, P. Stransky, J. Chavez, P. Cejnar and R. Filip.
The paper is available on arXiv.
The repository is currently being developed; please await its completion. Some data for the figures can be found in folder /data1.  The files contain Mathematica notebooks, necessary to recreate the Figures of the article.


## Usage

For using the code move to the folder /src and go to the corresponding file

For stationary states, manipulate parameters from input.dat and run as

```bash
julia main.jl
```

For the survival amplitude of quench dynamics, manipulate parameters from main_dqpts.jl and run as :

```bash
julia main_dqpts.jl
```

For the survival amplitude of quench dynamics, with initial thermal state, manipulate parameters from main_dqpts.jl and run as :

```bash
julia main_dqpts_thermal.jl
```
