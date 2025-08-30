# Fortran
# Computational Physics Assignments (Fortran Codes)

This repository contains my Fortran implementations of **six computational physics assignments** completed as part of coursework. Each assignment covers important topics in numerical methods, random processes, Monte Carlo methods, statistical mechanics, differential equations, Turing patterns, and molecular dynamics.

---

## 📂 Repository Structure

├── Assignment1/ # Random Numbers & Random Walks
├── Assignment2/ # Monte Carlo Integration & Random Distributions
├── Assignment3/ # 3D Ising Model Simulation
├── Assignment4/ # ODE & PDE Solvers
├── Assignment5/ # Turing Pattern Formation
├── Assignment6/ # Molecular Dynamics Simulation
└── README.md

---

## 📘 Assignment Summaries

### **Assignment 1: Random Numbers & Random Walks**
- Generate random numbers with different seeds.  
- Write results to files with proper formatting.  
- Calculate averages of random samples.  
- Study sums of random numbers and their distributions.  
- Simulate 1D random walks and analyze with Gaussian fitting.  

---

### **Assignment 2: Monte Carlo Integration**
- Numerical integration using trapezoidal rule with error analysis.  
- Verification using functions like \( \frac{1}{1+x^2}, \sin(x), \) and Gaussian.  
- Random number analysis: uniformity, correlations, standard deviation.  
- Generate random numbers with exponential & Gaussian distributions.  
- Monte Carlo integration with brute-force and importance sampling.  

---

### **Assignment 3: 3D Ising Model Simulation**
- Implement Ising model on an \( L \times L \times L \) cubic lattice with periodic boundary conditions.  
- Compute magnetization, energy fluctuations, susceptibility, and heat capacity.  
- Perform finite-size scaling across different lattice sizes.  
- Binder cumulant analysis.  

---

### **Assignment 4: ODE & PDE Solvers**
- **ODEs**:  
  - Solve \( \frac{dy}{dx} = y^2 + 1 \) using Euler, Modified Euler, Improved Euler, and RK4 methods.  
  - Compare results with the exact solution \( y = \tan(x) \).  
  - Solve nonlinear oscillator equations with different initial conditions.  
  - Simulate coupled oscillators in a circular lattice with RK4.  
- **PDEs**:  
  - Solve 2D Laplace’s equation on a \( 34 \times 34 \) plate.  
  - Implement Dirichlet and Neumann boundary conditions.  
  - Plot temperature distributions and analyze convergence.  

---

### **Assignment 5: Turing Pattern Formation**
- Implement reaction–diffusion equations (Turing model).  
- Study emergence of spatially periodic patterns.  
- Simulate and visualize stripe/spot-like structures.  
- Analyze role of diffusion coefficients and reaction rates.  

---

### **Assignment 6: Molecular Dynamics Simulation**
- **Q1–Q3:**  
  - Perform MD simulation with 2197 particles in a cubic box (20×20×20).  
  - Use Lennard-Jones potential with cutoff \( r_c = 2.5\sigma \).  
  - Implement periodic boundary conditions (PBC).  
  - Check energy and momentum conservation (NVE ensemble).  
  - Introduce thermostat (NVT ensemble) and study energy fluctuations.  
- **Q4–Q7:**  
  - Implement Verlet neighbour list to optimize pair interactions.  
  - Compute **radial distribution function \( g(r) \)** with binning.  
  - Study how \( g(r) \) evolves with particle number (1200, 2400, 3600).  
  - Identify peaks in \( g(r) \) and analyze liquid structure.  
- **Q8:**  
  - Verify Maxwell–Boltzmann velocity distribution of particles.  
  - Compare numerical histogram with exact theoretical distribution.  

---

## 🛠️ Tools & Language
- **Language**: Fortran 90/95  
- **Compiler**: `gfortran`  
- **Plotting**: External tools (Python/Matplotlib, gnuplot, or Excel) for visualization.  

---

## ▶️ How to Run
Compile and run any Fortran code using:

```bash
gfortran filename.f90 -o output
./output
📊 Results & Plots

Plots generated (distributions, random walk histograms, Ising model observables, ODE/PDE solutions, Turing patterns, MD radial distribution functions) are stored in respective assignment folders as .pdf or .png files.

📌 Notes

All codes are written in modular and well-commented Fortran for clarity.

Assignments demonstrate applications of computational physics methods such as Monte Carlo, ODE/PDE solvers, pattern formation, statistical mechanics, and molecular simulations.
