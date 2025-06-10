# Satellite Docking Mechanism Control and Velocity Regulation

This project presents a comparative analysis of advanced control strategies for autonomous satellite docking. The focus lies on implementing and evaluating:

- **Sliding Mode Control (SMC)**
- **Model Predictive Control (MPC)**
- **Nonlinear Model Predictive Control (NMPC)**

These methods were tested via **MATLAB** and validated further using **GMAT (General Mission Analysis Tool)** simulations to ensure real-world mission feasibility.

## 🚀 Key Features

- Simulation of Clohessy-Wiltshire dynamics
- Custom SMC implementation with real-time tuning
- Comparative analysis across fuel efficiency, thrust usage, and precision
- **GMAT-based trajectory validation** for mission-level accuracy
- Data visualizations of trajectory, velocity, and energy consumption

## 📁 Project Structure

- `Dissertation_Report.pdf`: Full technical report
- `MATLAB_Code/`: Simulation scripts and control implementations
- `GMAT_Simulations/`: Sample GMAT scripts and output logs
- `Figures/`: Visual results from both MATLAB and GMAT simulations

## 📊 Results Summary

| Controller | Position Error | Energy Used (J) | Max Thrust (N) |
|-----------|----------------|------------------|----------------|
| NMPC      | 0.21 m         | 102,971          | 11.2           |
| MPC       | 1.31 m         | 117,846          | 24.2           |
| SMC       | 4.94 m         | 316,785          | 2079 (chatter) |

## 🛠️ Technologies Used

- MATLAB R2023a
- GMAT R2020a

https://github.com/user-attachments/assets/afb779ab-a32b-46e8-b21c-ee6a48b5433e

https://github.com/user-attachments/assets/c86bbb83-57a5-468c-b148-86d074545bc6
