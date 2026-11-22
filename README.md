# Active Suspension Control — Half-Car Model (MATLAB / Simulink / Stateflow)

Modeling, structural analysis, and advanced control synthesis for an automotive active suspension system based on a **nonlinear 8-state Half-Car model**, developed in **MATLAB, Simulink, and Stateflow**.

📄 **Full Technical Report:** [Download PDF](doc/Half_Car_Control_Report.pdf)

---

## Project Overview

The goal of this project is to maximize ride comfort and vehicle stability (minimizing body heave, pitch angle, and their respective velocities) under realistic road disturbances:
- **Road Bump:** 4 m length, 15 cm height at 40 km/h
- **Uneven Road Profile:** ±10–15 cm random irregularity at 40 km/h

## Key Features & Control Strategies

1. **Modeling & Analysis**
   - Nonlinear 8-state Half-Car dynamic model (sprung mass heave & pitch, front/rear unsprung masses) and small-angle linearization.
   - Internal stability verification via eigenvalues and Lyapunov LMI feasibility; full reachability and observability analysis.
2. **Classical Control**
   - Static state-feedback **Pole Placement** validated on the nonlinear Simulink plant.
3. **LMI-Based Optimal Control (Static & Dynamic Output Feedback)**
   - Synthesis and comparison of **$L_1$**, **$H_2$**, and **$H_\infty$** controllers (full state feedback and partial output measurement feedback).
4. **Multi-Objective Control**
   - Mixed **$H_2 / H_\infty$** static control and **$H_\infty$ with regional pole placement ($\mathcal{R}$-stability)**.
5. **Robust Control (Polytopic LPV)**
   - Robust **$H_2$** controller synthesized over a 4-vertex polytopic representation accounting for **±15% parametric uncertainty** on front and rear suspension stiffness ($k_{sf}, k_{sr}$).
6. **Gain-Scheduling with Simulink Stateflow**
   - Finite State Machine implemented in **Stateflow** to switch dynamically between three driving setups (**Normal**, **Comfort**, **Sport**) using offline-synthesized $H_2 / H_\infty$ gains on the LPV Simulink plant.

## Repository Structure

- `src/` — MATLAB scripts (LMI synthesis, system matrices, plotting) and Simulink/Stateflow models (`.slx`)
- `docs/` — Complete technical report with mathematical derivations, LMI formulations, and simulation benchmarks
