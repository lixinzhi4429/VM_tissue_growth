# 2D Vertex Model with Cell Division, Mechanical Feedback, and Polarity Alignment

This repository contains a **MATLAB implementation of a two-dimensional vertex model for epithelial tissues**, incorporating **cell growth and division**, **mechanical feedback**, and **cell polarity alignment**.  
The code is designed for **theoretical and computational studies of tissue mechanics, morphogenesis, and active matter**, and is suitable for producing reproducible simulation data and figures for research publications.

---

## Overview

Vertex models represent epithelial tissues as polygonal tilings, where cell mechanics emerge from energy penalties associated with **area, perimeter, and interfacial tensions**.  
This implementation extends the standard 2D vertex model by including:

- **Cell growth and division**
- **Mechanically regulated feedback on growth and division**
- **Cell polarity fields with alignment interactions**
- **Topological rearrangements (T1 transitions)**

The framework is modular and intended for **systematic exploration of how mechanics, polarity, and division couple to shape tissue dynamics**.

---

## Model Features

### 🔹 Mechanical Energy
Each cell contributes to the total energy via:
- Area elasticity
- Perimeter elasticity
- Edge (line) tension

A typical energy functional has the form:
\[
E = 1/2 \sum_i \frac{K_A}{2}(A_i - A_{0,i})^2
+ 1/2 \sum_i \frac{K_P}{2}(P_i - P_{0,i})^2
+ 2 \sum_{e} \gamma_e L_e
\]

---

## Requirements

- MATLAB R2021a or newer  
- No external toolboxes required


## Citation

If you use this code in published work, please cite the corresponding paper.
https://doi.org/10.1101/2024.09.03.610990

---

## License

MIT License

---

## Contact

Dr. Xinzhi Li 
xli3267@gatech.edu
School of Physics, Georgia Institute of Technology
