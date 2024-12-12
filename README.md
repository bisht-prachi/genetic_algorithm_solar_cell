# Parameter Extraction of Solar Cells using Genetic Algorithm

## Overview
This project implements a **Genetic Algorithm (GA)** to extract key parameters of solar cells from their raw current-voltage (V-I) data. The model uses both the **Single Diode Equation (SDE)** and **Double Diode Equation (DDE)** to accurately represent solar cell behavior, including parasitic effects. This work provides insights into theoretical modeling and solar cell performance, aiding in the optimization of photovoltaic devices for space and terrestrial applications.

## Features
- Models solar cells as unbiased pn diodes under illumination using SDE and DDE.
- Incorporates effects of parasitic resistances: **Series Resistance (Rs)** and **Shunt Resistance (Rsh)**.
- Employs **Genetic Algorithm** to optimize parameters:
  - **SDE:** Isc, Io, Rs, Rsh, η.
  - **DDE:** Isc, Io1, Io2, Rs, Rsh, η1, η2.
- Efficiently handles non-linear, implicit, transcendental equations.
- Supports analysis of single-junction and multi-junction solar cells.

## Motivation
Parameters like Isc (short-circuit current), Io (reverse saturation current), and parasitic resistances significantly influence the performance of solar cells but cannot be measured directly. Traditional analytical methods fail to account for all parameters accurately. Heuristic approaches like GA provide a robust alternative for multi-parameter optimization.

## Genetic Algorithm (GA)
### Key Features
- **Initialization**: A random population of chromosomes is created.
- **Selection**: Fittest chromosomes are selected for reproduction.
- **Crossover**: Top 30% of the population undergoes crossover to produce the next generation.
- **Mutation**: Random bit-flipping introduces diversity.
- **Termination**: Process ends when the error falls below a set tolerance or a maximum number of generations is reached.

### Fitness Function
The GA minimizes the error between measured and calculated current values:

\[
f = \sqrt{\frac{1}{n} \sum_{j=0}^{n} \left(\frac{I_{j, \text{measured}} - I(V_j)}{I_{j, \text{measured}}}\right)^2}
\]

### Implementation Details
- Operating temperature: **300 K**.
- Population size: **1000 chromosomes**.
- Number of generations: ~**100**.
- Parameter bounds: Determined via visual inspection or preliminary direct search methods.

## Results

### Application to Single Junction Solar Cells
The GA was applied to GaInP, (In)GaAs, and Ge solar cells:

### V-I Curve Validation
The extracted parameters reproduce V-I characteristics closely matching the measured data:

![image](https://github.com/user-attachments/assets/48a110a6-87ee-4373-ba6b-02be1f447515)
![image](https://github.com/user-attachments/assets/1b0dc06d-5b00-48b3-87b2-ea030ab1c733)

### Application to Single Junction Solar Cells
The GA was applied to triple-junction solar cells stack:

![image](https://github.com/user-attachments/assets/84dd1137-8f81-4eca-bd56-d12ae0d3ee85)
![image](https://github.com/user-attachments/assets/839a6b1b-e6a3-41f6-aab4-a73856e10aa6)




## How to Use
### Prerequisites
- C++ compiler supporting C++11 or higher (e.g., GCC, Clang).

### Running the Code
1. Clone this repository:
   ```bash
   git clone https://github.com/your-repo/solar-cell-ga.git
   ```
2. Navigate to the project directory:
   ```bash
   cd solar-cell-ga
   ```
3. Compile the source code:
   ```bash
   g++ -o sd_ga sd_ga.cpp
   ```
4. Run the executable:
   ```bash
   ./sd_ga
   ```

## References
1. Kassap, S. *Principles of Electronic Materials and Devices*.
2. Goldberg, D. E. *Genetic Algorithms in Search, Optimization, and Machine Learning*.
3. Holland, J. H. *Adaptation in Natural and Artificial Systems*.

## Acknowledgments
- Project by **Prachi Bisht** under the guidance of **Shri Suresh E. P., Solar Panel Division, URSC, ISRO**.
- Special thanks to the URSC Solar Panel Division for their support and resources.
