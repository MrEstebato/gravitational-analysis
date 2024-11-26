# gravitational-analysis

This project simulates a dynamic Voronoi diagram where points orbit a central point, creating an ever-changing Voronoi tessellation.

## Animation

![Voronoi Animation](animacion_voronoi.gif)

---

## Getting Started

### Running the Simulation

1. Clone this repository:
   ```bash
   git clone https://github.com/MrEstebato/gravitational-analysis.git
   cd gravitational-analysis
   ```

2. Open `solarSystem.m` in MATLAB.

3. Run the script:
   ```matlab
   solarSystem
   ```

The simulation will display a dynamic Voronoi diagram with orbiting points.

---

## File Descriptions

- **`solarSystem.m`:** Main script for initializing the simulation and rendering the dynamic Voronoi diagram.
- **`voronoiDiagram.m`:** Function for computing and plotting Voronoi cells dynamically.

---

## Customization

Modify the following parameters in `solarSystem.m` to customize the simulation:


- **`padding`:** Space around the plot for better visualization.

## Execution
Give the following data once you excecute `solarSystem.m` as inputs required by the program:
- A float number that indicates the x position of the central point.
- A floar number that indicates the y position of the central point.
- A float number that indicates the mass of the central point.
- An integer number (n) that indicates the number of points to be provided.
- Then, introduce n times the following data, each of them for one point.
1. A float number that indicates the distance from the central point.
2. A float number that indicates the the angular speed of the point in grades by frame (°/frame). 
3. A float number that indicates the mass of the point.

## Example
Copy and paste the following information to see an example of how works the program:
0
0
10
8
1
0.05
5
1.5
0.15
3
2
0.065
5
2.5
0.03
1
3
0.07
7
3.5
0.04
2
4
0.08
5
4.5
0.01
7