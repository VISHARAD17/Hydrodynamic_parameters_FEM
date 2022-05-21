## Added Mass and Damping Coefficient of a typical ship (Barge) section in Heave, Sway and Roll motion
Implementation of **Finite element method** to calculate Added Mass and Damping Coefficient using `MATLAB` using radiation and diffraction pheomena.
Steps :
- Find out boundary conditions
- Discretize the domain region into rectagular nodes.
- Calculate general Stiffness Matrix by calculating element stiffness matrix.
- Incorporate Boundary condition into the Global Siffness Matrix.
- Calculate radiation and force matrix using stiffness Matrix.
- Calculate Added mass and Daming coefficient for each mode of motion
