# Electromagnetic Waves
 Solvers for Maxwells Equations

## Summary

This experiment was meant to develop a modular and easy to use class and methods for doing vector calculus.  I later experimented with makeing the entire simulation environment itself a class, with its evolution routines enclosed as internal methods.  A "C++ inspired" implementation, this makes the files rather light, reduces the number of functions a user would need to memorize, and reduces the amount of code required to generate a rudimentary script.  The vector calculus routines included here are by no means exhaustive, but enough to solve Maxwell's equations in a vaccum, and are readily adaptable to other vector field type PDEs such as Navier-Stokes (if indeed one really wanted to spend all that CPU time on a Navier-Stokes experiment), or other vector valued systems like the Landau-Lifshitz Equations.

## Components

### field.m
This class handles all the vector calculus, it contains the methods for initializing and plotting vector fields in 2D with three components.  Coordinate systems are also realized as vector fields.  This class also has the methods for taking derivatives with boundary conditions.

### shape.m
This class handles geometry and is utilized by the simulation environment MaxSystemDielectric.m to stand for solid dielectric materials.  Each shape has minimal properties: position, radius, and ep stands for the dielectric constant which is 1 for a vacuum.  And methods for drawing the shape, and detecting whether a point is inside or outside each shape.  

### MaxSystem.m
Class and methods for solving maxwells equations in a vacuum.

### MaxSystemDielectric.m
Class and methods for solving maxwells equations in a mixed vacuum/solid dielectrics system.  solid geometry is stored as a member array of the class shapes.  Unlike the simpler vacuum, this takes 3 steps to initialize, (1) set initial conditions, (2) append shapes, (3) set the spatial dependence of the dielectric permeability epsilon.  see the first three methods of this class.
