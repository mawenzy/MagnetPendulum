# MagnetPendulum
This project holds scripts to calculate the movement of a magnetic pendulum under the influence of magnets. Two example results can be seen here:

<p align="center">
  <img src="examples/3_magnets_2500x2500.png" width="350"/>
  <img src="examples/4_magnets_2500x2500.png" width="350"/>
</p>

## Setup
To simplify the calculations, the pendulum only moves in the x-y plane.
The forces which are applied to the pendulum are: *drag*, *spring force* and *magnetic force*. This results in the differential equation:

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\ddot{x}&space;&plus;&space;\gamma&space;\cdot&space;\dot{x}&space;&plus;&space;\kappa&space;\cdot&space;x&space;=&space;F_M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ddot{x}&space;&plus;&space;\gamma&space;\cdot&space;\dot{x}&space;&plus;&space;\kappa&space;\cdot&space;x&space;=&space;F_M" title="\ddot{x} + \gamma \cdot \dot{x} + \kappa \cdot x = F_M" /></a>
</p>

## Calculation
Using the differential equation, the movement of the pendulum in the two-dimensinal x-y space can be calculated numerically. 

To make those colourful pictures, for every point in a limited area, the differntial equation is solved numerically, until the pendulum ends up at one of the magnets. The starting point then gets a value assigned to it corresponding to the magnet. From this file you are able to create a picture.

<p align="center">
  <img src="examples/2_magnets_2500x2500.png" width="350"/>
  <img src="examples/8_magnets_2500x2500.png" width="350"/>
</p>

## Content

### [magneticPendulum.c](src/magneticPendulum.c)

Source code to create a file from which a picture can be created. This file contains a matrix of numbers, one for each pixel.

Usage for the compiled programm:
```
magnets.exe [-qhu] [-o outputFile] [-f force] [-m magnets] [-z zoom] [-y drag] [-k spring_constant]
```
Disclaimer:

The programs stop condition is solely based on the kinetic energy the pendulum has left. After the energy drops below a certain point, then the program stops and calculates the color/number corresponding to the circle sector the pendulum is in. This stop condition is far from physically correct. Given a setup without any magnets, the pendulum would behave as a regular harmonic pendulum. When it reaches the turning point, it\`s kinetic energy is zero, terminating the program, although the pendulum still has potential energy. Of course the pendulum would move further in this case, and if drag is taken into account, move closer and closer to the equilibrium position.

### [magnet.txt](src/magnet.txt)

Example script, which would allow you to convert the file containing numbers from [magneticPendulum.c](src/magneticPendulum.c) to an image using gnuplot.
