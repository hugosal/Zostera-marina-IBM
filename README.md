# Zostera-marina-IBM
Code for an individual based  model that simulates eelgrass (Zostera marina) growth

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3381611.svg)](https://doi.org/10.5281/zenodo.3381611)

This repository has the code necessary to make a simulation of the development of an eelgrass (Zostera marina) meadow.
The model uses instances of object classes to represent the individuals that form a seagrass  meadow. The member data of each Class represents some attribute of a real-world object, and the class function the processes that affect them. 
The simulated seagrass of the model develop under certain ambient conditions and are located spatially in a point of a “world.”
The code is in Python 3, and to perform a simulation the scipy, pandas, numpy, matplotlib, tqdm modules are required additionally to some other base modules.
To save an animation of the meadow growth ffmpeg is required.
The repository also has the necessary files and code to run the simulations used in Salinas (2019) to validate the model.
A complete description and discussion of the model is presented in:

Salinas, Hugo. (2019). Simulación dinámica de una pradera de Zostera marina. Master of Science degree in Marine Ecology thesis. Center for Scientific Research and Higher Education of Ensenada (CICESE), Baja. California, Mexico.
