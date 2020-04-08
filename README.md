# Zostera-marina-IBM

Code for an individual-based model that simulates eelgrass (*Zostera marina*) growth.
The model represents individuals as objects, which develop inside of a spatially explicit world while being exposed to environmental conditions.
Test simulations are available.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3381611.svg)](https://doi.org/10.5281/zenodo.3381611)

## Overview

Four components are needed to make a simulation with the model :

* The code for the model classes (*Zostera_model.py*)
* An initial set of Zostera class objects (**Initial**)
* A grid representing the world in which the objects are located (**World**)
* A time series of fortnightly environmental conditions (**Environment**)

The model simulates the development of the **Initial** inidividuals, developing inside the **World** when exposed to the **Environment** in a fortnightly time step. The output of a simulation is a database with morphological measures of a sample of individuals each time step.

## Prerequisites

Download or clone the necessary files and data sets from the [repository](https://github.com/hugosal/Zostera-marina-IBM).

To run a simulation, the directory *inputs* must exist along *zostera_model.py* and must contain the datasets corresponding to **World**, **Environment**, and **Initial**.

The **Environment** file must be a CSV file. The columns in the file must correspond in order to Sea surface temperature (C°), water temperature anomaly (C°), global horizontal irradiance (GHI; kW fortnight<sup>-1</sup> m<sup>-2</sup>), and air exposure at a depth of 0.5 m (h fortnight<sup>-1</sup>). Optionally, it can contain a fifth row, with the datetime module timestamp of the corresponding row of environmental conditions. The timestamp is only used to plot the corresponding date in an output animation. Each row corresponds to a fortnight of environmental conditions.

The **World** corresponds to a pickable list of three numpy mgrids; the first is the grid with x coordinates, the second with the y coordinates, and the third the depth at each cell.

The **Initial** file may be one of two types. It can be a meadow obtained from the output of a previous simulation. Or it can be a CSV file that will be used to create the corresponding objects. This csv file must contain n rows, corresponding to n new rhizomes to create. The first to the third column of the file contains the longitude coordinate of the first phytomer of the individual, the latitude, and the orientation of the rhizome (rad), respectively. The subsequent columns contain each of the internode lengths of the phytomers in the rhizome.

The test correspond to the simulation of the year 2000 in the Punta Banda Estuary in Baja California, Mexico.

The test simulation of the 2000 is carried out by using as inputs *sample_environment_2000.csv*, *sample_founding_rhizomes_2000.csv* and *cannal_200m_broad_4m_prof* as **Environment**, **Initial** and **World**, respectively. 

In the case of the World for this simulation the file can be obtained by running:
```
$ python sample_world.py
```
This will create the file cannal_200m_broad_4m_prof and print a "Success" message.

A script is available to prepare the files for the test simulation.
```
$ bash zostera_model_test_preparer.sh
```

## Running a simulation

Run the *Zostera_model.py* code and specify the names of the **Initial**, **Environment**, and **World** files to use, in that order:

```
$ python zostera_model.py sample_founding_rhizomes_2000.csv sample_environment_2000.csv cannal_200m_broad_4m_prof.dat
```
A progress bar should be displayed, and at the end, print a *Success* message.

The output of a simulation will be created in a directory named output. The output is comprised of a .dat file, an optional mp4 file, and a CSV file.  The .dat file contains the meadow at the end of a simulation, and this file can be used to run a new simulation. The mp4 file is an animation of the development of the meadow. The CSV file contains the internode length of the phytomer comprising a sample of 20 randomly sampled individuals of the meadow. The data can be analyzed with the supporting .R files.

## Notes
The model uses pseudo-random number generators. The number seed is set in the file *Zostera_model.py*.
The required python modules are;  **random**, **os**, **sys**, **pickle**, **csv**, **pandas**, **math**, **numpy** , **scipy**, **matplotlib**, and **tqdm**.
The file *zostera_plotter.py* can be used to make an animation of the development of the meadow. This function is called directly by the *Zostera_model.py*, but is disabled to reduce time and resource consumption. The plot obtained is a mp4 file, and the ffmpeg program is required to create it. To make an animation change the call from the main loop to:
```
plot_meadow(output, var_maps[0], [grid_x, grid_y, depth_x_y], save=True, show=True)
```
## Authors

**Hugo Salinas** - [hugosal](https://github.com/hugosal)

## License

This project is licensed under the GNU General Public License v3.0 - see the [https://github.com/hugosal/Zostera-marina-IBM/blob/master/LICENSE](LICENSE.md) file for details

## Acknowledgments

* This project was funded by CONACYT [grant number 836404/634358] and CICESE.

## References

Salinas, Hugo. (2019). Simulación dinámica de una pradera de Zostera marina. Master of Science degree in Marine Ecology thesis. Center for Scientific Research and Higher Education of Ensenada (CICESE), Baja. California, Mexico.
