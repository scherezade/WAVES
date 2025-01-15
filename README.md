# WAVES
WAVES (Wire Array for Variable Exhaust Scannig):   Algorithm for the reconstruction of the time-varying plume cross-section of unsteady electric thrusters. 

Scherezade Barquero, Mario Merino, Jaume Navarro Cavallé. Space and Plasmas Propulsion Team (EP2). Universidad Carlos III de Madrid.

---------------------------------------------------------------------------------------------------------------------------------------------------------

A novel diagnostic system to characterize the transient exhaust of electrical plasma thrusters is presented. The proposed approach time-reconstructs the
cross-sectional expansion of their plume. This experimental technique consists of an array of electrostatic wire probes biased at the ion saturation regime. 
The two-dimensional and time-dependent ion current distribution is reconstructed from the probe data using a variable separation algorithm. 

The inverse problem to solve constitutes a nonlinear system of equations and is described in the following journal paper:
"Reconstruction of the transient plume cross-section of a Pulsed Plasma Thruster", S. Barquero, M. Merino, J. Navarro-Cavallé. 
Plasma Sources Science and Technology, 2025. (Currently under review)

---------------------------------------------------------------------------------------------------------------------------------------------------------

This repository contains the MATLAB implementation of the proposed algorithm, along with a linearized version of the problem. Both are included in the file 
'transientplume_crosssection_scanner_algorithm.m'. To solve the problem and perform postprocessing, run the provided script.
Key considerations:
1) The script includes a USER section for customization.
2) The probe data is available at the following URL/DOI: 10.5281/zenodo.13820945

---------------------------------------------------------------------------------------------------------------------------------------------------------

The WAVES repository can be cited as follows:
BALSERA BARQUERO, M. S., Merino, M., & NAVARRO CAVALLE, J. (2025). scherezade/WAVES: WAVES (v1.0). Zenodo. https://doi.org/10.5281/zenodo.14653075

 


