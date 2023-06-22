# Codes_Manipulative_Cheat
This repository is for open access to source codes and data for the manuscript entitled **"The evolution of manipulative cheating"**.

**Codes** contains source codes of all models: *numerical* is the code for generating analytical predictions, and *simulative* is the code for running individual-based simulations. We have also uploaded some R codes for plotting and performing harmonic regression in *analysis* folder. For the C files in *codes/simulative*, **please** make sure the location of random number generators files match the description in C files (i.e., *#include "../../dSFMT-src-2.2.3/dSFMT.c"*) because those files are essential for running. Also, **please** read the notes at the beginning of each file as they provide descriptions of key parameters and instructions for running. The new file in *code/simulative/recombination* is the code used in section 4.8 of the Appendix, with continuous probability of exchanging trait values between individuals.

**Data** contains the raw data for making figures, all files are organised according to the figure numbers.

Feel free to email *ming.liu.ac[at]gmail.com* if you have any problem! (Ming Liu)
