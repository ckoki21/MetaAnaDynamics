# MetaAnaDynamics
This repository contains the trajectory data and Stan code for the paper "Data-Driven Modelling of Kinetochore Dynamics: Space-Time Organisation of the Human Metaphase Plate."

To gain insights into chromosome dynamics at the metaphase-anaphase transition, we developed a tracking algorithm that achieves near-complete tracking of fluorescently labeled kinetochores. This was achieved using an endogenous label of a kinetochore protein. We performed live-cell imaging of non-transformed human RPE1 cells using lattice light-sheet microscopy (LLSM) and generated tracks with our tracking-pairing pipeline. Data were collected at a high temporal resolution of 2.05 seconds per z-stack over extended timescales, typically tens of minutes, starting during prometaphase and continuing through to anaphase.

This repository includes data from 26 DMSO-treated cells (MetaAnaDynamics/Data/DMSO) and 33 nocodazole washout treated cells(MetaAnaDynamics/Data/Nocodazole_washout).

To run the models, go to R_files --> Choose one of the "Asymmetry", "MetaphaseAnaphaseAsymmetry", "TimeDependent" folders, depending on the preferred variation of the biophysical models.  -->choose "run_...._models.R" to run all the trajectories of a cell. Note: Within the R file, choose the cell that you want to run by commenting out all the others!



