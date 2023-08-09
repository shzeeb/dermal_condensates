# dermal_condensates

This readme contains description of the code and the data accompanying the article: Rapid mechanosensitive migration and dispersal of newly divided mesenchymal cells aid their recruitment into dermal condensates.

1. dermal_condensates_v10.m: codes up the model described in the appendix Section 1 (model in which upon division jump a mother daughter pair is produced which perform persistent motion for a period of 180 minutes. Saves data in .mat format after running. All parameters are explained at the beginning of the code.

2. dermal_condensates_v9.m: cells diffuse and divide and get adsorbed when hit a condensate. At division, cells perform a jump over 10 mins. No persistent motion here.

3. dermal_condensates_v11.m: Cell proliferation occurs but the parent cell does not produce a daughter cell. The parent cell perform the solo mitotic jump followed by persistent motion for 180 mins (3 hrs).

4. plot_props.m: uses datasets in (a), (b) and (c) to plot mean proportions over a sample of 8 repeats and plot the corresponding 95% confidence interval on top. User specifies which of the datasets (a), (b) or (c) they would like to plot proportions for. There is also the option to plot simulated or experimentally observed proportions. This can be elected by changing the values of the binary variable "simulated" at the top of the .m script where a 1 plots the simulated proportions and 0 plots the experimentally observed proportions.

5. msd_non_dividers.m: Estimate the diffusion coefficient by calculating the
mean squared displacement (msd) for all the tracks in each of the M=8 datasets (d).

6. persistence_angles_func.m: analyses the persistence angles of dividing cells.

7. persistence_VM_fit.m: finds an optimal parameter kappa_star of von-Mises (VM) which minimises the sum of squared residuals between the empirical distribution and the fitted curve.

8. division_jump_size.m: estimates the average size of the jump at mitosis.

9. mitosis_angle_analysis.m: analyses division angles of dividing cells and fits a VM distribution to the angles.

10: vmrand.m: generates random number from the VM distribution with given parameters mu and kappa. Ref: Dylan Muir (2023). vmrand(fMu, fKappa, varargin) (https://www.mathworks.com/matlabcentral/fileexchange/37241-vmrand-fmu-fkappa-varargin), MATLAB Central File Exchange

11. follicles.m: generates a triangular lattice of hair follicles with hexagonal packing. Used in (1), (2) and (3).

12. pdf_dist_2: analysis of the distances between cells and follicles for divided and removed cells and for cells which did not divide and got adsorbed.
End result: a probability distribution of starting a distance r away from the neatest follicle conditional on the fact that the cell divided (or did not divide) and got adsorbed later on.

13: distance_density_plot.m: normalises the distribution obtained in (12) are plots it.

14. msd_persistence.m: Calculates mean square displacement for cells after division to check for persistence. Track length up to 18 length corresponding to 180 minutes (or 3 hrs) is used as is the suggested duration of persistence according to data.

DATA SETS:

(a) dermal_condensates_with_persistence_1/2/3/4.mat (four files). Simulations for model with persistence.

(b) dermal_condensates_no_persistence_1/2/3/4.mat (four files). Simulations for model without persistence.

(c) dermal_condensates_no_progeny_1/2/3/4.mat (four files). Simulations for model where no daughter cell is produced upon division.

(d) Video 1/2/.../8 (8 folders). Empirical cell tracking files for dividers (mother and daughter) and non-dividers.

(e) Division_angle_data (1 folder). Contains 7 csv files containing division angle data of cells from experiments.


