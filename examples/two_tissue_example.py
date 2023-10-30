#Load packages
import numpy as np
import synmorph as sm

#Specify simulation parameters
grn_params = {"n_AVE_cells": 20, #Number of DVE cells to initialise with
              "AVE_alpha_dir": 0.2, #Migration orientation coupling strength \alpha
              "non_AVE_alpha_dir": 0., #Migration coupling strength for non-DVE cells. Set to 0
              "AVE_v0": 0.05, #DVE migration speed v_{0,DVE}
              "non_AVE_v0": 0., #speed for nonDVE cells
              "AVE_alpha0": 0., #Initial orientation of migration. Over-riden in grn_ave_couple_orientation
              "boundary_frac": 0.08, #\phi_{bound}
              "AVE_A0": 0.54, #Value of A_0 for DVE cells. Subsequently, non-DVE A_0s are recomputed to prevent compression/stretch in the bulk
              "exe_frac": 0.0, #\phi_{exVE}. If this is set to 0, no exVE cells are found.
              "AVE_p0": 3.4, # Effective descriptor of \Gamma_{DVE}, given rescaling defined in Supplementary Modelling
              "nonAVE_p0": 4.2,#Likewise for \Gamma_{emVE}
              "ExEVE_p0": 4.0} #and \Gamma_{exVE}
tissue_params = {"L": 17.0,
                 "A0": 1.,
                 "P0": 3.2,
                 "kappa_A": 1.,
                 "kappa_P": 0.2,
                 "W": (np.array(((0.0, 0, 0, 0.1), (0, 0, 0, 0.5), (0, 0, 0, 0.5),
                                 (0.1, 0.5, 0.5, 0.1))) * 1).astype(np.float32), #J_ij in the text
                 "a": 0., # a and k define "soft" repulsion among cell centres, set to 0. See Barton et al., 2017.
                 "k": 0.}
active_params = {"v0": 0, #active velocity. This is set to 0 in synmorph > grn > grn_ave_couple_orientation.py
                 "Dr": 5e-3} #Persistence timescale aka extent of noise in active migration.
init_params = {"init_noise": 0.1, #Amount of noise in the initial hexagonal tiling of cells
               "c_type_proportions": (1.0, 0.0)} #proportions of cell types. Over-riden in grn_ave_couple_orientation.py
run_options = {"equiangulate": True, #Determines if equiangulation protocol (similar to Barton et al., 2017) is used to speed up repeated Delaunay Triangulations.
               "equi_nkill": 10} #Max iterations of equiangulation protocol before a new triangulation is generated.
simulation_params = {"dt": 0.1, #time-step \delta_t
                     "tfin": 300, #final time t_{fin}
                     "tskip": 20, #Number of time steps skipped during saving
                     "dt_grn": 0.1,#time-step of the "grn" module, i.e. the module updating DVE flocking behaviour etc.
                     "grn_sim": "grn_ave_couple_orientation", #refers to the flocking module
                     "tinit": 10, #Total time to do passive initalisation
                     "random_seed": int(2023)} #Random seed
save_options = {"save": "hdf5",#Saving mode. HDF5 is the most compact.
                "result_dir": "../scan_results",#Saving file
                "name": "AVE_example_full",#Simulation name
                "compressed": True}#Whether to zip the saved file.

scan_dict = {"tissue_params": tissue_params, "active_params": active_params, "init_params": init_params,
             "run_options": run_options, "simulation_params": simulation_params, "grn_params": grn_params,
             "save_options": save_options} #Assemble the parameters.


sim = sm.simulation(tissue_params=tissue_params,
                    active_params=active_params,
                    init_params=init_params,
                    simulation_params=simulation_params,
                    grn_params=grn_params,
                    run_options=run_options,
                    save_options=save_options) #Instantiate the simulation class

sim.save_dir_plots = "results" #Folder to save the plots
sim.simulate(progress_bar=True) #Simulate

sim.animate_c_types(n_frames=15,
                    c_type_col_map=["#4bdb71", "#ffbb4d","#ffbb4d","white"],
                    file_name="two_tissue") #Animate the simulation.

