---LUTsyn model demonstration files---
by Duy-Tan Jonathan Pham (duytanph@usc.edu)
August 11, 2021

See paper: "Bridging hierarchies in multi-scale models 
		   of neural systems: Look-up tables enable computationally 
		   efficient simulations of nonlinear synaptic dynamics"
		   (Pham et al., 2021)
		   
Running Instructions:

- Download look-up table files at: https://senselab.med.yale.edu/modeldb/data/267103/LUTs.zip
	- ensure the two folders containing the LUT files are located in your working directoty
- Execute terminal command 'nrnivmodl' within the working directory that contains all files including .mod files
- Execute terminal command 'python LUTsyn_example_main.py' to run example code that compares multiple synapse model outputs
		   
The description for each file included here are as follows:

- LUTsyn_example_main.py	

		This file contains the main python script to run
		a demonstration of the LUTsyn model, as well as the
		kinetic and exponential models and plots a comparison
		between the models. Both AMPA and NMDA receptors are
		simulated and compared. This was tested using Python 3.5
		and NEURON v7.6.7.
		-	NOTE: Please be sure to run 'nrnivmodl' in terminal/command 
			prompt in your working directory that contains the .mod files 
			before running this script.
								
- LUTsyn_example_functions.py
		
		This file contains the supporting functions that are used in
		LUTsyn_example_main.py which take care of executing NEURON
		simulations and instantiating the synapse mechanisms.
		
- AMPA16v8_noNetCon.mod

		.mod file that implements the 16 state kinetic model for AMPA
		
- E3_NMDA_v2.mod

		.mod file that implements NMDA dynamics using a linear triple exponential
		
- exp2syn_v2.mod

		A modified version of NEURON's native exp2syn mechanism to represent AMPA dynamics
		using a linear double exponential function
		
- LUTsyn_AMPA_4th_E3_dtc.mod

		Implementation of the AMPA LUTsyn model.
		
- LUT_AMPA_kinetic_g_4th_M300ms_d1ms_NT0.1_DIFF_v3.npy

		This file represents the look-up table data structure that holds the amplitude values
		for the AMPA LUTsyn model.
		
- LUTsyn_NMDA_5th_E3_dtc.mod

		Implementation of the NMDA LUTsyn model.
		
- LUT_NMDA_kinetic_OSP_5th_M1000ms_d5ms_NT0.1_DIFF_v3.npy

		This file represents the look-up table data structure that holds the amplitude values
		for the NMDA LUTsyn model.
		
- negative_init.hoc

		File that implements an initial simulation that has reached steady-state
		
- NMDA_v6_3_opt.mod

		.mod file that implements the 15 state kinetic model for NMDA
		
- NTDiffusion.mod

		.mod file that implements presynaptic glutamate release (concentration).
		This mechanism runs in series with the kinetic models.
		
- vecstim.mod
		
		.mod file for creating vectors that contain stimulation patterns
		
		
Changelog
---------
2022-12: Updated MOD files to contain valid C++ and be compatible with
the upcoming versions 8.2 and 9.0 of NEURON.		