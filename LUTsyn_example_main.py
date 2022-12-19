# LUTsyn_example_main.py
# Code by Duy-Tan Jonathan Pham (duytanph@usc.edu)
# July 8, 2021

# This file is a Python script for demonstrating the LUTsyn synapse models described in
# (Pham, 2021). This file executes multiple NEURON simulations that compare the LUTsyn
# model to other synapse models and plots the results. Both AMPA and NMDA cases are
# simulated. This script has been tested using Python (v3.5) and NEURON (v7.6.7).
#
# NOTE: Please be sure to run 'nrnivmodl' in terminal/command prompt in your working
# directory that contains the .mod files before running this script.

######################################################################################
# This software is Copyright © 2021 The University of Southern
# California. All Rights Reserved.
#
# Permission to use, copy, modify, and distribute this software
# and its documentation for educational, research and non-profit
# purposes, without fee, and without a written agreement is
# hereby granted, provided that the above copyright notice, this
# paragraph and the following three paragraphs appear in all copies.
#
# Permission to make commercial use of this software may be obtained by contacting:
# USC Stevens Center for Innovation
# University of Southern California
# 1150 S. Olive Street, Suite 2300
# Los Angeles, CA 90115, USA

# This software program and documentation are copyrighted by The
# University of Southern California. The software program and
# documentation are supplied "as is", without any accompanying
# services from USC. USC does not warrant that the operation of the
# program will be uninterrupted or error-free. The end-user understands
# that the program was developed for research purposes and is advised
# not to rely exclusively on the program for any reason.

# IN NO EVENT SHALL THE UNIVERSITY OF SOUTHERN CALIFORNIA BE
# LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT
# OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
# UNIVERSITY OF SOUTHERN CALIFORNIA HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF SOUTHERN CALIFORNIA
# SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
# BASIS, AND THE UNIVERSITY OF SOUTHERN CALIFORNIA HAS NO OBLIGATIONS
# TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
######################################################################################


from LUTsyn_example_functions import *
import matplotlib.pyplot as plt
############################################# MAIN #############################################

# possible synapse models: 'LUTsyn_AMPA', 'E2_AMPA', 'Kinetic_AMPA'
#                          'LUTsyn_NMDA', 'E3_NMDA', 'Kinetic_NMDA'

# set stimulation parameters
freq = 10
stim_params = {}
stim_params['sim_time'] = 20000 # 20 second simulation
stim_params['freq'] = freq      # 10 Hz mean input firing rate
stim_params['tau_AHP'] = 0.035  # Poisson process parameter (seconds)

# run AMPA simulations
sim_LUTsyn_AMPA = nrn_sim(stim_params,'AMPA','LUTsyn_AMPA')  # LUTsyn AMPA
sim_E2 = nrn_sim(stim_params,'AMPA','E2_AMPA')               # double exponential AMPA
sim_kinetic_AMPA = nrn_sim(stim_params,'AMPA','Kinetic_AMPA')# Kinetic AMPA

# run NMDA simulations
sim_LUTsyn_NMDA = nrn_sim(stim_params,'NMDA','LUTsyn_NMDA')  # LUTsyn NMDA
sim_E3 = nrn_sim(stim_params,'NMDA','E3_NMDA')               # triple exponential NMDA
sim_kinetic_NMDA = nrn_sim(stim_params,'NMDA','Kinetic_NMDA')# Kinetic NMDA

# Calculate NRMSE values (round to 4 decimal spots)
NRMSE_LUTsyn_AMPA = round(calc_NRMSE(sim_kinetic_AMPA['g'], sim_LUTsyn_AMPA['g'] ), 4)
NRMSE_E2_AMPA = round(calc_NRMSE(sim_kinetic_AMPA['g'], sim_E2['g'] ), 4)

NRMSE_LUTsyn_NMDA = round(calc_NRMSE(sim_kinetic_NMDA['osp'], sim_LUTsyn_NMDA['osp'] ), 4)
NRMSE_E3_NMDA = round(calc_NRMSE(sim_kinetic_NMDA['osp'], sim_E3['osp'] ), 4)

#### PLOTTING ####
t = sim_LUTsyn_AMPA['t']

# AMPA plot comparison
plt.figure(1)
# multiply by 1000 to convert from nS to pS
plt.plot(t, sim_kinetic_AMPA['g']*1e3, 'b-', label = 'Kinetic AMPA', alpha = 0.6)
plt.plot(t, sim_LUTsyn_AMPA['g']*1e3, 'g--', label = 'LUTsyn_AMPA, NRMSE: ' + str(NRMSE_LUTsyn_AMPA), alpha = 0.6)
plt.plot(t, sim_E2['g']*1e3, 'r:', label = 'Double Exponential, NRMSE: ' + str(NRMSE_E2_AMPA), alpha = 0.6)
plt.title('AMPA Comparison at ' + str(freq) + ' Hz')
plt.ylabel('Conductance (pS)')
plt.xlabel('Time (ms)')
plt.legend()

# NMDA plot comparison
plt.figure(2)
plt.plot(t, sim_kinetic_NMDA['osp'], 'b-', label = 'Kinetic NMDA', alpha = 0.6)
plt.plot(t, sim_LUTsyn_NMDA['osp'], 'g--', label = 'LUTsyn_NMDA, NRMSE: ' + str(NRMSE_LUTsyn_NMDA), alpha = 0.6)
plt.plot(t, sim_E3['osp'], 'r:', label = 'Triple Exponential, NRMSE: ' + str(NRMSE_E3_NMDA), alpha = 0.6)
plt.title('NMDA Comparison at ' + str(freq) + ' Hz')
plt.ylabel('Open-state probability')
plt.xlabel('Time (ms)')
plt.legend()

plt.show()