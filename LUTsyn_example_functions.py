# LUTsyn_example_functions.py
# Code by Duy-Tan Jonathan Pham (duytanph@usc.edu)
# July 8, 2021

# This file contains supporting functions for LUTsyn_example_main.py

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


from neuron import h, numpy_element_ref
import numpy as np
import time as clock

###############################################################################################
# genSyn function instantiates a synapse mechanism (written in NMODL for NEURON)
# at a specified dendrite location.
#
# Inputs
# 	dend - NEURON section in which the synapse mechanism will be instantiated at
#	synmodeltype - a string that specifies the type of synapse to instantiate
#
# Outputs
#	syn - the NEURON mechanism that represents the newly created synapse
# 	synweight - the synaptic weight of the synapse (this is only used by the E2_AMPA mechanism
#	NT - Neurotransmitter mechanism, which is ONLY used in series with the kinetic models
###############################################################################################
def genSyn(dend, synmodeltype):

    NT = 0 # only used if using kinetic models
    synweight = 1

    ### AMPA models ###
    # double exponential synapse for AMPA
    if synmodeltype == "E2_AMPA":
        syn = h.Exp2Syn_v2(dend(0.5)) # create synapse

        # set synapse parameters
        syn.tau1 = 0.9
        syn.tau2 = 4.4
        synweight = 2.146142e-04   # scaling factor

    # LUTsyn model for AMPA
    elif synmodeltype == "LUTsyn_AMPA":
        # load look-up table data structure into memory
        gain_array_AMPA_E3 = np.load('./LUT_AMPA_kinetic_g_4th_M300ms_d1ms_NT0.1_DIFF_v3/LUT_AMPA_kinetic_g_4th_M300ms_d1ms_NT0.1_DIFF_v3.npy')
        LUT_ptr = numpy_element_ref(gain_array_AMPA_E3,0) # create pointer to LUT

        syn = h.LUTsyn_AMPA_4th_E3_dtc(dend(0.5)) # create synapse

        # set synapse parameters
        syn.basis_gain = 0.00021258154267359922
        h.setpointer(LUT_ptr, 'gain_array', syn) # link LUT pointer to LUTsyn mechanism (NOTE: this syntax might not
                                                 # work in versions of NEURON later than 7.6.7)
        syn.scalar = 1
        syn.tc1 = 1.05
        syn.tc2 = 3.8
        syn.tc3 = 16.0
        syn.wtc2 = 0.963468100127
        syn.wtc3 = 1 - syn.wtc2
        syn.factor = 2.2071045693250277

        # first order parameters
        syn.o1_tc1 = 1.05
        syn.o1_tc2 = 3.8
        syn.o1_tc3 = 16.0
        syn.factor1 = 2.2071045693250277

        # second order parameters
        syn.o2_tc1 = 1.05
        syn.o2_tc2 = 3.95
        syn.o2_tc3 = 19.0
        syn.factor2 = 2.150761928070988

        # third order parameters
        syn.o3_tc1 = 1.1
        syn.o3_tc2 = 3.85
        syn.o3_tc3 = 18.0
        syn.factor3 = 2.2546969539333945

        # fourth order parameters
        syn.o4_tc1 = 1.05
        syn.o4_tc2 = 4.2
        syn.o4_tc3 = 19.5
        syn.factor4 = 2.0721351225715905

    elif synmodeltype == "Kinetic_AMPA":
        # Note: kinetic model has 2 mechanisms in series with each other: NT diffusion and kinetic state model
        NT = h.NTDiffusion(dend(0.5))    # neurotransmitter mechanism
        NT.Radius = 0.060
        NT.k = 1.32

        syn = h.AMPA16v8_noNC(dend(0.5)) # create synapse
        syn.nbAMPAR = 5
        h.setpointer(NT._ref_NTConcentration, 'Glu', syn)

    ### NMDA models ###
    elif synmodeltype == "E3_NMDA":
        syn = h.E3_NMDA_v2(dend(0.5)) # create synapse

        # set synapse parameters
        syn.tau1 = 20.3877436154
        syn.tau2 = 26.6830234133
        syn.tau3 = 158.729359569
        syn.wtau2 = 0.963468100127
        syn.wtau3 = 1 - syn.wtau2
        syn.factor = 8.670516899305603
        syn.scalar = 2.4594359743378894e-05

    elif synmodeltype == "LUTsyn_NMDA":
        # load look-up table data structure into memory
        gain_array_NMDA_E3 = np.load('./LUT_NMDA_kinetic_OSP_5th_M1000ms_d5ms_NT0.1_DIFF_v3/LUT_NMDA_kinetic_OSP_5th_M1000ms_d5ms_NT0.1_DIFF_v3.npy')
        LUT_ptr = numpy_element_ref(gain_array_NMDA_E3,0) # create pointer to LUT

        syn = h.LUTsyn_NMDA_5th_E3_dtc(dend(0.5)) # create synapse

        # set synapse parameters
        syn.basis_gain = 2.4594359743378894e-05
        h.setpointer(LUT_ptr, 'gain_array', syn) # link LUT pointer to LUTsyn mechanism (NOTE: this syntax might not
                                                 # work in versions of NEURON later than 7.6.7)
        syn.scalar = 1
        syn.gran = 5
        syn.tc1 = 18
        syn.tc2 = 23
        syn.tc3 = 148
        syn.wtc2 = 0.963468100127
        syn.wtc3 = 1 - syn.wtc2
        syn.factor = 9.336247512713125

        # first order parameters
        syn.o1_tc1 = 18
        syn.o1_tc2 = 23
        syn.o1_tc3 = 148
        syn.factor1 = 9.336247512713125

        # second order parameters
        syn.o2_tc1 = 17
        syn.o2_tc2 = 20
        syn.o2_tc3 = 140
        syn.factor2 = 12.845340591656415

        # third order parameters
        syn.o3_tc1 = 18
        syn.o3_tc2 = 21
        syn.o3_tc3 = 144
        syn.factor3 = 13.378037918262585

        # fourth order parameters
        syn.o4_tc1 = 18
        syn.o4_tc2 = 22
        syn.o4_tc3 = 168
        syn.factor4 = 10.878447743899189

        # fifth order parameters
        syn.o5_tc1 = 18
        syn.o5_tc2 = 21
        syn.o5_tc3 = 140
        syn.factor5 = 13.402875584141828

    elif synmodeltype == "Kinetic_NMDA":
        # Note: kinetic model has 2 mechanisms in series with each other: NT diffusion and kinetic state model
        NT = h.NTDiffusion(dend(0.5))    # neurotransmitter mechanism
        NT.Radius = 0.060
        NT.k = 1.32

        syn = h.NMDA_v6_3_opt(dend(0.5)) # create synapse
        syn.nbNMDAR = 1.235
        h.setpointer(NT._ref_NTConcentration, 'Glu', syn)

    return syn, synweight, NT

###############################################################################################
# genStim function returns a NEURON Vector that contains the pulse times of a
# Poisson process with a specified mean firing rate, refractory period, and simulated time
#
# Inputs
# 	freq - specifies the mean input firing rate of stimulation (in Hz)
#	tau_AHP - the time constant of afterhyperpolarization (AHP) measured in seconds
#             which characterizes the refractory period
#	tstop - specifies the total time length of the stimulation
#
# Outputs
# 	event_times - NEURON h.Vector() that contains the pulse stimulation time series
###############################################################################################
def genStim(freq,tau_AHP,tstop):
    # Create Random Number Generator
    ranGen = np.random
    ranGen.seed(0)

    # create pulse stimulus pattern as a poisson process
    event_times = h.Vector()
    if freq > 0:
        tres = 0.00025 # seconds
        refract = -1e9
        time = 0
        while time < tstop:
            roll = ranGen.uniform(0,1)
            refract = 1-np.exp(-(time-refract)/tau_AHP)
            if (roll <= freq*tres*refract):
                event_times.append(time*1000)
                refract = time+tres

            time += tres

    else:
        event_times.append(tstop+10)

    return event_times

###############################################################################################
# nrn_sim function executes a NEURON simulation with specified stimulus parameters, receptor,
# and allows user to specify which synapse model to run
#
# Inputs
# 	stim_params - dictionary containing stimulation parameters such as simulation time, input
#                 firing rate, and refractory period
#	receptor - string that indicates what receptor is being simulated at the synapse
#              "AMPA" or "NMDA"
#	synmodeltype - string that specifies which synapse model to simulate
#   	possible values: "E2_AMPA", "LUTsyn_AMPA", "Kinetic_AMPA", "E3_NMDA", "LUTsyn_NMDA",
#						 "Kinetic_NMDA"
#
# Outputs
# 	outputs - dictionary containing multiple traces of interest such as membrane voltage, EPSC,
#             and receptor conductance
###############################################################################################
def nrn_sim(stim_params,receptor,synmodeltype):

    # Set up NEURON simulation
    h.load_file("stdrun.hoc")
    h.tstop = stim_params['sim_time']
    h.dt = 0.1  # timesteps
    h.steps_per_ms = 1.0/h.dt
    h.load_file('negative_init.hoc')  # allow simulation to reach resting state
    h.v_init = -73
    h.celsius = 35
    h.stdinit()

    # create test NEURON cell ("ball and stick")
        # note that this ball and stick cell is a simplification and is not representative
        # of the cell compartments reported in the LUTsyn paper (Pham, 2021)
    soma = h.Section(name='soma')
    soma.diam = 15   # microns
    dendrite = h.Section(name='dendrite')
    dendrite.diam = 5
    dendrite.L = 3000
    dendrite.connect(soma(1))

    # create synapse mechanism
    syn, synweight,NT = genSyn(dendrite, synmodeltype)

    # create stimulation pattern
    sim_time = stim_params['sim_time'] # ms
    freq = stim_params['freq']         # mean input firing rate (Hz)
    tau_AHP = stim_params['tau_AHP']   # seconds
    tstop = sim_time / 1000.		   # seconds
    event_times = genStim(freq, tau_AHP, tstop) # generate input pulse train using Poisson process
    input_stim = h.VecStim()     # turn pulse train into a NEURON stimulation vector
    input_stim.play(event_times)

    # Create NetCon between synapse and input stimulation
    # (if using kinetic model, then input connects to NTDiffusion Mechanism, not synapse)
    if synmodeltype == 'Kinetic_AMPA' or synmodeltype == 'Kinetic_NMDA':
        nc = h.NetCon(input_stim,NT)
    else:
        nc = h.NetCon(input_stim,syn)
    nc.weight[0] = synweight
    nc.delay = 1

    # Set up NEURON recording vectors
    t = h.Vector()      # time vector
    v = h.Vector()      # membrane voltage
    g = h.Vector()      # receptor conductance
    i = h.Vector()      # receptor mediated EPSC

    t.record(h._ref_t)
    v.record(syn._ref_v1)
    g.record(syn._ref_g)
    i.record(syn._ref_i)

    # if simulating NMDA receptor, record open-state probability too
    if receptor == "NMDA":
        osp = h.Vector()
        if synmodeltype == 'Kinetic_NMDA':
            osp.record(syn._ref_open_total)
        else:
            osp.record(syn._ref_open)

    # Execute NEURON simulation (print out runtimes)
    print('Simulating: ')
    print('            Receptor  : ', receptor)
    print('            Syn model : ', synmodeltype)
    wall_t0 = clock.time()

    h.run()  # RUN NEURON simulation

    wall_t1 = clock.time()
    runtime = wall_t1 - wall_t0
    print('            Run Time  : ', runtime)

    # create a dictionary that will return all useful information about the simulation
    outputs = {}

    # return output values as a dictionary
    outputs['runtime'] = runtime
    outputs['i'] = np.array(i)
    outputs['t'] = np.array(t)
    outputs['g'] = np.array(g)
    outputs['v'] = np.array(v)
    if receptor == 'NMDA':
        outputs['osp'] = np.array(osp)

    return outputs


###############################################################################################
# calc_NRMSE function calculates the Normalized Root Square Error between two time series of
# equal length
#
# Inputs
# 	ref - the reference time series. Both time series get normalized with respect to the maximum
#         value of this time series.
#   test - the time series that is being evaluated and assigned an error value
#
# Outputs
# 	rms_error_val - the NRMSE value between the two time series with ref used as reference
###############################################################################################
def calc_NRMSE(ref, test):
    data_max_val = np.max(np.absolute(ref))
    data_norm_val = np.divide(ref.T, data_max_val)

    estimate_norm_val = test / data_max_val
    squared_difference_val = np.power(np.subtract(data_norm_val, estimate_norm_val), 2)
    sum_data_norm_sq_val = np.sum(np.power(data_norm_val, 2))
    rms_error_val = np.sqrt(sum(squared_difference_val) / sum_data_norm_sq_val)
    return rms_error_val