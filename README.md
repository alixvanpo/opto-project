# Optogenetic Project
FitzHugh-Nagumo and Aliev-Panfilov models for single-cell and cardiac monolayer simulations with optogenetic control
## Table of Contents
* [Introduction](#introduction)
* [Single-Cell](#single-cell)
* [Grid](#grid)
* [Optogenetic Single-Cell](#optogenetic-single-cell)
* [Project Status](#project-status)

## Introduction
This optogenetic project runs in MATLAB (R2020b), and aims to simulate excitation patterns seen in single cardiomyocyte and 2D cardiac monolayer experiments. The two main models used were the FitzHugh-Nagumo and [Aliev-Panfilov](https://ieeexplore.ieee.org/document/9030412) from Iranian computational neuroscientist Zafarghandi and others.

The single-cell models were manipulated to include optogenetic control via channel rhodopsin current generation. It is described in [this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4407252/?report=reader#!po=4.05405) paper by colleagues Williams and Entcheva, including a Markov-style four-state model
for the channel.

Most codes follow the same structure:
* Parameters
* Integration Setup
* Steady State
* Initial Conditions
* Forward Euler
* Plotting

The parameters and differential equations were taken from papers, and most include comments describing them. Next is an integration setup, preparing for the forward euler integration near the end of the code. Then, the variables of interest are set to their steady state values in the absence of a stimulus. From these values, the initial conditions were set up for all variables that will be computed in the euler method (some of these are set to zero to indicate their dimensions). Finally, the integration for-loop appears along with a plotting code to create output (usually voltage against time). Most of the codes incorporate the stimulus here as a step function `Iapp`, turned *on* for a given time.

Below, you will find special details and instructions for each file.


## Single-Cell
**fhn_sc** is a simple FHN single-cell model. A 0.5 stimulus is applied for the first 2 seconds of integration.

**aliev_sc** is the AP single-cell model. The same stimulus as above is applied. Note that the plotting section has commented out a way to dimensionalise the system, which is an arbitrary unit conversion. This manipulation will resurface in various files to compare to experimental time series.

**steadyState_sc** is a helper function called in the main single-cell files, in order to run the model without any stimulus. The variables must be updated in each iteration, and should reach steady state values within 1000 time units.

## Grid
**fhn_grid** is the FHN 2D grid model of 100x100 cells. The output here will be a heat map of voltage, where each frame at each time step can be saved and turned into a video. Included in the code is various of ways of stimulating the grid of cells. The *position* section can zoom in on a cell in the grid and plot its time series. Here, `t_step` serves the purpose of chosing how often to display and save a frame. 

You can decide what actions your code will take in the *output* section by making variables equal to 1. `save_im` and `save_mov` will save the output in the same location as the matlab file, in the proper folder with the name `baseFileName` or `movieFileName`. `plot_ts` will display the time series of the cell you selected in *position*.

Types of stimuli:
* Left-bar (from `iapp_width` and `iapp_height`)
* Disk
* Import from csv (in same location as matlab file)
* Random grid of 0s and 1s

Decide in *forward euler* section how to apply your stimulus. Refer to *eulerStep* below to find how the variables are updated in time.

**ap_grid** is the AP 2D grid model of 100x100 cells. The only difference between this file and the FHN one described above is the naming convention. Here, the variable u describes the voltage and v the recovery. This was done to keep the variables consistent with what is found in the literature.

**steadyState** is a helper code that computes the steady state of the model in the absence of a stimulus. You may chose how many steps to take to reach a steady state, but it was found that 1000 time units was sufficient. This code can handle both the FHN and AP models.

**eulerStep** is a helper code that computes the laplacian (for the diffusion terms) and uses the result to compute the next integration in the forward euler computation. The variables are stored in *pages* or *sheets* of the alias variable S. This code can handle both the FHN and AP models.

**del2_noflux** is a helper code that overrides the built-in del2 function in MATLAB. The modifications here were made to handle the no-flux boundary conditions necessary for the models.

## Optogenetic Single-Cell
**fhn_sc_opto** is the single-cell FHN model with optogenetic stimulation. The formalism for the channel rhodopsin current was obtained from the paper mentioned in *Introduction*. Note that the Markov-style four-state model is self-contained, and the light-sensitive parameters were simplified to constants (refer to paper and commented section to get full model). The initial conditions were found to be 100% in the C1 state, and progress through the states to allow current across the membrane. 

**ap_sc_opto** is the single-cell AP model with optogenetic stimulation. The code is relatively the same as the `fhn_sc_opto` code.

**fhn_comp** is a helper function built to compare optogenetic versus electrical stimulation of the FHN single-cell optogenetic model. The two parameters controlling these events are:
* conductance of the ChR2 `g_chr`
* applied current strength `iapp_stim`

**ap_comp** is a helper function built to compare optogenetic versus electrical stimulation of the AP single cell optogenetic model. The same parameters as above were chosen to compare the events.

**plotfigs** has two sections to create different, meaningful outputs. The first makes use of the `fhn_comp` and `ap_comp` files to compare both stimulations on the same time series. You may change the input values of the functions to investigate the effects of various conductances and applied stimulus strengths. The second plots the evolution of the channels' states. This is where you can see that all channels start in the C1 state, and progress into different states as time passes.

## Project Status
The last planned changes were about the voltage-dependence in the channel rhodopsin current formalism. The term is known to be dimensionless, but the input values were being investigated. The next phase of the project is to move forward with the optogenetic implementation, and build a FHN and AP 2D grid model with light-sensitive stimulation.

