
This repository contains code associated with Oude Lohuis et al., 2022 “Multisensory task demands temporally extend the causal requirement for visual cortex in perception” published in Nature Communications. 

Code Files:
The code is organized in directories per main figure. Code associated with producing data presented in the supplementary information is found within the folder of the main figure that the supplementary figure belongs to. All files are named by their corresponding Figures in the paper:

	•	Behavior
	•	Neural
	•	GLM
	•	AUC
	•	Optogenetics
	•	VisuotactileTask
	•	Noise correlations

Within each directory the preprocessed data necessary for the reproduction of the figures is present:

	•	Data1_1: three example behavioral sessions
	•	Data1_2: all behavioral sessions
	•	Data2_1: recording sessions with V1 neurons sampled from NE, UST, MST mice (excluding optogenetic inactivation) 
	•	Data3_1: GLM results 
	•	Data4_1: AUC single neuron results 
	•	Data4_2: AUC population activity results 
	•	Data5_1: Raw LFP during photostimulation
	•	Data5_2: Firing rates during photostimulation
	•	Data5_3: Behavioral data from optogenetic inactivation sessions
	•	Data6_1: Behavioral data from visuotactile task. Raw data and model fits.
	•	Data6_2: Neuronal data from visuotactile task. Spike sorted data.
	•	Data6_3: Behavioral data from optogenetic inactivation sessions for visuotactile task. Raw data.
	•	Data7_1 V1 activity during photostimulation 

There are two additional directories: 

Utils: set of helper functions

CreateData scripts: original scripts that load the data, filter out the relevant conditions and save the data accordingly. To execute these scripts the total preprocessed dataset is necessary (which is not included). These scripts merely serve to illustrate how the datasets were created.

This repository was tested on MATLAB R2016a, MATLAB R2020a with the following toolboxes:
	•	Optimization toolbox
	•	Curve fitting toolbox
	•	Statistics and machine learning toolbox
	•	Signal processing toolbox
	•	Image processing toolbox
	•	GUI layout toolbox

How to cite this material: 
    Companion code to Oude Lohuis et al., 2022 “Multisensory task demands temporally extend the causal requirement for visual cortex in perception”, doi:10.5281/zenodo.6451623, available at Https://gitlab.com/csnlab/olcese-lab/modid-project/2nd-bump
