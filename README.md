# Sensory-Thalamic-Localizer
Analysis code for the Auditory and Visual Sensory Thalamic Localizer Task
Released under the GNU Public License version 3.

Code authors: 
John C. Williams, Srineil Nizambad, Philip N. Tubiolo, Yash Patel, and Jared X. Van Snellenberg
Department of Psychiatry and Behavioral Health
Department of Biomedcial Engineering
Renaissance School of Medicine
Stony Brook University

Task Presentation code additionally found here:
Neurobehavioral Systems Archives of Neurobehavioral Experiments and Stimuli
http://www.neurobs.com/ex_files/expt_view?id=302

Accompanies the following manuscript:
John C. Williams, Philip N. Tubiolo, Zu Jie Zheng, Eilon B. Silver-Frankel, Dathy T. Pham, Natalka K. Haubold, Sameera K. Abeykoon, Anissa Abi-Dargham, Guillermo Horga, and Jared X. Van Snellenberg. (2024).
Functional Localization of the Human Auditory and Visual Thalamus Using a Thalamic Localizer Functional Magnetic Resonance Imaging Task.

To use this code, see the script: main_TL_analysis.m
It sets up the required variables as described in the comments, and calls the main function that performs the analysis: internal_TL_analysis.m.

1. Path to the participant's Sensory Thalamic Localizer task log file from Presentation
2. Participant's ID number
3. Desired path for the MGN and LGN functionally-defined regions of interest (fROIs)
4. Desired path for intermediates used, such as outputs from first-level modeling and contrasts
5. Paths to each run of unsmoothed BOLD fMRI data, supplied in a cell array.
6. Paths to each run of Smoothed BOLD fMRI data, supplied in a cell array.
7. Paths to the motion parameters for each run, supplied in a cell array.
8. Repetition time for the BOLD data, (TR) in seconds
9. Path of the FreeSurfer segmentation/parcellation atlas (Desikan-Killiany) in MNI space, Atlas.wmparc.2.nii
10. Participant's FreeSurfer thalamic segmentation, warped to MNI space.
    See: https://freesurfer.net/fswiki/ThalamicNuclei
    Iglesias JE, Insausti R, Lerma-Usabiaga G, Bocchetta M, Van Leemput K,
    Greve DN, van der Kouwe A; Alzheimer's Disease Neuroimaging Initiative;
    Fischl B, Caballero-Gaudes C, Paz-Alonso PM.
    A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI
    and histology. Neuroimage. 2018 Dec;183:314-326.
    doi: 10.1016/j.neuroimage.2018.08.012. Epub 2018 Aug 17.
    PMID: 30121337; PMCID: PMC6215335.
