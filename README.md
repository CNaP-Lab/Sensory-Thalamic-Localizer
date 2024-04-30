# Sensory-Thalamic-Localizer
Analysis code for the Auditory and Visual Sensory Thalamic Localizer Task<br>
Released under the GNU Public License version 3.

Code authors: <br>
John C. Williams, Srineil Nizambad, Philip N. Tubiolo, Yash Patel, and Jared X. Van Snellenberg<br>
Department of Psychiatry and Behavioral Health<br>
Department of Biomedcial Engineering<br>
Renaissance School of Medicine<br>
Stony Brook University

Task Presentation code additionally found here:<br>
Neurobehavioral Systems Archives of Neurobehavioral Experiments and Stimuli<br>
http://www.neurobs.com/ex_files/expt_view?id=302

Accompanies the following manuscript:<br>
John C. Williams, Philip N. Tubiolo, Zu Jie Zheng, Eilon B. Silver-Frankel, Dathy T. Pham, Natalka K. Haubold, Sameera K. Abeykoon, Anissa Abi-Dargham, Guillermo Horga, and Jared X. Van Snellenberg. (2024).<br>
Functional Localization of the Human Auditory and Visual Thalamus Using a Thalamic Localizer Functional Magnetic Resonance Imaging Task.<br>
bioRxiv 2024.04.28.591516; doi: https://doi.org/10.1101/2024.04.28.591516<br>
https://www.biorxiv.org/content/10.1101/2024.04.28.591516

To use this code, see the script: main_TL_analysis.m.<br>
It sets up the required variables as described in the comments, and calls the main function that performs the analysis: internal_TL_analysis.m.

Requirements:<br>
1. Path to the participant's Sensory Thalamic Localizer task log file from Presentation
2. Participant's ID number
3. Desired path for the MGN and LGN functionally-defined regions of interest (fROIs)
4. Desired path for intermediates used, such as outputs from first-level modeling and contrasts
5. Paths to each run of unsmoothed BOLD fMRI data, supplied in a cell array.
6. Paths to each run of Smoothed BOLD fMRI data, supplied in a cell array.
7. Paths to the motion parameters for each run, supplied in a cell array.
8. Repetition time for the BOLD data, (TR) in seconds
9. Path of the FreeSurfer segmentation/parcellation atlas (Desikan-Killiany) in MNI space, Atlas.wmparc.2.nii
10. Participant's FreeSurfer thalamic segmentation, warped to MNI space.<br>
    See: https://freesurfer.net/fswiki/ThalamicNuclei<br>
    Iglesias JE, Insausti R, Lerma-Usabiaga G, Bocchetta M, Van Leemput K, Greve DN, van der Kouwe A; Alzheimer's Disease Neuroimaging Initiative; Fischl B, Caballero-Gaudes C, Paz-Alonso PM.<br>
    A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology. Neuroimage. 2018 Dec;183:314-326.<br>
    doi: 10.1016/j.neuroimage.2018.08.012. Epub 2018 Aug 17. PMID: 30121337; PMCID: PMC6215335.
