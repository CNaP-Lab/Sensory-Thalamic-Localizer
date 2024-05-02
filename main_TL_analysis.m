%% Analysis code for the Auditory and Visual Sensory Thalamic Localizer Task
% Released under the GNU Public License version 3.

%% Code authors: 
% John C. Williams, Srineil Nizambad, Philip N. Tubiolo, Yash Patel, 
% and Jared X. Van Snellenberg
% Department of Psychiatry and Behavioral Health
% Department of Biomedcial Engineering
% Renaissance School of Medicine
% Stony Brook University

%% Task Presentation code additionally found here:
% Neurobehavioral Systems Archives of Neurobehavioral 
% Experiments and Stimuli
% http://www.neurobs.com/ex_files/expt_view?id=302

%% Accompanies the following manuscript:
% John C. Williams, Philip N. Tubiolo, Zu Jie Zheng, 
% Eilon B. Silver-Frankel, Dathy T. Pham, Natalka K. Haubold, 
% Sameera K. Abeykoon, Anissa Abi-Dargham, Guillermo Horga, 
% and Jared X. Van Snellenberg. (2024).
% Functional Localization of the Human Auditory and Visual Thalamus Using 
% a Thalamic Localizer Functional Magnetic Resonance Imaging Task.

%% This sets up the path based on the location of the code
code_file_path = mfilename('fullpath');
addpath( fullfile( fileparts( code_file_path ), 'CNaP_Lab_Dependencies' ) );
addpath( fullfile( fileparts( code_file_path ), 'spm12' ) );
addpath(fullfile( fileparts(code_file_path ), 'CanlabCore','Index_image_manip_tools'));
addpath(fullfile( fileparts(code_file_path ), 'CanlabCore','Image_space_tools'));
addpath(fullfile( fileparts(code_file_path ), 'CanlabCore','Data_processing_tools'));

%% This gets the location of the supplied AC and VC search region masks
% Users can supply their own AC and VC search region masks if desired.
supplied_masks_directory = fullfile( fileparts( code_file_path ), 'CNaP_Lab_Dependencies', 'masks' ) ;
left_AC_search_region_mask = fullfile( supplied_masks_directory, 'temporal_left.nii' ) ;
right_AC_search_region_mask = fullfile( supplied_masks_directory, 'temporal_right.nii' ) ;
left_VC_search_region_mask = fullfile( supplied_masks_directory, 'occipital_left.nii' ) ;
right_VC_search_region_mask = fullfile( supplied_masks_directory, 'occipital_right.nii' ) ;

%% Path to the participant's Sensory Thalamic Localizer task log file from Presentation
participant_TL_task_log_path = {'/path/to/data/50023/THL/50023-thal.log'};%example

%% Participant's ID
participantID = '50023';

%% Desired path for the MGN and LGN functionally-defined regions of interest (fROIs)
fROIoutputDirectory = '/desired/path/to/outputs/MGN_LGN_fROIs'; %MGN and LGN fROI output directory

%% Desired path for intermediates used, such as outputs from first-level modeling and contrasts
intermediatesOutputDirectory = '/desired/path/to/outputs/intermediates'; %Output directory for intermediates

%% These are the unsmoothed BOLD fMRI data for each run, in a cell array.
BOLD_fMRI_unsmoothed_data{1} = '/path/to/data/50023/THL/tTHL_fMRI_1.nii'; %example
BOLD_fMRI_unsmoothed_data{2} = '/path/to/data/50023/THL/tTHL_fMRI_2.nii'; %example
BOLD_fMRI_unsmoothed_data{3} = '/path/to/data/50023/THL/tTHL_fMRI_3.nii'; %example
BOLD_fMRI_unsmoothed_data{4} = '/path/to/data/50023/THL/tTHL_fMRI_4.nii'; %example

%%OPTIONAL
%% These are the Smoothed BOLD fMRI data for each run, in a cell array.
% This is an optional input if you would like to provide your own smoothed time series volumes
% as opposed to computed by the pipeline. However, note that the results using
% unsmoothed data were not optimal for the tested dataset relative to when
% a 4mm FWHM Gaussian smoothing kernel was used.
% Thus, some testing and comparison using the dataset in question is
% important if a user desires to attempt using data that is unsmoothed, or
% smoothed with a lesser smoothing kernel.
BOLD_fMRI_Smoothed_data{1} = '/path/to/data/50023/THL/stTHL_fMRI_1.nii'; %example
BOLD_fMRI_Smoothed_data{2} = '/path/to/data/50023/THL/stTHL_fMRI_2.nii'; %example
BOLD_fMRI_Smoothed_data{3} = '/path/to/data/50023/THL/stTHL_fMRI_3.nii'; %example
BOLD_fMRI_Smoothed_data{4} = '/path/to/data/50023/THL/stTHL_fMRI_4.nii'; %example

%% These are the motion parameters for each run, in a cell array.
runMotionParameters{1} = '/path/to/data/50023/THL/trimmedMovement_Regressors1.txt'; %example
runMotionParameters{2} = '/path/to/data/50023/THL/trimmedMovement_Regressors2.txt'; %example
runMotionParameters{3} = '/path/to/data/50023/THL/trimmedMovement_Regressors3.txt'; %example
runMotionParameters{4} = '/path/to/data/50023/THL/trimmedMovement_Regressors4.txt'; %example

%% Repetition time for the BOLD data, (TR) in seconds
% used by spm.stats.fmri_spec.timing.TR
TR = 0.8;  % seconds

%% Path of the FreeSurfer segmentation/parcellation atlas (Desikan-Killiany) in MNI space, Atlas.wmparc.2.nii
FreeSurfer_segmentation_Atlas_path = '/path/to/data/50023/Atlas_wmparc.2.nii'; % Used to make a gray matter mask for 1st level modeling

%% Participant's FreeSurfer thalamic segmentation, warped to MNI space.
%See: https://freesurfer.net/fswiki/ThalamicNuclei
% Iglesias JE, Insausti R, Lerma-Usabiaga G, Bocchetta M, Van Leemput K,
% Greve DN, van der Kouwe A; Alzheimer's Disease Neuroimaging Initiative;
% Fischl B, Caballero-Gaudes C, Paz-Alonso PM.
% A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI
% and histology. Neuroimage. 2018 Dec;183:314-326.
% doi: 10.1016/j.neuroimage.2018.08.012. Epub 2018 Aug 17.
% PMID: 30121337; PMCID: PMC6215335.
FreeSurfer_Thalamic_Segmentation = '/path/to/data/50023/50023_ThalamicNuclei.v13.T1_WARPED.nii'; %example

% OPTIONAL
% This is an optional input needed to smoothen the unsmoothed time series 
% if you do not intend to use user fed smoothed volumes 
smoothingFWHM = 4; %example

% Notes:
%   - If both ‘BOLD_fMRI_Smoothed_data’ and ‘smoothingFWHM’ are provided, an error will be thrown.
%   - If neither ‘BOLD_fMRI_Smoothed_data’ and ‘smoothingFWHM’ are provided, an error will be thrown.


%% Start the analysis
% considering user fed smoothed timeseries volumes
% optional input 'BOLD_fMRI_Smoothed_data' and 'smoothingFWHM' provided as name value pair
[obj] = internal_TL_analysis( ...
    BOLD_fMRI_unsmoothed_data, ...
    runMotionParameters, ...
    FreeSurfer_segmentation_Atlas_path, ...
    FreeSurfer_Thalamic_Segmentation, ...
    participant_TL_task_log_path, ...
    participantID, ...
    TR, ...
    left_AC_search_region_mask, ...
    right_AC_search_region_mask, ...
    left_VC_search_region_mask, ...
    right_VC_search_region_mask, ...
    intermediatesOutputDirectory, ...
    fROIoutputDirectory, ...
    'BOLD_fMRI_Smoothed_data', BOLD_fMRI_Smoothed_data, ... % either this or smoothingFWHM
    'SmoothingFWHM', smoothingFWHM ... % either this or BOLD_fMRI_Smoothed_data
    );
