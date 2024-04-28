restoredefaultpath;
clear all;
code_file_path = mfilename('fullpath');
addpath( fullfile( fileparts( code_file_path ), 'CNaP_Lab_Dependencies' ) );
addpath( fullfile( fileparts( code_file_path ), 'spm12' ) );
addpath(fullfile( fileparts(code_file_path ), 'CanlabCore','Index_image_manip_tools'));
addpath(fullfile( fileparts(code_file_path ), 'CanlabCore','Image_space_tools'));
addpath(fullfile( fileparts(code_file_path ), 'CanlabCore','Data_processing_tools'));

supplied_masks_directory = fullfile( fileparts( code_file_path ), 'CNaP_Lab_Dependencies', 'masks' ) ;
left_AC_search_region_mask = fullfile( supplied_masks_directory, 'temporal_left.nii' ) ;
right_AC_search_region_mask = fullfile( supplied_masks_directory, 'temporal_right.nii' ) ;
left_VC_search_region_mask = fullfile( supplied_masks_directory, 'occipital_left.nii' ) ;
right_VC_search_region_mask = fullfile( supplied_masks_directory, 'occipital_right.nii' ) ;

participant_TL_task_log_path = {'/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/rdoc/50023/THL/50023-thal.log'};%example

subjectID = '50023';

fROIoutputDirectory = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/OneSubject_Home/50023/dataset/THL/MGN_LGN_fROIs'; %MGN and LGN ROI output directory
intermediatesOutputDirectory = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/OneSubject_Home/50023/dataset/THL/intermediates'; %Output directory for intermediates

BOLD_fMRI_unsmoothed_data{1} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/tTHL_fMRI_1.nii'; %example
BOLD_fMRI_unsmoothed_data{2} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/tTHL_fMRI_2.nii'; %example
BOLD_fMRI_unsmoothed_data{3} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/tTHL_fMRI_3.nii'; %example
BOLD_fMRI_unsmoothed_data{4} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/tTHL_fMRI_4.nii'; %example

BOLD_fMRI_Smoothed_data{1} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/stTHL_fMRI_1.nii'; %example
BOLD_fMRI_Smoothed_data{2} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/stTHL_fMRI_2.nii'; %example
BOLD_fMRI_Smoothed_data{3} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/stTHL_fMRI_3.nii'; %example
BOLD_fMRI_Smoothed_data{4} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/stTHL_fMRI_4.nii'; %example

runMotionParameters{1} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/trimmedMovement_Regressors1.txt'; %example
runMotionParameters{2} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/trimmedMovement_Regressors2.txt'; %example
runMotionParameters{3} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/trimmedMovement_Regressors3.txt'; %example
runMotionParameters{4} = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/neil/trimmedMovement_Regressors4.txt'; %example

% provide the repetition time (TR) in seconds
% used by spm.stats.fmri_spec.timing.TR
TR = 0.8;  % seconds

% provide the path of the FreeSurfer segmentationatlas in MNI space,
% containing MNI parcels, such as 'Atlas.wmparc.2.nii'

FreeSurfer_segmentation_Atlas_path = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/OneSubject_Home/50023/Atlas_wmparc.2.nii'; % Used to make a gray matter mask for 1st level modeling

%Participant's FreeSurfer thalamic segmentation, warped to MNI space.
%See: https://freesurfer.net/fswiki/ThalamicNuclei
% Iglesias JE, Insausti R, Lerma-Usabiaga G, Bocchetta M, Van Leemput K,
% Greve DN, van der Kouwe A; Alzheimer's Disease Neuroimaging Initiative;
% Fischl B, Caballero-Gaudes C, Paz-Alonso PM.
% A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI
% and histology. Neuroimage. 2018 Dec;183:314-326.
% doi: 10.1016/j.neuroimage.2018.08.012. Epub 2018 Aug 17.
% PMID: 30121337; PMCID: PMC6215335.
FreeSurfer_Thalamic_Segmentation = '/mnt/jxvs2_01/Thal_Loc_Data/RDoC_Analysis/QCplots_HCP_4.2_for_final/50023/MNINonLinear/50023_ThalamicNuclei.v13.T1_WARPED.nii'; %example

[obj] = internal_TL_analysis( ...
    BOLD_fMRI_unsmoothed_data, ...
    BOLD_fMRI_Smoothed_data, ...
    runMotionParameters, ...
    FreeSurfer_segmentation_Atlas_path, ...
    FreeSurfer_Thalamic_Segmentation, ...
    participant_TL_task_log_path, ...
    subjectID, ...
    TR, ...
    left_AC_search_region_mask, ...
    right_AC_search_region_mask, ...
    left_VC_search_region_mask, ...
    right_VC_search_region_mask, ...
    intermediatesOutputDirectory, ...
    fROIoutputDirectory ...
    );


function [obj] = internal_TL_analysis(...
        BOLD_fMRI_unsmoothed_data, ...
        BOLD_fMRI_Smoothed_data, ...
        runMotionParameters, ...
        FreeSurfer_segmentation_Atlas_path, ...
        FreeSurfer_Thalamic_Segmentation, ...
        participant_TL_task_log_path, ...
        subjectID, ...
        TR, ...
        left_AC_search_region_mask, ...
        right_AC_search_region_mask, ...
        left_VC_search_region_mask, ...
        right_VC_search_region_mask, ...
        intermediatesOutputDirectory, ...
        fROIoutputDirectory ...
        )

    if (~exist(intermediatesOutputDirectory,'dir'))
        mkdir(intermediatesOutputDirectory);
    end
    if (~exist(fROIoutputDirectory,'dir'))
        mkdir(fROIoutputDirectory);
    end

    obj.spmuse.TR = TR;
    obj.home = intermediatesOutputDirectory;
    obj.par.par.name = subjectID;
    obj.par.par.home = intermediatesOutputDirectory;
    obj.td.logs.source = participant_TL_task_log_path;

    obj.tMNI_run = BOLD_fMRI_unsmoothed_data;
    obj.stMNI_run = BOLD_fMRI_Smoothed_data;
    obj.MPfile_name_run = runMotionParameters;

    [Y_grayMatter_mask,V_grayMatter_mask, ...
        Y_whiteMatter_mask,Y_CSF_mask, ...
        grayMatterMask_eroded,whiteMatterMask_eroded,csfMask_eroded] = ...
        TLmakeMasks(FreeSurfer_segmentation_Atlas_path);

    % run all the main 6 functions of the pipeline.
    obj = TLreadLogs(obj , participant_TL_task_log_path);

    obj = TLprep1l(obj,Y_grayMatter_mask,V_grayMatter_mask);

    obj = TLrun1l(obj);

    obj = TLmakeSearchSpace(obj,FreeSurfer_Thalamic_Segmentation,FreeSurfer_segmentation_Atlas_path);

    obj = TLgetROIs(obj, ...
        left_AC_search_region_mask,right_AC_search_region_mask, ...
        left_VC_search_region_mask,right_VC_search_region_mask, ...
        Y_grayMatter_mask,Y_whiteMatter_mask,Y_CSF_mask,...
        grayMatterMask_eroded,whiteMatterMask_eroded,csfMask_eroded, ...
        intermediatesOutputDirectory, ...
        fROIoutputDirectory);

end
