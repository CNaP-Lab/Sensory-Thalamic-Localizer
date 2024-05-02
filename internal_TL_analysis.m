function obj = internal_TL_analysis(varargin)

    % Inputs:
    %   BOLD_fMRI_unsmoothed_data - Cell array of strings, paths to unsmoothed data files.
    %   runMotionParameters - Cell array of strings, paths to motion parameter files.
    %   FreeSurfer_segmentation_Atlas_path - String, path to the atlas for segmentation.
    %   FreeSurfer_Thalamic_Segmentation - String, path to thalamic segmentation data.
    %   participant_TL_task_log_path - Cell array of path to task log file. 
    %   subjectID - String, identifier for the subject.
    %   TR - Numeric, repetition time.
    %   left_AC_search_region_mask, right_AC_search_region_mask, left_VC_search_region_mask, right_VC_search_region_mask - String, paths to masks.
    %   intermediatesOutputDirectory, fROIoutputDirectory - String, output directories.
    %
    % Optional Name-Value Pairs:
    %   'BOLD_fMRI_Smoothed_data' - Cell array of strings, paths to pre-smoothed data files (Optional).
    %   'smoothingFWHM' - Numeric, smoothing Kernel (Optional).
    
    %
    % Notes:
    %   - If both 'BOLD_fMRI_Smoothed_data' and 'smoothingFWHM' are provided, an error will be thrown.
    %   - If neither 'BOLD_fMRI_Smoothed_data' and 'smoothingFWHM' are provided, an error will be thrown.
    %
        p = inputParser;
    
        % Required parameters
        addRequired(p, 'BOLD_fMRI_unsmoothed_data', @iscell);
        addRequired(p, 'runMotionParameters', @iscell);
        addRequired(p, 'FreeSurfer_segmentation_Atlas_path', @ischar);
        addRequired(p, 'FreeSurfer_Thalamic_Segmentation', @ischar);
        addRequired(p, 'participant_TL_task_log_path', @iscell);
        addRequired(p, 'subjectID', @ischar);
        addRequired(p, 'TR', @isnumeric);
        addRequired(p, 'left_AC_search_region_mask', @ischar);
        addRequired(p, 'right_AC_search_region_mask', @ischar);
        addRequired(p, 'left_VC_search_region_mask', @ischar);
        addRequired(p, 'right_VC_search_region_mask', @ischar);
        addRequired(p, 'intermediatesOutputDirectory', @ischar);
        addRequired(p, 'fROIoutputDirectory', @ischar);
    
        % Optional parameters
        addOptional(p, 'BOLD_fMRI_Smoothed_data', [], @(x) iscell(x));
        addParameter(p, 'smoothingFWHM', -1, @isnumeric);  % Default to -1 if not provided
    
        parse(p, varargin{:});
        BOLD_fMRI_unsmoothed_data = p.Results.BOLD_fMRI_unsmoothed_data;
        runMotionParameters = p.Results.runMotionParameters;
        FreeSurfer_segmentation_Atlas_path = p.Results.FreeSurfer_segmentation_Atlas_path;
        FreeSurfer_Thalamic_Segmentation = p.Results.FreeSurfer_Thalamic_Segmentation;
        participant_TL_task_log_path = p.Results.participant_TL_task_log_path;
        subjectID = p.Results.subjectID;
        TR = p.Results.TR;
        left_AC_search_region_mask = p.Results.left_AC_search_region_mask;
        right_AC_search_region_mask = p.Results.right_AC_search_region_mask;
        left_VC_search_region_mask = p.Results.left_VC_search_region_mask;
        right_VC_search_region_mask = p.Results.right_VC_search_region_mask;
        intermediatesOutputDirectory = p.Results.intermediatesOutputDirectory;
        fROIoutputDirectory = p.Results.fROIoutputDirectory;
        BOLD_fMRI_Smoothed_data = p.Results.BOLD_fMRI_Smoothed_data;
        smoothingFWHM = p.Results.smoothingFWHM;
    
    
        disp('BOLD_fMRI_unsmoothed_data:');
        disp(BOLD_fMRI_unsmoothed_data);
    
        disp('runMotionParameters:');
        disp(runMotionParameters);
    
        disp('FreeSurfer_segmentation_Atlas_path:');
        disp(FreeSurfer_segmentation_Atlas_path);
    
        disp('FreeSurfer_Thalamic_Segmentation:');
        disp(FreeSurfer_Thalamic_Segmentation);
    
        disp('participant_TL_task_log_path:');
        disp(participant_TL_task_log_path);
    
        disp('subjectID:');
        disp(subjectID);
    
        disp('TR:');
        disp(TR);
    
        disp('left_AC_search_region_mask:');
        disp(left_AC_search_region_mask);
    
        disp('right_AC_search_region_mask:');
        disp(right_AC_search_region_mask);
    
        disp('left_VC_search_region_mask:');
        disp(left_VC_search_region_mask);
    
        disp('right_VC_search_region_mask:');
        disp(right_VC_search_region_mask);
    
        disp('intermediatesOutputDirectory:');
        disp(intermediatesOutputDirectory);
    
        disp('fROIoutputDirectory:');
        disp(fROIoutputDirectory);
    
        disp('BOLD_fMRI_Smoothed_data:');
        disp(BOLD_fMRI_Smoothed_data);
    
        disp('smoothingFWHM:');
        disp(smoothingFWHM);
    
        if ~isempty(BOLD_fMRI_Smoothed_data) && smoothingFWHM ~= -1
            error('You cannot provide both pre-smoothed data and a smoothing FWHM. Choose one.');
        elseif isempty(BOLD_fMRI_Smoothed_data) && smoothingFWHM == -1
            error('You must provide either pre-smoothed data or a smoothing FWHM.');
        end
    
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
        obj.MPfile_name_run = runMotionParameters;
    
        if isempty(BOLD_fMRI_Smoothed_data) && smoothingFWHM ~= -1
            disp("Will smoothen and save to obj.stMNI_run")
            obj = TLsmoothVols(obj, smoothingFWHM);
            disp([obj.stMNI_run]);
        elseif ~isempty(BOLD_fMRI_Smoothed_data) && smoothingFWHM == -1
            obj.stMNI_run = BOLD_fMRI_Smoothed_data;
        else
            % Error if neither smoothed data nor smoothing parameter is provided
            % Ideally should not reach here as there are enough checks before
            error('Either smoothed data must be provided or a smoothing FWHM must be specified.');
        end
    
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
    
    
    
        
        