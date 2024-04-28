function [obj] = TLmakeSearchSpace(obj,thalsegfile,maskfilein)
    
    % Open TL image
    leftMGNnum = 8115;
    leftLGNnum = 8109;
    rightMGNnum = 8215;
    rightLGNnum = 8209;

    fn = thalsegfile;
    [thalamusHeader,thalamusData] = tippVol(fn);

    leftMGN = thalamusData == leftMGNnum;
    leftLGN = thalamusData == leftLGNnum;
    rightMGN = thalamusData == rightMGNnum;
    rightLGN = thalamusData == rightLGNnum;

    numToDilateGeniculi = 2;
    for i = 1:numToDilateGeniculi
        leftMGN = dilate3d(logical(leftMGN));
        leftLGN = dilate3d(logical(leftLGN));
        rightMGN = dilate3d(logical(rightMGN));
        rightLGN = dilate3d(logical(rightLGN));
    end

    numToDilateMGNalone = 1;
    for i = 1:numToDilateMGNalone
        leftMGN = dilate3d(logical(leftMGN));
        rightMGN = dilate3d(logical(rightMGN));
    end


    %Get pulvinar and mediodorsal from Iglesias et al. 2018 TL segmentation:
    % "The main image feature when subdiving the thalamus is the boundary
    % between the mediodorsal and pulvinar nuclei, and all other nuclei.
    % This boundary is faint, but visible in most T1 scans.
    % Other than this internal boundary and the external edge of the thalamus,
    % the segmentation needs to rely on prior knowledge encoded in the atlas."
    leftPulvinar = (thalamusData >= 8120) & (thalamusData <= 8123);
    rightPulvinar = (thalamusData >= 8220) & (thalamusData <= 8223);
    bilateralMediodorsal = ((thalamusData >= 8112) & (thalamusData <= 8113)) | ((thalamusData >= 8212) & (thalamusData <= 8213));
    bilateralMediodorsal = dilate3d(bilateralMediodorsal); %Dilate 3D 1 voxel
    % Get regular Freesurfer structures
    maskOutFreeSurfer = { ...
        1000, ... %    ctx-lh-unknown
        2000, ... %    ctx-rh-unknown
        1016, ... %    ctx-lh-parahippocampal
        2016, ... %    ctx-rh-parahippocampal
        3016, ... %    wm-lh-parahippocampal
        4016, ... %    wm-rh-parahippocampal
        };
    maskOutDilatePutamen = { ...
        12, ... %  Left-Putamen
        51, ... %  Right-Putamen
        };
    maskOutDilatePallidum = { ...
        13, ... %  Left-Pallidum
        52, ... %  Right-Pallidum
        };
    maskOutDilateHippocampus = { ...
        17, ... %  Left-Hippocampus
        53, ... %  Right-Hippocampus
        };
    maskOutDilateChoroid = { ...
        31, ... %  Left-choroid-plexus
        63, ... %  Right-choroid-plexus
        };

    maskOutInsula = { ...
        1035, ... % ctx-lh-insula
        2035, ... % ctx-rh-insula
        3035, ... % wm-lh-insula
        4035, ... % wm-rh-insula
        };
    
    
    [maskOutCortex,~] = TLmask(maskfilein,'custom',[1000 3000]);
    [freesurferMaskOut,~] = TLmask(maskfilein,'custom',maskOutFreeSurfer);
    [putamenMask,~] = TLmask(maskfilein,'custom',maskOutDilatePutamen);
    [pallidumMask,~] = TLmask(maskfilein,'custom',maskOutDilatePallidum);
    [hippocampusMask,~] = TLmask(maskfilein,'custom',maskOutDilateHippocampus);
    [choroidMask,~] = TLmask(maskfilein,'custom',maskOutDilateChoroid);
    [insulaMask,~] = TLmask(maskfilein,'custom',maskOutInsula);


    freesurferMaskOut = logical(freesurferMaskOut);

    putamenMask = dilate3d(logical(putamenMask));
    pallidumMask = dilate3d(logical(pallidumMask));
    choroidMask = dilate3d(logical(choroidMask));
    insulaMask = dilate3d(logical(insulaMask));
    subcorticalMaskOut = putamenMask | pallidumMask | hippocampusMask | choroidMask | insulaMask;


    bilateralPulvinar = leftPulvinar | rightPulvinar;

    % Mask
    bilateralMaskOut = bilateralPulvinar | bilateralMediodorsal;

    leftLGN = leftLGN & ~bilateralMaskOut;
    rightLGN = rightLGN & ~bilateralMaskOut;
    leftMGN = leftMGN & ~bilateralMaskOut;
    rightMGN = rightMGN & ~bilateralMaskOut;

    %Now mask other structures
    leftLGN = leftLGN & ~freesurferMaskOut & ~maskOutCortex & ~subcorticalMaskOut;
    rightLGN = rightLGN & ~freesurferMaskOut & ~maskOutCortex & ~subcorticalMaskOut;
    leftMGN = leftMGN & ~freesurferMaskOut & ~maskOutCortex & ~subcorticalMaskOut;
    rightMGN = rightMGN & ~freesurferMaskOut & ~maskOutCortex & ~subcorticalMaskOut;

    %
    numToShrink = numToDilateGeniculi + numToDilateMGNalone - 2 + 1; % **** added 2
    direction = 'posterior';
    leftMGN = removeSlice(leftMGN, direction, numToShrink);
    rightMGN = removeSlice(rightMGN, direction, numToShrink);


    [leftMGN] = fixMGN(leftMGN,leftPulvinar);
    [rightMGN] = fixMGN(rightMGN,rightPulvinar);

    % Save
    thisSubjHome = fullfile(obj.home);
    left_LGN_header = thalamusHeader;
    right_LGN_header = thalamusHeader;
    left_LGN_fileName = 'left_LGN_SearchSpace.nii';
    right_LGN_fileName = 'right_LGN_SearchSpace.nii';
    left_LGN_header.fname = fullfile(thisSubjHome,left_LGN_fileName);
    right_LGN_header.fname = fullfile(thisSubjHome,right_LGN_fileName);
    spm_write_vol(left_LGN_header,double(leftLGN)); pause(eps); drawnow;
    spm_write_vol(right_LGN_header,double(rightLGN)); pause(eps); drawnow;

    obj.l_LGN = left_LGN_header.fname;
    obj.r_LGN = right_LGN_header.fname;

    left_MGN_header = thalamusHeader;
    right_MGN_header = thalamusHeader;
    thisSubjHome = fullfile(obj.home);
    left_MGN_fileName = 'left_MGN_SearchSpace.nii';
    right_MGN_fileName = 'right_MGN_SearchSpace.nii';
    left_MGN_header.fname = fullfile(thisSubjHome,left_MGN_fileName);
    right_MGN_header.fname = fullfile(thisSubjHome,right_MGN_fileName);
    spm_write_vol(left_MGN_header,double(leftMGN)); pause(eps); drawnow;
    spm_write_vol(right_MGN_header,double(rightMGN)); pause(eps); drawnow;

    obj.l_MGN=left_MGN_header.fname;
    obj.r_MGN=right_MGN_header.fname;
    
end

function [MGN] = fixMGN(MGN,pulvinar)
    plane = false(size(MGN,1),size(MGN,2),1);

    anyOne = squeeze(any(MGN,[1,2]));
    mostSuperiorMGN = find(anyOne,1,'last');

    anyOne = squeeze(any(pulvinar,[1,2]));
    mostInferiorPulvinar = find(anyOne,1,'first');

    while (mostInferiorPulvinar < mostSuperiorMGN)
        MGN(:,:,mostSuperiorMGN) = plane;
        anyOne = squeeze(any(MGN,[1,2]));
        mostSuperiorMGN = find(anyOne,1,'last');

        anyOne = squeeze(any(pulvinar,[1,2]));
        mostInferiorPulvinar = find(anyOne,1,'first');
    end

end
