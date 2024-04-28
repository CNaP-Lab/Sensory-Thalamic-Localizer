function [obj] = TLgetROIs(obj, ...
        left_AC_search_region_mask,right_AC_search_region_mask, ...
        left_VC_search_region_mask,right_VC_search_region_mask, ...
        Y_grayMatter_mask,Y_whiteMatter_mask,Y_CSF_mask,...
        grayMatterMask_eroded,whiteMatterMask_eroded,csfMask_eroded, ...
        intermediateOutputDirectory, ...
        fROIoutputDirectory)
    if ~isfield(obj.td,'job1l')
        tstr = ['No 1st level SPM job found for ' obj.home];
        warning(tstr);
        return
    end
    if ~exist([obj.td.job1l{1}.spm.stats.fmri_spec.dir{1} filesep 'SPM.mat'],'file')
        tstr = ['No SPM.mat found in ' obj.td.job1l{1}.spm.stats.fmri_spec.dir{1}];
        warning(tstr);
        return
    end

    mdir = obj.td.job1l{1}.spm.stats.fmri_spec.dir{1};

    conjob{1}.spm.stats.con.consess{1}.tcon.name = 'AUD - VIS';
    conjob{1}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
    conjob{1}.spm.stats.con.consess{1}.tcon.sessrep = 'both';

    conjob{1}.spm.stats.con.delete = 1;

    conjob{1}.spm.stats.con.spmmat = {[mdir filesep 'SPM.mat']};
    spm_jobman('run',{conjob});
    disp(['Calculated basic AUD - VIS contrast for TL model in ' mdir]);

    %This gets TL segmentation images!

    leftLGNfn = obj.l_LGN;
    rightLGNfn = obj.r_LGN;
    [mV_geniculate_left_LGN,mY_geniculate_left_LGN] = tippVol(leftLGNfn);
    [mV_geniculate_right_LGN,mY_geniculate_right_LGN] = tippVol(rightLGNfn);
    mY_geniculate_left_LGN = logical(mY_geniculate_left_LGN);
    mY_geniculate_right_LGN = logical(mY_geniculate_right_LGN);


    leftMGNfn = obj.l_MGN;
    rightMGNfn = obj.r_MGN;
    [mV_geniculate_left_MGN,mY_geniculate_left_MGN] = tippVol(leftMGNfn);
    [mV_geniculate_right_MGN,mY_geniculate_right_MGN] = tippVol(rightMGNfn);
    mY_geniculate_left_MGN = logical(mY_geniculate_left_MGN);
    mY_geniculate_right_MGN = logical(mY_geniculate_right_MGN);

    coactivation_percentile_threshold_LGN = 20;
    coactivation_percentile_threshold_MGN = 32;

    numVoxelsInLeftMask_LGN = ( sum(mY_geniculate_left_LGN(:)) + sum(mY_geniculate_right_LGN(:)) ) / 2;
    geniculiClusterThreshold_LGN = 1 - (coactivation_percentile_threshold_LGN / numVoxelsInLeftMask_LGN);
    disp(['Subject : ' obj.par.par.name ]); pause(eps); drawnow;
    disp(['Voxels in LGN mask : ' num2str(numVoxelsInLeftMask_LGN) ]); pause(eps); drawnow;
    disp(['LGN cluster threshold : ' num2str(geniculiClusterThreshold_LGN) ]); pause(eps); drawnow

    numVoxelsInLeftMask_MGN = ( sum(mY_geniculate_left_MGN(:)) + sum(mY_geniculate_right_MGN(:)) ) / 2;
    geniculiClusterThreshold_MGN = 1 - (coactivation_percentile_threshold_MGN / numVoxelsInLeftMask_MGN); 
    disp(['Subject : ' obj.par.par.name ]); pause(eps); drawnow;
    disp(['Voxels in MGN mask : ' num2str(numVoxelsInLeftMask_MGN) ]); pause(eps); drawnow;
    disp(['MGN cluster threshold : ' num2str(geniculiClusterThreshold_MGN) ]); pause(eps); drawnow;


    mV_temporal_left = spm_vol(left_AC_search_region_mask);
    mY_temporal_left = spm_read_vols(mV_temporal_left);
    mY_temporal_left = logical(mY_temporal_left);

    mV_temporal_right = spm_vol(right_AC_search_region_mask);
    mY_temporal_right = spm_read_vols(mV_temporal_right);
    mY_temporal_right = logical(mY_temporal_right);

    mV_occipital_left = spm_vol(left_VC_search_region_mask);
    mY_occipital_left = spm_read_vols(mV_occipital_left);
    mY_occipital_left = logical(mY_occipital_left);

    mV_occipital_right = spm_vol(right_VC_search_region_mask);
    mY_occipital_right = spm_read_vols(mV_occipital_right);
    mY_occipital_right = logical(mY_occipital_right);

    [contrastVolumeInfo,contrastData] = iimg_read_img([mdir filesep 'con_000' num2str(length(obj.run)+1) '.nii'],2); %Changed .img to nii JCW 01/23/2020

    if contrastVolumeInfo.mat(1) ~= mV_geniculate_left_LGN.mat(1)
        if abs(contrastVolumeInfo.mat(1)) == abs(mV_geniculate_left_LGN.mat(1))
            mY_geniculate_left_LGN = mY_geniculate_left_LGN(end:-1:1,:,:);
            mV_geniculate_left_LGN.mat(1,:) = -mV_geniculate_left_LGN.mat(1,:);
        else
            error('image file sizes do not match');
        end
    end
    if contrastVolumeInfo.mat(1) ~= mV_geniculate_right_LGN.mat(1)
        if abs(contrastVolumeInfo.mat(1)) == abs(mV_geniculate_right_LGN.mat(1))
            mY_geniculate_right_LGN = mY_geniculate_right_LGN(end:-1:1,:,:);
            mV_geniculate_right_LGN.mat(1,:) = -mV_geniculate_right_LGN.mat(1,:);
        else
            error('image file sizes do not match');
        end
    end

    [dat,dat2,dat3,dat4] = getACVC(contrastVolumeInfo,contrastData,Y_grayMatter_mask,mY_temporal_left,mY_temporal_right,mY_occipital_left,mY_occipital_right);

    iimg_write_images(dat,contrastVolumeInfo,fullfile(intermediateOutputDirectory,'ACl.nii'));
    iimg_write_images(dat2,contrastVolumeInfo,fullfile(intermediateOutputDirectory, 'ACr.nii'));
    iimg_write_images(dat3,contrastVolumeInfo,fullfile(intermediateOutputDirectory, 'VCl.nii'));
    iimg_write_images(dat4,contrastVolumeInfo,fullfile(intermediateOutputDirectory, 'VCr.nii'));
    obj.par.td.ROIs.ACl = fullfile(intermediateOutputDirectory,'ACl.nii');
    obj.par.td.ROIs.ACr = fullfile(intermediateOutputDirectory, 'ACr.nii');
    obj.par.td.ROIs.VCl = fullfile(intermediateOutputDirectory, 'ACr.nii');
    obj.par.td.ROIs.VCr = fullfile(intermediateOutputDirectory, 'VCr.nii');
    %disp(['Created auditory and visual cortex ROIs in ' intermediateOutputDirectory]);

    Y_AC_left = spm_read_vols(spm_vol(obj.par.td.ROIs.ACl));
    Y_AC_right = spm_read_vols(spm_vol(obj.par.td.ROIs.ACr));
    Y_VC_left = spm_read_vols(spm_vol(obj.par.td.ROIs.VCl));
    Y_VC_right = spm_read_vols(spm_vol(obj.par.td.ROIs.VCr));


    scans=obj.stMNI_run';


    [Yfraw,Yfraw_unsmoothed] = deal(cell(length(length(scans)),1));
    [nv,nv_unsmoothed] = deal(nan(length(length(scans)),1));

    for i = 1:length(scans)
        Yfraw{i} = spm_read_vols(spm_vol(scans{i}));
        nv(i) = size(Yfraw{i},4);
    end
    Yf = cat(4,Yfraw{:});


    unsmoothed_scans = obj.tMNI_run';

    for i = 1:length(unsmoothed_scans)
        Yfraw_unsmoothed{i} = spm_read_vols(spm_vol(unsmoothed_scans{i}));
        nv_unsmoothed(i) = size(Yfraw_unsmoothed{i},4);
    end
    Yf_unsmoothed = cat(4,Yfraw_unsmoothed{:});

    mY_geniculate_left_D1 = dilate3d((mY_geniculate_left_LGN & mY_geniculate_left_MGN));
    mY_geniculate_right_D1 = dilate3d((mY_geniculate_right_LGN & mY_geniculate_right_MGN));
    mY_geniculate_left_D5 = dilate3d(dilate3d(dilate3d(dilate3d(mY_geniculate_left_D1))));
    mY_geniculate_right_D5 = dilate3d(dilate3d(dilate3d(dilate3d(mY_geniculate_right_D1))));
    leftGeniculate_dilated_regressor_mask = mY_geniculate_left_D5 & ~(mY_geniculate_left_D1) & ~(mY_geniculate_right_D1);
    rightGeniculate_dilated_regressor_mask = mY_geniculate_right_D5 & ~(mY_geniculate_right_D1) & ~(mY_geniculate_left_D1);


    [MGNdat_left,MGNdat_right,LGNdat_left,LGNdat_right] = getMGNLGN(Yf,Yf_unsmoothed,Y_AC_left,Y_AC_right,Y_VC_left,Y_VC_right,whiteMatterMask_eroded,csfMask_eroded, ...
        grayMatterMask_eroded,nv,Y_whiteMatter_mask,mY_geniculate_left_LGN,mY_geniculate_right_LGN,mY_geniculate_left_MGN,mY_geniculate_right_MGN,contrastVolumeInfo, ...
        leftGeniculate_dilated_regressor_mask,rightGeniculate_dilated_regressor_mask, geniculiClusterThreshold_LGN,geniculiClusterThreshold_MGN);

    structuringElement = strel('cube',3);
    mV_geniculate_left_MGN.fname = fullfile(fROIoutputDirectory , 'MGNl.nii');
    MGNdat_left = imerode(dilate3d(logical(MGNdat_left)),structuringElement);
    MGNdat_left = MGNdat_left & mY_geniculate_left_MGN;
    spm_write_vol(mV_geniculate_left_MGN,double(MGNdat_left));
    obj.par.td.ROIs.MGNl = fullfile(fROIoutputDirectory , 'MGNl.nii');

    mV_geniculate_right_MGN.fname = fullfile(fROIoutputDirectory , 'MGNr.nii');
    MGNdat_right = imerode(dilate3d(logical(MGNdat_right)),structuringElement);
    MGNdat_right = MGNdat_right & mY_geniculate_right_MGN;
    spm_write_vol(mV_geniculate_right_MGN,double(MGNdat_right));
    obj.par.td.ROIs.MGNr = fullfile(fROIoutputDirectory , 'MGNr.nii');

    mV_geniculate_left_LGN.fname = fullfile(fROIoutputDirectory , 'LGNl.nii');
    LGNdat_left = imerode(dilate3d(logical(LGNdat_left)),structuringElement);
    LGNdat_left = LGNdat_left & mY_geniculate_left_LGN;
    spm_write_vol(mV_geniculate_left_LGN,double(LGNdat_left));
    obj.par.td.ROIs.LGNl = fullfile(fROIoutputDirectory , 'LGNl.nii');

    mV_geniculate_right_LGN.fname = fullfile(fROIoutputDirectory , 'LGNr.nii');
    LGNdat_right = imerode(dilate3d(logical(LGNdat_right)),structuringElement);
    LGNdat_right = LGNdat_right & mY_geniculate_right_LGN;
    spm_write_vol(mV_geniculate_right_LGN,double(LGNdat_right));
    obj.par.td.ROIs.LGNl = fullfile(fROIoutputDirectory , 'LGNr.nii');
    disp(['Created MGN and LGN fROIs in ' fROIoutputDirectory]); pause(eps); drawnow;

end

function [acdat,acdat2,vcdat,vcdat2] = getACVC(vinfo,cdat,gmY,tmY,tmY2,omY,omY2)

    [tdat,tdat2,odat,odat2] = deal(cdat);
    tdat(~tmY(:)|~gmY(:)) = 0;
    tdat2(~tmY2(:)|~gmY(:)) = 0;
    odat(~omY(:)|~gmY(:)) = 0;
    odat2(~omY2(:)|~gmY(:)) = 0;
    cl = iimg_indx2clusters(tdat,vinfo,quantile(cdat(tmY(:)),.9));
    cl2 = iimg_indx2clusters(tdat2,vinfo,quantile(cdat(tmY2(:)),.9));
    cl3 = iimg_indx2clusters(-odat,vinfo,quantile(-cdat(omY(:)),.9));
    cl4 = iimg_indx2clusters(-odat2,vinfo,quantile(-cdat(omY2(:)),.9));

    sz = [cl.numVox];
    sz2 = [cl2.numVox];
    sz3 = [cl3.numVox];
    sz4 = [cl4.numVox];

    cl(sz<10) = [];
    cl2(sz2<10) = [];
    cl3(sz3<10) = [];
    cl4(sz4<10) = [];

    cl = getMaxPeak(cl);
    cl2 = getMaxPeak(cl2);
    cl3 = getMaxPeak(cl3);
    cl4 = getMaxPeak(cl4);

    acdat = iimg_clusters2indx(cl,vinfo);
    acdat2 = iimg_clusters2indx(cl2,vinfo);
    vcdat = iimg_clusters2indx(cl3,vinfo);
    vcdat2 = iimg_clusters2indx(cl4,vinfo);


end

function [MGNdat_left,MGNdat_right,LGNdat_left,LGNdat_right] = getMGNLGN(Yf,Yf_unsmoothed,Y_AC_left,Y_AC_right,Y_VC_left,Y_VC_right,whiteMatterMask_eroded,csfMask_eroded, ...
        grayMatterMask_eroded,nv,Y_whiteMatter_mask,mY_geniculate_left_LGN,mY_geniculate_right_LGN,mY_geniculate_left_MGN,mY_geniculate_right_MGN,contrastVolumeInfo, ...
        leftGeniculate_dilated_regressor_mask, rightGeniculate_dilated_regressor_mask, geniculiClusterThreshold_LGN, geniculiClusterThreshold_MGN)

    [ACts, AC2ts, VCts, VC2ts, WMts, CSFts, GMts, ...
        leftGeniculateWMts, rightGeniculateWMts, ...
        ] = deal(nan(1,size(Yf,4)));

    leftGeniculate_WM_Mask = ( logical(leftGeniculate_dilated_regressor_mask(:)) | logical(rightGeniculate_dilated_regressor_mask(:)) ) & logical(Y_whiteMatter_mask(:)); %logical(leftGeniculate_regressor_mask(:))&logical(Y_whiteMatter_tpm(:));
    rightGeniculate_WM_Mask = ( logical(leftGeniculate_dilated_regressor_mask(:)) | logical(rightGeniculate_dilated_regressor_mask(:)) ) & logical(Y_whiteMatter_mask(:)); %logical(rightGeniculate_regressor_mask(:))&logical(Y_whiteMatter_tpm(:));

    %Note: combined local WM masks
    numElements = numel(Yf(:,:,:,1)); % How many elements in one?
    Yfvec = nan(numElements,size(Yf,4)); % Preallocate
    for j = 1:size(Yf,4)
        tYf = squeeze(Yf(:,:,:,j));
        tYf_unsmoothed = squeeze(Yf_unsmoothed(:,:,:,j));
        ACts(j) = mean(tYf(logical(Y_AC_left(:))));
        AC2ts(j) = mean(tYf(logical(Y_AC_right(:))));
        VCts(j) = mean(tYf(logical(Y_VC_left(:))));
        VC2ts(j) = mean(tYf(logical(Y_VC_right(:))));
        WMts(j) = mean(tYf_unsmoothed(logical(whiteMatterMask_eroded(:)))); %UNSMOOTHED
        CSFts(j) = mean(tYf_unsmoothed(logical(csfMask_eroded(:)))); %UNSMOOTHED
        GMts(j) = mean(tYf_unsmoothed(logical(grayMatterMask_eroded(:)))); %UNSMOOTHED

        leftGeniculateWMts(j) = mean(tYf_unsmoothed( leftGeniculate_WM_Mask ));  %UNSMOOTHED
        rightGeniculateWMts(j) = mean(tYf_unsmoothed( rightGeniculate_WM_Mask ));

        Yfvec(:,j) = tYf(:);
    end

    dmy = zeros(length(ACts),2);
    dmy(2:3:end,1) = 1;
    dmy(3:3:end,2) = 1;

    for i = 2:length(nv)
        dmy(:,i+1) = zeros(sum(nv),1);
        dmy(sum(nv(1:(i-1)))+1:sum(nv(1:(i-1)))+nv(i),i+1) = ones(nv(i),1);
    end
    for i = 2:length(nv)
        dmy(:,end+1) = dmy(:,1) .* dmy(:,i+1);
    end
    for i = 2:length(nv)
        dmy(:,end+1) = dmy(:,2) .* dmy(:,i+1);
    end

    meanACts = mean([ACts' AC2ts'],2);
    meanVCts = mean([VCts' VC2ts'],2);

    LGNleftGeniculateVecMask = mY_geniculate_left_LGN(:); 
    LGNrightGeniculateVecMask = mY_geniculate_right_LGN(:); 
    MGNleftGeniculateVecMask = mY_geniculate_left_MGN(:);
    MGNrightGeniculateVecMask = mY_geniculate_right_MGN(:); 


    combinedLeftMaskVec = LGNleftGeniculateVecMask|MGNleftGeniculateVecMask;
    combinedRightMaskVec = LGNrightGeniculateVecMask|MGNrightGeniculateVecMask;

    Yf_geniculate_left_vec = Yfvec(combinedLeftMaskVec,:)';
    Yf_geniculate_right_vec = Yfvec(combinedRightMaskVec,:)';

    regressor_wbRs_left = [WMts' CSFts' GMts' dmy];
    regressor_wbRs_right = [WMts' CSFts' GMts' dmy];
    if(any(isnan(leftGeniculateWMts)))
        pause(eps); drawnow;
        warning('NaNs in left WM regressor.')
    end
    if(any(isnan(rightGeniculateWMts)))
        pause(eps); drawnow;
        warning('NaNs in right WM regressor.')
    end
    regressor_wbRs_left = [regressor_wbRs_left leftGeniculateWMts'];
    regressor_wbRs_right = [regressor_wbRs_right rightGeniculateWMts']; % :D

    wbRs_left = partialcorr(Yf_geniculate_left_vec,[meanACts meanVCts],regressor_wbRs_left); % :D
    wbRs_right = partialcorr(Yf_geniculate_right_vec,[meanACts meanVCts],regressor_wbRs_right); % :D

    MGNdat_left = wbRs_left(:,1);
    LGNdat_left = wbRs_left(:,2);
    MGNdat_right = wbRs_right(:,1);
    LGNdat_right = wbRs_right(:,2);

    MGNleftInCombined = MGNleftGeniculateVecMask(combinedLeftMaskVec);
    LGNleftInCombined = LGNleftGeniculateVecMask(combinedLeftMaskVec);
    MGNrightInCombined = MGNrightGeniculateVecMask(combinedRightMaskVec);
    LGNrightInCombined = LGNrightGeniculateVecMask(combinedRightMaskVec);

    MGNdat_left(~MGNleftInCombined) = 0;
    LGNdat_left(~LGNleftInCombined) = 0;
    MGNdat_right(~MGNrightInCombined) = 0;
    LGNdat_right(~LGNrightInCombined) = 0;

    QT_MGN_left = quantile(MGNdat_left,geniculiClusterThreshold_MGN);  %.925 !!!
    QT_MGN_right = quantile(MGNdat_right,geniculiClusterThreshold_MGN);
    QT_LGN_left = quantile(LGNdat_left,geniculiClusterThreshold_LGN);
    QT_LGN_right = quantile(LGNdat_right,geniculiClusterThreshold_LGN);

    MGNdat_left(~(MGNdat_left > QT_MGN_left)) = 0;
    MGNdat_right(~(MGNdat_right > QT_MGN_right)) = 0;
    LGNdat_left(~(LGNdat_left > QT_LGN_left )) = 0;
    LGNdat_right(~(LGNdat_right > QT_LGN_right)) = 0;

    %% end test

    [MGNdatImgVec_left,MGNdatImgVec_right,LGNdatImgVec_left,LGNdatImgVec_right] = ...
        deal(nan(size(Yfvec,1),1));

    MGNdatImgVec_left(combinedLeftMaskVec) = MGNdat_left;
    MGNdatImgVec_right(combinedRightMaskVec) = MGNdat_right;
    LGNdatImgVec_left(combinedLeftMaskVec) = LGNdat_left;
    LGNdatImgVec_right(combinedRightMaskVec) = LGNdat_right;

    MGNcluster_left = iimg_indx2clusters(MGNdatImgVec_left,contrastVolumeInfo);
    MGNcluster_right = iimg_indx2clusters(MGNdatImgVec_right,contrastVolumeInfo);
    LGNcluster_left = iimg_indx2clusters(LGNdatImgVec_left,contrastVolumeInfo);
    LGNcluster_right = iimg_indx2clusters(LGNdatImgVec_right,contrastVolumeInfo);

    sz = [MGNcluster_left.numVox];
    sz2 = [MGNcluster_right.numVox];
    sz3 = [LGNcluster_left.numVox];
    sz4 = [LGNcluster_right.numVox];

    MGNcluster_left(max(sz) > sz) = [];
    MGNcluster_right(max(sz2) > sz2) = [];
    LGNcluster_left(max(sz3) > sz3) = [];
    LGNcluster_right(max(sz4) > sz4) = [];

    MGNcluster_left = getMaxPeak(MGNcluster_left);
    MGNcluster_right = getMaxPeak(MGNcluster_right);
    LGNcluster_left = getMaxPeak(LGNcluster_left);
    LGNcluster_right = getMaxPeak(LGNcluster_right);

    MGNdat_left = iimg_clusters2indx(MGNcluster_left,contrastVolumeInfo);
    MGNdat_right = iimg_clusters2indx(MGNcluster_right,contrastVolumeInfo);
    LGNdat_left = iimg_clusters2indx(LGNcluster_left,contrastVolumeInfo);
    LGNdat_right = iimg_clusters2indx(LGNcluster_right,contrastVolumeInfo);

    LGNdat_left(logical(MGNdat_left)) = 0;
    LGNdat_right(logical(MGNdat_right)) = 0;
    MGNdat_left(logical(LGNdat_left)) = 0;
    MGNdat_right(logical(LGNdat_right)) = 0;

    MGNdat_left = reshape(MGNdat_left,size(mY_geniculate_left_LGN));
    MGNdat_right = reshape(MGNdat_right,size(mY_geniculate_right_LGN));
    LGNdat_left = reshape(LGNdat_left,size(mY_geniculate_left_LGN));
    LGNdat_right = reshape(LGNdat_right,size(mY_geniculate_right_LGN));

end

function [cout] = getMaxPeak(cin,varargin)
    cout = cin;
    cout(:) =[];

    if isempty(varargin)
        numr = 1;
    else
        numr = varargin{1};
    end

    for i = 1:numr
        mx = zeros(size(cin));
        for j = 1:length(mx)
            mx(j) = max(cin(j).Z);
        end
        cout(i) = cin(mx==max(mx));
        cin(mx==max(mx)) = [];
    end


end
