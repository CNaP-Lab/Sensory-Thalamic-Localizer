function [Y,V] = TLmask(maskfilein,varargin)

    % provide the name of the proper file containing MNI parcels,
    % such as Atlas.wmparc.2.nii

    %     Y = TLmask(maskfilein);
    %     returns a 3d logical of in-mask voxels for the grey matter mask of
    %     the subject associated with the call.
    %     [Y,V] = TLmask(maskfilein);
    %     also returns V, an spm_vol structure for the associated mask.
    %     Y = TLmask(maskfilein,'white');
    %     returns white matter mask instead of grey matter.
    %     Y = TLmask(maskfilein,'csf');
    %     returns CSF instead of grey matter.
    %     Y = TLmask(maskfilein,'custom',{idx1,idx2,...});
    %     returns the values equal to any of the values idx1, idx2, etc. in the
    %     Atlas_wmparc.2.nii HCP/freesurfer parcellation image



    r = [1 3000];
    if ~isempty(varargin)
        switch lower(varargin{1})
            case 'white'
                r = [3000 5000];
            case 'brain'
                r = [];
            case 'csf'
                r = {4,43,5,44,14,15,72};
            case 'custom'
                r = varargin{2};
            case 'gray'
                r = [1 3000];
        end
    end


    fn = maskfilein;
    %             V = spm_vol(fn);
    %             Y = spm_read_vols(V);
    [V,Y] = tippVol(fn);
    if ~iscell(r)
        Y = Y >= r(1) & Y < r(2);
    else
        Yout = Y == r{1};
        for i = 2:length(r)
            Yout = Yout | Y == r{i};
        end
        Y = Yout;
    end

end
