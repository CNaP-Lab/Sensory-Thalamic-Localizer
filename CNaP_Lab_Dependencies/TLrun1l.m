function [obj] = TLrun1l(obj,varargin)

    %force = false;
    force = true;
    tlvl = class(obj);

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'force'
                    force = true;
                case 'tlvl'
                    tlvl = varargin{i+1};
            end
        end
    end


    mdir = [obj.home filesep 'model1l']; %the model directory name should be setable via user input both here and in SOprep1l.m

    if ~exist(mdir)
        mkdir(mdir)
    elseif ~isempty(dir([mdir filesep 'SPM.mat']))
        if force
            disp(['Deleting contents of ' mdir])
            delete([mdir filesep 'SPM.mat']);
            delete([mdir filesep '*.nii']);
            disp(['Deleted contents of ' mdir]);

        else
            tstr = ['.mat files exist in ' mdir '. Not running 1st level model.'];
            disp(tstr);
            out = [];
            return
        end
    end

    if ~isfield(obj.td,'job1l')
        tstr = ['No SPM job prepared for ' obj.home];
        disp(tstr);

        out = [];
        return
    end

    out = obj.td.job1l;


    if isa(obj,tlvl)

        disp('Estimating the following 1st level models:');

        for i = 1:length(out)
            estimate_job{i}.spm.stats.fmri_est.method.Classical = 1;
            estimate_job{i}.spm.stats.fmri_est.spmmat = {[out{i}.spm.stats.fmri_spec.dir{1} filesep 'SPM.mat']};
            disp(out{i}.spm.stats.fmri_spec.dir{1});
        end

        for i = 1:length(out)
            spm_jobman('run',out(i));
        end
        disp('Model spec done.');

        for i = 1:length(estimate_job)
            spm_jobman('run',{estimate_job(i)});
        end
        disp('Estimation done.');

    end

end


function out = unpackcell(in)
    out = {};
    for i = 1:length(in)
        out(end+1:end+length(in{i})) = in{i};
    end
end
