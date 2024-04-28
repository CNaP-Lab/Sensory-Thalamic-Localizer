function obj = TLprep1l(obj,Y_grayMatter_mask,V_grayMatter_mask)

    for i=1:length(obj.MPfile_name_run)
        obj.run(i).td.MPs=readmatrix(obj.MPfile_name_run{i});
    end

    for i=1:length(obj.MPfile_name_run)
        obj.run(i).name=['run' num2str(i)]; %auto generate the run names
    end


    obj.run(1).p.def.TR=obj.spmuse.TR;

    modelDir = [obj.home filesep 'model1l'];

    if isempty(obj.run)
        str = 'No runs found';
        error(str);
    end

    if isfolder(modelDir)
        delete([modelDir filesep '*.*']);
    else
        mkdir(modelDir);
    end

    Y_grayMatter_mask = dilate3d(Y_grayMatter_mask); %Dilate gray matter mask by 1
    V_grayMatter_mask.fname = [obj.par.par.home filesep 'spmmask.nii'];
    spm_write_vol(V_grayMatter_mask,Y_grayMatter_mask);


    if ~isfield(obj.td,'conditions')
        return
    end

    nscans = length(obj.run);

    for i = 1:length(obj.run)
        %run name below
        rn(i) = str2num(obj.run(i).name(end));

    end

    if(length(obj.td.conditions) < nscans)

        warning(['Missing run in task log for subject : ' num2str(obj.par.par.name)]); pause(eps); drawnow;

        warning(['Imaging : ' num2str(length(rn)) ' runs.  Task log : ' num2str(length(obj.td.conditions)) ' runs.']); pause(eps); drawnow;

        difference = length(rn) - length(obj.td.conditions);

        rn(end-(difference-1)) = [];
        nscans = nscans - difference;
    end

    conditions = obj.td.conditions(rn);


    for i = 1:length(conditions)
        for j = 1:length(conditions{i})
            if conditions{i}{j}==1
                reg1{i}(3 + (j-1)*3 + 1:3 + (j-1)*3 + 3,:) = 1;
                reg2{i}(3 + (j-1)*3 + 1:3 + (j-1)*3 + 3,:) = 0;
            elseif conditions{i}{j}==2
                reg1{i}(3 + (j-1)*3 + 1:3 + (j-1)*3 + 3,:) = 0;
                reg2{i}(3 + (j-1)*3 + 1:3 + (j-1)*3 + 3,:) = 1;
            end
        end
    end

    dreg = zeros(length(conditions{1})*3+3,2);
    dreg(1:3:end,1) = 1;
    dreg(2:3:end,2) = 1;


    firstlevel_job{1}.spm.stats.fmri_spec.timing.units = 'secs';
    firstlevel_job{1}.spm.stats.fmri_spec.timing.RT=obj.run(1).p.def.TR;
    firstlevel_job{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    for i = 1:nscans
        firstlevel_job{1}.spm.stats.fmri_spec.sess(i).hpf = 10000;
    end
    firstlevel_job{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    firstlevel_job{1}.spm.stats.fmri_spec.volt = 1;
    firstlevel_job{1}.spm.stats.fmri_spec.global = 'None';

    firstlevel_job{1}.spm.stats.fmri_spec.cvi = 'none';

    firstlevel_job{1}.spm.stats.fmri_spec.mask = {''};

    firstlevel_job{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    firstlevel_job{1}.spm.stats.fmri_spec.dir = {modelDir};

    fns=obj.stMNI_run';

    for j = 1:nscans
        list = spm5_image_list(size(obj.run(j).td.MPs,1), fns(j));
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).scans = list{1}';

        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(1).name = 'auditory';
        %firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(1).val = reg1{mod(j+1,2)+1};  %%% CHECK THIS
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(1).val = reg1{j}; % Edited JCW 01/23/2020 - fixed
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(2).name = 'visual';
        %firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(2).val = reg2{mod(j+1,2)+1};  %%% CHECK THIS
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(2).val = reg2{j}; % Edited JCW 01/23/2020 - fixed
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(3).name = 'volIntercept_1';
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(3).val = dreg(:,1);
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(4).name = 'volIntercept_2';
        firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(4).val = dreg(:,2);


        M = obj.run(j).td.MPs(:,1:6);

        M = [M M(:,1:6).^2];

        for k = 1:size(M,2)
            firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(k+4).name = ['M' num2str(k)];
            firstlevel_job{1}.spm.stats.fmri_spec.sess(j).regress(k+4).val = M(:,k);
        end

    end
    %firstlevel_job{1}.spm.stats.fmri_spec.sess(j).cond.orth = 0;

    obj.td.job1l = firstlevel_job;

end
