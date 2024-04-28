function obj = TLsmoothVols(obj,smoothingFWHM)
    %
    %     optional inputs formatted as with a call to fnames.
    f=obj.tMNI_run;
    for i=1:length(f)
        nv(i)=length(tippVol(f{i}));
    end

    disp('This will smooth the following volumes:'); pause(eps); drawnow;


    for i = 1:length(f)
        disp(f{i}); pause(eps); drawnow;
    end

    %     smooth_job{1}.spm.spatial.smooth.fwhm = [4 4 4];
    smooth_job{1}.spm.spatial.smooth.fwhm = [smoothingFWHM smoothingFWHM smoothingFWHM];

    smooth_job{1}.spm.spatial.smooth.dtype = 0;
    smooth_job{1}.spm.spatial.smooth.im = 0;
    smooth_job{1}.spm.spatial.smooth.prefix = 's';
    if length(f) > 4
        startpar(length(f));

        parfor i = 1:length(f)
            pause(eps); drawnow;
            jb = jobedit(smooth_job,spm5_image_list(nv(i),f(i)));
            spm_jobman('run',{jb});
            pause(eps); drawnow;
        end
    else
        for i = 1:length(f)
            pause(eps); drawnow;
            jb = jobedit(smooth_job,spm5_image_list(nv(i),f(i)));
            spm_jobman('run',{jb});
            pause(eps); drawnow;
        end
    end
    for i=1:length(obj.tMNI_run)
        [fa,fb,fc]=fileparts(obj.tMNI_run{i});
        obj.stMNI_run{i}=fullfile(fa,['s' fb fc]);
    end


    for i = 1:length(f)
        disp(['Smoothed volume ' f{i}]); pause(eps); drawnow;

    end

end


function jb = jobedit(jb,list)

    jb{1}.spm.spatial.smooth.data = list{1}';

end
