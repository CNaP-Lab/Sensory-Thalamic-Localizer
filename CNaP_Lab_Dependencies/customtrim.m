function obj = customtrim(obj,nrmv,varargin)
    % trim both images and MP at the same time

    if isempty(varargin)
        prf = 't';
    else
        prf = varargin{1};
    end

    runloop=1;

    % only providing rawMNI data
    if isfield(obj,'rawMNI_run') && ~isfield(obj,'sMNI_run') && ~isfield(obj,'tMNI_run') && ~isfield(obj,'stMNI_run')
        ref=1;
        images=obj.rawMNI_run;
        %providing only rawMNI & sMNI data
    elseif isfield(obj,'rawMNI_run') && isfield(obj,'sMNI_run') && ~isfield(obj,'tMNI_run') && ~isfield(obj,'stMNI_run')
        ref=2;
        images=obj.rawMNI_run;
        imagesb=obj.sMNI_run;
        runloop=2;
        % providing only tMNI & sMNI data
    elseif isfield(obj,'tMNI_run') && isfield(obj,'sMNI_run') && ~isfield(obj,'stMNI_run')
        ref=3;
        images=obj.sMNI_run;
    else
        error('Error - Please check the MNI time series data type(s) you have provided. You must provide either rawMNI and/or sMNI data. Only provide required data.');
        % return;
    end




    for uu=1:runloop

        if uu==2
            images=imagesb;
        end

        for i = 1:numel(images)
            [P,F,X] = fileparts(images{i});
            V = spm_vol(images{i});
            Y = spm_read_vols(V);

            V(1:nrmv) = [];
            Y(:,:,:,1:nrmv) = [];
            images{i} = [P  filesep prf F X];
            for k = 1:length(V)

                V(k).fname = images{i};
                V(k).n(1) = k;
                spm_write_vol(V(k),squeeze(Y(:,:,:,k)));
            end
        end


        if ref==1
            obj.tMNI_run=images;
        elseif (ref==2 && runloop==1)
            obj.tMNI_run=images;
        elseif (ref==2 && runloop==2)
            obj.stMNI_run=images;
        elseif ref==3
            obj.stMNI_run=images;
        end

    end


    % also gather and trim MPs, if needed
    % if obj.MPfile_trim==false
    for i=1:length(obj.MPfileraw_name_run)
        obj.run(i).td.MPs=readmatrix(obj.MPfileraw_name_run{i});
        obj.run(i).td.MPs(1:nrmv,:) = [];
    end
    % end

end
