function obj = TLreadLogs(obj,participant_TL_task_log_path)


    %set up
    [a,b,c]=fileparts(participant_TL_task_log_path);
    obj.td.logs.fnames={[b c]};
    obj.td.logs.copied=true;

    %specify name of dataset with arbitrary name
    obj.par.name = 'dataset';

    if ~isempty(obj.td) && isfield(obj.td,'logs') && isfield(obj.td.logs,'fnames')

        if isfield(obj.td,'conditions')
            obj.td = rmfield(obj.td,'conditions');
        end


        if length(obj.td.logs.fnames) ~= 1
            tstr = [num2str(length(participant_TL_task_log_path)) ' log files found for ' obj.par.par.name '. Skipping.'];
            warning(tstr);
            return;
        end

        fid = fopen(participant_TL_task_log_path{1});

        T = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s');

        TMP = strfind(T{4},'Subject_1_Run_');
        for i = 1:length(TMP)
            if isempty(TMP{i})
                TMP{i} = 0;
            end
        end
        runs = unique(T{4}(logical([TMP{:}])));

        if length(runs) < 1
            TMP = strfind(T{4},['Subject_' obj.par.par.name '_Run_']);

            if ( ~all(cellfun(@isempty,TMP)) )
                logSubjID = obj.par.par.name;

            else
                % What is the subject ID?
                filtered = T{4}(contains(T{4},'Subject_'));
                logSubjID = strtok(filtered{1},'Subject_');
                TMP = strfind(T{4},['Subject_' logSubjID '_Run_']);
            end
            for i = 1:length(TMP)
                if isempty(TMP{i})
                    TMP{i} = 0;
                end
            end
            runs = unique(T{4}(logical([TMP{:}])));
        end

        sidx = zeros(1,length(runs)+1);
        for i = 1:length(runs)
            try
                sidx(i) = find(strcmp(T{4},['Subject_1_Run_' num2str(i)]),1);
            catch
                try
                    sidx(i) = find(strcmp(T{4},['Subject_' logSubjID '_Run_' num2str(i)]),1);
                catch
                    tstr = ['Failed to find Run_' num2str(i) ' data in ' obj.td.logs.fnames{1}];

                    warning(tstr);
                    return;
                end
            end
        end
        sidx(end) = length(T{4});

        [obj.td.timing.pulses,obj.td.timing.vis_chk,obj.td.timing.aud_clp] = deal(cell(size(runs)));

        for r = 1:length(runs)
            idx = sidx(r):sidx(r+1);
            times = T{5}(idx);
            pulses = times(strcmp(T{3}(idx),'Pulse'));
            presents = times(strcmp(T{4}(idx),'check1')|strcmp(T{4}(idx),'check2')|strcmp(T{3}(idx),'Sound'));
            vis = times(strcmp(T{4}(idx),'check1')|strcmp(T{4}(idx),'check2'));
            aud = times(strcmp(T{3}(idx),'Sound'));
            try
                sstim = str2num(presents{1});
            catch err
                continue;
            end

            for i = 1:length(pulses)
                pls = str2num(pulses{i});
                if pls > sstim
                    fpi = i - 3;
                    tstrt = str2num(pulses{fpi})/10;
                    break
                end
            end
            for i = fpi:length(pulses)
                obj.td.timing.pulses{r}(i-(fpi-1))  = str2num(pulses{i})/10 - tstrt;

            end
            for i = 1:length(vis)
                obj.td.timing.vis_chk{r}(i) = str2num(vis{i})/10 - tstrt;

            end
            for i = 1:length(aud)
                obj.td.timing.aud_clp{r}(i) = str2num(aud{i})/10 - tstrt;
            end

            auditoryTimings = obj.td.timing.aud_clp{r}';
            visualTimings = obj.td.timing.vis_chk{r}';

            auditoryDiff = diff(auditoryTimings);
            visualDiff = diff(visualTimings);
            %This is what gets used:
            auditoryOnsetIndex = [1; find(auditoryDiff>1200)+1];
            visualOnsetIndex = [1; find(visualDiff>1200)+1];
            %
            auditoryIndexDiff = diff(auditoryOnsetIndex);
            visualIndexDiff = diff(visualOnsetIndex);
            auditoryWrong = auditoryIndexDiff < 9;
            visualWrong = visualIndexDiff < 68;

            if(any(auditoryWrong(1:end-1)))
                firstAuditoryWrongIndex = find(auditoryWrong,1,'first');
                if ~auditoryWrong(firstAuditoryWrongIndex+1)
                    error('Tried to fix auditory stimulus timings, failed.');
                end
                auditoryOnsetIndex(firstAuditoryWrongIndex+1) = auditoryOnsetIndex(firstAuditoryWrongIndex+1) + auditoryOnsetIndex(firstAuditoryWrongIndex);
                auditoryOnsetIndex(firstAuditoryWrongIndex) = [];
            end

            if(any(visualWrong(1:end-1)))
                firstVisualWrongIndex = find(visualWrong,1,'first');
                if ~visualWrong(firstVisualWrongIndex+1)
                    error('Tried to fix visual stimulus timings, failed.');
                end
                visualOnsetIndex(firstVisualWrongIndex+1) = visualOnsetIndex(firstVisualWrongIndex+1) + visualOnsetIndex(firstVisualWrongIndex);
                visualOnsetIndex(firstVisualWrongIndex) = [];
            end


            auditoryOnsetVec = obj.td.timing.aud_clp{r}(auditoryOnsetIndex);
            visualOnsetVec = obj.td.timing.vis_chk{r}(visualOnsetIndex);

            %                 sorry for not commenting this next line.
            [~,IDX] = sort([auditoryOnsetVec visualOnsetVec]);

            for i = 1:length(IDX)
                obj.td.conditions{r}{i} = (IDX(i) > 8) + 1;

            end

        end
        fclose(fid);

    end

end
