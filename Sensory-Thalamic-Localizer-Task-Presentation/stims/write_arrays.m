
runs = 2;
runlength = 24;%runlength * .2 must be an integer
stimpertrial = 9;

% cond = [zeros(runlength*.4,1) + 3; zeros(runlength*.2,1) + 2; zeros(runlength*.2,1) + 1; zeros(runlength*.2,1)];
cond = [zeros(runlength*.5,1) + 2; zeros(runlength*.5,1) + 1;]; %changed to just 50% aud and 50% vis trials

fid = fopen('array_code.txt','w');

fprintf(fid,'array<int>condition_array[runs][runlength] = {');

for i = 1:runs
    rcond = cond(randperm(length(cond)));
    for j = 1:runlength
        if j==1
            fprintf(fid,'{');
        end
        fprintf(fid,num2str(rcond(j)));
        if j ~= runlength
            fprintf(fid,',');
        else
            fprintf(fid,'}');
        end
    end
    if i~=runs
        fprintf(fid,',');
    end
end
fprintf(fid,'};\n');

fprintf(fid,['array<int>stim_array[runs][runlength][' num2str(stimpertrial) '] = {']);
rstim = randperm(runs*runlength*stimpertrial);
count = 0;
for i = 1:runs
    for j = 1:runlength
        if j==1
            fprintf(fid,'{');
        end
        for k = 1:stimpertrial
            if k==1
                fprintf(fid,'{');
            end
            count = count + 1;
            fprintf(fid,num2str(rstim(count)));
            if k~=stimpertrial
                fprintf(fid,',');
            else
                fprintf(fid,'}');
            end
        end
        if j~=runlength
            fprintf(fid,',');
        else
            fprintf(fid,'}');
        end
    end
    if i~=runs
        fprintf(fid,',');
    end
end
fprintf(fid,'};\n');

fclose(fid);