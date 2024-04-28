

files = dir('*.wav');

fid = fopen('array_code.txt','w');

fprintf(fid,'array {\n');

for i = 1:length(files)
    fprintf(fid,['\tsound{ wavefile{ filename = "' files(i).name '"; }; };\n']);
end

fprintf(fid,'}aud_array;');
fclose(fid);
    