[X,Y] = meshgrid(-210:20:210,-170:20:170);

fid = fopen('pictureCode.txt','w');

fprintf(fid,'picture {\n');

switcher = false;
for i = 1:numel(X)
    if (mod(i,2) && ~switcher ) || (~mod(i,2) && switcher)
        fprintf(fid,'\tbox black;\n');
    else
        fprintf(fid,'\tbox white;\n');
    end
    fprintf(fid,['\tx = ' num2str(X(i)) '; y = ' num2str(Y(i)) ';\n']);
    
    if ~mod(i,size(X,1))
        if switcher
            switcher = false;
        else
            switcher = true;
        end
    end
end
fprintf(fid,'}check1;\n\n');

fprintf(fid,'picture {\n');

switcher = true;
for i = 1:numel(X)
    if (mod(i,2) && ~switcher ) || (~mod(i,2) && switcher)
        fprintf(fid,'\tbox black;\n');
    else
        fprintf(fid,'\tbox white;\n');
    end
    fprintf(fid,['\tx = ' num2str(X(i)) '; y = ' num2str(Y(i)) ';\n']);
    
    if ~mod(i,size(X,1))
        if switcher
            switcher = false;
        else
            switcher = true;
        end
    end
end
fprintf(fid,'}check2;\n');
            
fclose(fid);