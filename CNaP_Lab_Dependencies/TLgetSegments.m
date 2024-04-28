function obj = TLgetSegments(obj,thalsegfile)

    %TLGETSEGMENTS
    % pull thalamic segmentation from preprocessed directory 

    obj.subj.name=obj.par.par.name;

    thisSubjHome = fullfile(obj.subj.home);
    % path to the preprocessed data where segmentation file lives
    thisSubjThalSegSourceFile = thalsegfile;

    % copy file to home dir
    copyfile(thisSubjThalSegSourceFile,thisSubjHome);

    disp(['Got TL segmentation for : ' obj.subj.name '.']); pause(eps); drawnow;

end

