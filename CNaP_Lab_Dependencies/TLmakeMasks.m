function [Y_grayMatter_mask,V_grayMatter_mask,...
        Y_whiteMatter_mask,Y_CSF_mask,...
        grayMatterMask_eroded,whiteMatterMask_eroded,csfMask_eroded] ...
        = TLmakeMasks(FreeSurfer_segmentation_Atlas_path)

    [Y_grayMatter_mask,V_grayMatter_mask] = TLmask(FreeSurfer_segmentation_Atlas_path);
    Y_whiteMatter_mask = logical(TLmask(FreeSurfer_segmentation_Atlas_path,'white'));
    Y_CSF_mask = logical(TLmask(FreeSurfer_segmentation_Atlas_path,'csf'));

    structuringElement = strel('cube',3);
    grayMatterMask_eroded = imerode(TLmask(FreeSurfer_segmentation_Atlas_path) , structuringElement);

    [~,whiteMatterMask_eroded] = TLerodeMasks(Y_whiteMatter_mask); %second output
    [~,csfMask_eroded] = TLerodeMasks(Y_CSF_mask);
end
