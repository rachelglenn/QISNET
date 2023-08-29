

function saveNifti(outdir, info, niftiName,segnifit)
    filename = sprintf('%s/%s',outdir, niftiName);
    disp(filename);
    disp('bitsperPixel');
    disp(info.BitsPerPixel);
    if info.BitsPerPixel == 16
        info.BitsPerPixel = 16;
        niftiwrite(int16(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);
        %niftiwrite(single(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);

    elseif info.BitsPerPixel == 32
        info.BitsPerPixel = 32;
        niftiwrite(single(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);
    else
        info.BitsPerPixel = 64;
        niftiwrite(double(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);

    end
end
