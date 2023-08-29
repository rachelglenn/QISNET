


clc; clear;
% TODO, optimize for the cluster and lamblist better. Check out
% results/liver_dices.txt to an idea about the range
fid = fopen('results/LiverDice_lambda.txt','w');
fprintf(fid,'%s\t%s\t%s\n','patientID','lambda', 'Dice');
fclose(fid);

topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..
% Optional fun : Print folder names to command window.
patientDice = zeros(size(length(subFolderNames)));
disp(subFolderNames);





for k = 1 :length(subFolderNames)
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    info = niftiinfo(patient + "/Art.bc.nii.gz");
    fprintf('BitsPerPixel %d', info.BitsPerPixel);
    art = niftiread(info);

    
    mask = niftiread(patient + '/QIS.nii.gz');
    if (size(art,1) ~= size(mask,1)) || (size(art, 2) ~= size(mask, 2))
        disp("sizes are different");
        mask = imresize3(mask, size(art));
    end



    outdir = append('results/', subFolderNames{k});
    disp(outdir);
    mkdir(outdir);
   
    save_figures(art, mask, outdir)
 

end



%[bestDiceImage, bestDice, quant, threshold] = maximumDice(( out_matrix), imgtruth, max_levels);
%bestDice =  calculateDice(bestDiceImage, mat2gray(imgtruth));

% 
%  RGB = label2rgb(imquantize(out_matrix, multithresh(out_matrix,max_levels)),'winter'); 
% 
% subplot (1,3,1), imshow(bestDiceImage, []); title ('Segmented Image');
% subplot (1,3,2), imshow(imgtruth, []); title ('Truth Image');
% subplot (1,3,3), imshow(RGB, []); title ('Segmented Multiple Levels Image');



% 


% filename = sprintf('%s/liver_dice.txt',outdir);
% fid = fopen(filename,'a+');
% fprintf(fid,'%s\t%d\t%f\t%f\t\n',patient, i, lamb, bestDice);
% fclose(fid);

% filename = sprintf('%s/liver_cluster.txt',outdir);
% fid = fopen(filename,'a+');
% fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d \n',clusterListUsed(:, pos));
% fclose(fid);



function save_figures(image, mask, outdir)
    for n = 1 : length(image(1,1,:))
        img = image(:,:,n);
        mask_img = mask(:,:,n);
        rgbImage = imoverlay(mat2gray(img),mat2gray(mask_img), 'red');
        filename = sprintf('%s/Comb_%d.png',outdir, n);
        imwrite( rgbImage,filename);
        
        filename = sprintf('%s/Pred_%d.png',outdir, n);
        imwrite(mat2gray(mask_img),filename);
    end
end