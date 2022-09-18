
clc; clear;

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

filename = sprintf('results/liver_dices.txt');
fid = fopen(filename,'w'); 
fprintf(fid,'paitentID\tDice\t\n');
fclose(fid);

for k = 1 : length(subFolderNames)
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    info = niftiinfo(patient + "/Art.bc.nii.gz");
    fprintf('BitsPerPixel %d', info.BitsPerPixel);
    art = niftiread(info);

    
    truth = niftiread(patient + '/Truth.raw.nii.gz');
    if (size(art,1) ~= size(truth,1)) || (size(art, 2) ~= size(truth, 2))
        disp("sizes are different");
        truth = imresize3(truth, size(art));
    end


    quantList = zeros(size(length(art(1,1,:))));
    segnifit = zeros(size(art),'single');
    truth = single(truth);
%     if info.BitsPerPixel == 16
%         segnifit = zeros(size(art),'int16');
%         truth = int16(truth);
%     elseif info.BitsPerPixel == 32
%         segnifit = zeros(size(art),'single');
%         truth = single(truth);
%     else
%         segnifit = zeros(size(art),'int64');
%         truth = double(truth);   
%     end    

    outdir = append('results/', subFolderNames{k});
    mkdir(outdir);
    
    filename = append(outdir,'/anay.txt');

    fid = fopen(filename,'w');
    fprintf(fid,'%s\t%s\t%s\t%s\n','n', 'quant', 'threshold', 'dice');
    fclose(fid);


    for n = 1 : length(art(1,1,:))
        fprintf("scan num:%d\n",n);
        img = art(:,:,n);
 
       %% 
       imgtruth = truth(:,:,n);
        %Perform Algorithm
        img = sgolayfilt(img, 2, 9, [], 1);
        segimg = fuzzyimage(img, single(imgtruth), n, outdir, subFolderNames{k} );
        %segimg = img;
        segnifit(:,:,n) = single(segimg);
%        if info.BitsPerPixel == 16
%             segnifit(:,:,n) = int16(segimg);
%        elseif info.BitsPerPixel == 32
%             segnifit(:,:,n) = single(segimg);
%        else   
%             segnifit(:,:,n) = double(segimg);
%        end
    end
    % Calculate Dice
    dice = calculateDice(double(segnifit),double(truth));
   
    filename = sprintf('results/liver_dices.txt');
    fid = fopen(filename,'a+'); 
    fprintf(fid,'%s\t%f\t\n', subFolderNames{k}, dice);
    fclose(fid);
    info = niftiinfo(patient + "/Art.bc.nii.gz");
    if info.BitsPerPixel == 16
        filename = sprintf('%s/QIS.nii',outdir);
        disp(filename);
        info.BitsPerPixel = 16;
     
        niftiwrite(int16(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);

    elseif info.BitsPerPixel == 32
        filename = sprintf('%s/QIS.nii',outdir);
        disp(filename);
        info.BitsPerPixel = 32;

        niftiwrite(single(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);
    else
        filename = sprintf('%s/QIS.nii',outdir);
        disp(filename);
        info.BitsPerPixel = 64;
      
        niftiwrite(double(segnifit),filename, info, 'Version', 'NIfTI1',  'Compressed',true);

     end
end

function biggerIm = resizeImage(img, scaleFactor)
     
    biggerIm = zeros(size(img));   
    for i = 1 : length(img(1,1,:))
        n = 1;  
        I =(mat2gray(img(:,:,i)));
        %Idouble = im2double(I); 
        %avg = mean2(Idouble);
        %sigma = std2(Idouble);
        img_c = imadjust(I,[],[],0.95);
        %img_c = imadjust(I,[avg-n*sigma avg+n*sigma],[]);
        %img_c = int16(img(:,:,n));

        biggerIm(:,:,n) = img_c; 
    end
end


