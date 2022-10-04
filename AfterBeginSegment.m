
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

filename = sprintf('results/liver_dices.txt');
fid = fopen(filename,'w'); 
fprintf(fid,'paitentID\tDice\t\n');
fclose(fid);

% Get parameters
lambda=load('results/liver_dices_get_lambda.txt');
%disp(lambda);

% Parameters
max_levels = 5;

clusterlist = {[0 0.14 0.28 0.42 0.56 0.70 0.90 1]};
% [0 0.13 0.26 0.39 0.52 0.75 0.91 1],...
% [0 0.18 0.36 0.54 0.72 0.90 0.97 1],...
% [0 0.17 0.34 0.51 0.68 0.85 0.95 1],...
%  [0.15, 0.29, 0.43, 0.57, 0.71, 0.85, 1],...
%  [0.15, 0.30, 0.46, 0.61, 0.76, 0.91, 1],...
% [0 0.62 0.71 0.80 0.83 0.93 0.98 1],...
% [0 0.61 0.69 0.77 0.82 0.83 0.97 1],...
%  [0 0.60 0.72 0.82 0.94 0.96 0.98 1],...
%  [0 0.63 0.74 0.79 0.82 0.88 0.97 1],...
% [0 0.70 0.74 0.79 0.82 0.88 0.98 1]}; 
%lamblist = [lambda(k,3)];
lamblist = [0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17];
lamblist = [0.16];
lamblist = [-10, -5.0, -1.0, 0.1, 0.5, 5.0, 7.0, 10, 15, 20, 25, 50 ];
%lamblist = [6,7,8,9,10,11,12,13,14,15];
for k = 1 :length(subFolderNames)
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

    outdir = append('results/', subFolderNames{k});
    mkdir(outdir);
    
%     filename = append(outdir,'/anay.txt');
% 
%     fid = fopen(filename,'w');
%     fprintf(fid,'%s\t%s\t%s\t%s\n','n', 'quant', 'threshold', 'dice');
%     fclose(fid);

    filename = sprintf('%s/liver_dice.txt',outdir);
    fid = fopen(filename,'w');
    fprintf(fid,'%s\t%s\t%s\t%s\t\n','patient', 'i', 'lamba', 'dice');
    fclose(fid);


    diceList = zeros(size(lamblist));
    [row_I,col_I,temp] = size(art);
    images = zeros(row_I, col_I, temp, length(lamblist));
    
    for j = 1 : length(lamblist)
        % Get Data
        lamb = [lamblist(j)];
        cluster = clusterlist{1};
        for n = 1 : length(art(1,1,:))
           
            
            img = art(:,:,n);
            imgtruth = truth(:,:,n);
            
            % Get rid of small lines
            %img = sgolayfilt(double(img), 3, 7, [], 2);

            % Perform image segmentation
            segimg = fuzzyimage((img), (imgtruth), ...
                n, outdir, subFolderNames{k}, lamb, cluster, max_levels );
            %segimg = img;
           
            segnifit(:,:,n) = single(segimg);

        end
        [bestDiceImage, dice, quant, threshold] = maximumDice(( segnifit), truth, max_levels);
        %thresh = multithresh(segnifit, max_levels);
        %bestDiceImage = mat2gray(segnifit >= thresh(end));
        % Calculate Dice
        %dice = calculateDice(mat2gray(segnifit),mat2gray(truth));
        diceList(j) = dice;
        images(:,:,:, j) = bestDiceImage;
        fprintf("dice %f lamb %f\n", dice, lamb);
        
    end

    % Get maximum Dice value
    fprintf("DiceList:", diceList);
    [value, pos] = max(diceList);
    segnifit = images(:,:,:, pos);


    filename = sprintf('results/liver_dices.txt');
    fid = fopen(filename,'a+'); 
    %If running begin, change dice to bestDice
    fprintf(fid,'%s \t %f \t %f \n', subFolderNames{k}, diceList(pos), lamblist(pos));
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
