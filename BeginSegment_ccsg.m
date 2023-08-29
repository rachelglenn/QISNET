
clc; clear;
% TODO, optimize for the cluster and lamblist better. Check out
% results/liver_dices.txt to an idea about the range
fid = fopen('results/LiverDice_lambda.txt','w');
fprintf(fid,'%s\t%s\t%s\n','patientID','lambda', 'Dice');
fclose(fid);

topLevelFolder = '/rsrch1/ip/rglenn1/data/Processed';
topLevelFolder = '/rsrch1/ip/rglenn1/data/ccsg/nifti';
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
%lambda=load('results/10_4_2022_liver_dices.txt');
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
lamblist = [-10, -5.0, -1.0, 0.1, 0.5, 5.0, 7.0, 10, 15, 20, 25, 50 ];
%lamblist = lambda(:,3);
lamblist = [1.5];
for k = 1 :length(subFolderNames)
    
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    info = niftiinfo(patient + "/image.nii.gz");
    fprintf('BitsPerPixel %d', info.BitsPerPixel);
    art = niftiread(info);
    manualScale = @(I) 255.0*(I-min(I(:)))/(max(I(:))-min(I(:)));

    %for n = 1 : length(art(1,1,:))
    %    segimg = single(art(:,:,n));
        %segimg = manualScale(segimg);
    %    fprintf('%f \t %f \t %f \t %f \n', std(segimg(:)), mean(segimg(:)), max(segimg(:)), min(segimg(:)));
    %end
    art = mat2gray( double(art) , [min(double(art(:))) max(double(art(:)))] );
    %figure(1); 
    %for n = 1 : length(art(1,1,:))
    %  segimg = single(art(:,:,n));
    %  imshow(segimg);
      %pixel_values = impixel
    %end
    %art = double(art);
    cluster_matrix = manualScale(art); % changem( art , 0.0 , -1024.0 );
    min_value = min(cluster_matrix(:));
    max_value = max(cluster_matrix(:));
    mean_value = mean(cluster_matrix(:));
    std_value = std(cluster_matrix(:));
    cluster_min = 0; %mean_value - std_value;
    cluster_max = 3.14; %mean_value + std_value;
    num_clusters = 20;
    clusterlist = {linspace(cluster_min,cluster_max, num_clusters)};
    fprintf('\nclusterlist:');
    disp(clusterlist);

    
    %truth = niftiread(patient + '/Truth.raw.nii.gz');
    %if (size(art,1) ~= size(truth,1)) || (size(art, 2) ~= size(truth, 2))
    %    disp("sizes are different");
    %    truth = imresize3(truth, size(art));
    %end

    quantList = zeros(size(length(art(1,1,:))));
    segnifit = zeros(size(art),'single');
    %truth = single(truth);

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
    %lamb = [lamblist(k)];

    filename = string(outdir) + '/image_results' + '.txt';
    fid = fopen(filename,'a+'); 
    
    for j = 1 : length(lamblist)
        % Get Data
        lamb = [lamblist(j)];
        cluster = clusterlist{1};
        
        for n = 1 : length(art(1,1,:))
            fprintf("\nOn %i \n", n);
            
            img = art(:,:,n);
            %imgtruth = truth(:,:,n);
            
            % Get rid of small lines
            %img = sgolayfilt(double(img), 3, 7, [], 2);

            % Perform image segmentation
            segimg = fuzzyimage((img), (img), ...
                n, outdir, subFolderNames{k}, lamb, cluster, max_levels );
            %segimg = img;
           
            segnifit(:,:,n) = segimg;
            %for x = 1 : length(segimg(:,:,n))

            %imshow(I2(:,:,n));

            fprintf('%f \t %f \t %f \t %f \n', std(double(segimg(:))), mean(double(segimg(:))), max(double(segimg(:))), min(double(segimg(:))));
            %I2 = 255*(segimg - min(segimg(:))) ./ (max(segimg(:)) - min(segimg(:))); %scale values between 0 and 255
            I2 = mat2gray( double(segimg) , [min(double(segimg(:))) max(double(segimg(:)))] );
            %I2 = uint8(255*mat2gray(segimg));
            figure(2);
            imshow(I2, [])
            figure(1);
            imshow(img);
            imwrite(I2(:,:), string(outdir) + '/QIS_' + string(n) + '.jpeg','JPEG');
            segnifit(:,:,n) = I2;
            %end
            

        end
        %[bestDiceImage, dice, quant, threshold] = maximumDice(( segnifit), truth, max_levels);
        %thresh = multithresh(segnifit, max_levels);
        %bestDiceImage = mat2gray(segnifit >= thresh(end));
        % Calculate Dice
        %dice = calculateDice(mat2gray(segnifit),mat2gray(truth));
        %diceList(j) = dice;
        %images(:,:,:, j) = segimg;
        %fprintf("dice %f lamb %f\n", dice, lamb);
        
    end
    fclose(fid);
    % Get maximum Dice value
    %fprintf("DiceList:", diceList);
    %[value, pos] = max(diceList);
    %segnifit = images(:,:,:, pos);


    %filename = sprintf('results/liver_dices.txt');
    %fid = fopen(filename,'a+'); 
    %If running begin, change dice to bestDice
    %fprintf(fid,'%s \t %f \t %f \n', subFolderNames{k}, diceList(pos), lamblist(pos));
    %fclose(fid);
    %disp("Creating the figure");
    %figure(2);
    info = niftiinfo(patient + "/image.nii.gz");
    niftiName = 'QIS.nii';
    niftisave = segnifit;
    %I2 = (segnifit - min(segnifit(:))) ./ (max(segnifit(:)) - min(segnifit(:))); %scale values between 0 and 255
    %I2 = cast(I2,'int8');
    
   %for n = 1 : length(I2(1,1,:))
   %     imshow(I2(:,:,n));
   %     imwrite(uint8(I2(:,:,n)), string(outdir) + '/QIS_' + string(n) + '.jpeg','JPEG');
   %end
    
    I2 = manualScale(segnifit);
    saveNifti(outdir, info, niftiName, uint8(I2));
    %imdisp(segnifit);
    fprintf("std %f",  std(I2(:)));
    fprintf("mean %f",mean(I2(:)));
    fprintf("min %f",  min(I2(:)));
    fprintf("min %f",  max(I2(:)));
    
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


