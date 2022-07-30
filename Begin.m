% input_loc = 'Data_SkullStripped/SkullStripped/';
% output_loc='Output_SkullStripped/';
% for i=1:80
%     path = strcat(input_loc, strcat( int2str(i),'.png'));
%     respath = strcat(input_loc, strcat( int2str(i),'.png'));
%     
%     sprintf( path);
%     im=imread (path);
%     %imshow(im);
%     fuzzyimage(im,i, output_loc);
% end
clc; clear;

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
for k = 1 : length(subFolderNames)
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    art = niftiread(patient + '/Pre.raw.nii.gz');
    truth = niftiread(patient + '/Truth.raw.nii.gz');
    diceList = zeros(size(length(art(1,1,:))));
    outdir = append('results/', subFolderNames{k});
    mkdir(outdir)
    for n = 1 : length(art(1,1,:))
        img = art(:,:,n);
        imgtruth = truth(:,:,n);
        % img = resizeImage(img, 5/6);
        % img = adjustIm(img);
        %imgtruth = resizeImage(imgtruth, 5/6);
        % imgtruth = adjustIm(imgtruth);
        % disp(size(img));
        segimg = fuzzyimage(img, imgtruth, n, outdir, subFolderNames{k} );
        % subplot(1,2,1);
        % imshow(img, []);
        % subplot(1,2,2);
        % imshow(segimg, []);
        % filename = sprintf('LiverOutput/Truth_%d.png',n );
        %imwrite(uint8(imgtruth),filename);
        % before_after_img = imtile({img,segimg});
        %disp("Size");
        %disp(ndims(segimg));
        %imshow(imgtruth, []);
        imga = double(segimg);
        imgb = double(imgtruth);
        filename = sprintf('%s/Truth_%d_%s.png',outdir, n, subFolderNames{k});
        imwrite(double(imgtruth),filename);
        diceList(n) = generalizedDice(imga,imgb);
        %hold;
    end
    patientDice(k) = mean(diceList);
end
fid = fopen('results/LiverDice.txt','w');
fprintf(fid,'%f\t%f\n',subFolderNames{:},patientDice{:});
fclose(fid);
% output_loc = 'LiverOutput/';
% dicomlist = dir(fullfile(inputBrain,'*.dcm'));
% exit
% total_Im = numel(dicomlist);
% disp(total_Im); 
% for cnt = 20 : total_Im - 50
%     disp('cnt');
%     disp(cnt);
%     filename = fullfile(inputBrain,dicomlist(cnt).name);
%     img = dicomread(filename); 
%     img = resizeImage(img, 5/6);
%     img = adjustIm(img);
%     % disp(size(img));
%     segimg = fuzzyimage(img ,cnt, output_loc);
%     subplot(1,2,1);
%     imshow(img, []);
%     subplot(1,2,2);
%     imshow(segimg, []);
%     before_after_img = imtile({img,segimg});
%     imshow(before_after_img, []);
%     hold;

% ivals = 50:10:192;


function biggerIm = resizeImage(img, scaleFactor)
    sz = size(img);
    xg = 1:sz(1);
    yg = 1:sz(2);
    F = griddedInterpolant({xg,yg},double(img));
    xq = (0:scaleFactor:sz(1))';
    yq = (0:scaleFactor:sz(2))';
    biggerIm = uint8(F({xq, yq}));
   
end

function adjIm = adjustIm(img)
    adjIm = imadjust(img, [0.3, 1.0], [0.3, 1.0]);

    % imshow(adjIm,[]);
    %imcontrast
    %hold;

end


%im=imread ('/home/Debanjan_cse/MUSIG Function Modified/SkullStripped/72.png');
%im=rgb2gray(im);
% %image(inputImages[:, :, 1], C);
% output_loc='LiverOutput/';
% ivals = 50:10:192;
% ni = length(ivals);
% for K = 1 : ni
%     i = ivals(K);
%     im = V(:,:,i); %n is the slice number that you want to visualize.
%     %imshow(im,[])
%     
%     %imshow(im);
%     
% end
