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

phase = "/Art.raw.nii.gz";
filename = sprintf('results/liver_dices.txt');
fid = fopen(filename,'w'); 
fprintf(fid,'paitentID\tDice\t\n');
fclose(fid);
for k = 1 : length(subFolderNames)
    patient = topLevelFolder+ "/"+subFolderNames{k};
	fprintf('Sub folder #%d = %s\n', k, patient);
    info = niftiinfo(patient + phase');
    art = niftiread(info);
    disp("PixelDimen");
    disp(info.PixelDimensions);
    
    truth = niftiread(patient + '/Truth.raw.nii.gz');
    diceList = zeros(size(length(art(1,1,:))));
    if info.BitsPerPixel == 16
        segnifit = zeros(size(art),'int16');
        truth = int16(truth);
    else
        segnifit = zeros(size(art),'single');
        truth = single(truth);
    end    
    info = niftiinfo(patient + phase');
    art = niftiread(info);
    outdir = append('results/', subFolderNames{k});
    mkdir(outdir);
    filename = sprintf('%s/test.nii',outdir);
    niftiwrite(art,filename, info, 'Version', 'NIfTI1',  'Compressed',true);
    for n = 1 : length(art(1,1,:))
        img = art(:,:,n);
        imgtruth = truth(:,:,n);
        filename = sprintf('%s/Before_%d.png',outdir, n );
        imwrite(uint8(img),filename); 
        if sum(imgtruth,'all') ~= 0
            segimg = fuzzyimage(img, imgtruth, n, outdir, subFolderNames{k} );
            %segimg = img;
    
            imga = double(segimg);
            imgb = double(imgtruth);
            filename = sprintf('%s/Truth_%d.png',outdir, n);
            imwrite(double(imgtruth),filename);
            diceList(n) = generalizedDice(imga,imgb);
            
            if info.BitsPerPixel == 16
                disp(info.BitsPerPixel);
                segnifit(:,:,n) = int16(segimg);
            else
                segnifit(:,:,n) = single(segimg);
            end   
        else
            segimg = fuzzyimage(img, imgtruth, n, outdir, subFolderNames{k} );
            if info.BitsPerPixel == 16
                %segnifit(:,:,n) = int16(segimg);
                segnifit(:,:,n) =  zeros(size(art(:,:,1)),'int16');
            else
                %segnifit(:,:,n) = single(segimg);
                segnifit(:,:,n) =  zeros(size(art(:,:,1)),'single');
            end   
        end
    end
    dice = generalizedDice(truth,segnifit);
    %filename = sprintf('%s/liver_dice.txt',outdir);
    filename = sprintf('results/liver_dices.txt');
    fid = fopen(filename,'a+'); 
    fprintf(fid,'%s\t%f\t\n', subFolderNames{k}, dice);
    fclose(fid);
    


    disp(size(segnifit));
    disp(size(art));
    filename = sprintf('%s/QIS.nii',outdir);
    disp(filename);
    %niftiwrite(segnifit,filename,info, 'Compressed',true);
    %info = niftiinfo(patient + '/Pre.raw.nii.gz');
    %[a, b, c ] =size(segnifit);
    %infoinfo.Filesize = (a*b*c)*31;

%     if info.BitsPerPixel == 16
%         segnifit = zeros(size(art),'int16');
%     else
%         segnifit = zeros(size(art),'single');
%     end 
    info = niftiinfo(patient + phase);
    niftiwrite(segnifit,filename, info, 'Version', 'NIfTI1',  'Compressed',true);
    patientDice(k) = mean(diceList);
end



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
    n = 2;  
    Idouble = im2double(img); 
    avg = mean2(Idouble);
    sigma = std2(Idouble);
    %adjIm = imadjust(img);
    %adjIm = imadjust(img,[avg-n*sigma avg+n*sigma],[]);
    %adjIm = histeq(adjIm);
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
