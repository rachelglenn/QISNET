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


inputBrain = 'LiverInput/13.000000-t1vibeqfstrap2bhFIL-72776';
output_loc = 'LiverOutput/';
dicomlist = dir(fullfile(inputBrain,'*.dcm'));
total_Im = numel(dicomlist);
disp(total_Im); 
for cnt = 20 : total_Im - 50
    disp('cnt');
    disp(cnt);
    filename = fullfile(inputBrain,dicomlist(cnt).name);
    img = dicomread(filename); 
    sz = size(img);
    xg = 1:sz(1);
    yg = 1:sz(2);
    F = griddedInterpolant({xg,yg},double(img));
    xq = (0:3/6:sz(1))';
    yq = (0:3/6:sz(2))';
    img = uint8(F({xq, yq}));
    % disp(size(img));
    segimg = fuzzyimage(img ,cnt, output_loc);
    % subplot(1,2,1);
    % imshow(img, []);
    % subplot(1,2,2);
    % imshow(segimg, []);
    % before_after_img = imtile({img,segimg});
    % imshow(before_after_img, []);

% ivals = 50:10:192;

end

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
