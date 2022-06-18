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


inputBrain = 'DataBrainInput/A00028185/sub-A00028185_ses-NFB3_T1w_brain.nii.gz';
maskBrain = 'DataBrainInput/A00028185/sub-A00028185_ses-NFB3_T1w_brainmask.nii.gz';
output_loc = 'DataBrainOutput/A00028185';
V = niftiread(inputBrain);
Mask = niftiread(maskBrain);
% tool = imtool3D(V);
% tool.setMask(Mask);



%image(inputImages[:, :, 1], C);
output_loc='DataBrainOutput/';
ivals = 50:10:192;
ni = length(ivals);
imshow(V(:,:,50), []);
for K = 1 : ni
    i = ivals(K);
    im = V(:,:,i); %n is the slice number that you want to visualize.
    %imshow(im,[])
     segimg = fuzzyimage(im,i, output_loc);
     subplot(1,2,1);
     imshow(im, []);
     disp(max(im(:)));
     subplot(1,2,2);
     imshow(segimg, []);
     pause(12);
    %imshow(im);
    
end
