%function fuzzy_image
function out_matrix = fuzzyimage(im, imgtruth, i, outdir, patient, lamb,cluster, max_levels)

wt=size(im,2); % dimention of matrix
ht=size(im,1);  % dimention of matrix
%out_matrix = zeros(ht,wt);    % new matrix intermediate matrix

pi=22./7;
%input_image = im;

outdir = append('results/', patient);


%int_matrix = zeros(ht,wt); % new matrix intermediate matrix
%out_matrix = zeros(ht,wt);    % new matrix intermediate matrix
%y = zeros(ht,wt);           % cardianlity matrix


in_matrix = im;

%inp =in_matrix;  % takng input gray image not doing anything


in_matrix = double (in_matrix)./255; % fuzzy_set
%in_matrix = im;
l=in_matrix;
in_matrix = double((pi/2)*in_matrix);

x = zeros (ht,wt);
w1 = zeros (ht,wt);
w2 = zeros (ht,wt);

a = 0;

[int_matrix,x] = Transit(in_matrix,ht,wt,lamb, cluster);   % I1 -> I2
int_matrix(isnan(int_matrix))=0;
[out_matrix,x] = Transit(int_matrix,ht,wt,lamb, cluster);  % I2 -> I3
out_matrix(isnan(out_matrix))=0;
[int_matrix,x] = Transit(out_matrix,ht,wt,lamb, cluster);  % I3 -> I2
int_matrix(isnan(int_matrix))=0;
[out_matrix,w1] = Transit(int_matrix,ht,wt,lamb, cluster); % I2 -> I3
%fprintf("max1 %f lamb\t %f\n", max(out_matrix(:)), lamb);
count=0;

out_matrix(isnan(out_matrix))=0;
while(a==0)
    %out_matrix = normalize(out_matrix);
    [int_matrix, x] = Transit(out_matrix,ht,wt,lamb, cluster);
    int_matrix(isnan(int_matrix))=0;
    [out_matrix, w2] = Transit(int_matrix,ht,wt,lamb, cluster); 
    a = check(w1,w2,ht,wt);
    w1 = w2; 
    count = count+1;
    %disp(a); 
    
    if(count>40)
        disp("a");
        
        break;
    end
end
%fprintf("max %f lamb\t %f\n", max(out_matrix(:)), lamb);
%disp(a); 
%subplot (1,3,1), imshow(im); title ('Original Image');
%subplot (1,3,2), imshow(l); title ('Gray Image');
%subplot (1,3,3), imshow(out_matrix); title ('segmented Image')
%figure, imshow(out_matrix,[])

out_matrix(isnan(out_matrix))=0;

%quant = 'max'; threshold = 0; 
%thresh = multithresh(out_matrix, max_levels);
%bestDiceImage = mat2gray(out_matrix >= thresh(end));

%quant8_I_min = out_matrix >= thresh(end); 
%quant8_I_max = out_matrix <= thresh(end - 1);
%quant8_I_max = out_matrix.*quant8_I_max;
%bestDiceImage = quant8_I_min.*quant8_I_max;  


%[bestDiceImage, bestDice, quant, threshold] = maximumDice(( out_matrix), imgtruth, max_levels);
%bestDice =  calculateDice(bestDiceImage, mat2gray(imgtruth));

% 
% RGB = label2rgb(imquantize(out_matrix, multithresh(out_matrix,max_levels)),'winter'); 
% 
% subplot (1,3,1), imshow(bestDiceImage, []); title ('Segmented Image');
% subplot (1,3,2), imshow(imgtruth, []); title ('Truth Image');
% subplot (1,3,3), imshow(RGB, []); title ('Segmented Multiple Levels Image');

% filename = append(outdir,'/anay.txt');
% fid = fopen(filename,'a+'); 
% fprintf(fid,'%d \t %f \n', i, bestDice);
% fclose(fid);


% input_image = mat2gray(im);
% rgbImage = imoverlay(input_image,mat2gray(bestDiceImage), 'red');
% filename = sprintf('%s/Comb_%d.png',outdir, i);
% imwrite( rgbImage,filename);
% 
% filename = sprintf('%s/Pred_%d.png',outdir, i);
% 
% imwrite(RGB,filename);
% 
% filename = sprintf('%s/liver_dice.txt',outdir);
% fid = fopen(filename,'a+');
% fprintf(fid,'%s\t%d\t%f\t%f\t\n',patient, i, lamb, bestDice);
% fclose(fid);

% filename = sprintf('%s/liver_cluster.txt',outdir);
% fid = fopen(filename,'a+');
% fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d \n',clusterListUsed(:, pos));
% fclose(fid);


end
