%function fuzzy_image
function out_matrix = fuzzyimage(im, imgtruth, i, outdir, patient)
%im=imread ('/home/Debanjan_cse/MUSIG Function Modified/SkullStripped/72.png');
%im=rgb2gray(im);
%im =imread(str);
%q0=str2;
%for q0 = 0.1:0.9 step 0.1
%for q0 = 1.0:-0.01:0.01
%filename = sprintf('%s/liver_output_%s.txt',outdir, patient);
%fileID = fopen(filename,'a+');
 
%fprintf(fileID,'%6s   %12s\n','q0','corelation'); % this for naming the variables
wt=size(im,2); % dimention of matrix
ht=size(im,1);  % dimention of matrix
out_matrix = zeros(ht,wt);    % new matrix intermediate matrix
%for lamb = 1:10% for lambda = 0.23 to 0.24
%for lamb = 0.1:0.01:0.2
pi=22./7;
input_image = im;


%.............. otsu Method...................%

%cluster =[0 0.16 0.32 0.48 0.64 0.80 0.96 1]; % Set 1 for cluaster 8
% [0 1.0/8 2.0/8.0 3.0/8.0 4.0/8.0 5.0/8.0 6.0/8.0 7.0/8.0],...
% Used
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

lamblist = linspace(-0.1,0.5,2); % 0.25:-0.001:0.24
lamblist = [-0.5];
numtrys = length(clusterlist)*length(lamblist);

images = zeros(ht, wt,numtrys );
diceList = zeros(size(numtrys));
k= 1;
j = 1;
lamblistUsed = zeros(size(numtrys));
clusterListUsed = zeros(8, numtrys);
outdir = append('results/', patient);



for n = 1 : numtrys
   
    if k > length(lamblist)
        k = 1;
        j =+ 1;
    end
    lamb = lamblist(k);
    cluster = clusterlist{j};
    k =+ 1;
   

    %lamb = -0.1;

    %q0=0.2
    %for q0 = 0.1:0.01:0.2      
    

    %in_matrix = imresize (rgb2gray(im),[128,128]);  % just for resize
    %in_matrix =im;



    int_matrix = zeros(ht,wt); % new matrix intermediate matrix
    out_matrix = zeros(ht,wt);    % new matrix intermediate matrix
    y = zeros(ht,wt);           % cardianlity matrix

    %in_matrix = rgb2gray(im);
    in_matrix = im;

    inp =in_matrix;  % takng input gray image not doing anything


    in_matrix = double (in_matrix)./255; % fuzzy_set
    l=in_matrix;
    in_matrix = double((pi/2)*in_matrix);

    x = zeros (ht,wt);
    w1 = zeros (ht,wt);
    w2 = zeros (ht,wt);

    a = 0;
    [int_matrix,x] = Transit(in_matrix,ht,wt,lamb, cluster);   % I1 -> I2
    [out_matrix,x] = Transit(int_matrix,ht,wt,lamb, cluster);  % I2 -> I3
    [int_matrix,x] = Transit(out_matrix,ht,wt,lamb, cluster);  % I3 -> I2
    [out_matrix,w1] = Transit(int_matrix,ht,wt,lamb, cluster); % I2 -> I3

    count=0;
    
    while(a==0)
        [int_matrix, x] = Transit(out_matrix,ht,wt,lamb, cluster);
        [out_matrix, w2] = Transit(int_matrix,ht,wt,lamb, cluster); 
        a = check(w1,w2,ht,wt);
        w1 = w2; 
        count = count+1;
        %disp(count); 
       
        if(count>40)
            disp("a");
            
            break;
        end
    end

    %subplot (1,3,1), imshow(im); title ('Original Image');
    %subplot (1,3,2), imshow(l); title ('Gray Image');
    %subplot (1,3,3), imshow(out_matrix); title ('segmented Image')
    %figure, imshow(out_matrix,[])

    out_matrix(isnan(out_matrix))=0;
	 

    max_levels = 5;
    %figure(1);
    [bestDiceImage, bestDice, quant, threshold] = maximumDice(out_matrix, imgtruth, max_levels);
    
    diceList(n) = bestDice;


    images(:,:,n) = bestDiceImage;
    lamblistUsed(n) = lamb;
    clusterListUsed(:,n) = cluster;
    %disp(quant);
    filename = append(outdir,'/anay.txt');
    fid = fopen(filename,'a+'); 
    fprintf(fid,'%d \t %d \t %s \t %f \n', i, quant, threshold, bestDice);
    fclose(fid);
    
end

RGB = label2rgb(imquantize(out_matrix, multithresh(out_matrix,max_levels)),'winter'); 
filename = sprintf('%s/Pred_%d.png',outdir, i);
imwrite(RGB,filename);

[value, pos] = max(diceList);
out_matrix = images(:,:,pos);

input_image = mat2gray(input_image);
rgbImage = imoverlay(input_image,mat2gray(out_matrix), 'red');
filename = sprintf('%s/Comb_%d.png',outdir, i);
imwrite( rgbImage,filename);


filename = sprintf('%s/liver_dice.txt',outdir);
fid = fopen(filename,'a+');
fprintf(fid,'%s\t%d\t%f\t%f\t\n',patient, i, lamblistUsed(pos), diceList(pos));
fclose(fid);

filename = sprintf('%s/liver_cluster.txt',outdir);
fid = fopen(filename,'a+');
fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d \n',clusterListUsed(:, pos));
fclose(fid);


end

function dice = calculateDice(seg, truth)
    seg_I = uint16(seg);

    B = im2bw(double(truth), graythresh(double(truth)));
    idx_img = find(B);
    idx_ref = find(seg_I);
    idx_inter = find(seg_I&B);
    dice = 2*length(idx_inter)/(length(idx_ref)+length(idx_img));
end
