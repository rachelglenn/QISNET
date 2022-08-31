%function fuzzy_image
function out_matrix = fuzzyimage(im, imgtruth, i, outdir, patient)
%im=imread ('/home/Debanjan_cse/MUSIG Function Modified/SkullStripped/72.png');
%im=rgb2gray(im);
%im =imread(str);
%q0=str2;
%for q0 = 0.1:0.9 step 0.1
%for q0 = 1.0:-0.01:0.01
filename = sprintf('%s/liver_output_%s.txt',outdir, patient);
fileID = fopen(filename,'a+');
 
fprintf(fileID,'%6s   %12s\n','q0','corelation'); % this for naming the variables
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
%[0 0.13 0.26 0.39 0.52 0.75 0.91 1],...
%[0 0.18 0.36 0.54 0.72 0.90 0.97 1],...
% [0 0.17 0.34 0.51 0.68 0.85 0.95 1],...
%  [0.15, 0.29, 0.43, 0.57, 0.71, 0.85, 1],...
%  [0.15, 0.30, 0.46, 0.61, 0.76, 0.91, 1],...
%[0 0.62 0.71 0.80 0.83 0.93 0.98 1],...
% [0 0.61 0.69 0.77 0.82 0.83 0.97 1],...
%  [0 0.60 0.72 0.82 0.94 0.96 0.98 1],...
%  [0 0.63 0.74 0.79 0.82 0.88 0.97 1],...
% [0 0.70 0.74 0.79 0.82 0.88 0.98 1]};

lamblist = linspace(-0.1,0.5,2); % 0.25:-0.001:0.24
lamblist = [-0.5];
numtrys = length(clusterlist)*length(lamblist);
%lamblist = linspace(-0.1,0.1,10);


images = zeros(ht, wt,numtrys );
diceList = zeros(size(numtrys));
k= 1;
j = 1;
lamblistUsed = zeros(size(numtrys));
clusterListUsed = zeros(8, numtrys);
previous_N = 1000.0;
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

    %filename = sprintf('DataBrainOutput/A00028185/%d_C8_S2_%f.png',i, lamb);
    %filename = sprintf('LiverOutput/%d_C8_S2_%f.png',i, lamb);
    %out_matrix=out_matrix*(2/pi);
    %imwrite(out_matrix,strcat(output_loc,filename));

    
    % B = out_matrix;
    %B=B-min(B(:)); % shift data such that the smallest element of A is 0
    %B=B/max(B(:)); % normalize the shifted data to 1   before_after_img = imtile({A,out_matrix});

    % Brain
    % before_after_img = imtile({A*255,B});
    %Liver
    %before_after_img = imtile({im*255,out_matrix});
    % before_after_img = imtile({A,B});

    out = out_matrix;  %output matrix
    r = corr2(inp,out);

    ka = [lamb  r ];

    %fprintf(fileID,'%6.8f   %12.8f\n',ka);    
    %A = im2bw(double(out_matrix), graythresh(out_matrix));
    thresh = multithresh(out_matrix, 3);
    %disp("thresh");
    %disp(thresh);
    %disp(isnan( A ));
    out_matrix(isnan(out_matrix))=0;
    %seg_I = imquantize(A, thresh);
    %RGB = label2rgb(seg_I); 	 
    %figure;
    seg_I = (out_matrix >=thresh(3));
    N = double(nnz( seg_I ));
    %fprintf("on image %f and nonzero pixels %f(i) %f \n", i, N,previous_N);
    if N < 1000.0
        %fprintf("on image %f and nonzero pixels(2) %f \n", i, N);
        seg_I = (out_matrix >=thresh(2));
        N = double(nnz( seg_I));
        if N < 1000.0
            %fprintf("on image %f and nonzero pixels(1) %f \n", i, N);
            seg_I = (out_matrix >=thresh(1));
            N = double(nnz( seg_I));  
        end
    end
    %fprintf("on image %f and nonzero pixels %f(f) \n", i, N)
%     previous_N = N;
%     if previous_N < 1000.0
%         fprintf("on image %d and nonzero pixels %d(0) \n", i, N);
%         previous_N = 1000.0;
%     end
    seg_I = uint16(seg_I);

    B = im2bw(double(imgtruth), graythresh(double(imgtruth)));
    idx_img = find(B);
    idx_ref = find(seg_I);
    idx_inter = find(seg_I&B);
    dist = 2*length(idx_inter)/(length(idx_ref)+length(idx_img));

    %Z =generalizedDice(imga,imgb);
    diceList(n) = dist;
    %disp("DiceList");
    %disp(diceList(n));
    images(:,:,n) = seg_I;
    lamblistUsed(n) = lamb;
    clusterListUsed(:,n) = cluster;
    %imshow(A, []);
end
[value, pos] = max(diceList);
%disp("Max Dice:");
%disp(value);


A = images(:,:,pos);




imshow(A)
axis off
title('RGB Segmented Image')
B = im2bw(double(A), graythresh(double(A)));
%imshow(A, []);
%A=A-min(A(:)); % shift data such that the smallest element of A is 0
%A=A/max(A(:)); % normalize the shifted data to 1     
filename = sprintf('%s/Pred_%d.png',outdir, i);
imwrite(A,filename);
R   = 1;  % Value in range [0, 1]
G   = 0;
B1   = 0;
input_image = mat2gray(input_image);
RGB = cat(3, B * R, B * G, B * B1);
out = imfuse(RGB,input_image,'blend','Scaling','joint');
%rgbImage = cat(3, B, A, A);
input_image = imadjust(mat2gray(input_image));
rgbImage = imoverlay(input_image,mat2gray(A), 'red');
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

% this for printing values
fclose(fileID);

 
%end

end

    
