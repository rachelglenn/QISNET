%function fuzzy_image
function out_matrix = fuzzyimage(im, imgtruth, i, outdir, patient)
%im=imread ('/home/Debanjan_cse/MUSIG Function Modified/SkullStripped/72.png');
%im=rgb2gray(im);
%im =imread(str);
%q0=str2;
%for q0 = 0.1:0.9 step 0.1
%for q0 = 1.0:-0.01:0.01
 
fileID = fopen('output.txt','w');
 
fprintf(fileID,'%6s   %12s\n','q0','corelation'); % this for naming the variables
wt=size(im,2); % dimention of matrix
ht=size(im,1);  % dimention of matrix
out_matrix = zeros(ht,wt);    % new matrix intermediate matrix
%for lamb = 1:10% for lambda = 0.23 to 0.24
%for lamb = 0.1:0.01:0.2
pi=22./7;
% B = im;

%lamblist = linspace(-0.1,0.1,10);
lamblist = linspace(-0.1,0.1,10); % 0.25:-0.001:0.24

images = zeros(ht, wt, length(lamblist));
diceList = zeros(size(length(lamblist)));

for n = 1 : length(lamblist)
    lamb = lamblist(n);
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
    %l=in_matrix;
    in_matrix = double((pi/2)*in_matrix);

    x = zeros (ht,wt);
    w1 = zeros (ht,wt);
    w2 = zeros (ht,wt);

    a = 0;
    [int_matrix,x] = Transit(in_matrix,ht,wt,lamb);   % I1 -> I2
    [out_matrix,x] = Transit(int_matrix,ht,wt,lamb);  % I2 -> I3
    [int_matrix,x] = Transit(out_matrix,ht,wt,lamb);  % I3 -> I2
    [out_matrix,w1] = Transit(int_matrix,ht,wt,lamb); % I2 -> I3

    count=0;
    
    while(a==0)
        [int_matrix, x] = Transit(out_matrix,ht,wt,lamb);
        [out_matrix, w2] = Transit(int_matrix,ht,wt,lamb); 
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
    %figure, imshow(out_matrix,[]);

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

    fprintf(fileID,'%6.8f   %12.8f\n',ka);    
    A = im2bw(double(out_matrix), graythresh(out_matrix));
    B = im2bw(double(imgtruth), graythresh(double(imgtruth)));
    idx_img = find(B);
    idx_ref = find(A);
    idx_inter = find(A&B);
    dist = 2*length(idx_inter)/(length(idx_ref)+length(idx_img));

    %Z =generalizedDice(imga,imgb);
    diceList(n) = dist;
    disp("DiceList");
    disp(diceList(n));
    images(:,:,n) = A;
    %imshow(A, []);
end
[value, pos] = max(diceList);
disp("Max Dice:");
disp(value);


A = images(:,:,pos);
B = im2bw(double(imgtruth), graythresh(double(imgtruth)));
%imshow(A, []);
%A=A-min(A(:)); % shift data such that the smallest element of A is 0
%A=A/max(A(:)); % normalize the shifted data to 1     
filename = sprintf('%s/Pred_%d_%s.png',outdir, i,patient);
imwrite(double(A),filename);

%rgbImage = cat(3, B, A, A);
rgbImage = imoverlay(A, B, [1 0 0]);
filename = sprintf('%s/Comb_%d_%s.png',outdir, i,patient);
imwrite(rgbImage,filename);

fid = fopen('results/LiverDice_lambda.txt','w');
fprintf(fid,'%s\t%f\t%f\n',patient, lamb(pos), diceList(pos));
fclose(fid);

% this for printing values
fclose(fileID);

 
%end

end

    