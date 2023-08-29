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


%in_matrix = double (in_matrix)./255; % fuzzy_set (RG Removed)
%in_matrix = im;
%l=in_matrix;
%in_matrix = double((pi/2)*in_matrix);

x = zeros (ht,wt);
w1 = zeros (ht,wt);
w2 = zeros (ht,wt);

a = 0;

[int_matrix,x] = Transit(in_matrix,ht,wt,lamb, cluster);   % I1 -> I2
%int_matrix(isnan(int_matrix))=0;
[out_matrix,x] = Transit(int_matrix,ht,wt,lamb, cluster);  % I2 -> I3
%out_matrix(isnan(out_matrix))=0;
%[int_matrix,x] = Transit(out_matrix,ht,wt,lamb, cluster);  % I3 -> I2
%int_matrix(isnan(int_matrix))=0;
[out_matrix,w1] = Transit(int_matrix,ht,wt,lamb, cluster); % I2 -> I3
%fprintf("max1 %f lamb\t %f\n", max(out_matrix(:)), lamb);
count=0;

%out_matrix(isnan(out_matrix))=0;
%disp("RG insert to troubleshoot");
%a= 1;

while(a==0)
    %out_matrix = normalize(out_matrix);
    [int_matrix, x] = Transit(out_matrix,ht,wt,lamb, cluster);
    %int_matrix(isnan(int_matrix))=0;
    [out_matrix, w2] = Transit(int_matrix,ht,wt,lamb, cluster); 
    a = check(w1,w2,ht,wt);
    w1 = w2; 
    count = count+1;
    %disp(a); 
    
    if(count>5)
        disp("a");
        
        break;
    end
end
%fprintf("max %f lamb\t %f\n", max(out_matrix(:)), lamb);

%out_matrix(isnan(out_matrix))=0;




end
