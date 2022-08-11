function [Z, w1]= Transit(I,ht,wt,q0,cluster)   % I is a input matrix, Ht is height of a matrix , wt is width of a matrixand q0 is lamda
pi=22./7;
%lamda=double (1./q0);
lamda= q0;
B=cardinality(I,ht,wt);
%A = ones(ht,wt);
A=(pi*2)*B;
%I = ones(ht,wt);
%A=(1.76*pi)*I;
alphas = 0;


%.............. otsu Method...................%

%cluster =[0 0.16 0.32 0.48 0.64 0.80 0.96 1]; % Set 1 for cluaster 8
% Used
% cluster =[0 0.14 0.28 0.42 0.56 0.70 0.90 1]; % Set 2 for cluaster 8
%cluster =[0 0.13 0.26 0.39 0.52 0.75 0.91 1]; % Set 3 for cluaster 8
%cluster =[0 0.18 0.36 0.54 0.72 0.90 0.97 1]; % Set 4 for cluaster 8
%cluster =[0 0.17 0.34 0.51 0.68 0.85 0.95 1]; % Set 5 for cluaster 8
%cluster = [0.15, 0.29, 0.43, 0.57, 0.71, 0.85, 1]; %RG from paper
% cluster = [0.15, 0.30, 0.46, 0.61, 0.76, 0.91, 1];
%.........Li.................%

% cluster = [0 0.62 0.71 0.80 0.83 0.93 0.98 1]; % Set 1 for cluster 8
%cluster = [0 0.61 0.69 0.77 0.82 0.83 0.97 1]; % Set 2 for cluster 8
%cluster = [0 0.60 0.72 0.82 0.94 0.96 0.98 1]; % Set 3 for cluster 8
%cluster = [0 0.63 0.74 0.79 0.82 0.88 0.97 1]; % Set 4 for cluster 8
%cluster = [0 0.70 0.74 0.79 0.82 0.88 0.98 1]; % Set 5 for cluster 8

cluster=(pi/2)*cluster;
d= size(cluster,2); % size of cluster

for a=1:ht

    for b=1:wt
        p=-1;
        sum=0.00;

        for i=1:d-1
            if (I(a,b) >=cluster(i)) && (I(a,b)<=cluster (i+1))
                %cluster (i)=(pi/2)*cluster(i);
                %cluster (i+1)=(pi/2)*cluster(i+1);

                alphas =double(B(a,b)/(cluster(i+1)-cluster(i))); % this is for a alpha calculation
                break;
            else
            end
        end
        %}

        for m=1:3
            c=-1;
            for n=1:3
                if(a+p<=0)||(a+p>=ht+1)||(b+c<=0)||(b+c>=wt+1)
                    sum=sum+0.00;
                else
                    w =(pi*2)*(pi/2 - (I(a+p,b+c) - I(a,b)));
                    z=(I(a+p,b+c))*cos(w-A(a,b));

                    y= 1/(alphas + exp(-(lamda)*(double(z)-B(a,b))));    % this a

                    sum=sum+y;
                end
                c=c+1;
            end
            p=p+1;
        end
        Z(a,b)=sum;
        w1(a,b)=w;
    end
end