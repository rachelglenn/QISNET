
function y= cardinality(x,ht,wt)
%x=[1 2 5;2 45 7;23 5 7;5 2 9;1 2 45];
%wt=size(x,2); % dimention of matrix
%ht=size(x,1);  % dimention of matrix
for a=1:ht
    for b=1:wt
        p=-1;
        sum=0.0;
        %disp("\n\n\n");
        for m=1:3
            c=-1;
            for n=1:3 %RG added |(p==c==0)
                if(a+p<=0)||(a+p>=ht+1)||(b+c<=0)||(b+c>=wt+1) ||(p==0 && c==0)
                    sum=sum+0;
                else
                    sum=sum+x(a+p,b+c);
                    %%if sum > 1000000.00
                    %fprintf("x(%d,%d)=%f\n", a+p, b+c, x(a+p, b+c));
                    %fprintf("p=%d c=%d\n", p, c);
                    %end

                end
                c=c+1;
            end
            p=p+1;
        end
        y(a,b)=sum;

    end
end
end
