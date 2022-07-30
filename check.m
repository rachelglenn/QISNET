function a= check(w1,w2,ht,wt)
for p=1:ht
    for q=1:wt
        x=w2(p,q)-w1(p,q);
        %if(x < 0.00001)   % changes made
 	    %if(x ~= 0.1)   % changes made
        if(x < 0.1)   % changes made
           a=1;
           %printf('changed the points: %f %f', p, q);
        else
           a=0;
           break;
       end
    end
end
