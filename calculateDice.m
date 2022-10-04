function dice = calculateDice(seg, truth)
    seg_I = seg;

    %B = im2bw(double(truth), graythresh(double(truth)));
    idx_img = find(truth);
    idx_ref = find(seg_I);
    
    if(any(isnan(idx_ref(:))) || any(isnan(idx_img(:))))
        dice = 0.0;
        disp("issue with calculating the dice");
    else
        idx_inter = find(uint16(seg_I) & uint16(truth));
        dice = 2*length(idx_inter)/(length(idx_ref)+length(idx_img));
        %dice = jaccard(double(seg), double(truth));
        %fprintf("calculate Dice: %f \n", dice);
    end

    
end

