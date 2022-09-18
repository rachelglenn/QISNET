function dice = calculateDice(seg, truth)
    seg_I = uint16(seg);

    %B = im2bw(double(truth), graythresh(double(truth)));
    idx_img = find(truth);
    idx_ref = find(seg_I);
    idx_inter = find(seg_I&truth);
    dice = 2*length(idx_inter)/(length(idx_ref)+length(idx_img));
%     try
%         dice = generalizedDice(double(seg), double(truth));
%     catch
%         dice = 0;
%     end
    disp("dice function");
    disp(dice);
end

