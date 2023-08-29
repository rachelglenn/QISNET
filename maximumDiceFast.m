function [bestDiceImage, bestDice, quant, theshold ] = maximumDiceFast(QIM,imgtruth,max_levels)
    % Threshold the various levels to get the maximum dice 

    thresh = multithresh(QIM, max_levels);
    %quantize = imquantize(QIM, thresh);
    
    %imwrite(mat2gray(quantize), "test3.png");
    %thresh = multithresh(QIM,thresh);
    thesh_range = thresh;
    valuesMax = [thesh_range max(QIM(:))];
    [quant8_I_max, index] = imquantize(QIM,thesh_range,valuesMax);
    valuesMin = [min(QIM(:)) thesh_range];
    quant8_I_min = valuesMin(index);
    thresh = double(thresh);
    
    bestDiceImage = quant8_I_max > thresh(end);
    bestDice = calculateDice(im2single(mat2gray(quant8_I_max)), im2single(mat2gray(imgtruth))); 
    fprintf("dice:%f\n", bestDice);
    %imshow(bestDiceImage);
    if bestDice <0.1
         bestDiceImage = quant8_I_max > thresh(end-1);
         bestDice = calculateDice(im2single(mat2gray(quant8_I_max)), im2single(mat2gray(imgtruth)));
    end
    if bestDice <0.1
         bestDiceImage = quant8_I_max > thresh(end-2);
         bestDice = calculateDice(im2single(mat2gray(quant8_I_max)), im2single(mat2gray(imgtruth)));
    end
    %imshow(bestDiceImage);
    if bestDice <0.1
         bestDiceImage = quant8_I_max > thresh(end-3);
         bestDice = calculateDice(im2single(mat2gray(quant8_I_max)), im2single(mat2gray(imgtruth)));
    end
    imshow(quant8_I_max);
    quant = 0;
    theshold = "max";
    
end
 function dice = calculateDice(seg, truth)
    seg_I = uint16(seg);

    B = im2bw(double(truth), graythresh(double(truth)));
    idx_img = find(B);
    idx_ref = find(seg_I);
    idx_inter = find(seg_I&B);
    dice = 2*length(idx_inter)/(length(idx_ref)+length(idx_img));
end
