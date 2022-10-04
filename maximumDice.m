function [bestDiceImage, bestDice, quant, theshold ] = maximumDice(QIM,imgtruth,max_levels)
    % Threshold the various levels to get the maximum dice 
    %QIM = imcomplement( mat2gray(QIM));
    thresh = multithresh(QIM, max_levels);
    %quantize = imquantize(QIM, thresh);
    
    %imwrite(mat2gray(quantize), "test3.png");
    [row_I,col_I,temp] = size(imgtruth);
    images_tmp = zeros(row_I, col_I, temp, max_levels*3);
    distList = zeros(1,max_levels*3);
    thresholdlist = zeros(1,max_levels*3);
    maxminlist = zeros(1, max_levels*3);

    %figure(5);
    %QIM = sgolayfilt(QIM, 2, 11, [], 1);
    %QIM = imresize(QIM, [row_I, col_I], 'nearest');
    for i=0:1:max_levels-2
%         min_level = i;
%         thresh = multithresh(QIM,max_levels);
%         thresh(end +1 ) = max(QIM(:)) + 0.4;
%         % RG modified 9_14_21
%         
        quant8_I_min = QIM >= thresh(end - (i));
      
        quant8_I_max = QIM <= thresh(end - i);
        quant8_I_max = QIM.*quant8_I_max;
        quant8_I_max = quant8_I_min.*quant8_I_max;
        %max and min list
        thresholdlist(i +1) = 0;
        %threshold
        maxminlist(i +1) = i;

%       
%         RGB = label2rgb(imquantize(QIM, multithresh(QIM, 10))); 
    	    RGB = quant8_I_max;
%         RGB = imbinarize(RGB);
        %imshow(RGB);
        pause(0.2);
        quant8_I_max = QIM > thresh(end -i);
        images_tmp(:,:,:, i +1) = quant8_I_max; %*ones(size(quant8_I_max));
        dice =  calculateDice(mat2gray(quant8_I_max), mat2gray(imgtruth));


        distList(i+1) = dice;
    end
% 
%     for i=0:1:max_levels-2
%         min_level = i;
%         thresh = multithresh(QIM,max_levels);
%         % RG modified 9_14_21
%         
%         quant8_I_min = QIM >= thresh(end - (i+1));
%         %quant8_I_min = (QIM).*quant8_I_min;
%         quant8_I_max = QIM <= thresh(end - i);
%         %quant8_I_max = QIM.*quant8_I_max;
%         quant8_I_max = quant8_I_min.*quant8_I_max; %*ones(size(quant8_I_max));
%         %max and min list
%         thresholdlist(i + max_levels-2) = 0;
%         %threshold
%         maxminlist(i+ max_levels-2) = i;
% 
% %       
% %         RGB = label2rgb(imquantize(QIM, multithresh(QIM, 10))); 
% %    	    RGB = quant8_I_max;
% %         RGB = imbinarize(RGB);
%         %imshow(RGB);
%         images_tmp(:,:,i+ max_levels-2 ) = quant8_I_max;
%         try
%            distList(i+ max_levels-2) = jaccard(double(quant8_I_max), double(imgtruth));
%         catch
%              distList(i+ max_levels-2) = 0;
%              
%         end
%    end



    [bestDice, pos] = max(distList);
    
    if thresholdlist(pos) == 0
        theshold = 'max';
    else
        theshold = 'min';
    end
    quant = maxminlist(pos);
    %fprintf("Dice max:%f\n", bestDice);
    %fprintf("Threshold cut:%d\n", quant);
    bestDiceImage = mat2gray(images_tmp(:,:,:,pos));
    %disp("size");
    %disp(size(bestDiceImage));
end


