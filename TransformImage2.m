function [outputImage] = TransformImage2(grayImage, freq, cutoffTime)
        outputImage = zeros(size(grayImage));
        %grayImage = imadjustn(grayImage);
        LevelsImg = unique(grayImage(:));       % determine number of unique code values
        disp(['image K: ',num2str(length(LevelsImg)),' distinct levels']);
        disp(['max level = ' num2str( max(LevelsImg) )]);
        disp(['min level = ' num2str( min(LevelsImg) )]);
        %Kdouble = double(rgrayImageRef);                  % cast uint16 to double
        %kmult = 65535/(max(max(Kdouble(:)))); % full range multiplier
        
        %Ref = uint16(kmult*Kdouble);   % full range 16-bit reference image
        %N = length(LevelsImg);     % number of unique 16-bit code values in image A.
        %figure(10)
        for j = 1 :length(grayImage(1,1,:))
            %B16bitUniq = imhistmatch(grayImage(:,:,n),Ref,N);
            %B16bitUniq  = TransformHistogram(grayImage, Ref, false);
            %B16bitUniq = imcontrast(grayImage(:,:,n), 0.2, 0.7);
              
            Kdouble = double(grayImage(:,:,j));                  % cast uint16 to double
            kmult = 65535/(max(max(Kdouble(:)))); % full range multiplier
            
            Ref = uint16(kmult*Kdouble);   % full range 16-bit reference image
            B16bitUniq = imadjust(Ref,stretchlim(Ref),[]);
            refImage=B16bitUniq;
%                         
%             % Saving the size of the input_image in pixels-
%             % M : no of rows (height of the image)
%             % N : no of columns (width of the image)
%             [M, N] = size(refImage);
%               
%             % Getting Fourier Transform of the input_image
%             % using MATLAB library function fft2 (2D fast fourier transform)
%             FT_img = fft2(double(grayImage(:,:,j)));
%               
%             % Assign the order value
%             n = 5; % one can change this value accordingly
%               
%             % Assign Cut-off Frequency
%             D0 = freq; % one can change this value accordingly
%               
%             % Designing filter
%             u = 0:(M-1);
%             v = 0:(N-1);
%             idx = find(u > M/2);
%             u(idx) = u(idx) - M;
%             idy = find(v > N/2);
%             v(idy) = v(idy) - N;
%               
%             % MATLAB library function meshgrid(v, u) returns 
%             % 2D grid which contains the coordinates of vectors 
%             % v and u. Matrix V with each row is a copy of v 
%             % and matrix U with each column is a copy of u 
%             [V, U] = meshgrid(v, u);
%               
%             % Calculating Euclidean Distance
%             D = sqrt(U.^2 + V.^2);
%               
%             % determining the filtering mask
%             H = 1./(1 + (D./D0).^(2*n));
%               
%             % Convolution between the Fourier Transformed 
%             % image and the mask
%             G = H.*FT_img;
%               
%             % Getting the resultant image by Inverse Fourier Transform 
%             % of the convoluted image using MATLAB library function  
%             % ifft2 (2D inverse fast fourier transform)   
%             I2 = real(ifft2(double(G))); 
% %                 
% %             % Displaying Input Image and Output Image 
% %             subplot(3, 1, 1), imshow(mat2gray(grayImage(:,:,j))), 
% %             subplot(3, 1, 2), imshow(mat2gray(I2)),
% %             
% %            
% %             subplot(3, 1, 3),
                 %imshow(mat2gray(I2))
            %pause(3);


            Iw = dct2(refImage); 


            %Reduce peaky high-frequency components
            idx = abs(Iw) > freq;% freq;
            idx(1:cutoffTime,:) = false;
            idx(:,1:cutoffTime) = false;
            Iw(idx) = Iw(idx)/100;
            % Convert to image
            I2 = idct2(Iw);
            %imshow(mat2gray(Iw), []);
            %I2 = mat2gray(I2);
            
            %imshow(mat2gray(I2), []);
%             set(gcf,'Position',[100 100 5000 5000])
%             subplot(1,3,1)
%             imshow(mat2gray(B16bitUniq))
%             title('B16bit64: 64 bins')
%             subplot(1,3,2)
%             imshow(mat2gray(grayImage(:,:,n)))
%             title(['original: ',num2str(N),' bins']) 
%             subplot(1,3,3)
%             imshow(mat2gray(Ref))
%             title(['Reference'])
             
            %imshow(B16bitUniq);
             %hold on
             outputImage(:,:,j) = I2;
             %outputImage(:,:,n) = grayImage(:,:,n);

        end

end