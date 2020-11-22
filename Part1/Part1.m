rgbImage = imread('IMG_20201007_125030-WangCenterLeft.jpg'); %read image
grayImage = rgb2gray(rgbImage); %convert to grayscale
[M,N] = size(grayImage);%get rows and columns
hist = zeros(1, 256);%initalize histogram
disprob = zeros(1, 256);%initialize discrete probability
cumdf = zeros(1,256);%initialize cummulative distribution 
histeqImage = uint8(zeros(M,N));%initialize and cast output image
%get histogram
for i = 1:M
    for j= 1:N
        k = grayImage(i,j);
        k = k + 1;%map 0 to 255 to 1 to 256
        hist(k) = hist(k) + 1;
    end
end
%get discrete probability
disprob = (hist / (M*N));

%get cummulative distribution function
cumdf(1) = disprob(1);
for k = 2:256
    cumdf(k) = disprob(k) + cumdf(k-1);
end
%get histogram equilization output
for i = 1:M
    for j = 1:N
        histeqImage(i,j)= 255 * cumdf(grayImage(i,j)+1);
    end
end
figure, imshow(grayImage);
title ('Grayscale Image');
figure, imshow(histeqImage);
title ('Histogram Equalization Image');
imwrite(grayImage, 'grayscaleImage.tiff');
imwrite(histeqImage, 'histeqImage.tiff');