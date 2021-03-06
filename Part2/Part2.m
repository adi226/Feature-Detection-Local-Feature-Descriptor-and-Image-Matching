laplacian = readmatrix('Laplacian_filter.txt');
mean = readmatrix('Mean_filter.txt');
rgbImage = imread('IMG_20201007_125030-WangCenterLeft.jpg'); %read image
grayImage = rgb2gray(rgbImage); %convert to grayscale
[M,N] = size(grayImage);%get rows and columns
lapfiltImage = zeros(M,N);
meanfiltImage = zeros(M,N);
%find laplacian filtered and mean filtered image
for i = 3:M-2
    for j= 3:N-2
        window = grayImage(i-2:i+2, j-2:j+2);
        lapfiltImage(i,j) = sum(sum(double(window).*(laplacian)));
        meanfiltImage(i,j) = sum(sum(double(window).*(mean)));
    end
end
%lapfiltImage = uint8(lapfiltImage);
%meanfiltImage = uint8(meanfiltImage);
lapmin = min(min(lapfiltImage));
meanmin = min(min(meanfiltImage));
lapmax = max(max(lapfiltImage));
meanmax = max(max(meanfiltImage));

for i = 3:M-2
    for j = 3:N-2
        lapfiltImage(i,j) = (lapfiltImage(i,j) - lapmin) * (255/(lapmax - lapmin));
        meanfiltImage(i,j) = (meanfiltImage(i,j) - meanmin) * (255/(meanmax - meanmin));
    end
end
lapfiltImage = uint8(lapfiltImage);
meanfiltImage = uint8(meanfiltImage);
figure, imshow(lapfiltImage);
title ('Laplacian Filtered Image');
figure, imshow(meanfiltImage);
title ('Mean Filtered Image');
imwrite(lapfiltImage, 'lapfiltImage.tiff');
imwrite(meanfiltImage, 'meanfiltImage.tiff');
