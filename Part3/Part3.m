rgbImage = imread('tm_100.jpg'); %read image
grayImage = rgb2gray(rgbImage);
[M,N] = size(grayImage);
gaussian = zeros(1,9);
gm = zeros(M,N);
smoothImage = grayImage;
smoothImage = double(smoothImage);
sigma = 0.5;
threshold = 40;
%get gaussian
for i = 1:9
    ep = -(((i-5)^2)/(2*sigma^2));
    gaussian(i) = exp(ep);
end

gaussian = gaussian/sum(gaussian);
%smooth image by rows and then columns seperately
for i = 5:M-5
    for j = 5:N-5
        window = grayImage(i,j-4:j+4);
        smoothImage(i,j) = sum(double(window).*gaussian);
    end
end

for i = 5:M-5
    for j = 5:N-5
        window = grayImage(i-4:i+4,j);
        smoothImage(i,j) = sum(double(window).*gaussian');
    end
end
%get edges
for i = 2:M-1
    for j = 2:M-1
        gx = smoothImage(i,j+1) - smoothImage(i,j);
        gy = smoothImage(i+1,j) - smoothImage(i,j);
        gm(i,j) = sqrt((gx)^2 + (gy)^2);
        if(gm(i,j) >= threshold)
            gm(i,j) = 255;
        else
            gm(i,j) = 0;
        end
    end
end
gm = uint8(gm);
smoothImage = uint8(smoothImage);
figure;
subplot(1,3,1); imshow(grayImage);title('Grayscale Image');
subplot(1,3,2); imshow(smoothImage);title('Smoothed Image');
subplot(1,3,3); imshow(gm);title('Threshold Image');
imwrite(gm, 'edge.tiff')
%}