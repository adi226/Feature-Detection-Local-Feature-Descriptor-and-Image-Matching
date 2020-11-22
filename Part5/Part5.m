%read the images
rgbImage1 = imread('IMG_20201007_125027-WangCenterRight.jpg'); %read image
rgbImage2 = imread('IMG_20201007_125030-WangCenterLeft.jpg');

%get feature vectors
fv1 = getCorners(rgbImage1);
fv2 = getCorners(rgbImage2);

%find affine transformation matrix
Amat = [346,123,1,0,0,0;0,0,0,346,123,1;321,144,1,0,0,0;0,0,0,321,144,1;
        378,211,1,0,0,0;0,0,0,378,211,1;340,184,1,0,0,0;0,0,0,340,184,1;
        341,179,1,0,0,0;0,0,0,341,179,1;321,174,1,0,0,0;0,0,0,321,174,1;
        344,172,1,0,0,0;0,0,0,344,172,1;335,155,1,0,0,0;0,0,0,335,155,1;
        321,149,1,0,0,0;0,0,0,321,149,1;329,205,1,0,0,0;0,0,0,329,205,1;
        322,197,1,0,0,0;0,0,0,322,197,1;258,189,1,0,0,0;0,0,0,258,189,1;
        231,192,1,0,0,0;0,0,0,231,192,1;239,258,1,0,0,0;0,0,0,239,258,1;
        235,255,1,0,0,0;0,0,0,235,255,1;238,251,1,0,0,0;0,0,0,238,251,1];
bmat = [180;111;155;131;209;201;175;176;177;169;153;161;180;161;169;143;154;
     136;161;195;155;187;78;174;46;178;52;250;50;244;54;240];
xmat = Amat\bmat;
disp(xmat);

%Verification for affine transformation matrix, point (227,255) in left
%image should map to (38,248) in right image
c = [227;255;1];
x = [1.1426,-0.0610,-204.1115;0.0290,1.0526,-30.0016];
xp = round(x * c);
disp(xp);

function featurevectors = getCorners(rgbImage)
    grayImage = rgb2gray(rgbImage);
    [M,N] = size(grayImage);
    A = zeros(M,N);
    B = zeros(M,N);
    C = zeros(M,N);
    R = zeros(M,N);
    corner = zeros(M,N);
    superimpose = rgbImage;
    %get the gaussian
    gauss1 = getGaussian(9,2);
    gauss2 = getGaussian(11,5.5);
    smoothImage = double(grayImage);
    %smooth the image
    smoothImage = smooth(smoothImage,gauss1,M,N);
    threshold = 10000;
    %calculate gradients
    for i = 2:M-1
        for j = 2:N-1
            gx = smoothImage(i,j+1) - smoothImage(i,j);
            gy = smoothImage(i+1,j) - smoothImage(i,j);
            A(i,j) = gx^2;
            B(i,j) = gy^2;
            C(i,j) = gx * gy;
        end
    end
    A = smooth(A,gauss2,M,N);
    B = smooth(B,gauss2,M,N);
    C = smooth(C,gauss2,M,N);
    precorpix = [];
    corpix = [];
    for i = 6:M-6
        for j = 6:N-6
            Mmat = [A(i,j),C(i,j);C(i,j),B(i,j)];
            Rval = det(Mmat)-0.04*(trace(Mmat))^2;
            if(Rval > threshold)
                R(i,j) = Rval;
                precorpix = [precorpix;[i,j]];
            end
        end
    end
%get corner indices
    for i = 1:length(precorpix)
       pt = precorpix(i,:);
       if(pt(1)>6)&&(pt(1)<M-6)&&(pt(2)>6)&&(pt(2)<N-6) % can be removed
           window = R(pt(1)-1:pt(1)+1,pt(2)-1:pt(2)+1);
           if(R(pt(1),pt(2)) == max(window(:)))
               corpix = [corpix;[pt(1), pt(2)]];
               corner(pt(1)-1:pt(1)+1,pt(2)-1:pt(2)+1) = 255;
           end   
       end
    end
%superimpose corner on image
    for i=1:M
        for j=1:N
            if(corner(i,j) == 255)
                superimpose(i,j,2) = 0;
                superimpose(i,j,3) = 0;
            end
        end
    end
    corner = uint8(corner);
    %{
    figure;
    subplot(1,3,1); imshow(rgbImage);title('RGB Image');
    subplot(1,3,2); imshow(superimpose);title('Superimposed Image');
    subplot(1,3,3); imshow(corner);title('Corners Image');
%}
    figure;
    imshow(superimpose);title('Superimposed Image');
    featurevectors = [];
    %get feature vectors
    for i = 1:length(corpix)
        pt = corpix(i,:);
        if(pt(1)>4)&&(pt(1)<M-4)&&(pt(2)>4)&&(pt(2)<N-4)
            sift = zeros(1,8);
            winx = A(pt(1)-4:pt(1)+4, pt(2)-4:pt(2)+4);
            winy = B(pt(1)-4:pt(1)+4, pt(2)-4:pt(2)+4);
            gmag = sqrt(winx + winy);
            gangle = rad2deg(atan(winy./winx));       
            gangle = wrapTo360(gangle);
            gangle = round(gangle/45)*45;
            for j=1:9
                for k=1:9
                    ang = (gangle(j,k)/45)+1;
                    if (ang == 360)
                        ang = 1;
                    end
                    sift(ang) = sift(ang)+ gmag(j,k);
                end
            end
            mxid = find(sift == max(sift(:)));
            sift = circshift(sift, abs(5-mxid));
            sift = round(sift);
            featurevectors = [featurevectors;{pt,sift}];
        end
    end
    %celldisp(featurevectors);
end

%function to smooth the image
function I = smooth(inp, gauss, M, N)
    I = inp;
    n = round(size(gauss,2)/2);
    for i = n:M-n
        for j = n:N-n
            window = inp(i,j-(n-1):j+(n-1));
            I(i,j) = sum(double(window).*gauss);
        end
    end

    for i = n:M-n
        for j = n:N-n
            window = inp(i-(n-1):i+(n-1),j);
            I(i,j) = sum(double(window).*gauss');
        end
    end
end
%function to get the gaussian
function g = getGaussian(n, sigma)
    g = zeros(1,n);
    for i = 1:n
        ep = -(((i-round(n/2))^2)/(2*sigma^2));
        g(i) = exp(ep);
    end
    g = g/sum(g);
end
