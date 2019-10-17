function pixel = degree2pix(degree,dir)
% this function convert degree value to pixel value
% it is better been used to calculte the pixel length from the central point
% On X axis: dir = 1; On Y axis: dir = 2
global SCREEN

a=SCREEN.widthPix / SCREEN.widthCM;
b=SCREEN.heightPix / SCREEN.heightCM;

if nargin == 1
    if abs(a-b)/min(a,b) < 0.05
        length = tand(degree) * SCREEN.distance;
        pixel = length / SCREEN.widthCM * SCREEN.widthPix;
    else
        error('Error in screen parameter or screen config, or you should define it is horizontal(1) / vertical(2).')
    end
elseif nargin == 2
    pixel = tand(degree) * SCREEN.distance;
    if dir == 1
        pixel = pixel / SCREEN.widthCM * SCREEN.widthPix;
    elseif dir == 2
        pixel = pixel / SCREEN.heightCM * SCREEN.heightPix;
    else
        error('Invalid value for dir.')
    end
end