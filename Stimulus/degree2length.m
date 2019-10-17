function length = degree2length(degree)
% this function convert degree value to pixel value
% it is better been used to calculte the pixel length from the central point
% On X axis: dir = 1; On Y axis: dir = 2
global SCREEN

a=SCREEN.widthPix / SCREEN.widthCM;
b=SCREEN.heightPix / SCREEN.heightCM;

if nargin == 1
    if abs(a-b)/min(a,b) < 0.05
        length = tand(degree) * SCREEN.distance;
    else
        error('Error in screen parameter or screen config, or you should define it is horizontal(1) / vertical(2).')
    end
end