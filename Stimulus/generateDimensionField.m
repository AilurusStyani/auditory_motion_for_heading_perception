function [x,z] = generateDimensionField(distance,degree,frustumLeft,frustumRight,frustumDepth)
if iscell(distance)
    distance = cell2mat(distance);
end
if iscell(degree)
    degree = cell2mat(degree);
end
posX = max(distance)*max(sind(degree))+ frustumRight;
negX = max(distance)*min(sind(degree)) + frustumLeft;
x = max(abs([posX,negX]))*2;
z = max(distance)*max(cosd(degree))+frustumDepth;