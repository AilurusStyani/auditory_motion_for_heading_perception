function [x,y,z] = generateSourceMetrixLocation();
global AUDITORY
xMetrix = linspace(AUDITORY.sourceX(1),AUDITORY.sourceX(2),AUDITORY.sourceXnum);
yMetrix = linspace(AUDITORY.sourceY(1),AUDITORY.sourceY(2),AUDITORY.sourceYnum);
zMetrix = linspace(AUDITORY.sourceZ(1),AUDITORY.sourceZ(2),AUDITORY.sourceZnum);
x = repmat(xMetrix,1,AUDITORY.sourceYnum*AUDITORY.sourceZnum);
y = repmat(sort(repmat(yMetrix,1,AUDITORY.sourceXnum)),1,AUDITORY.sourceZnum);
z = sort(repmat(zMetrix,1,AUDITORY.sourceXnum*AUDITORY.sourceYnum));