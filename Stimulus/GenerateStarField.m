function GenerateStarField()
 global VISUAL;
 global STARDATA;
totalDots = round(VISUAL.dimensionX*VISUAL.dimensionY*VISUAL.dimensionZ*VISUAL.density);
baseX=rand(1,totalDots)*VISUAL.dimensionX-VISUAL.dimensionX/2.0;
baseY=rand(1,totalDots)*VISUAL.dimensionY-VISUAL.dimensionY/2.0;
baseZ=rand(1,totalDots)*VISUAL.dimensionZ-VISUAL.dimensionZ;
STARDATA.x=zeros(1,3*totalDots);
STARDATA.y=zeros(1,3*totalDots);
STARDATA.z=zeros(1,3*totalDots);
starSize = degree2length(VISUAL.starSize);

%Vertex1
Vertex1X=baseX - starSize/2.0;
Vertex1Y=baseY - starSize/2.0;
Vertex1Z=baseZ;

%Vertex2
Vertex2X=baseX;
Vertex2Y=baseY + starSize/2.0;
Vertex2Z=baseZ;

%Vertex3
Vertex3X=baseX + starSize/2.0;
Vertex3Y=baseY - starSize/2.0;
Vertex3Z=baseZ;
j=1;

for i=1:totalDots
    STARDATA.x(j)=Vertex1X(i);
    STARDATA.y(j)=Vertex1Y(i);
    STARDATA.z(j)=Vertex1Z(i);
    j=j+1;
    STARDATA.x(j)=Vertex2X(i);
    STARDATA.y(j)=Vertex2Y(i);
    STARDATA.z(j)=Vertex2Z(i);
    j=j+1;
    STARDATA.x(j)=Vertex3X(i);
    STARDATA.y(j)=Vertex3Y(i);
    STARDATA.z(j)=Vertex3Z(i);
    j=j+1;
end
end

