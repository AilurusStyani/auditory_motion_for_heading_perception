function modifyStarField()

global VISUAL;
global STARDATA;
dimensionX = VISUAL.dimensionX;
dimensionY = VISUAL.dimensionY;
dimensionZ = VISUAL.dimensionZ ;
size = VISUAL.starSize ;
density = VISUAL.density ;
probability = VISUAL.probability ;
totalDots = dimensionX*dimensionY*dimensionZ*density;
totalModifyDots = floor(totalDots*(1.0-probability));
j = 1;
if totalModifyDots>1
    for i=1:totalDots
        if probability < rand()
            baseX=rand()*dimensionX-dimensionX/2.0;
            baseY=rand()*dimensionY-dimensionY/2.0;
            baseZ=rand()*dimensionZ-dimensionZ/4.0*3.0;
      
            %Vertex1
            Vertex1X=baseX - size/2.0;
            Vertex1Y=baseY - size/2.0;
            Vertex1Z=baseZ;
            
            %Vertex2
            Vertex2X=baseX;
            Vertex2Y=baseY + size/2.0;
            Vertex2Z=baseZ;
            
            %Vertex3
            Vertex3X=baseX + size/2.0;
            Vertex3Y=baseY - size/2.0;
            Vertex3Z=baseZ;
            
            STARDATA.x(j)=Vertex1X;
            STARDATA.y(j)=Vertex1Y;
            STARDATA.z(j)=Vertex1Z;
            j=j+1;
            STARDATA.x(j)=Vertex2X;
            STARDATA.y(j)=Vertex2Y;
            STARDATA.z(j)=Vertex2Z;
            j=j+1;
            STARDATA.x(j)=Vertex3X;
            STARDATA.y(j)=Vertex3Y;
            STARDATA.z(j)=Vertex3Z;
            j=j+1;
        else
            j=j+3;
        end
    end
end




