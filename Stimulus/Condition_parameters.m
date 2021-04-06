% 一种声源数条件用“[]”分隔，一种声源数下的各个声源条件用“；”分隔，各声源条件内部用","，声源数条件之间用“，”
% sourceDisrance_0 = eval('0.1*coordinateMuilty,0.3*coordinateMuilty');
% sourceDegree_0  = eval('25,26');

% a1 = 21; a2 = 21.1;
% AUDITORY.sourceNum     = {3};
% AUDITORY.sourceHeading = {[180,180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-a2,-a1; -0.1,0.1; a1,a2]}; % degree for position [-55,-35;35,55] [-30,-10;10,30]
% AUDITORY.sourceLifeTimeSplit = 1;

% a1 = 21; a2 = 21.1; a3 = 40; a4 = 40.1;
% AUDITORY.sourceNum     = {4};
% AUDITORY.sourceHeading = {[180,180,180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-a4,-a3; -a2,-a1; a1,a2; a3,a4]}; % degree for position [-55,-35;35,55] [-30,-10;10,30]
% AUDITORY.sourceLifeTimeSplit = 1;

% a1 = 21; a2 = 21.1; a3 = 25; a4 = 25.1; a5 = 30; a6 = 30.1; a7 = 35; a8 = 35.1;
% AUDITORY.sourceNum     = {8};
% AUDITORY.sourceHeading = {[180,180,180,180,180,180,180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-a8,-a7; -a6,-a5; -a4,-a3; -a2,-a1; a1,a2; a3,a4; a5,a6; a7,a8]}; % degree for position
% AUDITORY.sourceLifeTimeSplit = 1;

% AUDITORY.sourceNum = {1 1};
% AUDITORY.sourceHeading = {180, 180}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.3*coordinateMuilty, 0.1*coordinateMuilty,0.3*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-25,-24], [24,25]}; % degree for position
% AUDITORY.sourceLifeTimeSplit = 1;

% AUDITORY.sourceNum     = {2 3};
% AUDITORY.sourceHeading = {[180,180], [180,180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.3*coordinateMuilty; 0.1*coordinateMuilty,0.3*coordinateMuilty, 0.1*coordinateMuilty,0.3*coordinateMuilty; 0.1*coordinateMuilty,0.3*coordinateMuilty; 0.1*coordinateMuilty,0.3*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-30,-10;10,30], [-30,-10;10,30;10,30]}; % degree for position
% AUDITORY.sourceLifeTimeSplit = 1;

%     t = linspace(0,2,1000);
%     if i == 1
%         alSourcef(sources(i), AL.GAIN, 10.^sin(frequency*pi*t) );
%     elseif i == 2
%         alSourcef(sources(i), AL.GAIN, 10.^cos(frequency*pi*t + pi/2) );
%     end

%     if i == 2
%         alSourcef(sources(i), AL.GAIN, 1);
%     elseif i == 3
%         alSourcef(sources(i), AL.GAIN, 1);
%     elseif i == 4
%         alSourcef(sources(i), AL.GAIN, 1);
%     elseif i == 1
%         alSourcef(sources(i), AL.GAIN, 1);
%     end

% 振幅变化
        frequency = 1/20;
%         for i = 1:2
%             if i == 1
%                 alSourcef(sources(i), AL.GAIN, 10.^sin(frequency*pi*framei) );
%             elseif i == 2
%                 alSourcef(sources(i), AL.GAIN, 10.^cos(frequency*pi*framei + pi/2) );
%             end
%         end
            
%         frequency = 1/5;
%         for i = 1:2
%             if i == 1
%                 alSourcef(sources(i), AL.GAIN, square(frequency*framei, 30) /2 + 0.5 );
%             elseif i == 2
%                 alSourcef(sources(i), AL.GAIN, square(frequency*framei + pi, 30) /2 + 0.5 );
%             end
%         end

%         frequency = 1/5;
%         for i = 1:4
%             if i == 1
%                 alSourcef(sources(i), AL.GAIN, square(frequency*framei, 20) /2 + 0.5 );
%             elseif i == 2
%                 alSourcef(sources(i), AL.GAIN, square(frequency*framei + pi/2, 20) /2 + 0.5 );
%             elseif i == 3
%                 alSourcef(sources(i), AL.GAIN, square(frequency*framei + pi, 20) /2 + 0.5 );
%             elseif i == 4
%                 alSourcef(sources(i), AL.GAIN, square(frequency*framei + 3*pi/2, 20) /2 + 0.5 );
%             end
%         end

%         frequency = 1/5;
%         s = Shuffle(1:auditorySourcei{1});
%         for i = 1:auditorySourcei{1}
%             alSourcef(sources(s(i)), AL.GAIN, square(frequency*framei + 2*pi/8*(i-1), 10) /2 + 0.5 );
%         end