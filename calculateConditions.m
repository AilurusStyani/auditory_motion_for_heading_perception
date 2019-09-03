function [trialCondition,trialIndex] = calculateConditions()
% {headDegree headDistance headTime, ...
% visualDegree visualDistance visualTime, ...
% auditoryDegree auditoryDistance auditoryTime sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}

global TRIALINFO
global AUDITORY
global VISUAL

TRIALINFO.stimulusType = [0 1 2]; % 0 for visual only, 1 for auditory only, 2 for both provided

% 0 for segregation
inteCondition0 = {};
index0 = [];
if sum(ismember(TRIALINFO.intergration,0))
    audiCon = cell(length(AUDITORY.sourceNum),4);
    for i=1:length(AUDITORY.sourceNum)
        audiCon(i,:) = {AUDITORY.sourceNum{i},AUDITORY.sourceDegree{i},AUDITORY.sourceDistance{i},AUDITORY.sourceHeading{i}};
    end
    int0AtS = [sort(repmat(AUDITORY.headingTime',size(audiCon,1),1)),...
        repmat(audiCon,length(AUDITORY.headingTime),1)];
    int0AdAtS = [sort(repmat(AUDITORY.headingDistance',size(int0AtS,1),1)),...
        repmat(int0AtS,length(AUDITORY.headingDistance),1)];
    int0AdAdAtS = [sort(repmat(AUDITORY.headingDegree',size(int0AdAtS,1),1)),...
        repmat(int0AdAtS,length(AUDITORY.headingDegree),1)];
end

% 1 for intergration
inteCondition1 = {};
index1 = [];
if sum(ismember(TRIALINFO.intergration,1))
    
end

trialCondition = [inteCondition0;inteCondition1];
trialIndex = [index0;index1];