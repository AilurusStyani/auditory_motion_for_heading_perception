function calculateConditions()
% {visualDegree visualDistance visualTime, ...
% auditoryDegree auditoryDistance auditoryTime sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}

global TRIALINFO
global AUDITORY
global VISUAL

% TRIALINFO.stimulusType = [0 1 2]; % 0 for visual only, 1 for auditory only, 2 for both provided

% 0 for segregation
inteCondition0 = {};
% index0 = {};
if sum(ismember(TRIALINFO.intergration,0))
    audiCon = cell(length(AUDITORY.sourceNum),4);
    for i=1:length(AUDITORY.sourceNum)
        audiCon(i,:) = {AUDITORY.sourceNum{i},AUDITORY.sourceDegree{i},AUDITORY.sourceDistance{i},AUDITORY.sourceHeading{i}};
    end
    con0AtS = [sortrows(repmat(AUDITORY.headingTime',size(audiCon,1),1)),...
        repmat(audiCon,length(AUDITORY.headingTime),1)];
    con0AdAtS = [sortrows(repmat(AUDITORY.headingDistance',size(con0AtS,1),1)),...
        repmat(con0AtS,length(AUDITORY.headingDistance),1)];
    con0AS = [sortrows(repmat(AUDITORY.headingDegree',size(con0AdAtS,1),1)),...
        repmat(con0AdAtS,length(AUDITORY.headingDegree),1)];
    
    con0VdVt = [sortrows(repmat(VISUAL.headingDistance',length(VISUAL.headingTime),1)),...
        repmat(VISUAL.headingTime',length(VISUAL.headingDistance),1)];
    con0V = [sortrows(repmat(VISUAL.headingDegree',size(con0VdVt,1),1)),...
        repmat(con0VdVt,length(VISUAL.headingDegree),1)];
    
    % visual only
    con0T0 = cat(2,con0V,num2cell(nan(size(con0V,1),size(con0AS,2))));
    
    % auditory only
    con0T1 = cat(2,num2cell(nan(size(con0AS,1),size(con0V,2))),con0AS);
    
    % both provided
    con0T2={};
    % reserved %
    
    if sum(ismember(TRIALINFO.stimulusType,0))
        inteCondition0 = cat(1,inteCondition0,con0T0);
    end
    if sum(ismember(TRIALINFO.stimulusType,1))
        inteCondition0 = cat(1,inteCondition0,con0T1);
    end
    if sum(ismember(TRIALINFO.stimulusType,2))
        inteCondition0 = cat(1,inteCondition0,con0T2);
    end
end

% 1 for intergration
inteCondition1 = {};
% index1 = [];
if sum(ismember(TRIALINFO.intergration,1))
    audiCon = cell(length(AUDITORY.sourceNum),4);
    for i=1:length(AUDITORY.sourceNum)
        audiCon(i,:) = {AUDITORY.sourceNum{i},AUDITORY.sourceDegree{i},AUDITORY.sourceDistance{i},AUDITORY.sourceHeading{i}};
    end
    con1AtS = [sortrows(repmat(AUDITORY.headingTime',size(audiCon,1),1)),...
        repmat(audiCon,length(AUDITORY.headingTime),1)];
    con1AdAtS = [sortrows(repmat(AUDITORY.headingDistance',size(con1AtS,1),1)),...
        repmat(con1AtS,length(AUDITORY.headingDistance),1)];
    con1AS = [sortrows(repmat(AUDITORY.headingDegree',size(con1AdAtS,1),1)),...
        repmat(con1AdAtS,length(AUDITORY.headingDegree),1)];
    
    con1VdVt = [sortrows(repmat(VISUAL.headingDistance',length(VISUAL.headingTime),1)),...
        repmat(VISUAL.headingTime',length(VISUAL.headingDistance),1)];
    con1V = [sortrows(repmat(VISUAL.headingDegree',size(con1VdVt,1),1)),...
        repmat(con1VdVt,length(VISUAL.headingDegree),1)];
    
    % visual only
    con1T0 = cat(2,con1V,num2cell(nan(size(con1V,1),size(con1AS,2))));
    
    % auditory only
    con1T1 = cat(2,num2cell(nan(size(con1AS,1),size(con1V,2))),con1AS);
    
    % both provided
    con1T2 = cat(2,con1AS(:,1:3),con1AS);
    
    if sum(ismember(TRIALINFO.stimulusType,0))
        inteCondition1 = cat(1,inteCondition1,con1T0);
    end
    if sum(ismember(TRIALINFO.stimulusType,1))
        inteCondition1 = cat(1,inteCondition1,con1T1);
    end
    if sum(ismember(TRIALINFO.stimulusType,2))
        inteCondition1 = cat(1,inteCondition1,con1T2);
    end
end

TRIALINFO.trialConditions = [inteCondition0;inteCondition1];
% trialIndex = [index0;index1];