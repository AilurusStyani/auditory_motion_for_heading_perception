close all;
clear all;

dataPath = 'D:\BYC\project\OpenAl\auditory_motion_for_heading_perception\Stimulus\data';
files = dir(fullfile(dataPath,'auditoryMotion_*.mat'));
figureNum = 1;
colorIndex = {[0.9,0.2,0.2],[0.2,0.9,0.2],[0.2,0.2,0.9]};
for fileI = 1:length(files)
    nameIndex = strfind(files(fileI).name,'_');
    subName = files(fileI).name(nameIndex(1)+1:nameIndex(2)-1);
    dateNum = files(fileI).name(nameIndex(2)+1:nameIndex(2)+11);
%     if contains(subName,'test') || isempty(subName)
%         continue
%     end
    
    data = load(fullfile(dataPath,files(fileI).name));
    
    % block check
    if any(isnan(data.choiceTime),true)
        fprintf(2,[files(fileI).name '\n']);
        fprintf(2,['This file is skipped cause the experiment was not completed.\n\n']);
        continue
    end
    
    if size(data.choice,2)==2
        audiIndex = cell2mat(data.conditionIndex(logical(sum(isnan(cell2mat(data.conditionIndex(:,1:3))),2)>0),end));
        visualIndex = cell2mat(data.conditionIndex(logical(sum(isnan(cell2mat(data.conditionIndex(:,4:6))),2)>0),end));
        bothIndex = cell2mat(data.conditionIndex(logical(sum(isnan(cell2mat(data.conditionIndex(:,3:4))),2)<1),end));
    else
        error('Invalid data size or not the same version.')
    end
    
    % auditory only
    aParameter = cell2mat(data.conditionIndex(audiIndex,[4:6,end]));
    aSource = data.conditionIndex(audiIndex,[7:10,end]);
    aChoice = data.choice(audiIndex,:);
    if ~isequal(aParameter(:,end), aChoice(:,end),cell2mat(aSource(:,end)))
        error('There have some error for data extraction.');
    end
    
    aHeadDeg = aParameter(:,1);
    aUniqueDeg = unique(aHeadDeg);
    aRight = zeros(size(aUniqueDeg));
    aChoiceTimes = zeros(size(aUniqueDeg));
    
    for i=1:length(aUniqueDeg)
        uniqueIndex = aHeadDeg == aUniqueDeg(i);
        aRight(i) = sum(aChoice(uniqueIndex,1)-1);
        aChoiceTimes(i) = sum(uniqueIndex);
    end
    aPR = aRight./aChoiceTimes;
    fitData = [aUniqueDeg,aPR,aChoiceTimes];
    
    [aBias,aThreshold] = cum_gaussfit_max1(fitData(1:end,:));
    xi = min(aUniqueDeg):0.1:max(aUniqueDeg);
    y_fit = cum_gaussfit([aBias,aThreshold],xi);
    
    if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
    plot([0,0],[0,1],'-.k');
    hold on
    plot(aUniqueDeg,aPR,'*r');
    plot(xi,y_fit,'-r');
    set(gca, 'xlim',[min(aUniqueDeg)-3,max(aUniqueDeg)+3],'ylim',[0 1])
    xlabel('Heading degree');
    ylabel('Proportion of "right" choice');
    title(['Participant ' subName ' ' dateNum]);
    text(5,0.9,sprintf('\\it\\mu_{Apsy} = \\rm%6.3g\\circ',aBias),'color','r')
    text(5,0.8,sprintf('\\it\\sigma_{Apsy} = \\rm%6.3g\\circ', aThreshold),'color','r');
    
    % visual only
    vParameter = cell2mat(data.conditionIndex(visualIndex,[1:3,end]));
    vChoice = data.choice(visualIndex,:);
    vHeadDeg = vParameter(:,1);
    vUniqueDeg = unique(vHeadDeg);
    vRight = zeros(size(vUniqueDeg));
    vChoiceTimes = zeros(size(vUniqueDeg));
    
    for i = 1:length(vUniqueDeg)
        uniqueIndex = vHeadDeg == vUniqueDeg(i);
        vRight(i) = sum(vChoice(uniqueIndex,1)-1);
        vChoiceTimes(i) = sum(uniqueIndex);
    end
    vPR = vRight ./ vChoiceTimes;
    fitData = [vUniqueDeg, vPR, vChoiceTimes];
    [vBias, vThreshold] = cum_gaussfit_max1(fitData(1:end,:));
    xi = min(vUniqueDeg):0.1:max(vUniqueDeg);
    y_fit = cum_gaussfit([vBias,vThreshold],xi);
    
    plot(vUniqueDeg,vPR,'*b');
    plot(xi,y_fit,'-b');
    text(5,0.7,sprintf('\\it\\mu_{Vpsy} = \\rm%6.3g\\circ',vBias),'color','b')
    text(5,0.6,sprintf('\\it\\sigma_{Vpsy} = \\rm%6.3g\\circ', vThreshold),'color','b');
    
    % both provided
    bParameter = cell2mat(data.conditionIndex(bothIndex,[1:5,end]));
    bChoice = data.choice(bothIndex,:);
    vH = bParameter(:,1); aH = bParameter(:,4);
    bHeadDeg = mean([vH,aH],2);
    bUniqueDeg = unique(bHeadDeg);
    bRight = zeros(size(bUniqueDeg));
    bChoiceTimes = zeros(size(bUniqueDeg));
    
    for i = 1:length(bUniqueDeg)
        uniqueIndex = bHeadDeg == bUniqueDeg(i);
        bRight(i) = sum(bChoice(uniqueIndex,1)-1);
        bChoiceTimes(i) = sum(uniqueIndex);
    end
    bPR = bRight ./ bChoiceTimes;
    fitData = [bUniqueDeg, bPR, bChoiceTimes];
    [bBias, bThreshold] = cum_gaussfit_max1(fitData(1:end,:));
    xi = min(bUniqueDeg):0.1:max(bUniqueDeg);
    y_fit = cum_gaussfit([bBias,bThreshold],xi);
    
    plot(bUniqueDeg,bPR,'*k');
    plot(xi,y_fit,'-k');
    text(5,0.5,sprintf('\\it\\mu_{Bpsy} = \\rm%6.3g\\circ',bBias),'color','k')
    text(5,0.4,sprintf('\\it\\sigma_{Bpsy} = \\rm%6.3g\\circ', bThreshold),'color','k');
    pred = sqrt((vThreshold^2*aThreshold^2)/(vThreshold^2+aThreshold^2));
    text(5,0.3,sprintf('\\it\\sigma_{prepsy} = \\rm%6.3g\\circ', pred),'color','k');
    figureNum = figureNum +1;
end