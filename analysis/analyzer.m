close all;
clear all;

dataPath = 'D:\BYC\project\OpenAl\auditory_motion_for_heading_perception\Stimulus\data';
files = dir(fullfile(dataPath,'auditoryMotion_*.mat'));
figureNum = 1;
colorIndex = {[0.9,0.2,0.2],[0.2,0.9,0.2],[0.2,0.2,0.9]};
for fileI = 1:length(files)
    nameIndex = strfind(files(fileI).name,'_');
    subName = files(fileI).name(nameIndex(1)+1:nameIndex(2)-1);
    
    data = load(fullfile(dataPath,files(fileI).name));
    
    if size(data.choice,2)==1
        % find trials did not provide visual cues
        audiIndex = find(sum(isnan(cell2mat(data.conditionIndex(:,1:3))),2)>0);
        % find trials did not provide auditory cues
        visualIndex = find(sum(isnan(cell2mat(data.conditionIndex(:,4:6))),2)>0);
        trialNCheck = 0;
    elseif size(data.choice,2)==2
        audiIndex = cell2mat(data.conditionIndex(logical(sum(isnan(cell2mat(data.conditionIndex(:,1:3))),2)>0),end));
        visualIndex = cell2mat(data.conditionIndex(logical(sum(isnan(cell2mat(data.conditionIndex(:,4:6))),2)>0),end));
        trialNCheck = 1;
    end
    
    % auditory only
    aParameter = cell2mat(data.conditionIndex(audiIndex,[4:6,end]));
    aSource = data.conditionIndex(audiIndex,[7:10,end]);
    aChoice = data.choice(audiIndex,:);
    if trialNCheck
        if isequal(aParameter(:,end), aChoice(:,end),cell2mat(aSource(:,end)))
            disp('Trial number check passed.');
        else
            error('Invalid trial number.');
        end
    end
    aHeadDeg = aParameter(:,1);
    aSourceDeg = aSource(:,1:2); % cell
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
    plot(aUniqueDeg,aPR,'*');
    plot(xi,y_fit,'-');
    set(gca, 'xlim',[min(aUniqueDeg)-3,max(aUniqueDeg)+3],'ylim',[0 1])
    xlabel('Heading degree');
    ylabel('Proportion of "right" choice');
    title(['Participant ' subName]);
    text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',aBias),'color','b')
    text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', aThreshold),'color','b');
    
    figureNum = figureNum+1;
    uniqueSourceDeg = unique(cell2mat(aSourceDeg(:,2)));
    uniqueSourceN = length(uniqueSourceDeg);
    if uniqueSourceN ~=1
        for sourcei = 1:uniqueSourceN
            sourceIndex = cell2mat(aSource(:,2)) == uniqueSourceDeg(sourcei);
            sChoice = aChoice(sourceIndex,:);
            sHeadDeg = aParameter(sourceIndex,1);
            sUniqueDeg = unique(sHeadDeg);
            sRight = zeros(size(sUniqueDeg));
            sChoiceTimes = zeros(size(sUniqueDeg));
            
            for i=1:length(sUniqueDeg)
                suniqueIndex = sHeadDeg == sUniqueDeg(i);
                sRight(i) = sum(sChoice(suniqueIndex,1)-1);
                sChoiceTimes(i) = sum(suniqueIndex);
            end
            sPR = sRight./sChoiceTimes;
            fitData = [sUniqueDeg,sPR,sChoiceTimes];
            
            [sBias,sThreshold] = cum_gaussfit_max1(fitData(1:end,:));
            xi = min(sUniqueDeg):0.1:max(sUniqueDeg);
            y_fit = cum_gaussfit([sBias,sThreshold],xi);
            
            if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
            plot([0,0],[0,1],'-.k');
            hold on
            plot([uniqueSourceDeg(sourcei),uniqueSourceDeg(sourcei)],[0,1],'--','color',colorIndex{sourcei});
            plot(sUniqueDeg,sPR,'*','color',colorIndex{sourcei});
            plot(xi,y_fit,'-','color',colorIndex{sourcei});
            set(gca, 'xlim',[min(sUniqueDeg)-3,max(sUniqueDeg)+3],'ylim',[0 1])
            xlabel('Heading degree');
            ylabel('Proportion of "right" choice');
            title(['Participant ' subName]);
            text(5,1-0.2*sourcei,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',sBias),'color',colorIndex{sourcei})
            text(5,0.9-0.2*sourcei,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', sThreshold),'color',colorIndex{sourcei});
        end
    end
end