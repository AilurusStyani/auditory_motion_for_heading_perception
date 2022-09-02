close all;
clear all;

FileName=('Z:\LQY Experiment\DOTS\tester5\auditoryMotion_GHR9_2207211833.mat');
[pathstr,name]=fileparts(FileName);
load(fullfile(pathstr,name));


%get number of parameters (X here) u need to fit in AUDITORY struct
%coh Analyzer: "coherence"; heading Analyzer: "headingDegree"
Parameter_need_to_fit='headingDegree';% "coherence" or "headingDegree"
switch Parameter_need_to_fit
    case 'coherence'
          Number_X=length(getfield(AUDITORY, Parameter_need_to_fit));
          X=sort(cell2mat(getfield(AUDITORY, Parameter_need_to_fit)));  
          choice_X_list=[choice(:,1),cell2mat(conditionIndex(:,11))];
    case 'headingDegree'
          Number_X=length(getfield(AUDITORY, Parameter_need_to_fit));
          X=sort(cell2mat(getfield(AUDITORY, Parameter_need_to_fit)));  
          choice_X_list=[choice(:,1),cell2mat(conditionIndex(:,4))];
    otherwise disp('other value');
end

% a matrix to count rightchoice times in each X
Counter=zeros(Number_X,1);
for i=1:size(choice_X_list,1)
    for j=1:Number_X
      if (choice_X_list(i,2)==X(j)) &&  (choice_X_list(i,1)==2)
          Counter(j)=Counter(j)+1;
      end
    end
end

% fit psychometric curve       
aChoiceTimes(1:Number_X)=size(choice_X_list,1)/Number_X;
aChoiceTimes=aChoiceTimes';
aPR = Counter./aChoiceTimes;
figureNum = 1;
aUniqueDeg=X';
fitData = [aUniqueDeg,aPR,aChoiceTimes];
    
[aBias,aThreshold] = cum_gaussfit_max1(fitData(1:end,:));
xi = min(aUniqueDeg):0.1:max(aUniqueDeg);
y_fit = cum_gaussfit([aBias,aThreshold],xi);
% plot   
if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
plot([0,0],[0,1],'-.k');
hold on
plot(aUniqueDeg,aPR,'*');
plot(xi,y_fit,'-');
set(gca, 'xlim',[min(aUniqueDeg),max(aUniqueDeg)],'ylim',[0 1])
xlabel('Coherence');
ylabel('Proportion of "right" choice');
title(['Participant ']);
text(6,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',aBias),'color','b')
text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', aThreshold),'color','b');
saveas(gcf,fullfile(pathstr,name),'jpg');
%csv file for DDM
for i=1:size(choice_X_list,1)
    if ((choice_X_list(i,2)<0.5) && (choice_X_list(i,1)==1)) || ((choice_X_list(i,2)>0.5) && (choice_X_list(i,1)==2))
        choice_X_list(i,3)=1;
    end

    if ((choice_X_list(i,2)<0.5) && (choice_X_list(i,1)==2)) || ((choice_X_list(i,2)>0.5) && (choice_X_list(i,1)==1))
        choice_X_list(i,3)=0;
    end
    
    if choice_X_list(i,2)==0.5
            seed=rand(1);
            if (seed<0.5 && choice_X_list(i,1)==1) || (seed>0.5 && choice_X_list(i,1)==2)
            choice_X_list(i,3)=1;
            else choice_X_list(i,3)=0;
            end
    end   
end

cohList=[choiceTime(:,1),choice_X_list(:,2),choice_X_list(:,3),choice_X_list(:,1)];
colNames={'rt','coh','correct','trgchoice'};
cohTable=array2table(cohList,'VariableNames',colNames);
writetable(cohTable,'Z:\LQY Experiment\Cohtest\csvfile\auditoryMotion_CohtestC_2208221444.csv');