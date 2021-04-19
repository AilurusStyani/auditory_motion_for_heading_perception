function lifetimeF = calculateAuditoryLifetime(duration,splitNum,refreshRate)
% this function coded for auditory heading experiment in lifetime version.
% this function calculated how many frames in visual going compare to the auditory lifetime
% duration.
if iscell(duration)
    duration = max(cell2mat(duration));
else
    duration = max(duration);
end
lifetimeF = ceil((duration*refreshRate)/splitNum);
end