function lifetimeF = calculateAuditoryLifetime(duration,splitNum,refreshRate)
if iscell(duration)
    duration = max(cell2mat(duration));
else
    duration = max(duration);
end
lifetimeF = ceil((duration*refreshRate)/splitNum);
end