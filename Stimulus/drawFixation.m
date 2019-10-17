function drawFixation(fixationPosition,sizeP,win)
fixation = [fixationPosition(1)-sizeP fixationPosition(2)-sizeP fixationPosition(1)+sizeP fixationPosition(2)+sizeP];
Screen('FillOval', win, [255 255 0 255], fixation);
end