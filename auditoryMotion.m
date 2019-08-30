subjectName = inputdlg({'Please input participant''s initials.'},'1',1);
fileName = ['auditoryMotion_' subjectName{1} '_' datestr(now,'yymmddHHMM')];
saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

% set keyboard
KbName('UnifyKeyNames'); 
skipKey = KbName('space');
escape = KbName('esc');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
upArror = KbName('UpArrow');
cKey = KbName('c'); % force calibration

pageUp = KbName('pageup'); % increase binocular deviation
pageDown = KbName('pagedown'); % decrease binocular deviation

testMode = 0; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC
feedback = 0; % in practice block, set 1 to provide feedback. otherwise set 0
feedbackDuration = 1; % unit s

% parameters
walkDegree = [-10 -5 0 5 10];
walkDistance = [5 10];

nsources = 1; 


% Initialize OpenAL subsystem at debuglevel 2 with the default output device:
InitializeMatlabOpenAL(2);

% Generate one sound buffer:
buffers = alGenBuffers(nsources);

% Query for errors:
alGetString(alGetError)
