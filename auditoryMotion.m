clear all
close all

global TRIALINFO
global SCREEN
global AUDITORY
global VISUAL

subjectName = inputdlg({'Please input participant''s initials.'},'1',1);
fileName = ['auditoryMotion_' subjectName{1} '_' datestr(now,'yymmddHHMM')];
saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

% set keyboard
KbName('UnifyKeyNames'); 
skipKey = KbName('space');
escape = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
upArror = KbName('UpArrow');
cKey = KbName('c'); % force calibration

pageUp = KbName('pageup'); % increase binocular deviation
pageDown = KbName('pagedown'); % decrease binocular deviation

testMode = 0; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC
feedback = 0; % in practice block, set 1 to provide feedback. otherwise set 0
feedbackDuration = 1; % unit s

%% parameters
coordinateMuilty = 1; % convert cm to coordinate system for moving distance etc.
TRIALINFO.repetition = 15;
TRIALINFO.headingDegree = {-10 -5 0 5 10};
TRIALINFO.headingDistance = {5 10};
TRIALINFO.headingTime = {1 2};
TRIALINFO.stimulusType = [0 1 2]; % 0 for visual only, 1 for auditory only, 2 for both provided

TRIALINFO.initialPeriod = 1; % second
TRIALINFO.choicePeriod = 2; % second
TRIALINFO.intertrialInterval = 1; % second

% 1 for intergration, both visual and auditory use the parameters in TRIALINFO,
% 0 for segregation, visual cue will use VISUAL, and auditory will use AUDITORY.
TRIALINFO.intergration = [0 1]; 

% for SCREEN
if testMode
%     SCREEN.widthCM = 34.5*coordinateMuilty; % cm, need to be measured on your own PC
%     SCREEN.heightCM = 19.7*coordinateMuilty; % cm, need to be measured on your own PC
    SCREEN.widthCM = 37.5*coordinateMuilty; % cm, need to be measured on your own PC
    SCREEN.heightCM = 30*coordinateMuilty; % cm, need to be measured on your own PC
else
    SCREEN.widthCM = 120*coordinateMuilty; % cm
    SCREEN.heightCM = 65*coordinateMuilty; % cm
end
SCREEN.distance = 60*coordinateMuilty;% cm

% parameters for visual cue
VISUAL.headingDegree = TRIALINFO.headingDegree; % cell
VISUAL.headingDistance = TRIALINFO.headingDistance; % cell
VISUAL.headingTime = TRIALINFO.headingTime; % cell

VISUAL.fixationSizeD = 0.25;  % degree
VISUAL.fixationWindow = 2; % degree

VISUAL.density = 1000/(100*coordinateMuilty)^3;    % convert num/m^3 to num/cm^3
VISUAL.coherence = 100; % in percent
VISUAL.lifeTime = 3; % frame number

VISUAL.dimensionX = 400*coordinateMuilty;  % cm
VISUAL.dimensionY = 400*coordinateMuilty;  % cm
VISUAL.dimensionZ = 700*coordinateMuilty;  % cm
VISUAL.starSize = 0.1;    % degree

% parameters for auditory cue
AUDITORY.height = 5*coordinateMuilty; % cm

AUDITORY.headingDegree = TRIALINFO.headingDegree; % cell
AUDITORY.headingDistance = TRIALINFO.headingDistance; % cell
AUDITORY.headingTime = TRIALINFO.headingTime; % cell

AUDITORY.sourceNum = {1,2}; 
AUDITORY.sourceHeading = {180,[180 0]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
AUDITORY.sourceDistance = {50*coordinateMuilty , [50*coordinateMuilty 60*coordinateMuilty]}; % cm
AUDITORY.sourceDegree = {0 , [-10 10]}; % degree

%% trial conditions and order
[TRIALINFO.trialConditions,conditionIndex] = calculateConditions();
trialIndex = repmat(conditionIndex, TRIALINFO.repetition,1);
trialNum = size(trialIndex,1);
trialOrder = randperm(trialNum);


% Initialize OpenAL subsystem at debuglevel 2 with the default output device:
InitializeMatlabOpenAL(2);

% Generate one sound buffer:
buffers = alGenBuffers(AUDITORY.sourceNum);

% Query for errors:
alGetString(alGetError)