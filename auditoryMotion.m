clear all
close all

global TRIALINFO
global SCREEN
global AUDITORY
global VISUAL
global GL

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
TRIALINFO.fixationPeriod = 0; % second

% 1 for intergration, both visual and auditory use the parameters in TRIALINFO,
% 0 for segregation, visual cue will use VISUAL, and auditory will use AUDITORY.
TRIALINFO.intergration = [1];

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

TRIALINFO.deviation = 1.2; % initial binocular deviation, cm
deviationAdjust = 0.2; % how fast to adjust the deviation by key pressing, cm

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
calculateConditions();
% TRIALINFO.trialConditions =
% {visualDegree visualDistance visualTime, ...
%       1               2               3
%
% auditoryDegree auditoryDistance auditoryTime sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}
%       4               5               6                7              8                9                  10

trialIndex = repmat(1:size(TRIALINFO.trialConditions,1),TRIALINFO.repetition,1);
trialNum = size(trialIndex,1);
trialOrder = randperm(trialNum);

disp(['This block has  ' num2str(trialNum) ' trials']);

timePredicted = (TRIALINFO.fixationPeriod + TRIALINFO.initialPeriod + mean(cell2mat(TRIALINFO.headingTime)) + TRIALINFO.choicePeriod + ...
    TRIALINFO.intertrialInterval ) * trialNum;
disp(['This block will cost  ' num2str(timePredicted/60) ' minutes']);
calibrationInterval = 600; % unit second, it is better to re-calibration every 10-15 minutes
automaticCalibration = timePredicted > 1.3*calibrationInterval; % make automatic calibration (every 10 min in default) if the block takes more than 15 min.


%% initial opengl
visualCue = sum(ismember(TRIALINFO.stimulusType,[0 2]));
if visualCue
    if testMode
        Screen('Preference', 'SkipSyncTests', 1); % for debug/test
    else
        Screen('Preference', 'SkipSyncTests', 0); % for recording
    end
    
    AssertOpenGL;
    InitializeMatlabOpenGL;
    
    SCREEN.screenId = max(Screen('Screens'));
    PsychImaging('PrepareConfiguration');
    
    % Define background color:
    whiteBackground = WhiteIndex(SCREEN.screenId);
    blackBackground = BlackIndex(SCREEN.screenId);
    
    % Open a double-buffered full-screen window on the main displays screen.
    [win , winRect] = PsychImaging('OpenWindow', SCREEN.screenId, whiteBackground);
    SCREEN.widthPix = winRect(3);
    SCREEN.heightPix = winRect(4);
    SCREEN.center = [SCREEN.widthPix/2, SCREEN.heightPix/2];
    
    TRIALINFO.fixationSizeP = degree2pix(TRIALINFO.fixationSizeD/2);
    TRIALINFO.fixationPosition = [SCREEN.widthPix/2,SCREEN.heightPix/2];
    
    SCREEN.refreshRate = Screen('NominalFrameRate', SCREEN.screenId);
    
    %% the configuration of the Frustum
    calculateFrustum(coordinateMuilty);
    
    Screen('BeginOpenGL', win);
    glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
    glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
    % glEnable(GL_BLEND);
    % glEnable(GL_ALPHA_BLEND_CORRECTLY);
    % glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('EndOpenGL', win);
    
    %% initial eyelink
    if ~testMode
        tempName = 'TEMP1'; % need temp name because Eyelink only know hows to save names with 8 chars or less. Will change name using matlab's moveFile later.
        dummymode=0;
        
        el=EyelinkInitDefaults(win);
        %     el.backgroundcolour = backgroundColor;
        %     el.foregroundcolour = BlackIndex(el.window);
        %     el.msgfontcolour    = BlackIndex(el.window);
        %     el.imgtitlecolour   = BlackIndex(el.window);
        
        if ~EyelinkInit(dummymode)
            fprintf('Eyelink Init aborted.\n');
            cleanup;  % cleanup function
            Eyelink('ShutDown');
            Screen('CloseAll');
            return
        end
        
        testi = Eyelink('Openfile', tempName);
        if testi~=0
            fprintf('Cannot create EDF file ''%s'' ', fileName);
            cleanup;
            Eyelink('ShutDown');
            Screen('CloseAll');
            return
        end
        
        %   SET UP TRACKER CONFIGURATION
        Eyelink('command', 'calibration_type = HV9');
        %	set parser (conservative saccade thresholds)
        Eyelink('command', 'saccade_velocity_threshold = 35');
        Eyelink('command', 'saccade_acceleration_threshold = 9500');
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,HREF,GAZERES,AREA,STATUS,INPUT,HTARGET');
        Eyelink('command', 'online_dcorr_refposn = %1d, %1d', SCREEN.center(1), SCREEN.center(2));
        Eyelink('command', 'online_dcorr_maxangle = %1d', 30.0);
        % you must call this function to apply the changes from above
        EyelinkUpdateDefaults(el);
        
        % Calibrate the eye tracker
        EyelinkDoTrackerSetup(el);
        
        % do a final check of calibration using driftcorrection
        EyelinkDoDriftCorrection(el);
        
        Eyelink('StartRecording');
        
        Eyelink('message', 'SYNCTIME');	 	 % zero-plot time for EDFVIEW
        eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
        if eye_used == el.BINOCULAR % if both eyes are tracked
            eye_used = el.LEFTEYE; % use left eye
        end
        errorCheck=Eyelink('checkrecording'); 		% Check recording status */
        if(errorCheck~=0)
            fprintf('Eyelink checked wrong status.\n');
            cleanup;  % cleanup function
            Eyelink('ShutDown');
            Screen('CloseAll');
        end
        
        calibrateCkeck = tic;
        WaitSecs(1); % wait a little bit, in case the key press during calibration influence the following keyboard check
    end
end
%% initial openal
auditoryCue = sum(ismember(TRIALINFO.stimulusType,[1 2]));
if auditoryCue
    % Initialize OpenAL subsystem at debuglevel 2 with the default output device: 
    InitializeMatlabOpenAL(2);
    
    % Generate one sound buffer:
    buffers = alGenBuffers(AUDITORY.sourceNum);
    
    % Query for errors:
    alGetString(alGetError)
    
    soundfiles = dir(fullfile(pwd,'*.wav'));
    
    % Create a sound source:
    sources = alGenSources(AUDITORY.sourceNum);
    
    perm = randperm(AUDITORY.sourceNum);
end
%%
HideCursor(SCREEN.screenId);
GenerateStarField();