CloseOpenAL;
clear all
close all

global TRIALINFO
global SCREEN
global AUDITORY
global VISUAL
global GL
global FRUSTUM
global STARDATA
global AL

subjectName = inputdlg({'Please input participant''s initials.'},'1',1);
fileName = ['auditoryMotion_' subjectName{1} '_' datestr(now,'yymmddHHMM')];
saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

% set keyboard
KbName('UnifyKeyNames');
skipKey   = KbName('space');
escape    = KbName('ESCAPE');
leftKey   = KbName('LeftArrow');
rightKey  = KbName('RightArrow');
upArror   = KbName('UpArrow');
cKey      = KbName('c'); % force calibration
enter     = KbName('Return');

pageUp = KbName('pageup'); % increase binocular deviation
pageDown = KbName('pagedown'); % decrease binocular deviation

testMode = 1; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC
feedback = 1; % in practice block, set 1 to provide feedback. otherwise set 0
feedbackDuration = 1; % unit s

%% parameters
coordinateMuilty = 1; % convert cm to coordinate system for moving distance etc.
TRIALINFO.repetition      = 15;
TRIALINFO.headingDegree   = {-30 -15 0 15 30};
TRIALINFO.headingDistance = {10 15};
TRIALINFO.headingTime      = {1 2};
TRIALINFO.stimulusType     = [1]; % 0 for visual only, 1 for auditory only, 2 for both provided

TRIALINFO.initialPeriod       = 1; % second
TRIALINFO.choicePeriod        = 2; % second
TRIALINFO.intertrialInterval = 1; % second
TRIALINFO.fixationPeriod     = 0; % second
TRIALINFO.fixationSizeD      = 0.25; % degree

% 1 for intergration, both visual and auditory use the parameters in TRIALINFO,
% 0 for segregation, visual cue will use VISUAL, and auditory will use AUDITORY.
TRIALINFO.intergration = [1];

% for SCREEN
if testMode
    %     %     SCREEN.widthCM = 34.5*coordinateMuilty; % cm, need to be measured on your own PC
    %     %     SCREEN.heightCM = 19.7*coordinateMuilty; % cm, need to be measured on your own PC
    SCREEN.widthCM = 37.5*coordinateMuilty; % cm, need to be measured on your own PC
    SCREEN.heightCM = 30*coordinateMuilty; % cm, need to be measured on your own PC
    %     SCREEN.widthCM = 52.4*coordinateMuilty; % cm, need to measure in your own PC
    %     SCREEN.heightCM = 29.4*coordinateMuilty; % cm, need to measure in your own PC
else
    SCREEN.widthCM = 120*coordinateMuilty; % cm
    SCREEN.heightCM = 65*coordinateMuilty; % cm
end
SCREEN.distance = 60*coordinateMuilty;% cm

TRIALINFO.deviation = 1.2; % initial binocular deviation, cm
deviationAdjust     = 0.2; % how fast to adjust the deviation by key pressing, cm

% parameters for visual cue
VISUAL.headingDegree = TRIALINFO.headingDegree; % cell
VISUAL.headingDistance = TRIALINFO.headingDistance; % cell
VISUAL.headingTime = TRIALINFO.headingTime; % cell

VISUAL.fixationSizeD  = 0.25;  % degree
VISUAL.fixationWindow = 2; % degree

VISUAL.density   = 1000/(100*coordinateMuilty)^3;    % convert num/m^3 to num/cm^3
VISUAL.coherence = 100; % in percent
VISUAL.lifeTime  = 3; % frame number

VISUAL.dimensionX = 400*coordinateMuilty;  % cm
VISUAL.dimensionY = 400*coordinateMuilty;  % cm
VISUAL.dimensionZ = 700*coordinateMuilty;  % cm
VISUAL.starSize = 0.1;    % degree

% parameters for auditory cue
AUDITORY.height = 5*coordinateMuilty; % cm

AUDITORY.headingDegree = TRIALINFO.headingDegree; % cell
AUDITORY.headingDistance = TRIALINFO.headingDistance; % cell
AUDITORY.headingTime = TRIALINFO.headingTime; % cell

% % sample currently not work for double sources.
% AUDITORY.sourceNum     = {1,2};
% AUDITORY.sourceHeading = {0,[-15 30]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {10*coordinateMuilty , [5*coordinateMuilty 13*coordinateMuilty]}; % cm
% AUDITORY.sourceDegree = {0 , [-10 10]}; % degree

AUDITORY.sourceNum     = {1};
AUDITORY.sourceHeading = {0}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
AUDITORY.sourceDistance = {10*coordinateMuilty}; % cm
AUDITORY.sourceDegree = {0}; % degree for position

%% trial conditions and order
calculateConditions();
% TRIALINFO.trialConditions =
% {visualDegree visualDistance visualTime, ...
%       1               2               3
%
% auditoryDegree auditoryDistance auditoryTime sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}
%       4               5               6                7              8                9                  10

trialIndex = repmat(1:size(TRIALINFO.trialConditions,1),1,TRIALINFO.repetition);
trialNum = size(trialIndex,2);
trialOrder = randperm(trialNum);

disp(['This block has  ' num2str(trialNum) ' trials']);

timePredicted = (TRIALINFO.fixationPeriod + TRIALINFO.initialPeriod + mean(cell2mat(TRIALINFO.headingTime)) + TRIALINFO.choicePeriod + ...
    TRIALINFO.intertrialInterval ) * trialNum;
fprintf(1,'This block will cost  ');
fprintf(2,[num2str(timePredicted/60) ' '] );
fprintf(1,'minutes \n');

calibrationInterval = 600; % unit second, it is better to re-calibration every 10-15 minutes
automaticCalibration = timePredicted > 1.3*calibrationInterval; % make automatic calibration (every 10 min in default) if the block takes more than 15 min.
disp('Continue?')

% terminate the block if you feel it is too long
tic
while toc<2 % unit second
    [~, ~, keyCode]=KbCheck;
    if keyCode(escape)
        return
    end
end

%% initial opengl
if testMode
    Screen('Preference', 'SkipSyncTests', 1); % for debug/test
else
    Screen('Preference', 'SkipSyncTests', 0); % for recording
end

AssertOpenGL;
InitializeMatlabOpenGL;

if max(Screen('Screens')) > 1
    SCREEN.screenId = max(Screen('Screens'))-1;
else
    SCREEN.screenId = max(Screen('Screens'));
end
PsychImaging('PrepareConfiguration');

% Define background color:
whiteBackground = WhiteIndex(SCREEN.screenId);
blackBackground = BlackIndex(SCREEN.screenId);

% Open a double-buffered full-screen window on the main displays screen.
[win , winRect] = PsychImaging('OpenWindow', SCREEN.screenId, blackBackground);
SCREEN.widthPix = winRect(3);
SCREEN.heightPix = winRect(4);
SCREEN.center = [SCREEN.widthPix/2, SCREEN.heightPix/2];

TRIALINFO.fixationSizeP = degree2pix(TRIALINFO.fixationSizeD/2);
TRIALINFO.fixationPosition = [SCREEN.widthPix/2,SCREEN.heightPix/2];

SCREEN.refreshRate = Screen('NominalFrameRate', SCREEN.screenId);
% SCREEN.frameRate = SCREEN.refreshRate;
%% the configuration of the Frustum
calculateFrustum(coordinateMuilty);

Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
% glEnable(GL_BLEND);
% glEnable(GL_ALPHA_BLEND_CORRECTLY);
% glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('EndOpenGL', win);

GenerateStarField();

%% initial eyelink
if ~testMode
    tempName = 'TEMP1'; % need temp name because Eyelink only know hows to save names with 8 chars or less. Will change name using matlab's moveFile later.
    dummymode=0;
    
    el=EyelinkInitDefaults(win);
    %     el.backgroundcolour = BlackIndex(el.window);
    %     el.foregroundcolour = GrayIndex(el.window);
    %     el.msgfontcolour    = WhiteIndex(el.window);
    %     el.imgtitlecolour   = WhiteIndex(el.window);
    el.calibrationtargetsize=1;  % size of calibration target as percentage of screen
    el.calibrationtargetwidth=0.5; % width of calibration target's border as percentage of screen
    
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

%% initial openal
% Initialize OpenAL subsystem at debuglevel 2 with the default output device:
InitializeMatlabOpenAL(2);

nsources = max(cell2mat(AUDITORY.sourceNum));

% Generate one sound buffer:
buffers = alGenBuffers(nsources);

% Query for errors:
alGetString(alGetError)

soundFiles = dir(fullfile(pwd,'*.wav'));

alListenerfv(AL.VELOCITY, [0, 0, -1]);
alListenerfv(AL.POSITION, [0, 0, 0]);

% no idea whats this code for in OSX, but just left it here
if IsOSX
    alcASASetListener(ALC.ASA_REVERB_ON, 1);
    alcASASetListener(ALC.ASA_REVERB_QUALITY, ALC.ASA_REVERB_QUALITY_Max);
    alcASASetListener(ALC.ASA_REVERB_ROOM_TYPE, ALC.ASA_REVERB_ROOM_TYPE_Cathedral);
end

% Create a sound source:
sources = alGenSources(nsources);

% if only one source, it will have some problem in matlab,
if buffers == 0
    buffers = buffers+1;
end
if sources==0
    sources=sources+2;
end

for i = 1:nsources
    filei = mod(i,length(soundFiles))+1;
    soundName = fullfile(pwd,soundFiles(filei).name);
    
    [myNoise,freq]= psychwavread(soundName);
    %         myNoise = myNoise(:, 1);
    
    % Convert it...
    myNoise = int16(myNoise * 32767);
    
    alBufferData( buffers(i), AL.FORMAT_MONO16, myNoise, length(myNoise)*2, freq);
    
    % Attach our buffer to it: The source will play the buffers sound data.
    alSourceQueueBuffers(sources(i), 1, buffers(i));
    
    % Switch source to looping playback: It will repeat playing the buffer until its stopped.
    alSourcei(sources(i), AL.LOOPING, AL.TRUE);
end

HideCursor(SCREEN.screenId);

%% trial start
trialI = 1;
while trialI < trialNum+1
    [~, ~, keyCode]=KbCheck;
    if keyCode(escape)
        break
    end
    
    % TRIALINFO.trialConditions =
    % {visualDegree visualDistance visualTime, ...
    %       1               2               3
    %
    % auditoryDegree auditoryDistance auditoryTime ...
    %       4               5               6
    %
    % sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}
    %       7              8                9                  10
    
    conditioni = TRIALINFO.trialConditions(trialIndex(trialOrder(trialI)),:);
    visualHeadingi = cell2mat(conditioni(1:3));
    auditoryHeadingi = cell2mat(conditioni(4:6));
    auditorySourcei = conditioni(7:10);
    
    if ~sum(isnan(visualHeadingi))
        [vx,vy,vz,vfx,vfy,vfz] = calMove(visualHeadingi,SCREEN.refreshRate);
    else
        clear vx vy vz vfx vfy vfz
    end
    if ~sum(isnan(auditoryHeadingi))
        [ax,ay,az,~,~,~] = calMove(auditoryHeadingi,SCREEN.refreshRate);
    else
        clear ax ay az
    end

    % set auditory source
    for i = 1:auditorySourcei{1}
        
        % Set emission volume to 100%, aka a gain of 1.0:
        alSourcef(sources(i), AL.GAIN, 1);
        
        alSourcef(sources(i), AL.CONE_INNER_ANGLE, 360);
        alSourcef(sources(i), AL.CONE_OUTER_ANGLE, 360);
        alSource3f(sources(i), AL.DIRECTION, sind(auditorySourcei{end}(i)), 0, -cosd(auditorySourcei{end}(i)));
        
        alSource3f(sources(i), AL.POSITION, auditorySourcei{2}(i)*sind(auditorySourcei{3}(i)), 0, -auditorySourcei{2}(i)*cosd(auditorySourcei{3}(i)));
        
        % Sources themselves remain static in space:
        alSource3f(sources(i), AL.VELOCITY, 0, 0, 0);
        
        if IsOSX
            % Source emits some sound that gets reverbrated in room:
            alcASASetSource(ALC.ASA_REVERB_SEND_LEVEL, sources(i), 0.0);
        end
    end

    if exist('vx','var')
        frameNum = length(vx)-1;
    end
    
    if exist('ax','var')
        if exist('frameNum','var')
            if length(ax)-1 ~= frameNum
                frameNum = min(frameNum,length(ax)-1);
                [~, ~, ~] = DrawFormattedText(win, 'Auditory and visual cues have different duration! ','center',SCREEN.center(2)/2,[200 20 20]);
                [~, ~, ~] = DrawFormattedText(win, 'The stimulus will be given based on shorter one!!!','center',SCREEN.center(2),[200 20 20]);
                warning('Auditory and visual cues have different duration! The stimulus will be given based on shorter one!!!')
                Screen('TextBackgroundColor',win, [0 0 0 0]);
                Screen('DrawingFinished',win);
                Screen('Flip',win,0,0);
                WaitSecs(1);
            end
        else
            frameNum = length(ax)-1;
        end
    end
    
    frameTime = nan(1,frameNum);
    frameTI = GetSecs;
    aCurT = tic;
    aSt = GetSecs;
    aPosition = [0 0 0];
    
    % start giving frames
    for framei = 1:frameNum
        [~,~,keyCode] = KbCheck;
        if keyCode(escape)
            break
        elseif keyCode(skipKey)
            break
        end
        if exist('vx','var')
            % for visual cue
            if keyCode(pageUp)
                TRIALINFO.deviation = TRIALINFO.deviation + deviationAdjust;
                disp(['binocular deviation: ' num2str(TRIALINFO.deviation)]);
                calculateFrustum(coordinateMuilty);
            end
            if keyCode(pageDown)
                if TRIALINFO.deviation > deviationAdjust
                    TRIALINFO.deviation = TRIALINFO.deviation - deviationAdjust;
                    disp(['binocular deviation: ' num2str(TRIALINFO.deviation)]);
                    calculateFrustum(coordinateMuilty);
                end
            end
            
           %% draw for left eye
            Screen('BeginOpenGL', win);
            glColorMask(GL.TRUE, GL.FALSE, GL.FALSE, GL.FALSE);
            glMatrixMode(GL.PROJECTION);
            glLoadIdentity;
            glFrustum( FRUSTUM.sinisterLeft,FRUSTUM.sinisterRight, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;
            gluLookAt(vx(framei)-TRIALINFO.deviation,vy(framei),vz(framei),vx(framei)-TRIALINFO.deviation+vfx(framei),vy(framei)+vfy(framei),vz(framei)+vfz(framei),0.0,1.0,0.0);
            glClearColor(0,0,0,0);
            glColor3f(1,1,0);
            
            % draw the fixation point and 3d dots
            DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
            
           %% draw for right eye
            glColorMask(GL.FALSE, GL.TRUE, GL.FALSE, GL.FALSE);
            glMatrixMode(GL.PROJECTION);
            glLoadIdentity;
            glFrustum( FRUSTUM.dexterLeft,FRUSTUM.dexterRight, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;
            gluLookAt(vx(framei)+TRIALINFO.deviation,vy(framei),vz(framei),vx(framei)+TRIALINFO.deviation+vfx(framei),vy(framei)+vfy(framei),vz(framei)+vfz(framei),0.0,1.0,0.0);
            glClearColor(0,0,0,0);
            glColor3f(1,1,0);
            
            
            % draw the fixation point and 3d dots for right eye
            DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
            Screen('EndOpenGL', win);
            drawFixation(TRIALINFO.fixationPosition,TRIALINFO.fixationSizeP,win);
            
            Screen('Flip', win);
        else
            WaitSecs(0.01)
            drawFixation(TRIALINFO.fixationPosition,TRIALINFO.fixationSizeP,win);
            Screen('Flip', win);
        end
        
        if ~sum(isnan(auditoryHeadingi))
            % for auditory cue
            if toc(aCurT) <= auditoryHeadingi(3)
                ati = GetSecs - aSt;
                aSt = GetSecs;
                
                va = [sind(auditoryHeadingi(1))*auditoryHeadingi(2)/auditoryHeadingi(3),0,cosd(auditoryHeadingi(1))*auditoryHeadingi(2)/auditoryHeadingi(3)];
                aPosition = aPosition + va*ati;
                
                alListenerfv(AL.VELOCITY, va);
                alListenerfv(AL.POSITION, aPosition);
                
                if IsOSX
                    alcASASetListener(ALC.ASA_REVERB_ON, 1);
                    alcASASetListener(ALC.ASA_REVERB_QUALITY, ALC.ASA_REVERB_QUALITY_Max);
                    alcASASetListener(ALC.ASA_REVERB_ROOM_TYPE, ALC.ASA_REVERB_ROOM_TYPE_Cathedral);
                end
                
                if framei == 1
                    % Start playback for these sources:
                    alSourcePlayv(auditorySourcei{1}, sources(1:auditorySourcei{1}));
                end
            end
        end
        
        frameTime(framei) = GetSecs - frameTI;
        frameTI = GetSecs;
    end
    % Stop playback of all sources:
    alSourceStopv(nsources, sources);

    SCREEN.frameRate = round(1/nanmean(frameTime));
    disp(['Frame rate for this trial is ' num2str(SCREEN.frameRate)]);
    if SCREEN.refreshRate > SCREEN.frameRate
        fprintf(2,'FPS drop!!!!\n');
    end
    
    trialI = trialI +1;
    WaitSecs(TRIALINFO.intertrialInterval);
end

Screen('Flip', win);

if ~testMode
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    try
        fprintf('Receiving data file ''%s''\n',fileName);
        status=Eyelink('ReceiveFile',tempName ,saveDir,1);
        if status > 0
            fprintf('ReceiveFile status %d\n ', status);
        end
        if exist(fileName, 'file')==2
            fprintf('Data file ''%s'' can be found in '' %s\n',fileName, pwd);
        end
    catch
        fprintf('Problem receiving data file ''%s''\n',fileName);
    end
    
    cd (saveDir);
    save(fullfile(saveDir,fileName));
    movefile([saveDir,'\',tempName,'.edf'],[saveDir,'\',fileName,'.edf']);
    
    % shut down the eyelink
    Eyelink('ShutDown');
end

for i=1:nsources
    % Unqueue sound buffer:
    alSourceUnqueueBuffers(sources(i), 1, buffers(i));
end

% Wait a bit:
WaitSecs(0.1);

% Delete buffer:
alDeleteBuffers(nsources, buffers);

% Wait a bit:
WaitSecs(0.1);

% Delete sources:
alDeleteSources(nsources, sources);

% Wait a bit:
WaitSecs(0.1);

% Shutdown OpenAL:
CloseOpenAL;

Screen('CloseAll');