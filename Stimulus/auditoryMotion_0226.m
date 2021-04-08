% simulate auditory-visual heading perception
% add a lifetime for both visual and auditory stimulus
% environment: Matlab R2017a, Psychtoolbox 3, OpenAl
%
% ToDO:
% 1. Eyelink
% 2. Audio-visual both provided segregation condition
CloseOpenAL;
clear all STARDATA
close all

global TRIALINFO
global SCREEN
global AUDITORY
global VISUAL
global GL
global FRUSTUM
global STARDATA
global AL

subjectName = inputdlg({'Please input participant''s initials.'},'Subject Name',1,{''},'on');
if isempty(subjectName)
    return
end
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
cKey      = KbName('c'); % force calibration, temporally not in use
enter     = KbName('Return');

pageUp = KbName('pageup'); % increase binocular deviation
pageDown = KbName('pagedown'); % decrease binocular deviation

eyelinkMode = false; % 1/ture: eyelink is in recording; 0/false: eyelink is not on call
feedback = 1; % in practice block, set 1 to provide feedback. otherwise set 0
feedbackDuration = 1; % unit s

%% parameters
coordinateMuilty = 1; % convert m to coordinate system for moving distance etc.
TRIALINFO.repetition      = 15;
TRIALINFO.headingDegree   = {-20 -10 -5 -1 1 5 10 20};%{-80 -60 -30 -15 -8 -3 3 8 15 30 60 80}{-15 -7 -3 -1 1 3 7 15}
TRIALINFO.headingDistance = {0.3*coordinateMuilty};
TRIALINFO.headingTime      = {2}; % second
TRIALINFO.stimulusType     = [1]; % 0 for visual only, 1 for auditory only, 2 for both provided

TRIALINFO.choicePeriod        = 2; % second
TRIALINFO.intertrialInterval = 1; % second
TRIALINFO.fixationPeriod     = 0; % second
TRIALINFO.fixationSizeD      = 0.25; % degree

% 1 for intergration, both visual and auditory use the parameters in TRIALINFO,
% 0 for segregation, visual cue will use VISUAL, and auditory will use AUDITORY.
TRIALINFO.intergration = [1];

% for SCREEN
SCREEN.distance = 0.75*coordinateMuilty;% m
TRIALINFO.deviation = 0.001; % initial binocular deviation, m
deviationAdjust     = 0.001; % how fast to adjust the deviation by key pressing, m

% parameters for visual cue
VISUAL.headingDegree = TRIALINFO.headingDegree; % cell
VISUAL.headingDistance = TRIALINFO.headingDistance; % cell
VISUAL.headingTime = TRIALINFO.headingTime; % cell

VISUAL.fixationSizeD  = 0.5;  % degree
VISUAL.fixationWindow = 2; % degree

VISUAL.density   = 300;    % num/m^3
VISUAL.coherence = 0.2; % in percent
VISUAL.probability = VISUAL.coherence;
VISUAL.lifeTime  = 10; % frame number

VISUAL.starSize = 0.1;    % degree

% parameters for auditory cue
AUDITORY.height = 0.05*coordinateMuilty; % m

AUDITORY.headingDegree = TRIALINFO.headingDegree; % cell
AUDITORY.headingDistance = TRIALINFO.headingDistance; % cell
AUDITORY.headingTime = TRIALINFO.headingTime; % cell

% % sample currently not work for double sources.

% sourceNum = 10;
% sourceHeading = 180;
% sourceDistance = [0.1*coordinateMuilty,0.3*coordinateMuilty];
% maxDeg = 60;
% maxHeadingDeg = 20;
%
% AUDITORY.sourceNum     = {sourceNum};
% AUDITORY.sourceHeading = {180}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% for i = 1:sourceNum-1
%     AUDITORY.sourceHeading{end+1} = sourceHeading;
% end
% AUDITORY.sourceDistance = sourceDistance;
% for i = 1:sourceNum-1
%     AUDITORY.sourceDistance = [AUDITORY.sourceDistance;sourceDistance];
% end
% AUDITORY.sourceDistance = {AUDITORY.sourceDistance};
%
% AUDITORY.sourceDegree = [-maxHeadingDeg-0.01,-maxHeadingDeg;maxHeadingDeg,maxHeadingDeg+0.01];
% for i = 1:sourceNum/2-1
%     baseDeg = maxHeadingDeg + i*(maxDeg - maxHeadingDeg)/(sourceNum/2);
%     addingDeg = [baseDeg,baseDeg+0.01];
%     AUDITORY.sourceDegree = [AUDITORY.sourceDegree;addingDeg];
%     addingDeg = [-baseDeg-0.01,baseDeg];
%     AUDITORY.sourceDegree = [AUDITORY.sourceDegree;addingDeg];
% end
% AUDITORY.sourceDegree = {AUDITORY.sourceDegree};
% AUDITORY.sourceLifeTimeSplit = 1;

% % a1 = -21.01; a2 = -21;
% a1 = 21; a2 = 21.01;
a1 = -0.01; a2 = 0.01;
AUDITORY.sourceNum = {1};
AUDITORY.sourceHeading = {180}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.11*coordinateMuilty]}; % m
AUDITORY.sourceDegree = {[a1,a2]}; % degree for position
AUDITORY.sourceLifeTimeSplit = 2;

% a1 = 21; a2 = 21.01;
% AUDITORY.sourceNum     = {2};
% AUDITORY.sourceHeading = {[180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-a2,-a1; a1,a2]}; % degree for position [-55,-35;35,55] [-30,-10;10,30]
% AUDITORY.sourceLifeTimeSplit = 2;

% a1 = 21; a2 = 21.1; a3 = 25; a4 = 25.1; a5 = 30; a6 = 30.1; a7 = 35; a8 = 35.1;
% AUDITORY.sourceNum     = {8};
% AUDITORY.sourceHeading = {[180,180,180,180,180,180,180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
% AUDITORY.sourceDistance = {[0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty; 0.1*coordinateMuilty,0.11*coordinateMuilty]}; % m
% AUDITORY.sourceDegree = {[-a8,-a7; -a6,-a5; -a4,-a3; -a2,-a1; a1,a2; a3,a4; a5,a6; a7,a8]}; % degree for position
% AUDITORY.sourceLifeTimeSplit = 2;

% random seed
seed = str2double(datestr(now,'yymmddHHMM'));
rng(seed,'twister');

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

timePredicted = (TRIALINFO.fixationPeriod + mean(cell2mat(TRIALINFO.headingTime)) + TRIALINFO.choicePeriod + ...
    TRIALINFO.intertrialInterval ) * trialNum;
fprintf(1,'This block will cost  ');
fprintf(2,[num2str(timePredicted/60) ' '] );
fprintf(1,'minutes \n');

% auto-calibrate for Eyelink, temporarily not used
% calibrationInterval = 600; % unit second, it is better to re-calibration every 10-15 minutes
% automaticCalibration = timePredicted > 1.3*calibrationInterval; % make automatic calibration (every 10 min in default) if the block takes more than 15 min.
disp('Continue? Or press any key to terminate.')

% terminate the block if you feel it is too long
tic
while toc<2 % unit second
    [~, ~, keyCode]=KbCheck;
    if keyCode(escape)
        return
    end
end

%% initial opengl
Screen('Preference', 'SkipSyncTests', 0); % for recording

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

[width, height] = Screen('DisplaySize', SCREEN.screenId);
SCREEN.widthM = width/1000; % mm to m
SCREEN.heightM = height/1000; % mm to m

TRIALINFO.fixationSizeP = degree2pix(TRIALINFO.fixationSizeD/2);
TRIALINFO.fixationPosition = [SCREEN.widthPix/2,SCREEN.heightPix/2];

SCREEN.refreshRate = Screen('NominalFrameRate', SCREEN.screenId);
% SCREEN.frameRate = SCREEN.refreshRate;
%% the configuration of the Frustum
calculateFrustum(coordinateMuilty);
VISUAL.dimensionY = SCREEN.heightM/SCREEN.distance*FRUSTUM.clipFar;

[VISUAL.dimensionX, VISUAL.dimensionZ] = generateDimensionField(AUDITORY.headingDistance,...
    VISUAL.headingDegree,FRUSTUM.checkLeft,FRUSTUM.checkRight,FRUSTUM.clipFar);
auditoryLifetimeF = calculateAuditoryLifetime(AUDITORY.headingTime,AUDITORY.sourceLifeTimeSplit,SCREEN.refreshRate);

Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
% glEnable(GL_BLEND);
% glEnable(GL_ALPHA_BLEND_CORRECTLY);
% glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('EndOpenGL', win);
Screen('FillRect', win ,blackBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);

GenerateStarField();

%% initial eyelink
if eyelinkMode
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
    pause(1); % wait a little bit, in case the key press during calibration influence the following keyboard check
end

%% initial openal
% Initialize OpenAL subsystem at debuglevel 2 with the default output device:
InitializeMatlabOpenAL(2);

% nsources = max(cell2mat(AUDITORY.sourceNum));
nsources = 5;

% Generate one sound buffer:
buffers = alGenBuffers(nsources);

% Query for errors:
alGetString(alGetError)

soundFiles = dir(fullfile(pwd,'*.wav'));

alListenerfv(AL.VELOCITY, [0, 0,-1]);
alListenerfv(AL.POSITION, [0, 0, 0]);
alListenerfv(AL.ORIENTATION,[0 0 -1 0 1 0]);

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
    filei = mod(i,length(soundFiles))+1; % filei = mod(i,length(soundFiles))+1;
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
    
    % Set emission volume to 100%, aka a gain of 1.0:
    alSourcef(sources(i), AL.GAIN, 1);
    
    alSourcef(sources(i), AL.CONE_INNER_ANGLE, 360);
    alSourcef(sources(i), AL.CONE_OUTER_ANGLE, 360);
end

% HideCursor(SCREEN.screenId);

choice = zeros(trialNum,2);
choiceTime = nan(trialNum,2);
conditionIndex = cell(trialNum,size(TRIALINFO.trialConditions,2)+1);
sourceLocation= cell(trialNum,AUDITORY.sourceLifeTimeSplit);

%% trial start
trialI = 1;
while trialI < trialNum+1
    [~, ~, keyCode]=KbCheck;
    if keyCode(escape)
        break
    end
    
    % TRIALINFO.trialConditions =
    % {visualDegree visualDistance visualTime, ...
    %       1                           2                        3
    %
    % auditoryDegree auditoryDistance auditoryTime ...
    %       4                                   5                           6
    %
    % sourceNum sourceDegree(:) sourceDistance(:) sourceHead(:)}
    %       7                      8                                  9                         10
    
    conditioni = TRIALINFO.trialConditions(trialIndex(trialOrder(trialI)),:);
    visualHeadingi = cell2mat(conditioni(1:3));
    auditoryHeadingi = cell2mat(conditioni(4:6));
    auditorySourcei = conditioni(7:10);
    visualPresent = ~any(isnan(visualHeadingi));
    soundPresent = any(~isnan(auditorySourcei{end}));
    lifeSources = auditorySourcei{1};
    if visualPresent
        [vx,vy,vz,vfx,vfy,vfz] = calMove(visualHeadingi,SCREEN.refreshRate);
    else
        clear vx vy vz vfx vfy vfz
    end
    if soundPresent
        [ax,ay,az,~,~,~] = calMove(auditoryHeadingi,SCREEN.refreshRate);
    else
        clear ax ay az
    end
    
    % set auditory source
    if soundPresent
        aLifetimei = 1;
        for i = 1:auditorySourcei{1}
            alSource3f(sources(i), AL.DIRECTION, sind(auditorySourcei{end}(i)), 0, -cosd(auditorySourcei{end}(i)));
            
            zPos = randi(sort(round((-auditoryHeadingi(2)*cosd(auditoryHeadingi(1))-auditorySourcei{3}{1}(i,:))*100)))/100;
            xPos = randi(sort(round((ax(1)+auditoryHeadingi(2)*sind(auditorySourcei{2}{1}(i,:)))*100)))/100;
            
            sourceLocation{trialI,aLifetimei} = cat(1,sourceLocation{trialI,aLifetimei},[xPos,0,zPos]);
            alSource3f(sources(i), AL.POSITION, xPos, 0, zPos);
            
            % Sources themselves remain static in space:
            alSource3f(sources(i), AL.VELOCITY, 0, 0, 0);
            
            if IsOSX
                % Source emits some sound that gets reverbrated in room:
                alcASASetSource(ALC.ASA_REVERB_SEND_LEVEL, sources(i), 0.0);
            end
        end
    end
    
    % delete the frameNum from last trial
    clear frameNum
    
    if visualPresent
        frameNum = length(vx)-1;
    end
    
    if soundPresent
        if exist('frameNum','var')
            if length(ax)-1 ~= frameNum
                frameNum = min(frameNum,length(ax)-1);
            end
        else
            frameNum = length(ax)-1;
        end
    end
    
    va = [sind(auditoryHeadingi(1))*auditoryHeadingi(2)/auditoryHeadingi(3),...
        0,cosd(auditoryHeadingi(1))*auditoryHeadingi(2)/auditoryHeadingi(3)];
    alListenerfv(AL.VELOCITY, va);
    alListenerfv(AL.ORIENTATION,[0 0 -1 0 1 0]);
    alListenerfv(AL.POSITION, [0 0 0]);
    if soundPresent
        alSourcePlayv(auditorySourcei{1}, sources(1:auditorySourcei{1}));
    end
    
    frameTime = nan(1,frameNum);
    frameTI = GetSecs;
    
    % start giving frames
    for framei = 1:frameNum
        if soundPresent
            if mod(framei,auditoryLifetimeF)==0 && framei~=frameNum
                alSourceStopv(auditorySourcei{1}, sources(1:auditorySourcei{1}));
                
                aLifetimei = aLifetimei+1;
                for i = 1:auditorySourcei{1}
                    if lifeSources+i>nsources
                        lifeSources = mod(lifeSources+i,nsources);
                    end
                    alSource3f(sources(i+lifeSources), AL.DIRECTION, sind(auditorySourcei{end}(i)), 0, -cosd(auditorySourcei{end}(i)));
                    
                    zPos = randi(sort(round((-auditoryHeadingi(2)*cosd(auditoryHeadingi(1))-auditorySourcei{3}{1}(i,:))*100)))/100;
                    xPos = randi(sort(round((ax(framei)+auditoryHeadingi(2)*sind(auditorySourcei{2}{1}(i,:)))*100)))/100;
                    sourceLocation{trialI,aLifetimei} = cat(1,sourceLocation{trialI,aLifetimei},[xPos,0,zPos]);
                    alSource3f(sources(i+lifeSources), AL.POSITION, xPos, 0, zPos);
                    
                    % Sources themselves remain static in space:
                    alSource3f(sources(i+lifeSources), AL.VELOCITY, 0, 0, 0);
                    
                    if IsOSX
                        % Source emits some sound that gets reverbrated in room:
                        alcASASetSource(ALC.ASA_REVERB_SEND_LEVEL, sources(i), 0.0);
                    end
                end
                alSourcePlayv(auditorySourcei{1}, sources(lifeSources+1:auditorySourcei{1}+lifeSources));
                lifeSources = lifeSources+auditorySourcei{1};
            end
        end
        [~,~,keyCode] = KbCheck;
        if keyCode(escape)
            break
        elseif keyCode(skipKey)
            break
        end
        
         %%
        disp(['source x ' num2str(xPos) ' source z ' num2str(zPos)]);
        disp(['subj x ' num2str(ax(framei)) ' subj z ' num2str(az(framei))]);
        
        if soundPresent
            % for auditory cue
            alListenerfv(AL.POSITION, [ax(framei),ay(framei),az(framei)]);
            
            if IsOSX
                alcASASetListener(ALC.ASA_REVERB_ON, 1);
                alcASASetListener(ALC.ASA_REVERB_QUALITY, ALC.ASA_REVERB_QUALITY_Max);
                alcASASetListener(ALC.ASA_REVERB_ROOM_TYPE, ALC.ASA_REVERB_ROOM_TYPE_Cathedral);
            end
        end
        
        if visualPresent
            % for visual cue
            modifyStarField();
            if mod(framei,VISUAL.lifeTime)==0
                GenerateStarField();
            end
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
            drawFixation(TRIALINFO.fixationPosition,TRIALINFO.fixationSizeP,win);
            Screen('Flip', win);
        end
        
        frameTime(framei) = GetSecs - frameTI;
        frameTI = GetSecs;
    end
    
    
    % Stop playback of all sources:
    if soundPresent
        if AUDITORY.sourceLifeTimeSplit==1
            alSourceStopv(auditorySourcei{1}, sources(1:auditorySourcei{1}));
        elseif AUDITORY.sourceLifeTimeSplit>1
            alSourceStopv(auditorySourcei{1}, sources(lifeSources-auditorySourcei{1}+1:lifeSources));
        end
    end
    
    SCREEN.frameRate = round(1/nanmean(frameTime));
    disp(['Frame rate for this trial is ' num2str(SCREEN.frameRate)]);
    if SCREEN.refreshRate*0.95 > SCREEN.frameRate
        fprintf(2,'FPS drop!!!!\n');
    end
    
    %% start choice
    if soundPresent
        correctAnswer = (auditoryHeadingi(1) >0)+1;
        if auditoryHeadingi(1) == 0
            correctAnswer = randi(2)-1;
        end
    else
        correctAnswer = (visualHeadingi(1) >0)+1;
        if visualHeadingi(1) == 0
            correctAnswer = randi(2);
        end
    end
    startChoice = tic;
    [~, ~, ~] = DrawFormattedText(win, 'What''s your heading direction?','center',SCREEN.center(2)/2,[200 200 200]);
    Screen('TextBackgroundColor',win, [0 0 0 0]);
    Screen('DrawingFinished',win);
    Screen('Flip',win,0,0);
    while toc(startChoice) <= TRIALINFO.choicePeriod
        [ ~, ~, keyCode ] = KbCheck;
        if keyCode(leftKey)
            choice(trialI,:) = [1,trialI];
            choiceTime(trialI,:) = [toc(startChoice),trialI];
        elseif keyCode(rightKey)
            choice(trialI,:) = [2,trialI];
            choiceTime(trialI,:) = [toc(startChoice),trialI];
        end
        if choice(trialI,1)
            break
        end
    end
    if feedback
        if choice(trialI,1) == correctAnswer
            % sound(0.2*sin(2*pi*25*(1:3000)/200)); % correct cue
            [~, ~, ~] = DrawFormattedText(win, 'You are right!','center',SCREEN.center(2)/2,[20 200 20]);
            if eyelinkMode
                Eyelink('message', ['Decision made ' num2str(trialI)]);
            end
        elseif choice(trialI,1)
            % sound(0.2*sin(2*pi*25*(1:3000)/600)); % wrong cue
            [~, ~, ~] = DrawFormattedText(win, 'Please try again.','center',SCREEN.center(2)/2,[200 20 20]);
            if eyelinkMode
                Eyelink('message', ['Decision made ' num2str(trialI)]);
            end
        else
            sound(0.2*sin(2*pi*25*(1:3000)/600)); % missing cue
            [~, ~, ~] = DrawFormattedText(win, 'Oops, you missed this trial.','center',SCREEN.center(2)/2,[200 20 20]);
            if eyelinkMode
                Eyelink('message', ['Missing ' num2str(trialI)]);
            end
        end
        Screen('TextBackgroundColor',win, [0 0 0 0]);
        Screen('DrawingFinished',win);
        Screen('Flip',win,0,0);
        pause(feedbackDuration);
    else
        if choice(trialI,1)
            sound(0.2*sin(2*pi*25*(1:3000)/200)); % response cue
            if eyelinkMode
                Eyelink('message', ['Decision made ' num2str(trialI)]);
            end
        else
            sound(0.2*sin(2*pi*25*(1:3000)/600)); % missing cue
            if eyelinkMode
                Eyelink('message', ['Missing ' num2str(trialI)]);
            end
        end
    end
    if choice(trialI,1)
        conditionIndex(trialI,:) = [conditioni,trialI];
        if eyelinkMode
            Eyelink('message', ['Trial complete ' num2str(trialI)]);
        end
        trialI = trialI +1;
    else
        trialOrder = [trialOrder trialOrder(trialI)];
        trialOrder(trialI) = [];
        if eyelinkMode
            Eyelink('message', ['Trial repeat ' num2str(trialI)]);
        end
    end
    pause(TRIALINFO.intertrialInterval);
end


Screen('Flip', win);

if eyelinkMode
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
    try
        alSourceUnqueueBuffers(sources(i), 1, buffers(i));
    catch
    end
end

% Wait a bit:
pause(0.1);

% Delete buffer:
try
    alDeleteBuffers(nsources, buffers);
catch
end

% Wait a bit:
pause(0.1);

% Delete sources:
try
    alDeleteSources(nsources, sources);
catch
end

% Wait a bit:
pause(0.1);

% Shutdown OpenAL:
CloseOpenAL;

% save result
save(fullfile(saveDir,fileName),'choice','choiceTime','conditionIndex','TRIALINFO','SCREEN','AUDITORY','VISUAL','seed','sourceLocation')
Screen('CloseAll');
cd(curdir);