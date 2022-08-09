% simulate auditory-visual heading perception
% add a lifetime for both visual and auditory stimulus
% environment: Matlab R2017a, Psychtoolbox 3, OpenAl
%
% ToDO:
% 1. Eyelink
% 2. Audio-visual both provided segregation condition
CloseOpenAL;
clear TRIALINFO SCREEN AUDITORY VISUAL GL FRUSTUM STARDATA AL
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
if strcmp(cell2mat(subjectName), 'test')
    testMode = true;
else
    testMode = false;
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
TRIALINFO.repetition      =10;
TRIALINFO.headingDegree   = {-15,15};
TRIALINFO.headingDistance = {0.3*coordinateMuilty};
TRIALINFO.headingTime      = {1}; % second
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
VISUAL.headingDegreeDelta = {0 20 -20 40 -40}; % delta degree for segregation condition

VISUAL.headingDistance = TRIALINFO.headingDistance; % cell
VISUAL.headingTime = TRIALINFO.headingTime; % cell

VISUAL.fixationSizeD  = 0.5;  % degree
VISUAL.fixationWindow = 2; % degree

VISUAL.density   = 300;    % num/m^3
VISUAL.coherence = 0.5; % in percent
VISUAL.probability = VISUAL.coherence;
VISUAL.lifeTime  = 3; % frame number

VISUAL.starSize = 0.1;    % degree

% parameters for auditory cue
AUDITORY.height = 0.05*coordinateMuilty; % m

AUDITORY.headingDegree = TRIALINFO.headingDegree; % cell
AUDITORY.headingDistance = TRIALINFO.headingDistance; % cell
AUDITORY.headingTime = TRIALINFO.headingTime; % cell

% % sample currently not work for double sources.
a1 = -10; a2=-10.1 ;a3=10;a4=10.1;a5=30;a6=30.01;a7=5;a8=5.01;a9=25;a10=25.01;
AUDITORY.sourceNum     = {20};
AUDITORY.sourceHeading = {[180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180]}; % degree, 0 for [0 0 -z], 90 for [x 0 0], -90 for [-x 0 0], 180 for [0 0 +z]
AUDITORY.sourceDistance = {[0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty;0.3*coordinateMuilty,0.4*coordinateMuilty]};%{[-2.7*coordinateMuilty,-2.7*coordinateMuilty;-1.7*coordinateMuilty,-1.7*coordinateMuilty;-0.7*coordinateMuilty,-0.7*coordinateMuilty;0.30*coordinateMuilty,0.30*coordinateMuilty]}; % m
AUDITORY.sourceDegree = {[-25,-24.9;-22.5,-22.4;-20,-19.9;-17.5,-17.4;-15,-14.9;-12.5,-12.4;-10,-9.9;-7.5,-7.4;-5,-4.9;-2.5,-2.4;2.4,2.5;4.9,5;7.4,7.5;9.9,10;12.4,12.5;14.9,15;17.4,17.5;19.9,20;22.4,22.5;24.9,25]}; % degree for position [-55,-35;35,55] [-30,-10;10,30]
AUDITORY.sourceLifeTimeSplit =1;%cell2mat(AUDITORY.sourceNum);

% parameter for coherence
AUDITORY.coherence =0.8; % the influenced sources number = round( (1-coherence) * courceNum )
AUDITORY.coherenceDirection = 0; % 0 random, 1 same as heading side in x axis
AUDITORY.coherenceVelocity = 1; % how many times of the heading velocity in x-axis component

% parameter for lift time
AUDITORY.sourceInitial = 0.03; % second
AUDITORY.sourceTerminal = 0.03; % second
AUDITORY.sourceDuration = max(cell2mat(AUDITORY.headingTime)); % second

if AUDITORY.sourceInitial+AUDITORY.sourceTerminal>AUDITORY.sourceDuration
    error('Invalid number for sourceInitial, sourceTerminal or sourceDuration.')
end

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
disp('Continue? Or press ESC to terminate.')

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
    Screen('Preference', 'SkipSyncTests', 1); % for test
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

[width, height] = Screen('DisplaySize', SCREEN.screenId);
SCREEN.widthM = width/1000; % mm to m
SCREEN.heightM = height/1000; % mm to m

TRIALINFO.fixationSizeP = degree2pix(TRIALINFO.fixationSizeD/2);
TRIALINFO.fixationPosition = [SCREEN.widthPix/2,SCREEN.heightPix/2];

SCREEN.refreshRate = Screen('NominalFrameRate', SCREEN.screenId);

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

nsources = max(cell2mat(AUDITORY.sourceNum));

muiltyInitialTime = linspace(0, max(cell2mat(AUDITORY.headingTime))-AUDITORY.sourceDuration,nsources);

HideCursor(SCREEN.screenId);

choice = zeros(trialNum,2);
choiceTime = nan(trialNum,2);
conditionIndex = cell(trialNum,size(TRIALINFO.trialConditions,2)+1);
sourceLocation= cell(trialNum,nsources);

%% trial start
trialI = 1;
sourceFileList = cell(nsources,trialNum);

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
    soundPresent = ~any(isnan(auditorySourcei{1}));
    
    if visualPresent
        [vx,vy,vz,vfx,vfy,vfz] = calMove(visualHeadingi,SCREEN.refreshRate);
    else
        clear vx vy vz vfx vfy vfz
    end
    
    % set auditory source
    if soundPresent
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
        
        joker = randi(length(soundFiles));
        for i = 1:nsources
            filei = mod(i+joker,length(soundFiles))+1;
            soundName = fullfile(pwd,soundFiles(filei).name);
            sourceFileList{i,trialI} = soundFiles(filei).name;
            [fileSample,freq]= psychwavread(soundName);
            
            initialAmp = linspace(0,1,freq*AUDITORY.sourceInitial);
            terminalAmp = linspace(1,0,freq*AUDITORY.sourceTerminal);
            fileSample(1:length(initialAmp),:) = fileSample(1:length(initialAmp),:).*initialAmp';
            fileSample(end-length(terminalAmp)+1:end,:) = fileSample(end-length(terminalAmp)+1:end,:).*terminalAmp';
            myNoise  = [zeros(round(muiltyInitialTime(i)*freq),1);fileSample;zeros(round((max(cell2mat(AUDITORY.headingTime))-AUDITORY.sourceDuration-muiltyInitialTime(i))*freq),1)];
            
            % Convert it...
            myNoise = int16(myNoise * 32767);
            
            alBufferData( buffers(i), AL.FORMAT_MONO16, myNoise, length(myNoise)*2, freq);
            
            % Attach our buffer to it: The source will play the buffers sound data.
            alSourceQueueBuffers(sources(i), 1, buffers(i));
            
            % Switch source to looping playback: It will repeat playing the buffer until its stopped.
            alSourcei(sources(i), AL.LOOPING, AL.FALSE);
            
            % Set emission volume to 100%, aka a gain of 1.0:
            %     if nsources == 1
            %         alSourcef(sources(i), AL.GAIN, 1);
            %     else
            %         if i == 1
            %             alSourcef(sources(i), AL.GAIN, 1);
            %         elseif i ==2
            %             alSourcef(sources(i), AL.GAIN, 0.5);
            %         end
            %     end
            alSourcef(sources(i), AL.GAIN, 0.1);
            alSourcef(sources(i), AL.CONE_INNER_ANGLE, 360);
            alSourcef(sources(i), AL.CONE_OUTER_ANGLE, 360);
        end
    end
    if soundPresent
        [ax,ay,az,~,~,~] = calMove(auditoryHeadingi,SCREEN.refreshRate);
        audioCoherenceNum = round(auditorySourcei{1}*(1-AUDITORY.coherence));
        audioCoherenceIndex = randperm(auditorySourcei{1},audioCoherenceNum);
        if AUDITORY.coherenceDirection == 0
            seta = rand;
            audioCoherenceDirection = [sin(seta*2*pi), 0, cos(seta*2*pi)];
        elseif AUDITORY.coherenceDirection == 1
            audioCoherenceDirection = [sind(auditoryHeadingi(1)), 0, 0];
        end
        sourceV = audioCoherenceDirection.*(auditoryHeadingi(2)/auditoryHeadingi(3)*AUDITORY.coherenceVelocity);
    else
        clear ax ay az
    end
    
    if soundPresent
        sourcePosition = cell(auditorySourcei{1},1);
        for i = 1:auditorySourcei{1}
            alSource3f(sources(i), AL.DIRECTION, sind(auditorySourcei{end}(i)), 0, -cosd(auditorySourcei{end}(i)));
            
            zPos = randi(sort(round((-auditoryHeadingi(2)*cosd(auditoryHeadingi(1))-auditorySourcei{3}{1}(i,:))*100)))/100;
            xPos = randi(sort(round((ax(1)+auditoryHeadingi(2)*sind(auditorySourcei{2}{1}(i,:)))*100)))/100;
            sourcePosition{i} = [xPos, 0, zPos];
            sourceLocation{trialI,i} = cat(1,sourceLocation{trialI,i},[sourcePosition{i}, 0]);
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
        sourceMovingInitialF = nan(size(muiltyInitialTime));
        sourceMovingTerminalF = nan(size(muiltyInitialTime));
        sourceMovingInitialF(audioCoherenceIndex) = round(muiltyInitialTime(audioCoherenceIndex).*SCREEN.refreshRate);
        sourceMovingTerminalF(audioCoherenceIndex) = sourceMovingInitialF(audioCoherenceIndex)+round(AUDITORY.sourceDuration*SCREEN.refreshRate);
        sourceMovingIndex = nan(size(muiltyInitialTime));
        sourceMovingIndex(sourceMovingInitialF==0)=1;
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
    
    startChoice = tic;
    % start giving frames
    for framei = 1:frameNum
        if visualPresent
            modifyStarField();
            if mod(framei,VISUAL.lifeTime)==0
                GenerateStarField();
            end
        end
        if soundPresent
            if mod(framei,auditoryLifetimeF)==0 && framei~=frameNum
                for i = 1:auditorySourcei{1}
                    alSource3f(sources(i), AL.DIRECTION, sind(auditorySourcei{end}(i)), 0, -cosd(auditorySourcei{end}(i)));
                    
                    zPos = randi(sort(round((-auditoryHeadingi(2)*cosd(auditoryHeadingi(1))-auditorySourcei{3}{1}(i,:))*100)))/100;
                    xPos = randi(sort(round((ax(framei)+auditoryHeadingi(2)*sind(auditorySourcei{2}{1}(i,:)))*100)))/100;
                    sourcePosition{i} = [xPos, 0 ,zPos];
                    sourceLocation{trialI,i} = cat(1,sourceLocation{trialI,i},[sourcePosition{i}, framei]);
                    alSource3f(sources(i), AL.POSITION, xPos, 0, zPos);
                    
                    % Sources themselves remain static in space:
                    alSource3f(sources(i), AL.VELOCITY, 0, 0, 0);
                    
                    if IsOSX
                        % Source emits some sound that gets reverbrated in room:
                        alcASASetSource(ALC.ASA_REVERB_SEND_LEVEL, sources(i), 0.0);
                    end
                end
            end
            if ismember(framei,sourceMovingInitialF)
                iList = find(sourceMovingInitialF == framei);
                sourceMovingIndex(iList) = 1;
            end
            if ismember(framei,sourceMovingTerminalF)
                tList = find(sourceMovingTerminalF == framei);
                sourceMovingIndex(tList) = 0;
            end
            for i = 1:auditorySourcei{1}
                if sourceMovingIndex(i) == 1
                    sourcePosition{i} = sourcePosition{i}+sourceV./SCREEN.refreshRate;
                    sourceLocation{trialI,i} = cat(1,sourceLocation{trialI,i},[sourcePosition{i}, framei]);
                    alSource3f(sources(i), AL.POSITION, sourcePosition{i}(1), sourcePosition{i}(2), sourcePosition{i}(3));
                    
                    % Sources themselves remain static in space:
                    alSource3f(sources(i), AL.VELOCITY, sourceV(1), sourceV(2), sourceV(3));
                elseif sourceMovingTerminalF(i) == framei
                    alSource3f(sources(i), AL.VELOCITY, 0, 0, 0);
                end
            end
        end
        [~,~,keyCode] = KbCheck;
        if keyCode(escape)
            break
        elseif keyCode(skipKey)
            break
        end
        
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
        
        % RT mode, 
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
        
        frameTime(framei) = GetSecs - frameTI;
        frameTI = GetSecs;
    end
    
    
    % Stop playback of all sources:
    if soundPresent
        alSourceStopv(auditorySourcei{1}, sources(1:auditorySourcei{1}));
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
    
    %     startChoice = tic;
    if ~choice(trialI,1)
        [~, ~, ~] = DrawFormattedText(win, 'What''s your heading direction?','center',SCREEN.center(2)/2,[200 200 200]);
        Screen('TextBackgroundColor',win, [0 0 0 0]);
        Screen('DrawingFinished',win);
        Screen('Flip',win,0,0);
        
        timei = [conditioni{3}, conditioni{6}];
        timei = max(timei(~isnan(timei)));
        while toc(startChoice) <= TRIALINFO.choicePeriod+timei
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
        if soundPresent
            for i=1:size(sourceLocation,2)
                sourceLocation{trialI,i} = [];
            end
        end
        if eyelinkMode
            Eyelink('message', ['Trial repeat ' num2str(trialI)]);
        end
    end
    tti = tic;
    if soundPresent
        if true%~bypassSourceGenerate
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
        end
    end
    pause(TRIALINFO.intertrialInterval-toc(tti));
    %pause(TRIALINFO.intertrialInterval);
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
save(fullfile(saveDir,fileName),'choice','choiceTime','conditionIndex','TRIALINFO','SCREEN','AUDITORY','VISUAL','seed','sourceLocation','muiltyInitialTime','sourceFileList')
Screen('CloseAll');
cd(curdir);