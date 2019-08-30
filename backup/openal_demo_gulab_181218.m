function openal_demo_gulab_181218
% test for openal and perception for the position of voice source
% BYC Dec/2018

CloseOpenAL;


walk_degree = [0 80 89 91 100 180]';
walk_speed = 3;
trial_repeats = 2;
trial_time = 3; % unit s
choice_duration = 2; % unit s

start_position = [0 0 -10];
axis_direction = [cosd(walk_degree) , zeros(length(walk_degree),1) , sind(walk_degree)];

sourcesdir = 'D:\MATLAB\R2017a\toolbox\Psychtoolbox\PsychDemos\SoundFiles\motor_a8.wav';

% Establish key mapping: ESCape aborts, Space toggles between auto-
% movement of sound source or user mouse controlled movement:
KbName('UnifyKeynames');
space = KbName('space');
esc = KbName('ESCAPE');
kb_left = KbName('LeftArrow');
kb_right = KbName('RightArrow');

groupCondition = repmat( axis_direction , trial_repeats ,1);
degreeCondition = repmat(walk_degree, trial_repeats , 1);
[trialTotalNum,~] = size(groupCondition);
order = randperm(trialTotalNum);

% Initialize OpenAL subsystem at debuglevel 2 with the default output device:
InitializeMatlabOpenAL(2);

% Generate one sound buffer:
buffers = alGenBuffers(length(sourcesdir));

% Query for errors:
alGetString(alGetError);

% Create a sound source:
sources = alGenSources(length(sourcesdir));

% Load it...
[mynoise,freq]= psychwavread(sourcesdir);

 % Convert it...
mynoise = int16(mynoise * 32767);

% Fill our sound buffer with the data from the sound vector. Tell AL that its
% a 16 bpc, mono format, with length(mynoise)*2 bytes total, to be played at
% a sampling rate of freq Hz. The AL will resample this to the native device
% sampling rate and format at buffer load time.
alBufferData( buffers, AL.FORMAT_MONO16, mynoise, length(mynoise)*2, freq);

% Attach our buffer to it: The source will play the buffers sound data.
alSourceQueueBuffers(sources, 1, buffers);

% Switch source to looping playback: It will repeat playing the buffer until
% its stopped.
alSourcei(sources, AL.LOOPING, AL.TRUE);

% Set emission volume to 100%, aka a gain of 1.0:
alSourcef(sources, AL.GAIN, 1);

alSourcef(sources, AL.CONE_INNER_ANGLE, 30);
alSourcef(sources, AL.CONE_OUTER_ANGLE, 270);
alSource3f(sources, AL.DIRECTION, 0, 0, -1);

choice = zeros(size(order));

for i = 1 : length(order)
    
% set source position
alSource3f(sources, AL.POSITION, 0, 2, 0);

% Sources themselves remain static in space:
alSource3f(sources, AL.VELOCITY, 0, 0, 0);

% Start playback for these sources:
alSourcePlayv(length(sourcesdir), sources);


    
% Velocity of listener:
alListenerfv(AL.VELOCITY, groupCondition(order(i),:));

% Start position of listener
alListenerfv(AL.POSITION, [0, 0, -50]);

curposition = start_position;
tstart = GetSecs;
cur_time = tic;
    while toc(cur_time) < trial_time
        
        t = GetSecs;
        telapsed = t - tstart;
        tstart = t;
        tdistance = walk_speed * telapsed;
        
        curposition = curposition + groupCondition(order(i),:) * tdistance;
        
        alListenerfv(AL.POSITION, curposition);
        
        % Pause for 10 milliseconds in order to yield the cpu to other processes:
        WaitSecs(0.01);
    end
    
    % Stop playback of all sources:
    alSourceStopv(length(sourcesdir), sources);
    
    choice_time = tic;
    disp('Where are you heading to?')
    while toc(choice_time) < choice_duration
        [ ~, ~, keyCode ] = KbCheck;
        if keyCode(esc)
            break;
        end
        if keyCode(kb_left)
            choice(i) = 1;
            disp('You choice: to my LEFT.')
            if groupCondition(order(i),1) <= 0
                disp(['Yes, you are correct. Degree: ' num2str( 90 - degreeCondition(order(i)))])
            else
                disp(['Nope, you are heading to your RIGHT. Degree: ' num2str( 90 - degreeCondition(order(i)))])
            end
            WaitSecs( choice_duration - toc(choice_time) );
            break
        elseif keyCode(kb_right)
            choice(i) = 2;
            disp('You choice: to my RIGHT.')
            if groupCondition(order(i),1) >0
                disp(['Yes, you are correct. Degree: ' num2str( 90 - degreeCondition(order(i)))])
            else
                disp(['Nope, you are heading to your LEFT. Degree: ' num2str( 90 - degreeCondition(order(i)))])
            end
            WaitSecs( choice_duration - toc(choice_time) );
            break
        end

    end
end

% Stop playback of all sources:
alSourceStopv(length(sourcesdir), sources);

for i=1
    % Unqueue sound buffer:
    alSourceUnqueueBuffers(sources(i), 1, buffers(i));
end

% Wait a bit:
WaitSecs(0.1);

% Delete buffer:
alDeleteBuffers(length(sourcesdir), buffers);

% Wait a bit:
WaitSecs(0.1);

% Delete sources:
alDeleteSources(length(sourcesdir), sources);

% Wait a bit:
WaitSecs(0.1);

% Shutdown OpenAL:
CloseOpenAL;

% Done. Bye.
return;