% This version shows checkerboard on two hemi-fields and sends LSL markers
% to the fNIRS laptop
sca; clc; close all; clear all; clearvars; 

% Break and issue an error message if the installed Psychtoolbox is not
% based on OpenGL or Screen() is not working properly.
AssertOpenGL;

% if ~IsLinux
%   error('Sorry, this demo currently only works on a Linux.');    
% end

% Set trial conditions ****************************************************
conditions = [1 2 3];
conditionNames={'Synchronization' 'Syncopation' 'Randomization'};  % conditionNames{1}
ConditionInstructions={'Press buttom when cross appear' 'Press buttom in between crosses' 'Press buttom whenever cross appear'}; 

% Block & Trial number of the experiment **********************************
% number of taps per trial/condition
numTaps=240;
% number of trials per block
numTrials=3;
% number of blocks
numBlock=1;
% total trial number
numtotal=numTrials*numBlock; 
% num of conditions in the experiment
numconditions=length(conditions);
% **********************************Block & Trial number of the experiment

% Randomization of the experiment #########################################
% set the random number seed as the date of today in formate such as 20210809
seed=input('enter the date in format YYYYMMDD:');
rng(seed);

% how many semirandom permutation set in the experiment 
numPerm=numtotal/numconditions;
% create the whole set of random conditions for the experiment
allPerm=[];
for p=1:numPerm
    allPerm=[allPerm randperm(numconditions)];
end
% ######################################### Randomization of the experiment 


% ########################################################################
try      
    
    % Here we call some default settings for setting up Psychtoolbox 
    PsychDefaultSetup(2);

    % Start with black screen
    % Removes the blue screen flash and minimize extraneous warnings.
    Screen('Preference', 'VisualDebugLevel', 1); 
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'SuppressAllWarnings', 1);

    % Get the screen numbers. This gives us a number for each of the screens
    % attached to our computer. For help see: Screen Screens?
    screens = Screen('Screens');
    % Draw we select the maximum of these numbers. So in a situation where we
    % have two screens attached to our monitor we will draw to the external
    % screen. When only one screen is attached to the monitor we will draw to
    % this. For help see: help max
    screenNumber = max(screens);

    % Define black and white (white will be 1 and black 0). This is because 
    % luminace values are (in general) defined between 0 and 1.0
    % For help see: help WhiteIndex and help BlackIndex
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    white = [white white white];
    black = [black black black];
    % Initialize some other colors
    red   = [1 0 0];
    blue  = [0 0 1];
    green = [0 1 0];
    grey  = [0.5 0.5 0.5];
    yellow  = [1 1 0];
    megenta = [1 0 1];

    % Open an on screen window and color it black
    % For help see: Screen Openwindow?
    % This will draw on a black backgroud with a size of [0 0 500 1000] and
    % return a window pointer windowPtr
%     [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black, [0 0 600 400]); 
    [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black); 

    % Get the size of the on screen windowPtr in pixels
    % For help see: Screen windowSize?
    [screenXpixels, screenYpixels] = Screen('windowSize', windowPtr);

    % Get the centre coordinate of the window in pixels
    % For help see: help RectCenter
    [xCenter, yCenter] = RectCenter(windowRect); 
    
    % Stimulus design ***********************************************
    % Just a dot
    
    % Enable alpha blending for anti-aliasing
    % For help see: Screen BlendFunction?
    % Also see: Chapter 6 of the OpenGL programming guide
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('ColorRange', windowPtr, 1, -1,1);

    % Set text display options for operating system other than Linux.
    % For help see: Screen TextFont?
    if ~IsLinux
        Screen('TextFont', windowPtr, 'Arial');
        Screen('TextSize', windowPtr, 18);
    end

    % Retreive the maximum priority number
    topPriorityLevel = MaxPriority(windowPtr) ; 
    % topPriorityLevel0 = MaxPriority(windowPtr0); 
    % set Priority once at the start of a script after setting up onscreen window.
    Priority(topPriorityLevel);
    
    % set 60-144 Hz Screen rate????????????????????????????????????????????
    % Measure the vertical refresh rate of the monitor ????????????????????
    ifi = Screen('GetFlipInterval', windowPtr);
    % Check if ifi=0.0167
    if round(1/ifi)~=60
      error('Error: Screen flash frequency is not set at 60Hz.');    
    end

    % Check if ifi=0.0083
%     if round(1/ifi)~=120
%       error('Error: Screen flash frequency is not set at 144Hz.');    
%     end

    % Check if ifi=0.0069
%     if round(1/ifi)~=144
%       error('Error: Screen flash frequency is not set at 144Hz.');    
%     end
    % ?????????????????????????????????????????????????????????????????????
    
    %  Randomization of the conditions ****************************************
    % Mean stimulus interval
    MeanTapInterval=2; % second
    NumFramesInterval=round(MeanTapInterval/ifi);  % on average 72 frames per stimulus 
    
    % condition 1
    Showframes1=[1:NumFramesInterval:NumFramesInterval*numTaps];
    
    % condition 2
    Showframes2=Showframes1;
    
    % condition 3
    RandomIntervals = round(NumFramesInterval + NumFramesInterval.*(rand(1,numTaps)-0.5)); % uniform distribution
    Showframes3=cumsum(RandomIntervals);
    
    % combine all 3 conditions
    Showframes=[Showframes1;Showframes2;Showframes3];
    % ***************************************** Randomization of the conditions 
    
    % Set size of the squares for photocell ###############################
    PhotosensorSize=30;
    % Positions of the four corners
    RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
    RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
    LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
    LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];
%     RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+230-PhotosensorSize screenXpixels screenYpixels/2+230+PhotosensorSize];
%     RightUpperSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+110-PhotosensorSize screenXpixels screenYpixels/2+110+PhotosensorSize];
%     LeftBottomSquare= [0 screenYpixels/2+230-PhotosensorSize PhotosensorSize*2 screenYpixels/2+230+PhotosensorSize];
%     LeftUpperSquare= [0 screenYpixels/2+110-PhotosensorSize PhotosensorSize*2 screenYpixels/2+110+PhotosensorSize];
    % ############################### Set size of the squares for photocell 

        
    % Numer of frames to wait when specifying good timing. Note: the use of
    % wait frames is to show a generalisable coding. For example, by using
    % waitframes = 2 one would flip on every other frame. See the PTB
    % documentation for details. In what follows we flip every frame. 
    % In order to reverse the checkerboard at 8Hz
    % We have to reverse the checkerboard every (1/8)/ifi frames
    waitframes = 1;
    
    % Setting time variables**********************************************
    % Length of one minute baseline
    BaselineLength = 2; 
    numFramesBaseline = round(BaselineLength / ifi / waitframes);

    % total number of frames per trial
%     numFrames=round(62.5/ifi/waitframes); % 125 buttom presses = about 62.5 seconds 
    % ********************************************** Setting time variables
    
    % Hide Mice:****************************************************************
    % Get handles for all virtual pointing devices, aka cursors:
    typeOnly='masterPointer'; 
    mice = GetMouseIndices(typeOnly);  
    % DuoMice: ***************************************************************
    % Hide the system-generated cursors. We do this, because only the
    % first mouse cursor is hardware-accelerated, HideCursorie., a GPU created
    % hardware cursor. All other cursors are software-cursors, created
    % by the Windowing system. These tend to flicker badly in our use
    % case. Therefore we disable all system cursor images and draw our
    % cursors ourselves for a more beautiful look:
    % Hide the cursor
    HideCursor(windowPtr,mice);
    %HideCursor(windowPtr,mice(1));

    % Baseline taking 60s  **********************************************
    instructionStart=['Hit any key and then look at the center of the screen for ' num2str(BaselineLength) ' seconds'];
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
    KbStrokeWait;
    % Create a fixation cross
    FixCrX=[xCenter-round(screenXpixels/100):xCenter+round(screenXpixels/100) repmat(xCenter,1,round(screenXpixels/50)+1)];
    FixCrY=[repmat(yCenter,1,round(screenXpixels/50)+1) yCenter-round(screenXpixels/100):yCenter+round(screenXpixels/100)];
    
    % LSL markers to the local network
    % send data into the outlet, sample by sample
%         outlet1.push_sample(1); % send data number 1
    % send markers into the outlet
%     mrk = markers{1};
%     outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array
    
%     tic
%     startTime = now;
%     numberOfSecondsElapsed = 0;
%     while numberOfSecondsElapsed < 60

    % get a timestamp and begin the baseline 
    vbl = Screen('Flip', windowPtr);
    n=1;
    while n < numFramesBaseline
        % If esc is press, break out of the while loop and close the screen
        [keyIsDown, keysecs, keyCode] = KbCheck;
        if keyCode(KbName('escape'))
            Screen('CloseAll');
            break;
        end

        % Update the while loop with time
%         numberOfSecondsElapsed = (now - startTime) * 10 ^ 5;
        
        % Show the fixation cross
        Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenXpixels/300, white, [0 0], 2);
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

        % LSL marker to check screen flip frequency
%         % send markers into the outlet
%         mrk = markers{5};
%         outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

        n=n+1;
    end
    % send markers into the outlet
%     mrk = markers{2};
%     outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array
%     BaselineDuration=toc
    

    pause(5);% to separate baseline markers and trial markers
    instructionStart=['OK. Press a key to start!']; % Tell subject to open eye and start
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
    KbStrokeWait;
    % ************************************************ Baseline taking 60s 

    % initialize some variables
    TrialDurations=[]

    %############################### Loop through block
    %*******************************Loop through block
    for block=1:numBlock

        
        %******************** Within block loop
        % how many trials 
        for t=1:numTrials 

            % Breakout by hitting esc
            [keyIsDown, keysecs, keyCode] = KbCheck;  
            if keyCode(KbName('escape'))
            Screen('CloseAll');
            break;
            end

            % Show trial and block number at the bottom
            Showtrial=['Beginning trial ' num2str(t) ' / ' num2str(numTrials) ', in block ' num2str(block) ' / ' num2str(numBlock)];
            DrawFormattedText2(Showtrial,'win',windowPtr,...
            'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);

            % flip to screen and pause for 1 sec
            vbl=Screen('Flip', windowPtr);
            pause(2); % or WaitSecs(1);

            % pick a condition from randomized set allPerm
            conditionSelected = allPerm(numTrials*(block-1)+t);
            Showframeselected=Showframes(conditionSelected,:);

            % show instruction for each trial / condition
            DrawFormattedText2(conditionNames{conditionSelected},'win',windowPtr,...
            'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
            Screen('Flip', windowPtr);
            pause(3); % for subject to read instruction

            % LSL markers to the local network
            % send data into the outlet, sample by sample
        %     outlet1.push_sample(1); % send data number 1
            % send markers into the outlet
%             mrk = markers{3};
%             outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

            % initial or reset trial frame number
            n=1;
            % get a timestamp and begin the trial
            vbl = Screen('Flip', windowPtr);
            
            % Run the trial #################
            run trial_1_f.m

            % send markers into the outlet
%             mrk = markers{4};
%             outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

            if block ~= numBlock | t ~= numTrials % only bypass the last trial
                % Show Resting
                between_block_rest= 5; % in seconds
                Resting = ['Rest for at least' num2str(between_block_rest) ' s. \n Then hit a key to continue.'];
                DrawFormattedText2(Resting,'win',windowPtr,...
                    'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
                vbl=Screen('Flip',windowPtr);
                % Rest 20 sec
                pause(between_block_rest);
            end
            % hit a key to continue
            KbStrokeWait;  

        end
        
    end
% initials trial frame number
% Show The End
TheEnd = ['The End'];
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
vbl=Screen('Flip',windowPtr);
WaitSecs(3)
% hit a key to continue
KbStrokeWait;     

%*************************************
Priority(0);   
sca;

catch
    sca;
    psychrethrow(psychlasterror);
end  
