% time for resting EEG (EO=eye open; EC= eye close)
TimeRestingEEG=150; % 2.5 min = 150 seconds 
numFramesRestEye=round (TimeRestingEEG/ifi); 

%*************************************************************************
% Start taking eye open and eye close resting stage EEG
instructionStart=['Hit any key and then look at the center of the screen for ' num2str(TimeRestingEEG/60) ' min']
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
Screen('Flip',windowPtr);
% hit a key to continue
KbStrokeWait;
pause(2); % given time to return the position of arm after key press 

% Flash photocells once at the start and end of eye open period 
% Flash photocells once at the start of eye open period ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, RightBottomSquare);
Screen('Flip',windowPtr);



% get a timestamp and begin taking resting EEG ****************************
vbl = Screen('Flip', windowPtr);
i=1;
while i<numFramesRestEye
% If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    % Show the fixation cross
    Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenXpixels/400, white, [0 0], 2);
    % Show the central dot
    % Screen('DrawDots', windowPtr, [xCenter;yCenter], Thickness, white, [0 0], 2);
    % Flip the black screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    % update the while loop
    i=i+1;
end

% Flash photocells once to mark the end of open eye resting EEG collection ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, RightBottomSquare);
Screen('Flip',windowPtr);

%****************************************************************************
instructionStart=['Hit a key and then close your eyes to rest for ' num2str(TimeRestingEEG/60) ...
    ' min. %\n Open your eye when you hear a beep.']
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
Screen('Flip',windowPtr);
% hit a key to continue
KbStrokeWait;
pause(2); % given time to return the position of arm after key press 

% Flash photocells to mark the start and end of eye closed period 
% Flash once to start ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, RightBottomSquare);
Screen('Flip',windowPtr);

% get a timestamp and begin taking resting EEG
vbl = Screen('Flip', windowPtr);
i=1;
while i<numFramesRestEye
% If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    % Flip the black screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    % update the while loop
    i=i+1;
end
    

% Flash photocells again to end ###### (upper left and bottom right)
i=1;
while i < 2;
    Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    Screen('Flip',windowPtr);
    i=i+1;
end


% ************************************************************************
% play a beep sound to signal subject to open eyes
Beeper(400,0.9,1.5);
% show the end
instructionStart=['OK. Press a key to start!'] % Tell subject to open eye and start
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
Screen('Flip',windowPtr);
% hit a key to continue
KbStrokeWait;
pause(2); % given time to return the position of arm after key press 
%*************************************************************************
