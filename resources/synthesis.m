function [] = synthesis(f,v)
    %%  Code to synthesize music
    %   Author: Eeyi Oon
    %   ECE 485 Final Project
    %   f is the filename of the file to read from 
    %   v is the variable name of the variable from the .mat file
    %%  Constants
    global figuresAllowed;
    figuresAllowed = 0;
    global fs; global beatToSec; global freqMultiplier; global frequencies;
    global X; global cs;
    cs = 340; % speed of sound (m/s)
    fs = 32000; % sampling frequency (Hz)
    beatToSec = 2; % (inverse seconds)
    freqMultiplier = 110*2^(29/12); % frequency of A  
    frequencies=freqMultiplier*2.^((1:20)/12); %octave of frequencies
    X = cs/fs; % sampling interval: speed of sound/fs
    
    %%  Receive matrix of data
        % extract from .mat file
        s = load(f,v);
        % extract matrix from struct
        songMatrix = getfield(s,v);
    %   first row amps, second row notes (not actual frequency), third row durations
    %songMatrix = [8,8,10,12,10,8,10,10,8,5,3,5,8,8,8,12,15,17,15,17,17,15,12,12,15,13,12,10,8,5,3,5,8,12,15,17,15,12,10,8;
    %1,0.5,0.5,1,0.75,0.25,1,0.5,0.5,1,0.5,0.5,1,1,1,0.5,0.5,2,2,1,0.5,0.5,1,0.5,0.5,0.5,0.5,0.5,0.5,1,0.5,0.5,1,0.5,0.5,1,0.5,0.5,2,2];

    %%  Creating the melody note by note 
    melody = 0;
    NoteLengths= songMatrix(2,:)./beatToSec;
    for i = 1:length(NoteLengths);
        %t1 = (0: 1/fs : NoteLengths(i) -1/fs); 
        %t1 is the time taken up by an individual note at fs, the given
        %sampling frequency
        % fprintf('hey %d',i);
        % switch to different functions to create each note
        n1 = guitar(songMatrix(1,i),NoteLengths(i));
        %TODO: don't forget to scale by amplitude 
        melody = [melody, n1];
    end 
    %keyboard
    %melody = resample(melody,fs,fs.*2);
    soundsc(melody,fs)
end

function [note] = guitar(pitch,durvec)
    %%  Constants
    global figuresAllowed;
    global fs; global beatToSec; global freqMultiplier; global frequencies;
    global X; global cs; global frequencies;
%    rel = 1.0595; %relative shift in pitch for each note
%    strlen = (cs/frequencies(1))/2; % length in meters of a guitar string
    %%  Initial condition (position)
        pluckHeight = 1; 
        %       Calculate location of pluck
            % TODO: make this work with frequencies: how to map frequency
            % to position
        %% TODO: question here: how to get a specific frequency from this 
        n =1.095;
%         waveL = cs/frequencies(pitch);
%         len = waveL/2; % location of fret in meters
        len = fs*X/frequencies(pitch);
        %len = strlen/(rel^pitch); % location of fret in meters
        xp = (3/4)*len;% location of pluck at three 
        %       assuming the initial position is a simple triangle function
        slope1 = pluckHeight/xp;
        slope2 = -pluckHeight/(len-xp);
        %       Create the simple functions
        pInit1 =@(x) slope1.*x; pInit2 =@(x) slope2.*(x-xp)+pluckHeight;
         x1   = 0:X:xp; x2 = xp:X:len;
        p0  = [pInit1(x1) pInit2(x2)];
        if (figuresAllowed); figure(1); plot(p0); end;
     %% Delay line
      i = 1; reflection = durvec*fs;
      ur = p0./2; ul = p0./2; 
      xs = round(length(p0)/2); % location of sampling
     %keyboard
     output= zeros(1,reflection);
     while (i < reflection)
         %% TODO: apply filter here
         %        if (i < 8); disp('Before Shift'); keyboard; end
        urOld = ur; ulOld = ul;
        ur = [-0.99*ulOld(1) urOld(1:length(urOld)-1)];
        %if (ur(i)==0); keyboard; end; 
        %temp1 = -0.99*ur(length(ur))
        ul = [ulOld(2:length(ulOld)) -0.99*urOld(length(urOld))];
        %if ((i>length(output))||(xs>length(ur))); keyboard; end;
        %temp3 = size(ul)
        output(i) = ur(xs)+ul(xs);
        %if (i < 8); disp('After Shift'); keyboard; end;
        %if (sum(ur) == 0); keyboard; end;
        % if (mod(i,10) == 0); keyboard; end;
        i = i+1;
        
     end
%      appFreq = fs/(len/X)
%         soundsc(sin(appFreq.*(0:reflection)),8000)
%         pause(1)
%         clear sound
        
     if (figuresAllowed); figure(3);
         hold on;
         plot(output); 
         %plot((1:length(ur)),ur,'.b',(1:length(ul)),ul,'-r');
     end;
    %keyboard
        %% TODO: end filter here 
     %% Filter 
     %keyboard
       note = output; 
       
        if (figuresAllowed); figure(2); plot(note); end;

end
