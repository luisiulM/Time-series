function xx = TouchTone(keyName,fs,fig);
% TouchTone: Create a signal vector of tones 
% which correspond to the keyName.
%
% USAGE: [xx]=TouchTone(keyName,fs);
%
% INPUT
% keyName : Valid character / key name
% fs      : Sampling frequency
%
% OUTPUT
% xx      : Signal vector that is the 
%           DTMF tone.
% 

% Define the valid key namesx
dtmf.keys = ...
['1','2','3','Y';
    '4','5','6','T';
    '7','8','9','S';
    '*','0','#','F'];

%save dtmf.keys.mat dtmf.keys;

% Check that keyName is valid!
if (str2num(keyName) > 9) | isempty(find(keyName == dtmf.keys))==1; 
disp('You blew it! - Press valid key next time...');
xx = 0;
else;

% Define tables of DTMF frequencies
dtmf.colTones = ones(4,1)*[1209,1336,1477,1633];
dtmf.rowTones = [697;770;852;941]*ones(1,4);

% Determine DTMF frequencies for keyName
[i,j] = find(keyName == dtmf.keys);

% Define duration of xx signal, T = 0.20 sec
T = 0.2;

% T equals N = T*fs + 1 samples
N = T*fs;

% discrete time index vector
n = (0:N-1);

dtmf.colTones(i,j)
dtmf.rowTones(i,j)

% Create the DTMF tone corresponding to keyName
what_1 = 2*pi*dtmf.colTones(i,j)/fs;
what_2 = 2*pi*dtmf.rowTones(i,j)/fs; 
xx = cos(what_1*n) + cos(what_2*n); 

% Play signal xx
%soundsc(xx,fs);

figure(fig);
plot(n,xx,'-*');

end;
