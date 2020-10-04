function [RiffAndEcho,riff_D]=CreateRiffAndEcho(riff,a,e,td)

D = round(td*8192);

% Create "dry" sound
riff_D = riff(1:end-e*D);

% Create echoes
echoes = zeros(e, length(riff_D));
for i = 2 : e+1
   echoes(i-1,:) = riff(((i-1)*D)+1:end-((e-(i-1))*D));
end

% Initialize output
RiffAndEcho = riff_D;

% Add "dry" sound and echoes
for j = 1 : e
   RiffAndEcho = RiffAndEcho + a^j*echoes(j,:)';
end

% Play "dry" sound followed by processed
% sound.
soundsc(riff_D); 
soundsc(RiffAndEcho);
end

