% osebx.txt is available in Fronter. The data is downloaded from finance.yahoo.com.

% Read the data. osebx.txt is assumed to reside in the same directory as the script.
% Note that MATLAB does not support mixed data types in numeric arrays, 
% but in cell arrays.
% The first two rows and the first column of osebx.txt are strings, so we need to
% offsets when reading the file.

% Note that it is also possible to read the data directly using the
% function csvread, but then the time vector must be specified manually
% since csvread only reads numeric data
% data = csvread('osebx.txt', 2, 1); 
% t = linspace(1996 + 2/365, 2016 + 14/365, n);

%% Open file with read premission (r)
fid = fopen('osebx.txt', 'r'); % fid = File ID

%% Skip 2 (header) lines, read comma seperated values with lines consisting
% of 1 date (%D) and 6 floating point numbers (%f)
cell_data = textscan(fid, '%D %f %f %f %f %f %f', 'Headerlines', 2, 'Delimiter', ','); 

% Close file
fclose(fid); % Semicolon suppresses the output

%% Convert from cell format to numeric format by adding the cells to the 
% columns of a matrix called data  
data = [cell_data{2}, cell_data{3}, cell_data{4}, cell_data{5}, cell_data{6}, cell_data{7}];

data(1:2, :)

% Reverse the data order (last observation first)
data = flipud(data); % flip Up Down

% Using closing prices. This is the 4th column in the data (see osebx.txt)
P = data(:, 4);

% # observations
n = length(P);
n

%% Convert time column from cell to vector and ensure it is in MATLABs
% datenum format
t = datenum(cell_data{1});
% Flip to ascending order to match data
t = flipud(t);

% Start and end time, convert date number to date string (str), print
datestr(t(1))
datestr(t(end))

%% Chunks of code can be specified by by using two percentage symbols.
%---------------------------------------------------------------------
% Code within this chunk can be run individually


% Plot. 'k-' specifies that the plot should be a black line.
figure;
plot(P, 'k-');
%---------------------------------------------------------------------

%% Time on the first axis.
% 'hold on' allows for several lines in the same plot
% note that the semicolon terminates a statement and allows multiple
% statements on a single line, the following two lines are equivalent to:
% figure;
% plot(t, P, 'k-'); 
% hold on;
figure;
plot(t, P, 'k-'); 
hold on;
line([datenum(2008,9,1) datenum(2008,9,1)], [0, 700], 'Color', 'k')

% Mark 'x' axis as time axis
datetick('x')
% Specify the axes
% Format: [ymin, ymax]
ylim([0, 700]);

%% Logarithmic prices
logP = log(P);
figure;
plot(logP, 'k-');

%% Log returns
r = diff(logP);

figure;
plot(r, 'k-');

% Specify axis labels
xlabel('Index');
ylabel('log returns');

%% Financial crisis
figure;
plot(r(3000:4000), 'k-');
xlabel('Index');
ylabel('log returns');

% Specify plot title
title('Financial crisis');

%% Use subplot to show two plots in the same window
figure;
subplot(2, 1, 1); % Arguments: number of subplots, col. position, row pos.
plot(logP, 'k-');
subplot(2, 1, 2);
plot(cumsum(normrnd(0, 1, 1, n)), 'k-');

figure;
subplot(2, 1, 1);
plot(r, 'k-');
subplot(2, 1, 2);
plot(normrnd(0, 1, 1, n), 'k-');


%%
x = linspace(0.5, 1.5, 100);
figure;
plot(x, log(x), 'k-'); hold on;
line(x, x-1, 'Color', 'k');
