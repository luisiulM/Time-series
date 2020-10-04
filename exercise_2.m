%% 1.5a ???

t_s = 1;
n = 200;
step = 1;

time = t_s:step:n;
w_t = transpose(wgn(n,1000,0));
s_t = zeros(1,n);

for j = 1:n
    if j < 101
       s_t(j) = 0;
    else
       s_t(j) = 10*exp(-(j - 100)/20)*cos(2*pi*j/4); 
    end
end

x_t = repmat(s_t,1000,1) + w_t;

mean_x = @(t) mean(x_t(t,:));

plot(time,mean(x_t, 1))
hold on
%plot(time,mean_x(time),'r')

%% 1.5b ???

%y = Autocov(x_t,n)

Cov(h)=0;

for h=1:n
    
for i=1:(n-h+1)    
   Cov(h)=Cov(h)+( x_t(i+h-1)-mean_x(i+h-1) ).*(( x_t(i)-mean_x(i) ));
end
end

y = Cov/n;

plot(time,y)

%% 1.20a ???



%% 1.20b ???



%% 3

fid = fopen('osebx.txt', 'r'); % fid = File ID

cell_data = textscan(fid, '%D %f %f %f %f %f %f', 'Headerlines', 2, 'Delimiter', ',');
fclose(fid); % Semicolon suppresses the output

% columns of a matrix called data  
data = [cell_data{2}, cell_data{3}, cell_data{4}, cell_data{5}, cell_data{6}, cell_data{7}];
data(1:2, :)

% Reverse the data order (last observation first)
data = flipud(data); % flip Up Down

% Using closing prices. This is the 4th column in the data (see osebx.txt)
P = data(:, 4);
