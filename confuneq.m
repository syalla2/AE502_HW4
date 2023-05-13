%confuneq3b.m 
%  constraint equation for orbit fitting problem 
function [c,ceq]=confuneq3b(x)
global mu;
deg = pi/180;
mu = 398600.8;
Re = 6378.137;
f = 1/298.26;
H = .233;
phi = 40.1164*deg;
lon = -88.2434; %deg

r = x(1:3);
v = x(4:6);
deg = pi/180;
coe = coe_from_sv(r,v,mu);
%% load csv file
obs_table = table2array(readtable('AE502_HW4.csv'));

%% Time
time_table = obs_table(:, 5)';
t1_0 = time_table(1) * 24 * 3600;
t2_0 = time_table(2) * 24 * 3600;
t3_0 = time_table(3) * 24 * 3600;

t1 = t1_0 - t1_0;
t2 = t2_0 - t1_0;
t3 = t3_0 - t1_0;

ra_table  = obs_table(:, 1)'; %deg
dec_table = obs_table(:, 2)'; %deg

partial_day = time_table - fix(time_table);
day         = fix(time_table);
day_obsJ0 = day + 2400000.5;
J0 = day_obsJ0;
JC = (J0 - 2451545.0)./36525;
GST0 = 100.4606184 + 36000.77004*JC + 0.000387933*JC.^2 - 2.583e-8*JC.^3; %[deg]
GST0 = mod(GST0, 360);  % GST0 range [0..360]

GST = GST0 + 360.98564724.*partial_day;

% Correction from Greenwich to Urbana-Champaign in degrees
LT = GST   - 88.2434;

%% Extracting Avg Bstar from Celestrak TLE for OneWeb
oneweb_text = fopen('oneweb_tle.txt');
onewebtle   = textscan(oneweb_text, '%s %s %s %s %s %s %s %s %s');
temp = onewebtle{1, 7};
temp1 = temp(1:2:end, :);
temp2 = char(temp1);
temp3 = sum(isspace(temp2), 2);
temp4 = circshift(temp2(boolean(temp3), :), 1, 2);
temp5 = temp2;
temp5(boolean(temp3), :) = circshift(temp2(boolean(temp3), :), 1, 2);
temp6 = [temp5(:, 1:end-2), repmat('e', length(temp5), 1), temp5(:, end-1:end)];
temp7 = str2num(temp6);
temp8 = temp7./100000;
temp9 = mean(temp8);

%% SGP4 init values
sgp4init_bstar = temp9;
sgp4init_ecco  = coe(2);

initsgp4 = time_table(2) + 2400000.5 - 2433281.5;
sgp4init_epoch = initsgp4;

sgp4init_argpo = coe(5);
sgp4init_inclo = coe(4);


%% Actual Code
c=[-coe(7)];
ceq = [];

end