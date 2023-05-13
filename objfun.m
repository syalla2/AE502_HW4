%objfun3.m 
%  Objective function for orbit fitting 
function g=objfun3(x) 
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

% If the inputs to Mean Anomoly are no long 
if(~isreal(-sqrt(1 -coe(2)^2) * sin(coe(6))) | ~isreal(-coe(2) - cos(coe(6))) | ~isreal( pi - ...
    coe(2) * (sqrt(1 -coe(2)^2) * sin(coe(6))) / (1 + coe(2) * cos(coe(6))) ))
    g = 0;
    return
end


M = atan2(-sqrt(1 -coe(2)^2) * sin(coe(6)), -coe(2) - cos(coe(6))) + pi - ...
    coe(2) * (sqrt(1 -coe(2)^2) * sin(coe(6))) / (1 + coe(2) * cos(coe(6)));

sgp4init_mo    = M;

sgp4init_no    = sqrt(mu / (coe(7)^3));
sgp4init_nodeo = coe(3);

%%

% Constant from Vallado
whichconst = 72; opsmode = 'a'; satrec.classification = 'U';
            satrec.intldesg = '44057U';
            satrec.ephtype = 0;
            satrec.elnum   = 0;
            satrec.revnum  = 0;
epoch     = sgp4init_epoch;
xbstar    = sgp4init_bstar;

% Setting below two values arbitrarily to 0;
xndot     = 0;
xnddot    = 0;

xecco     = sgp4init_ecco;
xargpo   =  sgp4init_argpo;
xinclo    = sgp4init_inclo;
xmo       = sgp4init_mo;
xno_kozai = sgp4init_no * 60;
xnodeo    = sgp4init_nodeo;
satrec = sgp4init(whichconst, opsmode, satrec, epoch, xbstar, xndot, xnddot, ...
         xecco, xargpo, xinclo, xmo, xno_kozai, xnodeo);

%% Finding r & v for all future times
time_table = obs_table(:, 5)';
time_table_minutes = time_table*24*60;
t2_minutes0 = time_table_minutes(2);
t_from_t2 = time_table_minutes - t2_minutes0;

[~, r_1predict, v_1predict] = sgp4(satrec, t_from_t2(1));
[~, r_2predict, v_2predict] = sgp4(satrec, t_from_t2(2));
[~, r_3predict, v_3predict] = sgp4(satrec, t_from_t2(3));
[~, r_4predict, v_4predict] = sgp4(satrec, t_from_t2(4));
[~, r_5predict, v_5predict] = sgp4(satrec, t_from_t2(5));
[~, r_6predict, v_6predict] = sgp4(satrec, t_from_t2(6));
[~, r_7predict, v_7predict] = sgp4(satrec, t_from_t2(7));
[~, r_8predict, v_8predict] = sgp4(satrec, t_from_t2(8));
[~, r_9predict, v_9predict] = sgp4(satrec, t_from_t2(9));


%% Turning positions into topocentric Right Ascension and Declination
LT = GST   - 88.2434;
theta = LT*deg;
fac1 = Re/sqrt(1-(2*f - f*f)*sin(phi)^2);
fac2 = (Re*(1-f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*sin(phi);
for i = 1:length(theta)
    R(i,1) = (fac1 + H)*cos(phi)*cos(theta(i));
    R(i,2) = (fac1 + H)*cos(phi)*sin(theta(i));
    R(i,3) = fac2;
end


R_predict   = [r_1predict; r_2predict; r_3predict; r_4predict; r_5predict; r_6predict; ...
    r_7predict; r_8predict; r_9predict]; 

rho_predict = R_predict - R;
mag_rho = vecnorm(rho_predict')';

l = rho_predict(:, 1)./mag_rho; m = rho_predict(:, 2)./mag_rho; 
n = rho_predict(:, 3)./mag_rho;
local_dec = asin(n);
local_ra = acos(l./cos(local_dec));
local_ra(m <= 0) = 2*pi - local_ra(m <= 0);

%% Find Errors
ra_predict_deg = local_ra*180/pi;
dec_predict_deg = local_dec*180/pi;

re_error = ra_predict_deg - ra_table';
dec_error = dec_predict_deg - dec_table';

g = sqrt(re_error'*re_error + dec_error'*dec_error);

end