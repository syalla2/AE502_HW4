clear all; clc

% This script owes its bones to Example 5_11 from Curtis
global mu
deg = pi/180;
mu = 398600.8;
Re = 6378.137;
f = 1/298.26;
H = .233;
phi = 40.1164*deg;
lon = -88.2434; %deg

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

% The code below for LT was modified from the following website:
% https://smallsats.org/2013/04/11/greenwich-sidereal-time/
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

t = [t1 t2 t3];
ra = [ra_table(1) ra_table(2) ra_table(3)]*deg;
dec = [dec_table(1) dec_table(2) dec_table(3)]*deg;
theta = [LT(1) LT(2) LT(3)]*deg;

%% The below optimization is modified from Curtis Example Code 5.11
%...
%...Equations 5.64, 5.76 and 5.79:
fac1 = Re/sqrt(1-(2*f - f*f)*sin(phi)^2);
fac2 = (Re*(1-f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*sin(phi);
for i = 1:3
    R(i,1) = (fac1 + H)*cos(phi)*cos(theta(i));
    R(i,2) = (fac1 + H)*cos(phi)*sin(theta(i));
    R(i,3) = fac2;
    rho(i,1) = cos(dec(i))*cos(ra(i));
    rho(i,2) = cos(dec(i))*sin(ra(i));
    rho(i,3) = sin(dec(i));
end

%...Algorithms 5.5 and 5.6:
[r, v, r_old, v_old, f1, f3, g1, g3, ffdot1, ffdot3, ggdot1, ggdot3] = ...
    gauss(rho(1,:), rho(2,:), rho(3,:), ...
R(1,:), R(2,:), R(3,:), ...
t(1), t(2), t(3));


%...Algorithm 4.2 for the initial estimate of the state vector
% and for the iteratively improved one:
coe_old = coe_from_sv(r_old,v_old,mu);
coe = coe_from_sv(r,v,mu);



%...Echo the input data and output the solution to
% the command window:
fprintf('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
fprintf('\n Example 5.11: Orbit determination by the Gauss method\n')
fprintf('\n Radius of earth (km) = %g', Re)
fprintf('\n Flattening factor= %g', f)
fprintf('\n Gravitational parameter (km^3/s^2) = %g', mu)
fprintf('\n\n Input data:\n/');
fprintf('\n Latitude (deg)= %g', phi/deg);
fprintf('\n Altitude above sea level (km) = %g/', H);
fprintf('\n\n Observations:')
fprintf('\n Right')
fprintf('Local')
fprintf('\n Time (s) Ascension (deg) Declination (deg)')
fprintf('Sidereal time (deg)')
for i = 1:3
    fprintf('\n %9.4g %11.4f %19.4f %20.4f', ...
    t(i), ra(i)/deg, dec(i)/deg, theta(i)/deg);
end
fprintf('\n\n Solution:\n')
fprintf('\n Without iterative improvement...\n')
fprintf('\n');
fprintf('\n r (km) = [%g, %g, %g]' , ...
r_old(1), r_old(2), r_old(3))
fprintf('\n v (km/s) = [%g, %g, %g]', ...
v_old(1), v_old(2), v_old(3))
fprintf('\n');
fprintf('\n Angular momentum (km^2/s) = %g', coe_old(1))
fprintf('\n Eccentricity = %g ', coe_old(2))
fprintf('\n RA of ascending node (deg) = %g', coe_old(3)/deg)
fprintf('\n Inclination (deg) = %g', coe_old(4)/deg)
fprintf('\n Argument of perigee (deg) = %g', coe_old(5)/deg)
fprintf('\n True anomaly (deg) = %g', coe_old(6)/deg)
fprintf('\n Semimajor axis (km) = %g' , coe_old(7))
fprintf('\n Periapse radius (km) = %g', coe_old(1)^2 ...
/mu/(1 + coe_old(2)))
%...If the orbit is an ellipse, output the period:
if coe_old(2)<1
    T = 2*pi/sqrt(mu)*coe_old(7)^1.5;
    fprintf('\n Period:')
    fprintf('\n Seconds = %g', T)
    fprintf('\n Minutes = %g', T/60)
    fprintf('\n Hours = %g', T/3600)
    fprintf('\n Days = %g', T/24/3600)
end
fprintf('\n\n With iterative improvement...\n')
fprintf('\n');
fprintf('\n r (km) = [%g, %g, %g]', ...
r(1), r(2), r(3))
fprintf('\n v (km/s) = [%g, %g, %g]', ...
v(1), v(2), v(3))
fprintf('\n');
fprintf('\n Angular momentum (km^2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n RA of ascending node (deg) = %g' , coe(3)/deg)
fprintf('\n Inclination (deg) = %g', coe(4)/deg)
fprintf('\n Argument of perigee (deg) = %g' , coe(5)/deg)
fprintf('\n True anomaly (deg) = %g', coe(6)/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
fprintf('\n Periapse radius (km) = %g', coe(1)^2 ...
/mu/(1 + coe(2)))
%...If the orbit is an ellipse, output the period:
if coe(2)<1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf('\n Period:')
    fprintf('\n Seconds = %g', T)
    fprintf('\n Minutes = %g', T/60)
    fprintf('\n Hours = %g', T/3600)
    fprintf('\n Days = %g', T/24/3600)
end

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
temp9 = median(temp8);

%% SGP4 init values
sgp4init_bstar = temp9;
sgp4init_ecco  = coe(2);

initsgp4 = time_table(2) + 2400000.5 - 2433281.5;
sgp4init_epoch = initsgp4;

sgp4init_argpo = coe(5);
sgp4init_inclo = coe(4);

M = atan2(-sqrt(1 -coe(2)^2) * sin(coe(6)), -coe(2) - cos(coe(6))) + pi - ...
    coe(2) * (sqrt(1 -coe(2)^2) * sin(coe(6))) / (1 + coe(2) * cos(coe(6)));
sgp4init_mo    = M;

sgp4init_no    = sqrt(mu / (coe(7)^3));
sgp4init_nodeo = coe(3);

%% Constant from Vallado
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

error = sqrt(re_error'*re_error + dec_error'*dec_error);


%% Minimization Problem
x0= [r, v];
x_b=fmincon(@objfun,x0,[],[],[],[],[],[], @confuneq);
coe_min = coe_from_sv(x_b(1:3), x_b(4:6), mu);

% DIsplay format from Curtis
fprintf('\n\n After minimization ...\n')
fprintf('\n');
fprintf('\n r (km) = [%g, %g, %g]', ...
x_b(1), x_b(2), x_b(3))
fprintf('\n v (km/s) = [%g, %g, %g]', ...
x_b(4), x_b(5), x_b(6))
fprintf('\n');
fprintf('\n Angular momentum (km^2/s) = %g', coe_min(1))
fprintf('\n Eccentricity = %g', coe_min(2))
fprintf('\n RA of ascending node (deg) = %g' , coe_min(3)/deg)
fprintf('\n Inclination (deg) = %g', coe_min(4)/deg)
fprintf('\n Argument of perigee (deg) = %g' , coe_min(5)/deg)
fprintf('\n True anomaly (deg) = %g', coe_min(6)/deg)
fprintf('\n Semimajor axis (km) = %g', coe_min(7))
fprintf('\n Periapse radius (km) = %g', coe_min(1)^2 ...
/mu/(1 + coe_min(2)))
%...If the orbit is an ellipse, output the period:
if coe_min(2)<1
    T = 2*pi/sqrt(mu)*coe_min(7)^1.5;
    fprintf('\n Period:')
    fprintf('\n Seconds = %g', T)
    fprintf('\n Minutes = %g', T/60)
    fprintf('\n Hours = %g', T/3600)
    fprintf('\n Days = %g', T/24/3600)
end



fprintf('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')

% �����������������������������������