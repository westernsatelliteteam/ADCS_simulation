Ts = 1;

%% Set mass properties
vehicle.mass = 1;                      % CubeSat mass [kg]
vehicle.inertia = eye(3);              % CubeSat Moments of Inertia [kg*m^2]

%% Set controller properties
gains.Kp = 0.00001;                    % Controller Proportional Gain
gains.Ki = 0.0000000001;               % Controller Integral Gain
gains.Kd = 0.01;                       % Controller Derivative Gain

%% Set default initial orbital state and attitude
sim_t0Vec = [2019 12 31 12 0 0];
mjd = mjuliandate(sim_t0Vec);
sim_t0 = juliandate(sim_t0Vec);            % Simulation start date
epochVec = [2000 1 1 12 0 0];
epoch = 2451545.0;                     % Epoch [Julian date] (J2000)
dAT = 37;                              % Difference between TAI and UTC [s]
dUT1 = deltaUT1(mjd,'action','none');  % Difference between UTC and UT1 [s]
% Conpute transformation between J2000 and MOD ECI frames
R_ijk2j2000 = dcmIJK2J2000(epochVec, dAT)';

a = 6786233.13;                        % Semi-major Axis [m]
ecc = 0.0010537;                       % Eccentricity
incl = 51.7519;                        % Inclination [deg]
RAAN = 95.2562;                        % Right ascension of the ascending node [deg]
argp = 93.4872;                        % Argument of perigee [deg]
nu  = 302.9234;                        % True anomoly [deg]
euler = [0 0 0];                       % Euler angles [deg];
pqr = [0 0 -0.05168];                  % Body angular rates [deg/s];
pm = polarMotion(mjd,'action','none')*180/pi; % Polar displacement in x and y axis [deg deg]
dCIP = deltaCIP(mjd,'action','none')*180/pi;  % Adjustment to CIP in x and y axes [deg deg]
lod = 0;                               % Excess length of day [s]

% Kepler to RV_eci
small = 1e-12;
if ( ecc < small )
    if incl < small || abs(incl-pi)< small % circular equatorial orbit
        [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, 0, 0,...
            0, 'truelon', truelon);
    else % circular inclined
        [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, RAAN, 0,...
            0, 'arglat', arglat);
    end
else
    if ( ( incl<small) || (abs(incl-pi)<small) ) % elliptical equatorial
        [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, 0, 0,...
            nu, 'lonper', lonper);
    else % elliptical inclined
        [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, RAAN, argp,...
            nu);
    end
end

% RV_tod to RV_j2000
r_j2000 = R_ijk2j2000*r_ijk;
v_j2000 = R_ijk2j2000*v_ijk;

% RV_j2000 to RV_ecef
[r_ecef, v_ecef] = eci2ecef(sim_t0Vec, r_j2000, v_j2000, 'dAT', dAT,...
    'dUT1', dUT1, 'pm', pm, 'dCIP', dCIP, 'lod', lod);

% R_ecef to R_lla
lla = ecef2lla(r_ecef(:)');

% V_ecef to V_ned
v_ned = dcmecef2ned(lla(1), lla(2))*v_ecef(:);

% V_ned to V_body
uvw = angle2dcm(euler(3)*pi/180, euler(2)*pi/180, euler(1)*pi/180, 'ZYX')*v_ned;

[~, thGAST] = greenwichSRT(sim_t0, dUT1, dAT);

% Time Parameters
initCond.simStartDate.JD = sim_t0;
initCond.simStartDate.dateVector = sim_t0Vec;
initCond.CoordEpoch.JD = epoch;
initCond.CoordEpoch.dateVector = epochVec;

% Earth Properties
initCond.EarthProps.LG = thGAST;

% Vehicle State
initCond.lla = lla(:)';
initCond.v_ned = v_ned(:)';
initCond.uvw = uvw(:)';
initCond.euler = euler(:)';
initCond.pqr = pqr(:)';

clearvars -except initCond gains Ts vehicle

load('buses.mat');

setupCubeSatVisualization();

function dcm = dcmIJK2J2000(epochVec, dAT)
%Compute tranformation from ECI with mean equinox at epoch to J2000

% Seconds for UTC
ssTT = epochVec(end) + dAT + 32.184;
% Julian date for terrestrial time
jdTT = mjuliandate(epochVec(1),epochVec(2),epochVec(3),epochVec(4),epochVec(5),ssTT);
% Number of Julian centuries since J2000 for terrestrial time.
tTT = (jdTT - 51544.5)/36525;
tTT2 = tTT.*tTT;
tTT3 = tTT2.*tTT;
% Zeta, theta and z represent the combined effects of general precession
zeta = convang((2306.2181*tTT + 0.30188*tTT2 + 0.017998*tTT3)/3600,'deg','rad');
theta = convang((2004.3109*tTT - 0.42665*tTT2 - 0.041833*tTT3)/3600,'deg','rad');
z = convang((2306.2181*tTT + 1.09468*tTT2 + 0.018203*tTT3)/3600,'deg','rad');
% ECI vector with mean equinox at epoch to ECI vector in J2000
dcm = angle2dcm(-zeta,theta,-z,'ZYZ')';
end