% This file contains some basic constants for GRACE:
%  - physical
%  - geodetic (WGS84)
%  - sun, moon, planets
%
% Matthias Weigelt

% Some numbers (1986 recommended values, see Physics Today, 08/94, part 2)
clight    =   2.99792458e8;          % speed of light [m/s]
G         =   6.67259e-11;           % gravitational constant [m^3 /(kg s^2)]
au        = 149.597870691e9;         % astronomical unit [m]

% WGS84 defining constants:
ae        =  6378136.3;              % semi-major axis of ellipsoid [m]
GM        =  3.986004415e14;         % geocentric grav. constant [m^3 / s^2]
J2        =  1.08262982131e-3;       % earth's dyn. form factor (= -C20 unnormalized)
Omega     =  7.292115085e-5;         % mean ang. velocity [rad/s]

% WGS84 derived constants:
flat      =  1/298.25765;            % flattening
J4        = -2.37091120053e-6;       % -C4,0  unnormalized
J6        =  6.08346498882e-9;       % -C6,0  unnormalized
J8        = -1.42681087920e-11;      % -C8,0  unnormalized
J10       =  1.21439275882e-14;      % -C10,0 unnormalized

% other constants:
Omega_dot = -4.5e-22;                % ang. accleration [rad/s^-2]
g_earth   =  9.780327e0;             % mean Earth gravity in [m/s^-2]
C20_dot   =  1.162755e-11;           % C2,0 normalized 
C21_dot   = -0.337e-11;              % C2,1 normalized
S21_dot   =  1.606e-11;              % S2,1 normalized

% Gravitational constants of the planets, sun and moon 
GMmoon    =  4902818954300;          % grav. constant of the Moon  [m^3 / s^2]
GMsun     =  1.326801149665e+020;    % grav. constant of the Sun   [m^3 / s^2]
GMmercury =  22032892180000;         % grav. constant of Mercury   [m^3 / s^2]
GMvenus   =  324855044150000;        % grav. constant of Venus     [m^3 / s^2]
GMearth   =  3.986004415e14;         % grav. constant of Earth     [m^3 / s^2]
GMmars    =  42828018915000;         % grav. constant of Mars      [m^3 / s^2]
GMjupiter =  1.267124841e+017;       % grav. constant of Juptier   [m^3 / s^2]
GMsaturn  =  3.7931005114e+016;      % grav. constant of Saturn    [m^3 / s^2] 
GMuranus  =  5.7939433488e+015;      % grav. constant of Saturn    [m^3 / s^2] 
GMneptune =  6.834733937e+015;       % grav. constant of Saturn    [m^3 / s^2] 
GMpluto   =  870772995000;           % grav. constant of Saturn    [m^3 / s^2] 

% calculation constants
rho     = 180/pi;                    % convert from rad to deg   
sec2deg = 15/3600;                   % convert from time seconds to hours to degree (where 24h = 360deg).

