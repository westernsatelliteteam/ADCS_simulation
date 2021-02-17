clear;
clc;

%% Time Step
Ts = 1;

%% Set mass properties
vehicle.mass = 1;                      % CubeSat mass [kg]
vehicle.inertia = eye(3);              % CubeSat Moments of Inertia [kg*m^2]

%% Set controller properties
gains.Kp = 0.00001;                    % Controller Proportional Gain
gains.Ki = 0.0000000001;               % Controller Integral Gain
gains.Kd = 0.01;                       % Controller Derivative Gain

%% Set Visualization
variantVisualization = 0;
visOff = Simulink.Variant('variantVisualization == 0');
visSL3D = Simulink.Variant('variantVisualization == 1');

%% Startup
calcInitialConditions();

load('buses.mat');