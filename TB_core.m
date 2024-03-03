% Initialization
tic
format long

% Physical constants and initial parameters under modern climate
rho0 = 1027; % Reference density (kg/m^3)
S0 = 35; % Reference salinity (psu)
cp = 4000 * rho0; % (J/m^3/K)
betat = 1.8 * 10^(-4); betatL = 1.8 * 10^(-4); betatH = 1 * 10^(-4); % Thermal expansion coefficient (K^-1)
betas = 7.6 * 10^(-4); % Saline expansion coefficient (psu^-1)
qe = 16 * 10^6; % Mean state AMOC across OSNAP section under modern climate (m^3/s)
ht = 0.37 * 10^15; % Mean state AMOC-related heat transport across OSNAP section under modern climate (W)
Hse = 7.5 * 10^6; % % Mean state AMOC-related salt transport across OSNAP section under modern climate Hse = S0 * Fe (psu*m^3/s)
dTe = ht / (cp * qe); % (K or ^oC)
dSe = Hse / qe; % from equilibrium (psu)
V1e = 3.2 * 10^16; % Low latitude box volume (m^3)
V2e = 3.8 * 10^15; % High latitude box volume (m^3)
A1 = 1.064 * 10^13; % Low latitude box area (m^2)
A2 = 0.618 * 10^13; % High latitude box area (m^2)

% Derived parameters
restor_time2 = V2e / A2 / 0.7 / 365; % High latitude box thermal restoring time 1/gamaH (year)
gamat = V2e / A2 / restor_time2 / (3600 * 24 * 365); % (m/s)
restor_time1 = 1 / gamat / (3600 * 24 * 365) * V1e / A1; % Low latitude box thermal restoring time 1/gamaL (year)
dT_a = (A1 + A2) / (A1 * A2) * ht / cp / gamat + dTe; % Prescribed atmosphere temperature difference between two boxes (K or ^oC)
T_a_tot = 12.5; % Prescribed atmosphere temperature sum of two boxes (^oC)
T1_a = 0.5 * (dT_a + T_a_tot); % Prescribed atmosphere temperature in low latitude box (^oC)
T2_a = 0.5 * (T_a_tot - dT_a); % Prescribed atmosphere temperature in high latitude box (^oC)
T1e = T1_a - qe * dTe / (A1 * gamat); % Ocean temperature in low latitude box (^oC)
T2e = T2_a + qe * dTe / (A2 * gamat); % Ocean temperature in high latitude box (^oC)
S1e = S0 + V2e / (V1e + V2e) * dSe; % Ocean salinity in low latitude box (psu)
S2e = S0 - V1e / (V1e + V2e) * dSe; % Ocean salinity in high latitude box (psu)
drhoe = rho0 * (betas * dSe - betat * dTe); % Density difference between two boxes (kg/^m3)
rho1e = rho0;  % Ocean density in low latitude box (kg/^m3)
rho2e = rho1e - drhoe; % Ocean density in high latitude box (kg/^m3)
k0 = -qe / drhoe; % Coefficient between AMOC and density differences ((m^3*s-1)/(kg*m^-3))
OSTe = qe .* dSe; % Ocean salt transport (m^3/s*psu)
FWe = Hse / S0; % Ocean freshwater transport (m^3/s)
ratio = V2e / V1e; % volume ratio
V1 = V1e; % %Low latitude box volume, same with V1e
V2 = V1 * ratio; % High latitude box volume, same with V2e

% Time-stepping and RK-4 scheme initialization
dt = 0.1; % time step (year)
T = 3000; % simulation length (year)
numSteps = T / dt + 1; % Number of time steps
% Pre-allocate arrays
T1 = zeros(numSteps, 1); T2 = zeros(numSteps, 1); 
S1 = zeros(numSteps, 1); S2 = zeros(numSteps, 1); 
q = zeros(numSteps, 1);
T1d = zeros(T / dt, 4); T2d = zeros(T / dt, 4);
S1d = zeros(T / dt, 4); S2d = zeros(T / dt, 4);

% Advection delays and feedback parameters
% You can change the parameters in this part
delayL = 5; % Advective delay from low latitude to OSNAP section (year)
delayH = 20; % Advective delay from high latitude to OSNAP section (year)
delayC = 5; % Coupled time delay F'=c*q'(t-tauc) (year)
c = 0.005; % Coupled feedback strength
a = 1.4 * 10^(-7); % Damping coefficient ((m^3/kg)^2*s-1)

tauL = round(delayL / dt); 
tauH = round(delayH / dt);
tauC = round(delayC / dt);

% Struct for parameters
p = struct('T1_a', T1_a, 'T2_a', T2_a, 'Hse', Hse, 'V1', V1, 'V2', V2, 'A1', A1, 'A2', A2,...
           'betat', betat, 'betatH', betatH, 'betatL', betatL, 'betas', betas, 'k', k0, 'gamat', gamat,...
           'S0', S0, 'S1e', S1e, 'T1e', T1e, 'S2e', S2e, 'T2e', T2e, 'rho0', rho0, 'rho2e', rho2e, 'rho1e', rho1e,...
           'qe', qe, 'FWe', FWe, 'tauL', tauL, 'tauH', tauH, 'tauC', tauC, 'c', c, 'a', a);

% Determine initial condition starting index based on delays and feedback
init = max(tauL, tauH) + tauC + 1;

% Load initial conditions from file or set initial conditions
%{ 
% Load initial conditions from file
load('ocn_initial.mat'); 
T1(1:init) = ocn_ctrl.T1(end-init+1:end); 
T2(1:init) = ocn_ctrl.T2(end-init+1:end); 
S1(1:init) = ocn_ctrl.S1(end-init+1:end); 
S2(1:init) = ocn_ctrl.S2(end-init+1:end);
q(1:init) = ocn_ctrl.R(end-init+1:end);
%}
% Alternative simpler initial condition for testing or first run. You may
% give a perturbation around the mean state as the initial condition. For example:
T1(1:init) = T1e; T2(1:init) = T2e; S1(1:init) = S1e; S2(1:init) = 1.001 * S2e; q(1:init) = qe;

T1s(1:init) = T1(1:init); T2s(1:init) = T2(1:init); S1s(1:init) = S1(1:init); S2s(1:init) = S2(1:init);
T1ss(1:init) = T1(1:init); T2ss(1:init) = T2(1:init); S1ss(1:init) = S1(1:init); S2ss(1:init) = S2(1:init);
T1sss(1:init) = T1(1:init); T2sss(1:init) = T2(1:init); S1sss(1:init) = S1(1:init); S2sss(1:init) = S2(1:init);
% Initial density calculations
rho1(1:init) = rho1e * (1 + betas * (S1(1:init) - S1e) - betatL * (T1(1:init) - T1e));
rho2(1:init) = rho2e * (1 + betas * (S2(1:init) - S2e) - betatH * (T2(1:init) - T2e));
FW_new(1:init) = NaN;
Hs_new(1:init) = NaN;

% Main simulation loop using RK-4 scheme for numerical integration
for n = init:T / dt
    % Ensure physical bounds on temperatures
    T2(n) = max(min(T2(n), T1(n)), -1.8);
    % RK-4 Integration
    [T1d1(n), T2d1(n), S1d1(n), S2d1(n), Hs_new(n + 1)] = TB_delay(T1, T2, S1, S2, q, n, p, dt);
    T1s(n) = T1(n) + 0.5 * dt * T1d1(n);
    T2s(n) = T2(n) + 0.5 * dt * T2d1(n);
    S1s(n) = S1(n) + 0.5 * dt * S1d1(n);
    S2s(n) = S2(n) + 0.5 * dt * S2d1(n);
    
    [T1d2(n), T2d2(n), S1d2(n), S2d2(n), Hs_new(n + 1)] = TB_delay(T1s, T2s, S1s, S2s, q, n, p, dt);
    T1ss(n) = T1(n) + 0.5 * dt * T1d2(n);
    T2ss(n) = T2(n) + 0.5 * dt * T2d2(n);
    S1ss(n) = S1(n) + 0.5 * dt * S1d2(n);
    S2ss(n) = S2(n) + 0.5 * dt * S2d2(n);
    
    [T1d3(n), T2d3(n), S1d3(n), S2d3(n), Hs_new(n + 1)] = TB_delay(T1ss, T2ss, S1ss, S2ss, q, n, p, dt);
    T1sss(n) = T1(n) + dt * T1d3(n);
    T2sss(n) = T2(n) + dt * T2d3(n);
    S1sss(n) = S1(n) + dt * S1d3(n);
    S2sss(n) = S2(n) + dt * S2d3(n);
    
    [T1d4(n), T2d4(n), S1d4(n), S2d4(n), Hs_new(n + 1)] = TB_delay(T1sss, T2sss, S1sss, S2sss, q, n, p, dt);
    % Update state variables with RK-4 formula
    T1(n + 1) = T1(n) + (1 / 6) * dt * (T1d1(n) + 2 * T1d2(n) + 2 * T1d3(n) + T1d4(n));
    T2(n + 1) = T2(n) + (1 / 6) * dt * (T2d1(n) + 2 * T2d2(n) + 2 * T2d3(n) + T2d4(n));
    S1(n + 1) = S1(n) + (1 / 6) * dt * (S1d1(n) + 2 * S1d2(n) + 2 * S1d3(n) + S1d4(n));
    S2(n + 1) = S2(n) + (1 / 6) * dt * (S2d1(n) + 2 * S2d2(n) + 2 * S2d3(n) + S2d4(n));
     
    rho1(n+1) = rho1e * (1 + betas * (S1(n+1) - S1e) - betatL * (T1(n+1) - T1e)); 
    rho2(n+1) = rho2e * (1 + betas * (S2(n+1) - S2e) - betatH * (T2(n+1) - T2e)); 
    q(n+1) = k0 * (rho2(n+1 - tauH) - rho1(n+1 - tauL));
    
    T_tot(n + 1) = (V1 * T1(n + 1) + V2 * T2(n + 1)) / (V1 + V2);
    S_tot(n + 1) = (V1 * S1(n + 1) + V2 * S2(n + 1)) / (V1 + V2);
    OST(n + 1) = abs(q(n + 1)) * (S1(n + 1 - tauL) - S2(n + 1 - tauH));
    OHT(n + 1) = cp * abs(q(n + 1)) * (T1(n + 1 - tauL) - T2(n + 1 - tauH));
    rho2a(n + 1) = rho2e * (betas * (S2(n + 1) - S2e) - betatH * (T2(n + 1) - T2e));
    dmp(n + 1) = -a * rho2a(n + 1)^2 * (S2(n + 1) - S2e);
end

% Example plot of AMOC strength over time
figure(1);
plot(0:dt:T, q / 10^6, 'LineWidth', 1);
grid on;
ylabel('AMOC Strength (Sv)');
xlabel('Time (years)');
title('AMOC Oscillation Over Time');

% Save the results. The results can be used as the initial condition at your next run.
% ocn_ctrl.q=q; ocn_ctrl.T1=T1; ocn_ctrl.T2=T2; ocn_ctrl.S2=S2; ocn_ctrl.S1=S1;
% save('ocn_initial.mat','ocn_ctrl');
toc

% Supporting functions (TB_delay)
function [T1d, T2d, S1d, S2d, Hs_new] = TB_delay(T1, T2, S1, S2, q, n, p, dt)
    % TB_delay computes the derivatives for T1, T2, S1, S2, and Hs_new based on current states and parameters
    % Inputs:
    % - T1, T2, S1, S2: Arrays of temperature and salinity for low and high latitude boxes
    % - q: Array of AMOC strength
    % - n: Current timestep index
    % - p: Parameters struct containing physical and model parameters
    % - dt: Time step size
    % Outputs:
    % - T1d, T2d: Derivatives of temperature in low and high latitude boxes
    % - S1d, S2d: Derivatives of salinity in low and high latitude boxes
    % - Hs_new: New value for the salt transport

    % Select appropriate indices for delays based on RK-4 step
    tau_hi = n - p.tauH;
    tau_lo = n - p.tauL;
    rho1a = p.rho1e * (p.betas * (S1(n) - p.S1e)  -p.betatL * (T1(n) - p.T1e)); %density anomaly
    rho2a = p.rho2e * (p.betas * (S2(n) - p.S2e) - p.betatH * (T2(n) - p.T2e));
    
    % Calculate new salt transport
    Hs_new = p.S0 * (p.FWe - p.c * p.qe + p.c * abs(q(n - p.tauC)));
    
    % Calculate temperature and salinity derivatives
    T1d = (p.gamat / (p.V1 / p.A1) * (p.T1_a - T1(n)) + abs(q(n)) / p.V1 * (T2(tau_hi) - T1(tau_lo)) - p.a * rho1a^2 * (T1(n) - p.T1e)) * 365 * 86400;
    T2d = (p.gamat / (p.V2 / p.A2) * (p.T2_a - T2(n)) + abs(q(n)) / p.V2 * (T1(tau_lo) - T2(tau_hi)) - p.a * rho2a^2 * (T2(n) - p.T2e)) * 365 * 86400;
    S1d = (Hs_new / p.V1 + abs(q(n)) / p.V1 * (S2(tau_hi) - S1(tau_lo)) - p.a * rho1a^2 * (S1(n) - p.S1e)) * 365 * 86400;
    S2d = (-Hs_new / p.V2 + abs(q(n)) / p.V2 * (S1(tau_lo) - S2(tau_hi)) - p.a * rho2a^2 * (S2(n) - p.S2e)) * 365 * 86400;    
    end

