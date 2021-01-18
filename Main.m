%% Main script 
% for running simulations in "Real-time Minimization of Mechanical Specific 
% Energy with Multivariable Extremum Seeking"
% This script requires acces to: 
% formation1.m, formation2.m, sat.m, squareWave.m and smoothstep.m

% Magnus Nystad, 18 Jan 2021, magnus.nystad@ntnu.no

tic
%% Detournay model parameters
r = 12.25*25.4/2;           % r = bit radius [mm], 12.25 inch diameter bit.
A_bit = pi*r^2;             % Bit area, mm^2
g = 9.81;
rho = 0;                    % coring bit constant. 0 value = full bit.
%% Load formation
formation1

%% conversion functions
W_to_w = @(W)(W*g/r/(1-rho));           % takes in W [kg], returns w [N/mm]
w_to_W = @(w)(w*r*(1-rho)/g);           % takes in w [N/mm], returns W [kg]
RPM_to_omega = @(RPM)(RPM*2*pi/60);     % takes in RPM, returns omega [rad/s]
omega_to_RPM = @(om)(om*60/2/pi);       % takes in omega [rad/s], returns RPM
d_to_V = @(d,RPM)(d*RPM/60);            % takes in d [mm/rev] % RPM, returns V [mm/s]
V_to_ROP = @(V)(V*3600/1000);           % takes in V [mm/s], returns ROP (m/hr]
t_to_T = @(t)(t*r^2*(1-rho^2)/2/1000);  % takes in t [N/mm], returns T [Nm]

%% get parameters based on the Detournay drilling model & penalty functions
getTorque = @(W)(t_to_T(detournay_normalized_t(W_to_w(W))));                    % function that takes in WOB [kg] and converts to Torque [Nm]
getROP = @(W,RPM)(V_to_ROP(d_to_V(detournay_normalized_d(W_to_w(W)),RPM)));     % function that takes in WOB [kg] and RPM and calculates the ROP [m/hr] 
getMSE1 = @(W)(1/A_bit*W*g);                                                    % first term in MSE formulation, [MPa]
getMSE2 = @(W,RPM)(1/A_bit*getTorque(W)*RPM_to_omega(RPM)*3600/getROP(W,RPM));  % second term in MSE formulation, [MPa]
getMSE = @(W,RPM)(getMSE1(W) + getMSE2(W,RPM));                                 % total MSE, dysfunction free, [MPa]
getROP_penalty = @(W,RPM,W_n,RPM_n)(V_to_ROP(d_to_V(detournay_penalty_normalized_d(W_to_w(W),W_n,RPM_n),RPM)));     % function that takes in WOB [kg] and RPM and calculates the ROP [m/hr]
getMSE2_penalty = @(W,RPM,W_n,RPM_n)(1/A_bit*getTorque(W)*RPM_to_omega(RPM)*3600/getROP_penalty(W,RPM,W_n,RPM_n));  % Calculate second MSE term with dysfunctions, [MPa]
getMSE_penalty = @(W,RPM,W_n,RPM_n)(getMSE1(W) + getMSE2_penalty(W,RPM,W_n,RPM_n));                                 % Calculate full MSE value with dysfunctions, [MPa]

%% Simulation parameters
time_max = 5000;            % [s]
samplingFreq = 1;           % [1/s]
dt = 1/samplingFreq;        % [s]
time = 0:dt:time_max-dt;    % time vector for plots

% Dither signal(s)
period = 120;               % period of WOB dither signal, [s] (120-240)
A_WOB = 200;                % amplitude of WOB dither signal, [kg] (200)
A_RPM = 2;                  % amplitude of RPM dither signal, [rpm] (2)
t_start_WOB = 60;           % time before WOB dither signal is added, [s]
t_start_RPM = t_start_WOB;  % Time before RPM dither signal is added, [s], possibility of 90 deg phase shift from start of WOB dither

% Least-squares gradient estimation
t0 = t_start_WOB + period;  % Initialization period before adaptation, [s]
K_esc_WOB = 2.5;            % ESC gain for WOB adaptation, suggested: 1-3
K_esc_RPM = 0.02;           % ESC gain for RPM adaptation, suggested: 0.01-0.02
k_WOB = 0.001;              % 0.001 for clean data
k_RPM = 0.05;               % 0.05 for clean data
t_LSE = period;             % sliding window "size" for LS gradient estimation, [s]
g1 = 0;                     % Initial gradient 1, dMSE/dWOB
g2 = 0;                     % initial gradient 2, dROP/dWOB
g3 = 0;                     % Initial gradient 3, dMSE/dRPM
g3_avg = 0;                 % initial avg gradient 3
g4 = 0;                     % initial gradient 4, dROP/dRPM
g5 = 0;                     % initial gradient 5, dT/dWOB
adapt_WOB = 0;              % initial WOB adaptation
adapt_RPM = 0;              % initial RPM adaptation

% Constraint handling
ROP_limit = 200;            % [m/hr] high value --> no constraint handling
T_limit = 100000;           % [Nm] high value --> no constraint handling
SF1 = 2;                    % Safety factor in predictive constraint handling
K_p = 0.05;                 % Proporional gain, reactive constraint handling
K_i = 0.001;                % Integral gain, reactive constraint handling
delta_WOB_penalty = 0;
WOB_penalty_last = 0;
e_sum = 0;
rho_ROP = 10;               % Rho parameter in modified objective function

% Preallocate vectors
WOB_buffer = zeros(1,t_LSE);
ROP_buffer = zeros(1,t_LSE);
T_buffer = zeros(1,t_LSE);
RPM_buffer = zeros(1,t_LSE);
MSE_buffer = zeros(1,t_LSE);
dMSEdRPM_buffer = zeros(1,t_LSE);
MSE_per_WOB_buffer = zeros(1,t_LSE);
WOB_store = zeros(1,time_max);
WOB_0_store = zeros(1,time_max);
WOB_SP_store = zeros(1,time_max);
ROP_store = zeros(1,time_max);
T_store = zeros(1,time_max);
RPM_store = zeros(1,time_max);
RPM_0_store = zeros(1,time_max);
RPM_SP_store = zeros(1,time_max);
MSE_store = zeros(1,time_max);
MSE_mavg_store = zeros(1,time_max);
dd_store = zeros(1,time_max);
adapt_WOB_store = zeros(1,time_max);
adapt_RPM_store = zeros(1,time_max);
dMSEdWOB_store = zeros(1,time_max);
dROPdWOB_store = zeros(1,time_max);
dMSEdRPM_store = zeros(1,time_max);
dMSEdRPM_mavg_store = zeros(1,time_max);
dROPdRPM_store = zeros(1,time_max);

% Initial parameters
WOB_SP = 13000;             % initial WOB setpoint, [kg]
WOB = 0;                    % actual WOB at t=0
WOB_0 = WOB_SP;             % initial "base" WOB, [kg]
RPM_SP = 90;                % initial RPM setpoint, [rpm]
RPM = 0;                    % actual RPM at t=0
RPM_0 = RPM_SP;             % initial "base" RPM
dd = 0;
torque_jump = 0;            % initiate torque jump val = 0
dd_fm_switch1 = 5;          % drilled depth [m] for formation shift 1
dd_fm_switch2 = 23;         % drilled depth [m] for formation shift 2
torque_jump_val = 0;        % ad hoc implementation of sudden torque increase, 1000-2000 Nm
torque_jump_time = 2000;    % time of torque increase

%% Simulation
for i = 1:time_max
    
    % formation shift code
    if dd > dd_fm_switch1 && dd_store(i-2) <= dd_fm_switch1
        formation2;
        disp('fm2')
        shift1 = i;     % time index @ formation shift 1
    elseif dd > dd_fm_switch2 && dd_store(i-2) <= dd_fm_switch2
        formation1;
        disp('fm1')
        shift2 = i;     % time index @ formation shift 2
    end
    
        
    % Ad hoc torque jump implementation
    if i == torque_jump_time
        torque_jump = torque_jump_val;
    end
    
    % get current parameters
    WOB = WOB + 1/4*(WOB_SP - WOB);                             % current WOB [kg], smooth dynamics with 4-5 sec rise time
    RPM = RPM + 1/3*(RPM_SP - RPM);                             % current RPM, smooth dynamics with 3 sec rise time
    %RPM = min(abs(RPM_SP-RPM),1)*sign(RPM_SP - RPM) + RPM;     % limear RPM dynamics, not used                      
    
    WOB_n = WOB/WOB_norm;                                       % Normalized WOB for Detournay calculation
    RPM_n = RPM/RPM_norm;                                       % Normalized RPM for Detournay calculation
    ROP = getROP_penalty(WOB,RPM,WOB_n,RPM_n);                  % Instantaneous ROP, [m/hr]
    T = getTorque(WOB) + torque_jump;                           % Current torque on bit, [Nm] ad hoc torque jump
    addMSE_jump = torque_jump*RPM_to_omega(RPM)*3600/ROP/A_bit; % Increase in MSE from simulated torque jump, ad hoc implementation
    MSE = getMSE_penalty(WOB,RPM,WOB_n,RPM_n) + addMSE_jump;    % Current MSE, [MPa]
    dd = dd + ROP/3600*dt;                                      % drilled depth, m
    
    %modified obj func for constraint handling
    ROP_penalty = (max(ROP,ROP_limit) - ROP_limit)/(ROP_limit*rho_ROP);
    MSE = MSE*(1+ROP_penalty);                                  % Modified objective function
    
    % Update buffers
    WOB_buffer = [WOB_buffer(2:end),WOB];
    ROP_buffer = [ROP_buffer(2:end),ROP];
    T_buffer = [T_buffer(2:end),T];
    RPM_buffer = [RPM_buffer(2:end),RPM];
    MSE_buffer = [MSE_buffer(2:end),MSE];
    T_avg = mean(T_buffer);
    MSE_avg = mean(MSE_buffer);
    
    if i > t_start_WOB                  % start WOB dither signal
        
        if i > t0                       % start adaptation, after 1 period of dither signal
            % Estimate gradient(s)
            if A_WOB > 0                % Test if WOB-dither is active, avoid polyfit issues
                p = polyfit(WOB_buffer,MSE_buffer,1);
                g1 = p(1);                                          % dMSE/dWOB [MPa/kg]
                p = polyfit(WOB_buffer,ROP_buffer,1);
                g2 = p(1);                                          % dROP/dWOB [m/hr/kg]
                adapt_WOB = sat(g1-k_WOB,2*k_WOB)*K_esc_WOB*dt;
                % Predictive torque constraint
                p = polyfit(WOB_buffer,T_buffer,1);                 %dT/dWOB
                g5 = p(1);
            end
             
            if A_RPM > 0        % after 1 adaptation-free period of RPM signal, and RPM dither activated
                p = polyfit(RPM_buffer,MSE_buffer,1);
                g3 = p(1);                                      % dMSE/dRPM [MPa/rpm]
                p = polyfit(RPM_buffer,ROP_buffer,1);
                g4 = p(1);                                      % dROP/dRPM [m/hr/rpm]
                
                dMSEdRPM_buffer = [dMSEdRPM_buffer(2:end),g3];
                g3_avg = mean(dMSEdRPM_buffer);                 % average gradient
                adapt_RPM = sat(g3-k_RPM,2*k_RPM)*K_esc_RPM*dt; % RPM adaptation
            end
            
        end
        
        % Predictive torque constraint
        if (T_avg + A_WOB*g5*SF1) > T_limit
            adapt_WOB = 0;
        end
        
        % reactive torque constraint
        if T > T_limit
            adapt_WOB = 0;
            adapt_RPM = 0;
            e_current = T-T_limit;
            e_sum = e_sum + e_current;
            WOB_penalty_current = K_p*e_current + K_i*e_sum;
            delta_WOB_penalty = WOB_penalty_current - WOB_penalty_last;
            WOB_penalty_last = WOB_penalty_current;
        else
            WOB_penalty_current = 0;
        end
       
        % WOB dither signal and adaptation
        WOB_0 = WOB_0 - adapt_WOB - WOB_penalty_current;                        % update WOB_0 value with ESC adaptation
        WOB_SP = WOB_0 + A_WOB*squareWave(i,t_start_WOB,period);                % update WOB_SP [kg]
        RPM_0 = RPM_0 - adapt_RPM;                                              % update RPM_0 value with ESC adaptation
        RPM_SP = RPM_0 + A_RPM*squareWave(i,t_start_RPM,period/2);              % update RPM_SP with oscillating square wave
    end
        
    % Parameter storage
    WOB_store(i) = WOB;
    WOB_0_store(i) = WOB_0;
    WOB_SP_store(i) = WOB_SP;
    ROP_store(i) = ROP;
    RPM_store(i) = RPM;
    RPM_0_store(i) = RPM_0;
    RPM_SP_store(i) = RPM_SP;
    T_store(i) = T;
    MSE_store(i) = MSE;
    MSE_mavg_store(i) = MSE_avg;
    dd_store(i) = dd;
    adapt_WOB_store(i) = adapt_WOB;
    adapt_RPM_store(i) = adapt_RPM;
    dMSEdWOB_store(i) = g1;
    dROPdWOB_store(i) = g2;
    dMSEdRPM_store(i) = g3;
    dMSEdRPM_mavg_store(i) = g3_avg;
    dROPdRPM_store(i) = g4;
end

toc
%% Detournay functions
function d = detournay_normalized_d(w)
% function that takes in normalized weight, w [N/mm], and calculates the
% corresponding depth of cut per rev, d [mm/rev]
global d_I d_II w_star

    if w < w_star
        d = d_I(w);
    else 
        d = d_II(w);
    end
end

function t = detournay_normalized_t(w)
% Function that takes in normalized weight, w [N/mm], and calculates the
% corresponding normalized torque, t [Nmm]
global t_I t_II w_star 

    if w < w_star
        t = t_I(w);
    else
        t = t_II(w);
    end
end

function d = detournay_penalty_normalized_d(w,WOB_n,RPM_n)
% Function that takes in normalized weight, w [N/mm], and calculates the
% (possibly) penalized depth of cut, d [mm]
global d_I d_II w_star bwhirl fwhirl stickslip c_vec a1 a2 a3 b1 b2 b3
    % Calculate "ideal" depth of cut at current w
    if w < w_star
        d = d_I(w);
    else
        d = d_II(w);
    end
    
    % Calculate penalty
    L1 = 0;
    L2 = 0;
    L3 = 0;
    
    if WOB_n > bwhirl(RPM_n)
        x = (WOB_n+RPM_n/a1-b1)/(a1+1/a1);          % find closest point on dysfunction curve
        y = bwhirl(RPM_n);                          % find corresponding WOB_norm value
        L1 = sqrt((WOB_n - y)^2 + (RPM_n - x)^2);   % normalized distance from point to dysfunction curve
    end
    if WOB_n < fwhirl(RPM_n)
        x = (WOB_n+RPM_n/a2-b2)/(a2+1/a2);          % find closest point on dysfunction curve
        y = fwhirl(RPM_n);                          % find corresponding WOB_norm value
        L2 = sqrt((WOB_n - y)^2 + (RPM_n - x)^2);   % normalized distance from point to dysfunction curve
    end
    if WOB_n > stickslip(RPM_n)
        x = (WOB_n+RPM_n/a3-b3)/(a3+1/a3);          % find closest point on dysfunction curve
        y = stickslip(RPM_n);                       % find corresponding WOB_norm value
        L3 = sqrt((WOB_n - y)^2 + (RPM_n - x)^2);   % normalized distance from point to dysfunction curve 
    end
    
    L = [L1,L2,L3];
    d = d*(1-smoothstep(L*c_vec',1));
  
end

