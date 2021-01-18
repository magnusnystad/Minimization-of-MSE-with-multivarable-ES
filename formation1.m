%% Detournay rock parameters

% Formation 1 - hard
UCS = 60;               % [MPa]
eta = 0.35;             % dimensionless, peak efficiency
epsilon = UCS/eta;      % [MPa], intrinsic specific energy
zeta = 0.7;             % bit constant in the range [0.5,0.8]. zeta = tan(theta + psi)
kappa = 5;              % number in the range [1,10]. Represents the rate of change of the contact length with depth of cut, d
sigma = UCS;            % [MPa], contact strength
my = 0.3;               % Coefficient of friction
gamma = 1;              % Bit constant that "encapsulates the influence of the orientation & distribution of contact forces on the bit.
l = 2;                  % [mm], length of cutter wear flat. "Objective measure of the bit bluntness"
                                            
% Detournay et al (2008) - drilling response of drag bits: theory and experiment                                        
E_0 = (1-my*gamma*zeta)*epsilon;            % eq. (12)
S_star = zeta*epsilon + kappa*sigma;        % eq. (25)
E_star = epsilon + my*gamma*sigma*kappa;    % eq. (28)
wf_star = sigma*l;                          % eq. (33) 

%% Detournay ROP model functions
global d_I d_II t_I t_II w_star
% Phase I drilling
d_I = @(w)(w/S_star);
t_I = @(w)(w*E_star/S_star);

% Phase II drilling
w_star = wf_star*(E_star-E_0)/(E_star - epsilon);
t_star = w_star*E_star/S_star;
d_star = w_star/S_star;

d_II = @(w)(d_star + (w - w_star)/zeta/epsilon);
t_II = @(w)(t_star + (w - w_star)/zeta);

%% Penalty parameters & functions
WOB_norm = 20000;       % 15000-20000
RPM_norm = 200;         % 150-200

% Normalized dysfunc functions
global bwhirl fwhirl stickslip c_vec a1 a2 a3 b1 b2 b3
a1 = -1.1;                          % backwhirl slope
b1 = 1.135;                          % backwhirl intercept
c1 = 1;%1.3;                           % penalty coefficient
bwhirl = @(RPM_n)(b1 +a1*RPM_n);

a2 = 0.5;                           % forward whirl slope, 0.8
b2 = -0.1;                          % forward whirl intercept, -0.2
c2 = 1;                           % penalty coefficient
fwhirl = @(RPM_n)(b2 + a2*RPM_n);

a3 = 1.2;                           % stick slip slope
b3 = 0.1;                          % stick slip intercept
c3 = 0.5; %0.8;                             % penalty coefficient in sstep approach
stickslip = @(RPM_n)(b3 + a3*RPM_n);

c_vec = [c1,c2,c3];                 % penalty coefficient vector: [c1, c2, c3]

