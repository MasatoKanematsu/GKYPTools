%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CP_type --- control performance type, 
%                        for example, minimize  gamma_1 or something...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization 
close all; clear all;
LMI = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Plant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ap = [   0,   1,   0;
         0,   0,   1;
       -10, -11,  -2]; 

%sys_L.analysis_type = 'bspline';
sys_L.analysis_type = 'fix';

np = 3;

Bp = [0;
      0;
      1];

Cp = [10, 0, 0];
Dp = 0;
sys_P.A = Ap; sys_P.B = Bp; sys_P.C = Cp; sys_P.D = Dp; sys_P.nx = np;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Controller 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = 0.05;
s = tf([1,0], 1);
kp = 1; ki = 1; kd = 1;
K = kp + ki/s + kd * s/(1+ tau*s);

% This expression is useful to define sdpvar controller.
temp = balreal(K);
[Ak, Bk, Ck, Dk] = ssdata(temp);
Bk = sdpvar(size(Bk, 1), 1, 'full'); Dk = sdpvar( 1, 1, 'full');
nk = size(Ak, 1);

% Warning ! numerical stability depends on the choice of ss expression relatively. If your calculation seems not to be stable, please use another ss expression.

%{
%K = k_p + k_i/s + k_d * s/(1+ tau*s);
Bk = sdpvar(2,1, 'full');
Dk = sdpvar(1,1, 'full');
Ak = [0,       0; 
      1, -1/tau];
nk = size(Ak, 1);

%Bk = [ Bk_11;
%       Bk_21];
%Bk = [ 21.4869;
%      -60.4783];

Ck = [0, 1];
%Dk = [ Dk_11];
%Dk = [ 3.2895 ];
%}

sys_K.A = Ak; sys_K.B = Bk; sys_K.C = Ck; sys_K.D = Dk; sys_K.nx = nk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Open Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Open Loop Characteristics
%[AL,BL,CL,DL] = series(Ak, Bk, Ck, Dk,  Ap, Bp, Cp, Dp)
AL = [   Ak, zeros(nk, np);
      Bp*Ck,           Ap];

BL = [   Bk; Bp*Dk];
CL = [Dp*Ck,    Cp];
DL = [Dp*Dk];

% Define other parameters
nu = 1;% the freedom of the input
ny = 1;% the freedom of the output
nl = nk + np;

sys_L.A = AL; sys_L.B = BL; sys_L.C = CL; sys_L.D = DL;
sys_L.nu = nu; sys_L.ny = ny; sys_L.nx = nl;

sys_L.time_domain = 'continuous'; % or 'discrete'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set LMI Condition
%[LMI_1, LMI_2, LMI_3, obj] = Load_LMI_Data_para(sys_L, Basis, CP_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI 1  --Line Property, Bandwidth--
% Phi for the specification for Continuous and Discrete
% WARNING! : This source code is not passed debug yet.
% w1 < w < w2
LMI_Cond_set1.freq_info.freq_range = 'MF';  % 'LF', 'MF', and 'HF'
LMI_Cond_set1.freq_info.w1 = 0.05; %[Hz]
LMI_Cond_set1.freq_info.w2 = 0.45; %[Hz]

LMI_Cond_set1.FDI_info.conic_type  = 'line';%'line', 'circle', and 'gen_circle'
LMI_Cond_set1.FDI_info.eq =  'ax + by < c'; %
LMI_Cond_set1.FDI_info.a  =  0; %
LMI_Cond_set1.FDI_info.b  =  1; %
LMI_Cond_set1.FDI_info.c  = -2; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI 2  --Line Property, Stability Margin--
% Phi for the specification for Continuous and Discrete
% WARNING! : This source code is not passed debug yet.

LMI_Cond_set2.freq_info.freq_range = 'MF';
LMI_Cond_set2.freq_info.w1 = 0.05;%[Hz]
LMI_Cond_set2.freq_info.w2 = 100; %[Hz]

LMI_Cond_set2.FDI_info.conic_type  = 'line'; %'line', 'circle', and 'gen_circle'
LMI_Cond_set2.FDI_info.eq =  'ax + by < c'; %
LMI_Cond_set2.FDI_info.a  =  -3; %'line', 'circle', and 'gen_circle'
LMI_Cond_set2.FDI_info.b  =   1; %'line', 'circle', and 'gen_circle'
LMI_Cond_set2.FDI_info.c  = sdpvar( 1, 1, 'full'); %'line', 'circle', and 'gen_circle'
obj = LMI_Cond_set2.FDI_info.c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI 3  --Circle Property, Robust Stability--
% Phi for the specification for Continuous and Discrete
LMI_Cond_set3.freq_info.freq_range = 'MF';
LMI_Cond_set3.freq_info.w1 = 7;  %[Hz]
LMI_Cond_set3.freq_info.w2 = 100;%[Hz]

%LMI_3.Pi.gamma = gamma_3; % This may cause some error... not a good manner.
LMI_Cond_set3.FDI_info.conic_type = 'circle';%'line', 'circle', and 'gen_circle'
r = 0.1;  r_sq = (0.1)^2;  
LMI_Cond_set3.FDI_info.r_sq   = r_sq; 
LMI_Cond_set3.FDI_info.r      = r; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective Function  --Set Objective Function--

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Objective Function 

LMI_Cond = {LMI_Cond_set1, LMI_Cond_set2, LMI_Cond_set3};
analysis_type = 'fix';
SV_type = 'BD';
knots_info = 0;

for i=1:1:size(LMI_Cond, 2)
%	[LMI_set] = Set_LMI_for_parametric_GKYP_synthesis(sys_L, LMI_Cond{i}, SV_type, analysis_type, knots_info);
	[LMI_set] = Set_LMI_for_GKYP_synthesis(sys_L, LMI_Cond{i}, SV_type, analysis_type, knots_info);
	LMI = [LMI, LMI_set];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Optimization Problem
%sol = optimize(LMI, obj.integral);
sol = optimize(LMI, obj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Post Prosessing                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI_Condのdouble変換
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI_Condの中に、sdpvarがあると後半のFDIの描写が出来ないため、
% 構造体の中のsdpvarのみ、doubleに変換する必要がある。
% かなりこの操作は高度で興味がある人は、文六さんのmatlabメモ参照

for i=1:1:size(LMI_Cond, 2) %全てのLMI_Condに対して、

  %構造体entry nameの所得
  fn_list = fieldnames(LMI_Cond{i}.FDI_info);

  %構造体エントリーの中から、sdpvar型を見つけて、doubleに変換する。
  for k=1:1:size(fn_list, 1) %全ての構造体エントリーに対して、
    %もし構造体エントリーのデータ型がsdpvarであるならば、doubleに変換
    if isa( LMI_Cond{i}.FDI_info.(fn_list{k}), 'sdpvar')
      LMI_Cond{i}.FDI_info.(fn_list{k}) = value( LMI_Cond{i}.FDI_info.(fn_list{k}) );
    end
  end

end

%save('LMI_Cond_Data.mat', 'LMI_Cond');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objのdouble変換
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All parameters are calculated here.
obj = value(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_Kのdouble変換
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set designed open loop characteristics(ss form)
sys_K.B = value(sys_K.B); sys_K.D = value(sys_K.D);
sys_K.ss = ss(sys_K.A, sys_K.B, sys_K.C, sys_K.D);%, sys_K.Ts);

sys_K0 = sys_K;
save('fixed_PID_controller.mat', 'sys_K0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_Lのdouble変換
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set designed open loop characteristics(ss form)
sys_L.B = value(sys_L.B); sys_L.D = value(sys_L.D);
sys_L.ss = ss(sys_L.A, sys_L.B, sys_L.C, sys_L.D);%, sys_L.Ts);

% set designed open loop characteristics(tf form)
%[num, den] = ss2tf(sys_L.A, sys_L.B, sys_L.C, sys_L.D);
%sys_L.tf = tf(num, den, sys_L.Ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Show Graph                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ここまでで既に伝達関数が定義できておりますので、各自好きなように
% 図示して、開ループ整形の結果をご確認ください。
% 自分は下記のプログラムを利用しております。
Ts = 20e-5;
%                          ↓システム構造体     x_axis range, y_axis range
%Draw_Nyquist_OpenLoop_ver6(sys_L, LMI_Cond, [-3, 1], [-2, 2] ,Ts);
system_cell = {sys_L.ss};
save_filename = 'nyquist_L.eps';
FDI_ON = 1;
Draw_Nyquist_OpenLoop_ver8(system_cell, LMI_Cond, [-3, 1], [-3, 1] ,save_filename, FDI_ON);




%Draw_Nyquist_OpenLoop_ver7(sys_L.ss, LMI_Cond, [-3, 1], [-3, 1] ,save_filename, FDI_ON);

%Draw_Nyquist_Plant_ver2(sys_Pd.ss, Ts);

%Plot_Plant_Bode
%Draw_Bode(sys_L.ss, Ts, [1, 10000], [-40, 20], 'bode_prop_L');



% グラフの表示
%Ts = 1/(20*10^3);
%Draw_Bode_ver7({sys_L.ss}, (Ts), [0.1, 10], [-50, 30], 'bode_L');







%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             パラメータ変動時のLMI_Condを計算                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = -0.3;
Ap = [               0,           1,           0;
                     0,           0,           1;
                 (-10-10*sigma),  (-11-sigma),  (-2-sigma)]; %２乗なども簡単に計算可能。

%sys_L.analysis_type = 'bspline';
sys_L.analysis_type = 'fix';
np = 3;
Bp = [0;
      0;
      1];
Cp = [10, 0, 0];
Dp = 0;

sys_P.ss = ss(Ap, Bp, Cp, Dp);
sys_L.ss = sys_P.ss * sys_K.ss;
[sys_L.A, sys_L.B, sys_L.C, sys_L.D] = ssdata(sys_L.ss);
sys_L.nx = size(sys_L.A, 1); 
sys_L.nu = 1; sys_L.ny = 1;

LMI_Cond_set2.freq_info.freq_range = 'MF';
LMI_Cond_set2.freq_info.w1 = 0.05;%[Hz]
LMI_Cond_set2.freq_info.w2 = 100; %[Hz]

LMI_Cond_set2.FDI_info.conic_type  = 'line'; %'line', 'circle', and 'gen_circle'
LMI_Cond_set2.FDI_info.eq     =  'ax + by < c'; %

LMI_Cond_set2.FDI_info.a  =  -3; %'line', 'circle', and 'gen_circle'
LMI_Cond_set2.FDI_info.b  =   1; %'line', 'circle', and 'gen_circle'
LMI_Cond_set2.FDI_info.c  = sdpvar( 1, 1, 'full'); %'line', 'circle', and 'gen_circle'
obj = LMI_Cond_set2.FDI_info.c;


% for LMI 
LMI = [];

knots_info = 0;
analysis_type = 'fix';
[LMI_set] = Set_LMI_for_parametric_GKYP_synthesis(sys_L, LMI_Cond_set2, knots_info, analysis_type);
LMI = [LMI, LMI_set];

sol = optimize(LMI, obj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI_Condのdouble変換
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI_Condの中に、sdpvarがあると後半のFDIの描写が出来ないため、
% 構造体の中のsdpvarのみ、doubleに変換する必要がある。
% かなりこの操作は高度で興味がある人は、文六さんのmatlabメモ参照

LMI_Cond = {LMI_Cond_set1, LMI_Cond_set2, LMI_Cond_set3};

for i=1:1:size(LMI_Cond, 2) %全てのLMI_Condに対して、

  %構造体entry nameの所得
  fn_list = fieldnames(LMI_Cond{i}.FDI_info);

  %構造体エントリーの中から、sdpvar型を見つけて、doubleに変換する。
  for k=1:1:size(fn_list, 1) %全ての構造体エントリーに対して、
    %もし構造体エントリーのデータ型がsdpvarであるならば、doubleに変換
    if isa( LMI_Cond{i}.FDI_info.(fn_list{k}), 'sdpvar')
      LMI_Cond{i}.FDI_info.(fn_list{k}) = value( LMI_Cond{i}.FDI_info.(fn_list{k}) );
    end
  end

end

LMI_Cond0 = LMI_Cond;
save('LMI_Cond_Data_PID_Plant.mat', 'LMI_Cond0');
obj = value(obj)


%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        END OF THIS PROGRAM                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
