function [LMI] = Set_LMI_for_GKYP_synthesis(...
                        sys, ...           % ss or tf system
                        LMI_Cond, ...      % please see below.
			SV_type, ...       % Searching Variable ('BD', or 'CD')
			analysis_type, ... % 'fix', or 'bspline'
			knots_info )       % please see below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Set_LMI_Specific_Conic_Analysis_4" が前身
% For copy paste 
%[LMI] = Set_LMI_for_GKYP_synthesis(sys,LMI_Cond,SV_type, analysis_type, knots_info ) 
%LMI = [LMI, LMI_set];
%
%SV_type = 'BD';
%analysis_type = 'fix';
%knots_info = 0;
%
%
% All variables are defined as complex.
% In the last line, complex LMI is converted to real LMI.
% ( This code can't solve B-Spline GKYP synthesis problem. )
%
% 複素変換を一番最後に一括で行うバージョン
% updated at 2016/08/02
% 
% sys構造体の要素
% A = sys.A; B = sys.B; C = sys.C; D = sys.D;
% nx = sys.nx;   nu = sys.nu;   ny = sys.ny;
% time_domain = sys.time_domain; % 'continuous' or 'discrete'
%
% LMI_Cond 
%   |
%   |---FDI_info  => Piの導出に用いる
%   |---freq_info => Phi, Psiの導出に用いる
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
LMI = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data

time_domain = sys.time_domain; % 'continuous' or 'discrete'
freq_range  = LMI_Cond.freq_info.freq_range;  % 'LF', 'MF', and 'HF'
conic_type  = LMI_Cond.FDI_info.conic_type; %'line', 'circle', and 'gen_circle'

A = sys.A; B = sys.B; C = sys.C; D = sys.D;
nx = sys.nx;   nu = sys.nu;   ny = sys.ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Some Matrix

In =     eye(nx, nx);
Jn = j * eye(nx, nx);
On =   zeros(nx, nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define SDP-Variable 
%WARNING : Actually these values are doubled 
%          for the transformation from complex SDP to real SDP

P = sdpvar(nx, nx, 'hermitian', 'complex');
Q = sdpvar(nx, nx, 'hermitian', 'complex');

LMI_set = [  real(Q), imag(Q);
            -imag(Q), real(Q)];

LMI = [LMI, LMI_set >= 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Set Phi, Psi, Pi                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are defined as nested function below.

Phi_Mat = Define_Phi_Mat(time_domain);
Psi_Mat = Define_Psi_Mat(time_domain, freq_range);
Pi_Mat  = Define_Pi_Mat(  conic_type,    SV_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define main LMI Condition 
%   最終的なLMIのセット
%
%       For the reduction of the computational cost, 
%       it is desirable to reduce matrix size here.
%       BSplineを使う場合は、
%       計算時間削減のためにも、LMIのサイズを下げれる場合は積極的に下げる。
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(SV_type, 'BD')
	% 1. Define W Matrix
	% 1. Wの定義
	set1 = [A,   eye(nx, nx); C, zeros(ny, nx)];
	
	PhiT_X_P = [ Phi_Mat{1,1} * P, Phi_Mat{2,1} * P;
	             Phi_Mat{1,2} * P, Phi_Mat{2,2} * P];
	
	PsiT_X_Q = [ Psi_Mat{1,1} * Q, Psi_Mat{2,1} * Q;
 	             Psi_Mat{1,2} * Q, Psi_Mat{2,2} * Q];

	W = set1 * (PhiT_X_P + PsiT_X_Q) * set1';

	% 2. Define V Matrix
	% 2. Vの定義
	V11 = zeros(nx, nx);
	V12 =  B * Pi_Mat{1,2}; %ここが要注意。ABIOとAICOで違う。
	V21 = (B * Pi_Mat{1,2})';
	V22 = D*Pi_Mat{1,2} + Pi_Mat{1,2}'*D' + Pi_Mat{2,2};

	V = [V11, V12
	     V21, V22];


	% 3. Define the LMI from W and V matrix.
	% 3. 全体の定義
	Mat_Data = [                   W + V,    [B; D] * Pi_Mat{1,1};
	                Pi_Mat{1,1} * [B; D]',           -Pi_Mat{1,1}];
%			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(SV_type, 'CD')
	% 1. Wの定義
	set1 = [A,  B;  eye(nx, nx), zeros(nx, ny)];
	
	Phi_X_P = [ Phi_Mat{1,1} * P, Phi_Mat{1,2} * P;
	            Phi_Mat{2,1} * P, Phi_Mat{2,2} * P];
	
	Psi_X_Q = [ Psi_Mat{1,1} * Q, Psi_Mat{1,2} * Q;
	            Psi_Mat{2,1} * Q, Psi_Mat{2,2} * Q];

	W = set1' * (Phi_X_P + Psi_X_Q) * set1;


	% 2. Vの定義
	V11 = zeros(nx, nx);
	V12 =  C' * Pi_Mat{1,2}; %ここが要注意。ABIOとAICOで違う。
	V21 = (C' * Pi_Mat{1,2})';
	V22 = D'*Pi_Mat{1,2} + Pi_Mat{1,2}'*D + Pi_Mat{2,2};
	V = [V11, V12
	     V21, V22];


	% 3. 全体の定義
	Mat_Data = [                   W + V,  [C, D]' * Pi_Mat{1,1};
	                Pi_Mat{1,1} * [C, D],           -Pi_Mat{1,1}];
%			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

LMI_set = [  real(Mat_Data), imag(Mat_Data);
            -imag(Mat_Data), real(Mat_Data)];
LMI = [LMI, LMI_set <= 0];
obj = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Define Nested Function(入れ子関数)                          %
% Please see 
% http://jp.mathworks.com/help/matlab/matlab_prog/nested-functions.html
% comment:この入れ子関数はスコープが異様に分かりづらいため、
%         基本的に非推奨だが、今回はソースコードの可視性のため用いている。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Define Phi Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi_Mat] = Define_Phi_Mat(time_domain);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if     strcmp(time_domain, 'continuous');
		Phi_Mat{1,1} =  On; Phi_Mat{1,2} =  In;
		Phi_Mat{2,1} =  In; Phi_Mat{2,2} =  On; 
	
	elseif strcmp(time_domain, 'discrete');
		Phi_Mat{1,1} =  In; Phi_Mat{1,2} =  On;
		Phi_Mat{2,1} =  On; Phi_Mat{2,2} = -In; 
	else
		error('MK Error : Please Check the source code in time_domain.');
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Define Psi Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Psi_Mat] = Define_Psi_Mat(time_domain, freq_range);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if strcmp(freq_range, 'LF')
		w = LMI_Cond.freq_info.w;
	
		if strcmp(time_domain, '')
	%		error('MK Error : Please Check the source code in freq_range');
			Psi_Mat{1,1} =     On;
			Psi_Mat{1,2} =    -Jn;
			Psi_Mat{2,1} =     Jn;
			Psi_Mat{2,2} = 2*w*In;
	
		elseif strcmp(time_domain, 'discrete')
			theta = sys.Ts * w;
	
			Psi_Mat{1,1} =                           On;
			Psi_Mat{1,2} =   -j*exp(  j* theta/2 ) * In;
			Psi_Mat{2,1} =    j*exp( -j* theta/2 ) * In;
			Psi_Mat{2,2} =    2*sin(     theta/2 ) * In;
	
		else 
			error('MK Error : Please Check the source code in freq_range');
	
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(freq_range, 'MF')
		w1 = LMI_Cond.freq_info.w1;
		w2 = LMI_Cond.freq_info.w2;
	
		if strcmp(time_domain, 'continuous')
			Psi_Mat{1,1} = -In;
			Psi_Mat{1,2} =  (w1+w2)/2*Jn;
			Psi_Mat{2,1} = -(w1+w2)/2*Jn;
			Psi_Mat{2,2} = -(w1*w2)*In;
	
		elseif strcmp(time_domain, 'discrete')
			theta_1 = sys.Ts * w1;
			theta_2 = sys.Ts * w2;
	
			Psi_Mat{1,1} =  On;
			Psi_Mat{1,2} =   exp(  j*(theta_1 + theta_2)/2 ) * In;
			Psi_Mat{2,1} =   exp( -j*(theta_1 + theta_2)/2 ) * In;
			Psi_Mat{2,2} =  -2*cos(  (theta_2 - theta_1)/2 ) * In;
		else 
			error('MK Error : Please Check the source code in freq_range');
	
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(freq_range, 'HF')
		w = LMI_Cond.freq_info.w;
		
		if strcmp(time_domain, 'continuous')
			Psi_Mat{1,1} =      On;
			Psi_Mat{1,2} =      Jn;
			Psi_Mat{2,1} =     -Jn;
			Psi_Mat{2,2} = -2*w*In;
	
		elseif strcmp(time_domain, 'discrete')
			theta = sys.Ts * w;
	
			Psi_Mat{1,1} =                          On;
			Psi_Mat{1,2} =   j*exp(  j* theta/2 ) * In;
			Psi_Mat{2,1} =  -j*exp( -j* theta/2 ) * In;
			Psi_Mat{2,2} =  -2*sin(     theta/2 ) * In;
	
		else 
			error('MK Error : Please Check the source code in freq_range');
	
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(freq_range, 'HF_1')
		w1 = LMI_Cond.freq_info.w1;
		
		Psi_Mat{1,1} =     0*In;
		Psi_Mat{1,2} =     1*Jn;
		Psi_Mat{2,1} =    -1*Jn;
		Psi_Mat{2,2} = -2*w1*In;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	else
		error('MK Error : Please Check the source code in freq_range.');
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end	

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Define Pi Matrix -FDI characteristics-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pi_Mat] = Define_Pi_Mat(conic_type, SV_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define some matrix
	
	Im = 1;
	Om = 0;
	Jm = j;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if strcmp(conic_type, 'line')
		% line condition
		% FDI_info.eq = 'ax + by < c'
		a = LMI_Cond.FDI_info.a;
		b = LMI_Cond.FDI_info.b;
		c = LMI_Cond.FDI_info.c;
	
		if strcmp(SV_type, 'BD')
			Pi_Mat{1,1} = 0 * Im;
			Pi_Mat{1,2} = a * Im - b * Jm;
			Pi_Mat{2,1} = a * Im + b * Jm;
			Pi_Mat{2,2} = -2 * c * Im;
	
		elseif strcmp(SV_type, 'CD')
			Pi_Mat{1,1} = 0 * Im;
			Pi_Mat{1,2} = a * Im + b * Jm;
			Pi_Mat{2,1} = a * Im - b * Jm;
			Pi_Mat{2,2} = -2 * c * Im;
		else 
			error('MK Error : Please Check the source code in freq_range');
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(conic_type, 'circle')
		% normal CIRCLE condition
		% FDI_info.eq = x^2 + y^2 <= r_sq
		r_sq = LMI_Cond.FDI_info.r_sq;
	
		if strcmp(SV_type, 'BD')
			Pi_Mat{1,1} = 1 * Im;
			Pi_Mat{1,2} = 0 * Om;
			Pi_Mat{2,1} = 0 * Om;
			Pi_Mat{2,2} = -r_sq * Im;
	
		elseif strcmp(SV_type, 'CD')
			Pi_Mat{1,1} = 1 * Im;
			Pi_Mat{1,2} = 0 * Om;
			Pi_Mat{2,1} = 0 * Om;
			Pi_Mat{2,2} = -r_sq * Im;
		else 
			error('MK Error : Please Check the source code in freq_range');
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(conic_type, 'gen_circle')
		% General CIRCLE condition
		% FDI_info.eq = (x-p)^2 + (y-q)^2 <= r_sq
		p = LMI_Cond.FDI_info.p;
		q = LMI_Cond.FDI_info.q;
		r_sq = LMI_Cond.FDI_info.r_sq;
	
		if strcmp(SV_type, 'BD')
			Pi_Mat{1,1} = 1 * Im;
			Pi_Mat{1,2} = -( p * Im - q * Jm );
			Pi_Mat{2,1} = -( p * Im + q * Jm );
			Pi_Mat{2,2} = ( (p^2+q^2) - r_sq ) * Im;
	
		elseif strcmp(SV_type, 'CD')
			Pi_Mat{1,1} = 1 * Im;
			Pi_Mat{1,2} = -( p * Im + q * Jm );
			Pi_Mat{2,1} = -( p * Im - q * Jm );
			Pi_Mat{2,2} = ( (p^2+q^2) - r_sq ) * Im;
	
		else 
			error('MK Error : Please Check the source code in freq_range');
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % For the main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               END OF THIS PROGRAM                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
