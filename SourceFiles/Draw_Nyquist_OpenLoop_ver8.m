function [] = Draw_Nyquist_OpenLoop_ver8(...
                 system_cell, ...
		 LMI_Cond, ...
		 XLim, YLim, ...
		 save_filename, ...
		 FDI_ON);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw Nyquist Diagram and Check the constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultAxesFontSize',20)
set(0,'defaultTextFontSize',100)
set(0,'defaultAxesFontName','Helvetica')
set(0,'defaultTextFontName','Helvetica')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wL =   1*10^(-2) * 2*pi;
wU =  20*10^3    * 2*pi;
%Ts =  1/(20*10^3); % Sampling Period
%FDI_ON = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logspace(a, b, N) 
% 10^a <--> 10^b, N-sample
w_para = logspace(log10(wL), log10(wU), 100001);

% Make (x,y) data for Nyquist Diagram

%{
system_1 = system_cell{1};
system_2 = system_cell{2};
system_3 = system_cell{3};

[re_1, im_1] = nyquist(system_1, w_para);%
re_1 = squeeze(re_1); im_1 = squeeze(im_1);

[re_2, im_2] = nyquist(system_2, w_para);%
re_2 = squeeze(re_2); im_2 = squeeze(im_2);

[re_3, im_3] = nyquist(system_3, w_para);%
re_3 = squeeze(re_3); im_3 = squeeze(im_3);
%}

figure()
set(0, 'DefaultlineLineWidth', 2); % change the default line width
axis square; % 
hold on;


bline_options = {'b-', 'b--', 'b:'}; %i番目の条件 base line options
marker_options = {'pr', 'r*', 'dm', 'pr', 'r*', 'dm'}; %i番目の条件
cline_options = {'pr', 'r*', 'dm', 'pr', 'r*', 'dm'}; % constraint line options

if ( size(LMI_Cond, 1) >= 6 )
	error('MK error : too much LMI_Cond_set !');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            周波数特性の表示               %
%             Plot Base Line                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:size(system_cell, 2)
	% Plot Base Line
	[re, im] = nyquist(system_cell{k}, w_para);%
	re = squeeze(re); im = squeeze(im);
	plot(re, im, bline_options{k}); % Plot Base line
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            有限周波数領域の可視化         %
% Plot Marker to show the upper and lower bound of the frequency region. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:size(system_cell, 2)
	% Plot Finite Frequency Boundary
	for i=1:1:size(LMI_Cond, 2)
		LMI_Cond_set = LMI_Cond{i};
		if ( strcmp(LMI_Cond_set.freq_info.freq_range, 'LF') | ...
		     strcmp(LMI_Cond_set.freq_info.freq_range, 'HF') )
			w1 = LMI_Cond_set.freq_info.w;
			w2 = LMI_Cond_set.freq_info.w;
		else 
			w1 = LMI_Cond_set.freq_info.w1;
			w2 = LMI_Cond_set.freq_info.w2;
		end
		[re_temp, im_temp] = nyquist( system_cell{k}, [w1, w2] );%
		re_i{i} = squeeze(re_temp); im_i{i} = squeeze(im_temp);

		plot(re_i{i},im_i{i}, marker_options{i}, 'MarkerSize', 10); 
		
	end

end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             FDIの可視化                   %
% Plot and Show the contraint 1, 2 and 3    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (FDI_ON == 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:size(LMI_Cond, 2)
	LMI_Cond_set = LMI_Cond{i};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[1] Line type
	if strcmp(LMI_Cond_set.FDI_info.conic_type, 'line')
		
		% for Line
		a = LMI_Cond_set.FDI_info.a;
		b = LMI_Cond_set.FDI_info.b;
		c = LMI_Cond_set.FDI_info.c;
		
		if ( not(b == 0) ) 
			x = linspace(-8, 8, 1001);
%			y = -a/b * x - c/b; %constraint 1
			y = -a/b * x + c/b; %constraint 1
		else
			y = linspace(-8, 8, 1001);
%			x = -b/a * y - c/a;
			x = -b/a * y + c/a;
		end
		
		plot(x,y, 'r--'); % Plot the constraint 1
%		legend_string = [legend_string, 'const. 1(high gain in LF)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[2] circle type
	elseif strcmp(LMI_Cond_set.FDI_info.conic_type, 'circle')
		r = LMI_Cond_set.FDI_info.r;
		t = linspace(0,2*pi,1001);
		x = r * sin(t);   y = r * cos(t);
		
		plot(x,y, 'm--'); % Plot the constraint i
%		legend_string = [legend_string, 'const. 1(high gain in LF)'];
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[3] general circle type
	elseif strcmp(LMI_Cond_set.FDI_info.conic_type, 'gen_circle')
		r = LMI_Cond_set.FDI_info.r;
		p = LMI_Cond_set.FDI_info.p;
		q = LMI_Cond_set.FDI_info.q;
		t = linspace(0,2*pi,1001);
		x = r * sin(t) + p;   y = r * cos(t) + q;
		
		plot(x,y, 'm--'); % Plot the constraint i
%		legend_string = [legend_string, 'const. 1(high gain in LF)'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	else
		error('MK error : Unknown conic_type');
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            その他のラインの表示           %
% Plot and Show other lines.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot stability circle 
t = linspace(0,2*pi,1001);
x = sin(t); y = cos(t);
%legend_string = [legend_string, 'const. 2(stab. margin in MF)'];
plot(x,y, 'k:'); % Plot stability circle 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             グラフの装飾                  %
% Formatting and annotation                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%legend(legend_string, 'Location', 'SouthEast', 'FontSize', 14);
grid on;
xlabel('Real');
ylabel('Imag');

xlim(XLim);
ylim(YLim);
%xlim([-1.5, 0.5]);
%ylim([-1.5, 0.5]);
%ax = gca;
%set(ax,'XTick',-0.2:0.1:0.2)
%set(ax,'YTick',-0.2:0.1:0.2)


grid on;
axis square


% realization of string concatenation
%legend_string = '''Loop transfer function L(s)''';
%legend_string = [legend_string, ',''constraint 1(high gain in LF)'''];
%eval(['legend(', legend_string, ',', '''Location''', ',' '''SouthEast''', ')']);

%legend_string = {'Loop transfer function L(s)'};
legend_strings = {'\sigma = -0.3', '\sigma = 0', '\sigma = 0.3'};
legend(legend_strings, 'Location', 'NorthWest', 'FontSize', 15);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        グラフ及びデータの保存             %
% Save Graph and other data                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Graph  
filename = save_filename;
%filename = 'Nyquist_OpenLoop.eps'
%graph_type = '-depsc2';
graph_type = '-dpng';

if strfind(filename, 'png')
	graph_type = '-dpng';
elseif strfind(filename, 'eps')
	graph_type = '-depsc2';
else
	error('MK error : Unknown save filename extension ( eps, or png can be used )');

end


%filename = sprintf('Nyquist_robust_PID_%d.png', Iter_counter);
print(filename, graph_type);

%graphというフォルダーがカレントディレクトリに存在すれば、
if exist('graph', 'dir')  	
	movefile( filename,  './graph/'); %move the file 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        END OF THIS PROGRAM                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
