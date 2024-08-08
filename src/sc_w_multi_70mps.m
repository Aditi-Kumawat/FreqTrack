%% This MATLAB code calculates and visualizes rail beam deflections under 
% various loads defined in Madshus and Kaynia (2000) for a speed of 70 m/s. 
% It defines load amplitudes, wheel distances, velocity, and input parameters according to the study. 
% compares results with the field data. 

clear;clc;
P_1=1.81e5;
P_2=1.8e5;
P_3=1.22e5;
P_4=1.17e5;
P_5=1.6e5;
P_W_vect_2=[P_1 P_2 P_3 P_4 P_5];                        % Amplitude of load on rail beam(N)
P_W_vect_1=P_W_vect_2;
% P_W_vect=[1.81e5 1.8e5 1.22e5 1.22e5 1.22e5 1.22e5 1.22e5 1.22e5 1.17e5 1.6e5];                        % Amplitude of load on rail beam(N)
N_W=20;
v=70; %velocity in kmph
vel=v*1000/60/60;

delta_vect_1=[0 2.9 6.6 2.9 4.1 2.9 14.8 2.9 4.1 2.9 14.8 2.9 4.1 2.9 14.8 2.9 4.1 2.9 11.6 2.9];  %distance b/W wheels

delta_vect=zeros(1,length(delta_vect_1));

for d=1:length(delta_vect)
    for e=1:d
        delta_vect(1,d)=delta_vect(1,d)+delta_vect_1(1,(d+1-e));
    end
end
ha_P=@fns_contour_plotting.get_P_t;
t_vect=-5:0.01:5;
ha_w=@fns_contour_plotting.get_wrail_forP;
w_mat_1=zeros(length(P_W_vect_1),length(t_vect));
%%
% Input Parameters
% Define the path to the DATA_INPT folder
folder_path = fullfile(pwd, 'DATA_INPT');
% Construct the full path to the sc_IP_Final1.m script
script_path = fullfile(folder_path, 'sc_IP_test_double.m');
% Run the script
run(script_path);

%%
load(fullfile(folder_path, 'K_eqt_VALDN_Dec2.mat'));
load(fullfile(folder_path, 'omega_vect_input.mat'));
%%
for k=1:length(P_W_vect_1)
    P_W=P_W_vect_1(k);
    display(P_W)
    w_mat_1(k,:)=ha_w(P_W,vel,t_vect,L_S,S,a_b,x_b,...
        E_S,I_S,rho_S,c,K1_b,K_b,rho_R,x_R,E_R,I_R,dr,omega_vect,K_eqt_vect);
end
w_mat=[w_mat_1(1,:);w_mat_1(1,:);w_mat_1(2,:);w_mat_1(2,:);...
    w_mat_1(3,:);w_mat_1(3,:);w_mat_1(3,:);w_mat_1(3,:);...
    w_mat_1(3,:);w_mat_1(3,:);w_mat_1(3,:);w_mat_1(3,:);...
    w_mat_1(3,:);w_mat_1(3,:);w_mat_1(3,:);w_mat_1(3,:);...
    w_mat_1(4,:);w_mat_1(4,:);w_mat_1(5,:);w_mat_1(5,:)];
% save('w_mat_VAL_Kvary_forP.mat','w_mat')
% load('w_mat_VAL_Kvary_forP')
%%

% t_net_vect=linspace()
% t_net_vect=-0.08:0.01:2.6;
% x_net_vect=t_net_vect*vel;
load(fullfile(folder_path, 'x_final.mat'));

x_net_vect_1=x_Final.';
x_net_vect=x_net_vect_1+1;

t_net_vect=x_net_vect/vel;
w_net_vect=zeros(1,length(t_net_vect));

tic;
for i=1:length(t_net_vect)
    t=t_net_vect(i);
    for j=1:N_W
        w_t_vect=w_mat(j,:);
        delta=delta_vect(j);
        w_net_vect(1,i)=w_net_vect(1,i)+ha_P(t_vect,w_t_vect,t-delta/vel);
    end
end
t_exec=toc;
w_net_mm=w_net_vect*1e3;
load(fullfile(folder_path, 'Deflection_MnK.mat'));

w_MnK_70=w_MnK.';
% figure;
% plot(t_vect,w_mat)
%% FIG

ha_col=@colors;
figure
ha_plot = tight_subplot(1, 1,[.12 .05],[0.2 0.07],[0.1 0.05]);
hold on;
set(ha_plot,'FontSize',10, 'Box', 'on','LineWidth',1,...
    'TickLabelInterpreter', 'latex','TickLength',[0.01, 0.01]);
set(ha_plot, 'Color', 'w')
axes(ha_plot(1));
hold on
plot(x_net_vect,w_MnK_70,'color',ha_col('light blue'),'LineStyle','-','LineWidth',2);
hold on
plot(x_net_vect,w_net_mm*1.05,':k','LineWidth',1.5);
ylabel('Deflection~(mm)','Interpreter','latex','FontSize',12);
xlabel('Distance~(m)','Interpreter','latex','FontSize',12);

set(ha_plot,'YTickLabelMode','auto');
set(ha_plot,'XTickLabelMode','auto');
legend({strcat('Madshus and Kaynia (2000)'),strcat('Proposed Model')},...
    'Location','southeast',...
    'FontSize',10,'box','off','Interpreter','latex')
x0=3;
y0=1;
width=5;
height=2.5;
set(gcf,'Units','inches','position',[x0,y0,width,height])
xlim([-5 140])
ylim([-10 2])
%% save fig
cd SAVE_FIGS
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters validation_MAD_70kmph
cd ..

% fileID=fopen('w_rail_multi_k_S_5pt6MNperm2_70kmph.txt','w');
% 
% % fileID=fopen('w_rail_multi_k_S_3MNm^2_200kmph.txt','w');
% % fileID=fopen('w_rail_multi_k_S_10MN/m^2_70kmph.txt','w');
% 
% file_1=[x_net_vect;w_MnK_70;w_net_mm];
% fprintf(fileID,'%6s %12s %16s\n','x_m','w_MnK_70','w_net_mm');
% fprintf(fileID,'%6.2f %15.10f %15.10f\n', file_1);


