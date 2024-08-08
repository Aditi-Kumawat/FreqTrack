% This code plots the normalized rail deflection in frequency domain 
% with velocity ratio in a contour plot for the Winkler model.

clear;clc;
% Input Parameters
% Define the path to the DATA_INPT folder
folder_path = fullfile(pwd, 'DATA_INPT');
% Construct the full path to the sc_IP_Final1.m script
script_path = fullfile(folder_path, 'sc_IP_Final1.m');
% Run the script
run(script_path);

% Axle Load
P_W=1.25e5;                        % Amplitude of load on rail beam(N)

% K_static Value
ha_k_st=@fns_contour_plotting.get_K_eqt_0_0n;
omega_0=0;
K_eqt_0=real(ha_k_st(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b,K_P,c_P));

% Parameters defined for normalisation
v_cr_con=(4*E_R*I_R*K_eqt_0/(rho_R^2))^0.25;
c_cr_con=2*sqrt(K_eqt_0*rho_R);
lambda_1=(K_eqt_0/4/E_R/I_R)^0.25;

% Input velocity/damping parameters
v_R_vect=linspace(0,2,1e3);
% v_R_vect=0.5; %ALPHA
vel_vect=v_R_vect.*v_cr_con;
% dr_vect=[0.05 0.3 1.1];

% Damping Ratio
dr_vect=0.05;

% Input damping parameters
c_R_vect=dr_vect*c_cr_con;

% Input normalised distance/time
nd_vect=-6:0.1:6;                   % Normalised distance in Mallik's graphs
x_vect=nd_vect/lambda_1;

% Input Frequecny
omega_vect_1=linspace(0,1.5e3,2.2e3);

% Function for calculating w(omega) for Winkler Model
ha_1=@fns_contour_plotting.get_w_omega_contour_Wink;

int_2D=zeros(length(vel_vect),length(omega_vect_1));
int_3D=zeros(length(vel_vect),length(omega_vect_1),length(c_R_vect));

tic;
for k=1:length(c_R_vect)
    c_R=c_R_vect(1,k);
    for i=1:length(vel_vect)
        vel=vel_vect(1,i);
        v_R=v_R_vect(1,i);
        display(v_R)
        int_3D(i,:,k)=ha_1(omega_vect_1,vel,x_R,E_R,I_R,rho_R,c_R,K_eqt_0);
    end
end

disp_3D=P_W*(int_3D);
disp_3D=abs(disp_3D);

for i_v=1:length(vel_vect)
    disp_3D(i_v,:)=disp_3D(i_v,:)/max(disp_3D(i_v,:));
end

% Script file for plotting Contour
plot_contour_hot;

% saveas(gca, fullfile(figuresdir, 'Contour_Kcnst'), 'jpg');
% saveas(gca, fullfile(figuresdir, 'Contour_Kcnst'), 'fig');
t_exec=toc;
