close all;clear;clc;
% Input Parameters
% Define the path to the DATA_INPT folder
folder_path = fullfile(pwd, 'DATA_INPT');
% Construct the full path to the sc_IP_Final1.m script
script_path = fullfile(folder_path, 'sc_IP_test2.m');
% Run the script
run(script_path);

% Find K_eqt
omega_0=-1e3:0.02:1e3;
K_eqt_1=real(fns_contour_plotting.get_K_eqt_0_val1(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b));
