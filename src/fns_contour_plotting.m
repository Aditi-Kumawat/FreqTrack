%% All functions needed for the repo
classdef fns_contour_plotting
    methods (Static)
        function [K_eqt_0,fn_phi]=get_K_eqt_0_0n(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho,c,K1_b,K_b,K_P,c_P)
            %% This function file calculates the value of K_eqt_0 for the given input parameters
            % omega_vect is in rad/s

            K_RP_vect=K_P+1i*c_P.*omega_0;

            fn_phi=zeros(1,length(omega_0));
            for k=1:length(omega_0)

                omega=omega_0(1,k);
                %     Disp=['omega = ',num2str(omega)];
                %     disp(Disp)
                rho_b=(rho*omega^2*L_S^4)/(E_S*I_S);   %non-dimensionalised inertial term
                c_b=(c*omega*L_S^4)/(E_S*I_S);         %non-dimensionalised damping term

                a1_b=K1_b;
                a2_b=(K_b-rho_b+1i*c_b);

                m1=sqrt((a1_b+sqrt(a1_b^2-4*a2_b))*0.5);
                m2=sqrt((a1_b-sqrt(a1_b^2-4*a2_b))*0.5);
                m3=-m1;
                m4=-m2;

                M=zeros(8,8);
                M(1,:)=[m1^2 m2^2 m3^2 m4^2 0 0 0 0];
                M(2,:)=[m1^3 m2^3 m3^3 m4^3 0 0 0 0];
                M(3,:)=[exp(a_b*m1) exp(a_b*m2) exp(a_b*m3) exp(a_b*m4) -exp(a_b*m1) -exp(a_b*m2) -exp(a_b*m3) -exp(a_b*m4)];
                M(4,:)=[m1*exp(a_b*m1) m2*exp(a_b*m2) m3*exp(a_b*m3) m4*exp(a_b*m4) -m1*exp(a_b*m1) -m2*exp(a_b*m2) -m3*exp(a_b*m3) -m4*exp(a_b*m4)];
                M(5,:)=[m1^3*exp(a_b*m1) m2^3*exp(a_b*m2) m3^3*exp(a_b*m3) m4^3*exp(a_b*m4) -m1^3*exp(a_b*m1) -m2^3*exp(a_b*m2) -m3^3*exp(a_b*m3) -m4^3*exp(a_b*m4)];
                M(6,:)=[m1^2*exp(a_b*m1) m2^2*exp(a_b*m2) m3^2*exp(a_b*m3) m4^2*exp(a_b*m4) -m1^2*exp(a_b*m1) -m2^2*exp(a_b*m2) -m3^2*exp(a_b*m3) -m4^2*exp(a_b*m4)];
                M(7,:)=[0 0 0 0 m1^3*exp(m1/2) m2^3*exp(m2/2) m3^3*exp(m3/2) m4^3*exp(m4/2)];
                M(8,:)=[0 0 0 0 m1*exp(m1/2) m2*exp(m2/2) m3*exp(m3/2) m4*exp(m4/2)];

                %Right-hand column matrix
                % % % % % % %     R=[0;0;0;0;(P*L^3/E/I);0;0;0]; %%%%%%%%%%
                phi_c_vect=[0;0;0;0;(L_S^3/E_S/I_S);0;0;0];  %phi_c_vect :column vector

                % A=[A1;A2;A3;A4;B1;B2;B3;B4];
                A=M\phi_c_vect;
                A1=A(1,1);A2=A(2,1);A3=A(3,1);A4=A(4,1);
                % B1=A(5,1);B2=A(6,1);B3=A(7,1);B4=A(8,1);
                fn_phi(1,k)=A1*exp(m1*x_b)+A2*exp(m2*x_b)+A3*exp(m3*x_b)+A4*exp(m4*x_b);   %x_b=a_b, I don't remember why I wrote it as x_b
            end

            K_sub_vect=-(1./fn_phi)./S;
            K_eqt_0=(K_RP_vect.*K_sub_vect)./(K_RP_vect+K_sub_vect);
            % K_eqt_vect=K_sub_vect;

        end


        function Int=get_w_omega_contour(omega_vect_1,vel,x_R,E_R,I_R,rho_R,c_R,omega_vect,K_eqt_vect)
            %% Just w(omega) for Proposed Model
            ha_k=@fns_contour_plotting.get_Komega_vect_intrpol;
            K_omega_vect=ha_k(omega_vect,K_eqt_vect,omega_vect_1);
            Int=((vel^3).*exp(-1i*(omega_vect_1).*(x_R/vel)))./(E_R*I_R*(omega_vect_1).^4+vel^4.*K_omega_vect-rho_R.*(omega_vect_1).^2.*vel^4+1i*(omega_vect_1)*c_R*vel^4);
        end

        function K_omega_vect=get_Komega_vect_intrpol(omega_vect,K_eqt_vect,omega_vect_1)
            %% This function file calculates the value of K(omega) at any given value of
            % omega through linear interpolation.

            o_lim_1=omega_vect(1,1);
            o_lim_2=omega_vect(1,end);
            o_0=omega_vect(1,2)-omega_vect(1,1);
            K_omega_vect=zeros(1,length(omega_vect_1));

            for j=1:length(omega_vect_1)
                omega=omega_vect_1(1,j);
                if omega>=o_lim_1 && omega<o_lim_2
                    i_num=((omega-o_lim_1)/o_0)+1;
                    i=floor(i_num);
                    K_omega_vect(1,j)=K_eqt_vect(i)+((K_eqt_vect(i+1)-K_eqt_vect(i))/(omega_vect(i+1)-omega_vect(i)))*(omega-omega_vect(i));
                else
                    K_omega_vect(1,j)=0;
                end
            end
        end


        function [w_st_max,w_vect_st]=get_static_max(x_vect,K_eqt_0,P_W,E_R,I_R)
            %% This function file calculates the magnitude of the static deflection (w_vect_st) and
            % maximum static deflection (w_st_max) for the rail beam on winkler
            % foundation(K_eqt_0=K(omega=0))

            w_vect_st=zeros(1,length(x_vect));
            % BM_vect_st=zeros(1,length(x_vect));
            K_b=K_eqt_0/E_R/I_R;

            m_1=exp(pi*1i/4)*K_b^0.25;
            m_2=exp(-pi*1i/4)*K_b^0.25;

            M=zeros(2,2);
            M(1,:)=[m_1 m_2];
            M(2,:)=[m_1^3 m_2^3];

            %Right-hand column matrix
            R=[0;P_W/(2*E_R*I_R)];

            % A=[A2;A4];
            A=M\R;
            A1=A(1,1);A2=A(2,1);

            for k=1:length(x_vect)
                x_b=x_vect(1,k);
                if x_b<0
                    w_vect_st(1,k)=A1*exp(m_1*x_b)+A2*exp(m_2*x_b);
                    %         BM_vect_st(1,k)=E_R*I_R*(m_1^2*A1*exp(m_1*x_b)+m_2^2*A2*exp(m_2*x_b));
                else
                    w_vect_st(1,k)=A1*exp(-m_1*x_b)+A2*exp(-m_2*x_b);
                    %         BM_vect_st(1,k)=E_R*I_R*(m_1^2*A1*exp(-m_1*x_b)+m_2^2*A2*exp(-m_2*x_b));
                end
            end
            w_st_max=max(abs(w_vect_st));
            % BM_st_max=max(abs(BM_vect_st));
        end

        function int_4=get_w_rail_omega(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,omega_vect,K_eqt_vect)

            %% This function file calculates the integrand (int_4) for the inverse fourier
            % transform of the deflection function i.e rail deflection as a function of omega
            % (W(omega)), when K is function of
            % omega and K(omega) has been calculated using linear interpolation.

            ha_k=@fns_contour_plotting.get_Komega_vect_intrpol;
            K_omega_vect=ha_k(omega_vect,K_eqt_vect,omega_vect_1);
            % I_1=((vel^3).*exp(-1i*(omega_vect-1i*o_I).*(x_R/vel)))./(E_R*I_R*(omega_vect-1i*o_I).^4+vel^4.*K_omega_vect-rho_R.*(omega_vect-1i*o_I).^2.*vel^4);
            I_1=((vel^3).*exp(-1i*(omega_vect_1-1i*o_I).*(x_R/vel)))./(E_R*I_R*(omega_vect_1-1i*o_I).^4+vel^4.*K_omega_vect-rho_R.*(omega_vect_1-1i*o_I).^2.*vel^4+1i*(omega_vect_1-1i*o_I)*c_R*vel^4);
            int_4=I_1.*exp(1i*omega_vect_1*t_vect(1:end));
        end


        function int_4=get_w_rail_omega_Kcnst(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,K_eqt_0)
            %% This function file calculates the integrand for the inverse fourier
            % transform of the deflection function (W(omega)), when K is constant equal
            % to K_eqt_0

            % sc_rail_3n_Kcnst.m
            I_1=((vel^3).*exp(-1i*(omega_vect_1-1i*o_I).*(x_R/vel)))./(E_R*I_R*(omega_vect_1-1i*o_I).^4+vel^4.*K_eqt_0-rho_R.*(omega_vect_1-1i*o_I).^2.*vel^4+1i*(omega_vect_1-1i*o_I)*c_R*vel^4);
            int_4=I_1.*exp(1i*omega_vect_1*t_vect(1:end));
        end

        function w_t=get_P_t(t_vect,w_t_vect,t_ny)
            %% This function file finds the value of sleeper load (P_t) at any given value
            % time (t) using linear interpolation.
            % Uses the "P_t_vect" and "t_vect" generated for any given velocity from
            % the code sc_P_t_vect.m

            t_lim_1=t_vect(1,1);
            t_lim_2=t_vect(1,end);
            d_t0=t_vect(1,2)-t_vect(1,1);

            if t_ny>=t_lim_1 && t_ny<t_lim_2
                i_num=((t_ny-t_lim_1)/d_t0)+1;
                i=floor(i_num);
                w_t= w_t_vect(i)+ ((w_t_vect(i+1)-w_t_vect(i))/(t_vect(i+1)-t_vect(i)))*(t_ny-t_vect(i));
            else
                w_t=0;
            end

        end
        function w_t_vect=get_wrail_forP(P_W,vel,t_vect,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b,rho_R,x_R,E_R,I_R,dr,omega_vect,K_eqt_vect)
            %% This script file calculates the deflection of rail beam in time domain (w_t_vect_4),
            % at a given values of x_R at a given velocity(vel)/velocity ratio, as the wheel moves over the rail beam.
            % and the max value of deflection (W_ABS_MAT_3n_Kcnst) at all the input values of v_R_vect.

            % Input parameters: sc_IPn_3 --> K_S=1e7 N/m^2
            % Rail beam UIC60
            % K_eqt is constant= K_eqt_0 --> Winkler Model

            %%%%%%%%%%%%%%%%% K(OMEGA=0) %%%%%%%%%%%%%%%%%%%%%
            ha_k_st=@fns_contour_plotting.get_K_eqt_0_val1;
            omega_0=0;
            K_eqt_0=real(ha_k_st(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b));

            v_cr_con_3n=(4*E_R*I_R*K_eqt_0/(rho_R^2))^0.25;
            c_cr_con=2*sqrt(K_eqt_0*rho_R);
            c_R=dr*c_cr_con;
            % lambda_1=(K_eqt_0/4/E_R/I_R)^0.25;
            % v_R_vect=linspace(0,2,1e3);
            % v_R_vect=[1 1.1];
            % dr_vect=[0.05 0.3 1.1];

            %------------------ STATIC-ANALYSIS ------------------%
            % nd_vect=-10:0.1:10;                   % Normalised distance in Mallik's graphs
            %------------------------------------------------------%
            % t_vect=-1:0.01:1;
            % x_vect=nd_vect/lambda_1;

            v_R=vel/v_cr_con_3n;
            x_vect=t_vect*vel;
            ha_st=@fns_contour_plotting.get_static_max;
            [w_st_max_val1,w_vect_st]=ha_st(x_vect,K_eqt_0,P_W,E_R,I_R);

            lim_1=0;
            lim_2=1e4;

            %%%%%%%%%% NOTE: OMEGA_VECT_INPUT IS IN RADIANS/SECOND %%%%%%%%%%%%
            o_I=0.0001;
            ha_1=@fns_contour_plotting.get_w_rail_omega;

            if v_R>=0 && v_R<=0.000800160032006401
                w_t_vect=w_vect_st;
            else
                In_t_vect=real(integral(@(omega_vect_1)ha_1(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,omega_vect,K_eqt_vect),lim_1,lim_2,'ArrayValued',true));

                In_t_vect_1=In_t_vect.*exp(o_I*t_vect);

                w_t_vect=fliplr(-(P_W/pi).*In_t_vect_1);    % Deflection when K is function of omega + No Damping
            end
        end

        function [K_eqt_0,fn_phi]=get_K_eqt_0_val1(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho,c,K1_b,K_b)
            %% This function file calculates the value of K_eqt_0 for the given input
            % parameters

            % omega_vect is in rad/s

            % K_RP_vect=K_P+1i*c_P.*omega_0;

            fn_phi=zeros(1,length(omega_0));
            for k=1:length(omega_0)

                omega=omega_0(1,k);
                %     Disp=['omega = ',num2str(omega)];
                %     disp(Disp)
                rho_b=(rho*omega^2*L_S^4)/(E_S*I_S);   %non-dimensionalised inertial term
                c_b=(c*omega*L_S^4)/(E_S*I_S);         %non-dimensionalised damping term

                a1_b=K1_b;
                a2_b=(K_b-rho_b+1i*c_b);

                m1=sqrt((a1_b+sqrt(a1_b^2-4*a2_b))*0.5);
                m2=sqrt((a1_b-sqrt(a1_b^2-4*a2_b))*0.5);
                m3=-m1;
                m4=-m2;

                M=zeros(8,8);
                M(1,:)=[m1^2 m2^2 m3^2 m4^2 0 0 0 0];
                M(2,:)=[m1^3 m2^3 m3^3 m4^3 0 0 0 0];
                M(3,:)=[exp(a_b*m1) exp(a_b*m2) exp(a_b*m3) exp(a_b*m4) -exp(a_b*m1) -exp(a_b*m2) -exp(a_b*m3) -exp(a_b*m4)];
                M(4,:)=[m1*exp(a_b*m1) m2*exp(a_b*m2) m3*exp(a_b*m3) m4*exp(a_b*m4) -m1*exp(a_b*m1) -m2*exp(a_b*m2) -m3*exp(a_b*m3) -m4*exp(a_b*m4)];
                M(5,:)=[m1^3*exp(a_b*m1) m2^3*exp(a_b*m2) m3^3*exp(a_b*m3) m4^3*exp(a_b*m4) -m1^3*exp(a_b*m1) -m2^3*exp(a_b*m2) -m3^3*exp(a_b*m3) -m4^3*exp(a_b*m4)];
                M(6,:)=[m1^2*exp(a_b*m1) m2^2*exp(a_b*m2) m3^2*exp(a_b*m3) m4^2*exp(a_b*m4) -m1^2*exp(a_b*m1) -m2^2*exp(a_b*m2) -m3^2*exp(a_b*m3) -m4^2*exp(a_b*m4)];
                M(7,:)=[0 0 0 0 m1^3*exp(m1/2) m2^3*exp(m2/2) m3^3*exp(m3/2) m4^3*exp(m4/2)];
                M(8,:)=[0 0 0 0 m1*exp(m1/2) m2*exp(m2/2) m3*exp(m3/2) m4*exp(m4/2)];

                %Right-hand column matrix
                % % % % % % %     R=[0;0;0;0;(P*L^3/E/I);0;0;0]; %%%%%%%%%%
                phi_c_vect=[0;0;0;0;(L_S^3/E_S/I_S);0;0;0];  %phi_c_vect :column vector

                % A=[A1;A2;A3;A4;B1;B2;B3;B4];
                A=M\phi_c_vect;
                A1=A(1,1);A2=A(2,1);A3=A(3,1);A4=A(4,1);
                % B1=A(5,1);B2=A(6,1);B3=A(7,1);B4=A(8,1);
                fn_phi(1,k)=A1*exp(m1*x_b)+A2*exp(m2*x_b)+A3*exp(m3*x_b)+A4*exp(m4*x_b);   %x_b=a_b, I don't remember why I wrote it as x_b
            end

            K_eqt_0=-(1./fn_phi)./S;
            % K_eqt_0=(K_RP_vect.*K_sub_vect)./(K_RP_vect+K_sub_vect);
            % K_eqt_vect=K_sub_vect;

        end

        function w_omega=get_w_omega_contour_Wink(omega_vect_1,vel,x_R,E_R,I_R,rho_R,c_R,K_eqt_0)
            %% Just w(omega) for unit loading for Winkler Model
            % sc_W_Kconst.m

            w_omega=((vel^3).*exp(-1i*(omega_vect_1).*(x_R/vel)))./(E_R*I_R*(omega_vect_1).^4+vel^4.*K_eqt_0-rho_R.*(omega_vect_1).^2.*vel^4+1i*(omega_vect_1)*c_R*vel^4);

        end
    end
end