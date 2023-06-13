clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fig / Video setting
nfig=0; 
video_steps = 50; %200
take_video = 0; %1 to save video of individual cases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters
L = 3; %length of network in m
a = 0.2; 
khat = 1; %stiffness between agents
f1 = 0.08; f2 = 0.1; f3 = 10; %Hz source frequency
disp(['Simulating frequencies: ' num2str(f1) ' Hz, ' num2str(f2) ' Hz, ' num2str(f3) ' Hz.'])
% wave properties

T1 = 1/f1; % time period of propagated in seconds
T2 = 1/f2;
T3 = 1/f3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DSR parameters
kmy = 0;

for gamma_dsr = 1:5:51
    kmy = kmy+1;

    beta2_dsr = 1;   %=1 for DSR
    delta_t = 10^(-4); 
    dt = delta_t;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% constrained parameters 
    D = 1; %number of spatial dimensions
    n = 1*round(L/a) %number of agents

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% network definition
    % matrix A
    connection_struct = [-0.5 1 -0.5];
    A = zeros(n,n);
    A(1,1) = 1; A(1, 2) = -0.5; 
    A(n,n) = 1; A(n,n-1) = -1;
    for i=2:1:n-1
        A(i,i-1:1:i+1) = connection_struct;
    end

    % A(1,n) = -0.5;
    % A(n,1) = -0.5;
    % A(n,n-1) = -0.5;

    B = [0.5; zeros(n-1,1)];

    lambda_A = eig(A);

    video_name = ['Simulation gamma 1.mp4'];

    c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
    beta1_dsr = 0.9*((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A)

    zeta_dsr = (1-beta2_dsr)*L/(pi*c*beta2_dsr*dt) + beta1_dsr*pi*c/(4*gamma_dsr*L)
    omega_0 = pi*c/(2*L);
    predicted_settling_time = 6/(zeta_dsr*omega_0)

    %Momentum parameter
    gamma_m = gamma_dsr;
    beta2_m_1 = 0.999;
    beta2_m_2 = 1;
    beta1_m = 0;

    %Nesterov parameter
    gamma_n = gamma_dsr;
    beta2_n = 0.999;
    beta1_n = beta2_n*gamma_n*dt*(1+beta2_n);
    beta2_n_2 = 1;
    beta1_n_2 = gamma_n*dt*beta2_n_2*(1+beta2_n_2);

    v = sqrt(gamma_dsr*a^2/(2*1*delta_t*beta2_dsr)) %wave velocity in m/s

    tend = max(0.9*L/v, 6) %s - duration of simulation time -- BEFORE REFLECTION 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    predicted_settling_time = tend;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %time variable
    ts = 0:delta_t:T2/2; %time period of initial BC creating pulse
    t = [ts (ts(end)+delta_t):delta_t:tend]; %time vector

    %topological radius
    threshold_dist = 2*a;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% define input source
    input1 = sin(2*pi*ts/T2);
    input2 = 0.5*sin(2*pi*ts/T3);
    input = input1 + input2;
    initialzeros = 500; %2
    Is1 = [zeros(1,initialzeros) input1 zeros(1, length(t)-length(ts)-initialzeros)];
    Is2 = [zeros(1,initialzeros) input2 zeros(1, length(t)-length(ts)-initialzeros)];
    Is = [zeros(1,initialzeros) input zeros(1, length(t)-length(ts)-initialzeros)];

    Is_unitstep = ones(size(t));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %number of discrete x points
    num_of_x = n; %equal to number of agents in 1d network

    %number of time steps 
    num_of_t = length(t);

    Idsr1 = zeros(num_of_x, 3);  
    Im_1 = zeros(num_of_x,3);
    Im_2 = zeros(num_of_x, 3);
    Im_3 = zeros(num_of_x, 3);

    %DSR (Internal damping)
    thetadsr = 0*pi/2*ones(num_of_x, 3); %initially the angle is pi/2
    xdsr = linspace(a, num_of_x*a, num_of_x);
    ydsr = 0*xdsr;
    % for k = 0:14
    %     xdsr(1,k+1) = 1.5/pi*sin(k*24*pi/180);
    %     ydsr(1,k+1) = 1.5/pi*cos(k*24*pi/180);
    % end
    vel = 1; 
    leader_xytheta_dsr = [xdsr(1) ydsr(1) thetadsr(1,1);
                   xdsr(1) ydsr(1) thetadsr(1,2)];
    follower_xytheta_dsr = [xdsr(end) ydsr(end) thetadsr(end,1);
                   xdsr(end) ydsr(end) thetadsr(end,2)];
    time_vec = [t(1:2)];

    %M1 (viscous damping)
    thetam1 = thetadsr; xm1 = xdsr; ym1 = ydsr;
    leader_xytheta_m1 = leader_xytheta_dsr;
    follower_xytheta_m1 = follower_xytheta_dsr;

    %M2 (undamped)
    thetam2 = thetadsr; xm2 = xdsr; ym2 = ydsr;
    leader_xytheta_m2 = leader_xytheta_dsr;
    follower_xytheta_m2 = follower_xytheta_dsr

    %detect formation breaking
    formationbroken_dsr = 0;
    formationbroken_m1 = 0;
    formationbroken_m2 = 0;

    %defining variables to take snaps of wave at two instants
    snap1_taken = 0;
    snap2_taken = 0;
    snap3_taken = 0;
    snap4_taken = 0;
    snap5_taken = 0;
    snap6_taken = 0; 
    for k = [3:1:num_of_t] 
        t(k)/t(end)

        iii=3;

        %frequency f1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        bigA = get_topological_A_matrix(xdsr, ydsr, threshold_dist, n);

        Delta_kminus1 = bigA*[Idsr1(:,iii-1); Is(k-1)];
        Delta_kminus2 = bigA*[Idsr1(:,iii-2); Is(k-2)];

        Idsr1(:,iii) = Idsr1(:,iii-1) - gamma_dsr*delta_t*Delta_kminus1 - beta1_dsr*(Delta_kminus1-Delta_kminus2) ...
                            + beta2_dsr*(Idsr1(:,iii-1)-Idsr1(:,iii-2));% + gamma_dsr*delta_t*B(i)*Is(k-1); 


        thetadsr = Idsr1*pi/2; %thetadsr + pi/2*Idsr1;
        xdsr = xdsr + vel*delta_t*cos(thetadsr(:,3)');
        ydsr = ydsr + vel*delta_t*sin(thetadsr(:,3)');

        Idsr = Idsr1;
        f=f2;

        Idsr1(:,1:2) = Idsr1(:,2:3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        bigA = get_topological_A_matrix(xm1, ym1, threshold_dist, n);

        Delta_kminus1 = bigA*[Im_1(:,iii-1); Is(k-1)];
        Delta_kminus2 = bigA*[Im_1(:,iii-2); Is(k-2)];

        Im_1(:,iii) = Im_1(:,iii-1) - gamma_m*delta_t*Delta_kminus1 - beta1_m*(Delta_kminus1-Delta_kminus2) ...
                        + beta2_m_1*(Im_1(:,iii-1)-Im_1(:,iii-2));

        diag_A_m1 = diag(bigA(1:n,1:n));

        for i=1:1:n
           if (diag_A_m1(i) == 0)
               xm1(i) = xm1(i);
               ym1(i) = ym1(i);
           else
               thetam1(i,:) = Im_1(i,:)*pi/2; %thetadsr + pi/2*Idsr1;
               xm1(i) = xm1(i) + vel*delta_t*cos(thetam1(i,3)');
               ym1(i) = ym1(i) + vel*delta_t*sin(thetam1(i,3)');
           end
        end

        Im_1(:,1:2) = Im_1(:,2:3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        bigA = get_topological_A_matrix(xm2, ym2, threshold_dist, n);

        Delta_kminus1 = bigA*[Im_2(:,iii-1); Is(k-1)];
        Delta_kminus2 = bigA*[Im_2(:,iii-2); Is(k-2)];

        Im_2(:,iii) = Im_2(:,iii-1) - gamma_m*delta_t*Delta_kminus1 - beta1_m*(Delta_kminus1-Delta_kminus2) ...
            + beta2_m_2*(Im_2(:,iii-1)-Im_2(:,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 

        diag_A_m2 = diag(bigA(1:n,1:n));

        for i=1:1:n
            if (diag_A_m2(i) == 0)
                xm2(i) = xm2(i);
                ym2(i) = ym2(i);
            else
                thetam2(i,:) = Im_2(i,:)*pi/2; 
                xm2(i) = xm2(i) + vel*delta_t*cos(thetam2(i,3)');
                ym2(i) = ym2(i) + vel*delta_t*sin(thetam2(i,3)');
            end
        end
        Im_2(:,1:2) = Im_2(:,2:3);         
        if (t(k) >= 7 && snap6_taken == 0)
            k_6 = k;
            snap6_taken = 1;
            signal_6_undamped = Im_2(:,iii-1);
            signal_6_viscous = Im_1(:,iii-1);
            signal_6_internal = Idsr1(:,iii-1);
            signal_6_ideal = Im_3(:,iii-1);
        end
    end
    dy(kmy) = ydsr(1)-ydsr(15);    
end

gamma=1:5:51;
plot(gamma,dy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function
function bigA = get_topological_A_matrix(xdsr, ydsr, radius, N)
%returns graph lapacian (N rows and columns where last N+1 column is for the source node
%and connected only to the leader agent) based on current x y location of agents with
%topological dist of given radius 
    bigA = zeros(N, N+1);

    pos_of_agents = sqrt(xdsr.^2+ydsr.^2);

    for i=1:1:N
            for j=1:1:N
                if ( j ~= i && abs(pos_of_agents(j)-pos_of_agents(i)) < radius )
                    bigA(i,j) = bigA(i,j) - 1;                   
                end    
            end
    end   
    %leader agent
    bigA(1,N+1) = bigA(1, N+1) - 1;
    for i=1:1:N           
            sum_of_row = - sum(bigA(i,:));           
            if (sum_of_row > 0)
                bigA(i,:) = bigA(i,:)/sum_of_row;         
                bigA(i,i) = 1;
            end
    end
end

function vf = filtered_derivative(t, x, wf, n)
% FILTERED_dERIVATIVE  take derivative of a filtered signal.
%   vf = filtered_derivative(t, x, wf, n) filters x to xf and takes its
%   derivative vf

    A = -wf*eye(n); B = wf*eye(n);
    C = -wf*eye(n); D = wf*eye(n);
    
    sys = ss(A, B, C, D);
    
    vf = lsim(sys, x, t);

end

function [t, I, Idot] = second_order_simulation(t, Is, n, c_s, c_sd, cd)

        %number of discrete x points
            num_of_x = n; %equal to number of agents in 1d network

            %number of time steps 
            num_of_t = length(t);

            I = zeros(num_of_x, num_of_t);
            Idot = zeros(num_of_x, num_of_t);

            %adding the initial condition on left end with sin input
            Idsr(1,:) = Is;
            for k=3:1:num_of_t
               for i=1:1:(num_of_x) 

                   if i==1 
                       Deltai_kminus1 =  Idsr(i,k-1) - 0.5*Idsr(i+1,k-1) -0.5*Is(k-1) ; 
                       Deltai_kminus2 =  Idsr(i,k-2) - 0.5*Idsr(i+1,k-2) -0.5*Is(k-2) ; 
                   elseif i == num_of_x        
                       Deltai_kminus1 = -1*Idsr(i-1,k-1) + Idsr(i,k-1); 
                       Deltai_kminus2 = -1*Idsr(i-1,k-2) + Idsr(i,k-2);         
                   else
                       Deltai_kminus1 = -0.5*Idsr(i-1,k-1) + Idsr(i,k-1) - 0.5*Idsr(i+1,k-1);  
                       Deltai_kminus2 = -0.5*Idsr(i-1,k-2) + Idsr(i,k-2) - 0.5*Idsr(i+1,k-2); 
                   end       
                       Idsr(i,k) = Idsr(i,k-1) - gamma_dsr*delta_t*Deltai_kminus1 - beta1_dsr*(Deltai_kminus1-Deltai_kminus2) ...
                                + beta2_dsr*(Idsr(i,k-1)-Idsr(i,k-2));
               end
            end
end

