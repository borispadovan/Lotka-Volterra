% =========================================================================
% Matlab program to draw a phase portrait for the LotkaVolterra Predetor
% Pey model. The results are shown at each time step. 
% In addition, the user has the option of plotting a time series graph for x or y. 
% Set the parameter choice = 1 for a time series plot for x.
% Set choice = 2 for a time series plot for y. 
% Equation parameters alpha, beta, gamma and detla can be changed by the
% user.
% Equations are solved using a numerical ODE solver. 
%
% James Adams 3/4/14
% =========================================================================

function LotkaVolterra_JA
clear  % Clears command window
clc    % Clears command history
clf    % Removes anything in the figure window before simulation. 


iterations = 1;  % Sets initial interation count to 1;
pausetime = 0;%0.1;  % Shows solutions at each time step. 
runtime = 580;    % Duration time of simulation.

% ================ Equation parameter values ==============================
N1  = 0.15;     % preys growth factor
L12 = 0.0001;   % loss of preys when preys meet predators
C11 =10^(-8);  % competition between preys
L13 = 0.00012;

M2  = 1.5;   % predators decrease factor
G21 = 0.00001;  % growth of predators due to 
C22 = 10^(-8);      % competition between predators
C23 = 10^(-5);

M3  = 1.46;
G31 = 10^(-5);
C33 = 0;
C32 =10^(-7);     % competition with native predator

% =============== Initial conditions for x and y ==========================
initialx = M2/G21; % Balance solution
initialy = N1/L12;  % Balance solution
initialz = 8; % Alien predator
%initialx = 150000; % initial number of preys 
%initialy = 1600;   % initial number of predators 

fprintf('----------------------------------\nLotka-Volterra Predetor Prey model \n\n----------------------------------')
%fprintf('\n\nParameter values set,')
%fprintf('\n\nalpha = %2.6f \nbeta = %2.6f \ngamma = %2.6f \ndelta = %2.6f ',alpha,beta,gamma,delta)


% Solves equations using numerical ODE solver 45 (nonstiff runge kutta)
deq1=@(t,x) [...
     x(1)*N1 - L12*x(1)*x(2) - C11*x(1)*x(1)- L13*x(1)*x(3) ; ...
    -x(2)*M2 + G21*x(1)*x(2) - C22*x(2)*x(2)- C23*x(2)*x(3) ;...
    -x(3)*M3 + G31*x(1)*x(3) - C33*x(3)*x(3)- C32*x(3)*x(2) ];
[t,sol] = ode45(deq1,[0 runtime],[initialx initialy initialz]);

arraysize = size(t);  % Sets time array size for the for loop.

%============ Solutions are plotted at each time step =====================

for i = 1 : max(arraysize) 
    subplot(4,1,1)
    plot(sol(iterations,1),sol(iterations,2),'.','color',[rand; rand; rand],'markersize',14,'MarkerFaceColor','b');                               
    %plot(sol(iterations,1),sol(iterations,2),'.','color','k','markersize',14,'MarkerFaceColor','b');                               
    hold on
    title(['Lotka-Volterra Equations       t = ' num2str(t(iterations)) ' months'],'fontsize',12)
    xlabel('x','fontsize',12)
    ylabel('y','fontsize',12)
    axis([min(sol(:,1)) max(sol(:,1)) min(sol(:,2)) max(sol(:,2))])
    
    subplot(4,1,2)
    %text(0.1,0.5,'Time Series graph will be shown at the end of the simulation')
    
    iterations = iterations + 1;   % Adds 1 to the iteration count. 
    pause(pausetime)
end

% ==== Plots time series of x or y graph depending on choice ============== 

    subplot(4,1,2)
    plot(t(:,1),sol(:,1),'b.','markersize',10,'MarkerFaceColor','b')
    title(['Time series for x'],'fontsize',12)
    xlabel('t [months]')
    ylabel('x')
    axis([min(t(:,1)) max(t(:,1)) min(sol(:,1)) max(sol(:,1))])

    %title(['Time series for y - run time = ' num2str(max(t)) ' months '],'fontsize',12)
    subplot(4,1,3)
    plot(t(:,1),sol(:,2),'r.','markersize',10,'MarkerFaceColor','b')
    title(['Time series for y'],'fontsize',12)
    xlabel('t [months]')
    ylabel('y')
    axis([min(t(:,1)) max(t(:,1)) min(sol(:,2)) max(sol(:,2))])
    
    subplot(4,1,4)
    plot(t(:,1),sol(:,3),'g.','markersize',10,'MarkerFaceColor','b')
    title(['Time series for z'],'fontsize',12)
    xlabel('t [months]')
    ylabel('z')
    if(initialz>0)
        axis([min(t(:,1)) max(t(:,1)) min(sol(:,3)) max(sol(:,3))])
    end
    
end

%====================== End of program ====================================
