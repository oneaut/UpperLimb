function x = RK4Step( F,x0,u,tau,p )
%RK4Step(F,x0,u,tau) takes an anonymous function (F), the previous time 
%step (x0), and the input at that time step (u) and solves for the next 
%iteration using the fourth order Runge-Kutta method. The anonymous 
%function that is given in the input is the differential equation that 
%you are trying to solve.

%% Solve Using Runge-Kutta 

if nargin == 4;
    % Solve for coefficients
    k1 = F(x0,u);
    k2 = F(x0+k1*tau/2,u);
    k3 = F(x0+k2*tau/2,u);
    k4 = F(x0+k3*tau,u);
    
    % Solve for next time step
    x = x0 + (1/6)*tau*(k1+2*(k2+k3)+k4);
else
    % Solve for coefficients
    k1 = F(x0,u,p);
    k2 = F(x0+k1*tau/2,u,p);
    k3 = F(x0+k2*tau/2,u,p);
    k4 = F(x0+k3*tau,u,p);
    
    % Solve for next time step
    x = x0 + (1/6)*tau*(k1+2*(k2+k3)+k4);
end


end