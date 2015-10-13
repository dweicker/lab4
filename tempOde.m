function [timeStep,cpuTime,hMax] = tempOde(efficiency)
%TEMPODE This function compares two built-in ode solvers (ode23 and ode23s)
%for the solution of the temperature problem
N = [10 20 40];
timeStep = zeros(3,2);
cpuTime = zeros(3,2);
hMax = zeros(3,2);

for i=1:3
    u0 = zeros(N(i),1);
    %We set, if needded, the options for ode23 and ode23s
    switch efficiency
        case 'sparse'
            e = ones(N(i),1);
            S = spdiags([e e e],-1:1,N(i),N(i));
            A = spdiags([e -2*e e],-1:1,N(i),N(i));
            A(N(i),N(i)-1) = 2;
            options = odeset('Jacobian',N(i)*N(i)*A,'JPattern',S);
        otherwise
            options = [];
    end
    %We use ode23
    tstart = tic;
    [t,~] = ode23(@(t,u) tempe(N(i),t,u),[0 2],u0,options);
    cpuTime(i,1) = toc(tstart); 
    %We use ode23s
    tstart = tic;
    [tStiff,~] = ode23s(@(t,u) tempe(N(i),t,u),[0 2],u0,options);
    cpuTime(i,2) = toc(tstart); 
    %We compute the time steps
    h = diff(t);
    timeStep(i,1) = length(h);
    hMax(i,1) = max(h);
    hStiff = diff(tStiff);
    timeStep(i,2) = length(hStiff);
    hMax(i,2) = max(hStiff);
end

end

function dudt = tempe(Nx,t,u)
hx = 1/Nx;

e = ones(Nx,1);
A = spdiags([e -2*e e],-1:1,Nx,Nx);
A(Nx,Nx-1) = 2;
b = sparse(Nx,1);
b(1) = (t<=1);

dudt = (A*u+b)./(hx*hx);
end
