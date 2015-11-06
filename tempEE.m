function [u,x,t] = tempEE(Nx,ht,tend)
%TEMPEE Solve the temperature problem with Euler Explicit 
% This function solve the problem given in LAB4 with the FDM for the
% spatial derivatives and with Euler Explicit method for the time
% derivative
% INPUT : hx is the spatial stepsize
%         ht is the time stepsize
%         tend is the final time  
% OUTPUT : u is the numerical solution
%          x is the number of points in the x-axis
%          t is the discretized t-axis (from 0 to tend)

hx= 1/Nx;
coeff = ht/(hx*hx);
x = 0:hx:1; %length = Nx+1
t = 0:ht:tend;

u = zeros(Nx+1,floor(tend/ht+1)); 

e = ones(Nx,1);
A = spdiags([e -2*e e],-1:1,Nx,Nx);
A(Nx,Nx-1) = 2;
b = sparse(Nx,1);

i = 1;
for ti = ht:ht:tend
    b(1,1) = boundCondLeft(ti);
    u(2:Nx+1,i+1) = (speye(Nx)+coeff*A)*u(2:Nx+1,i)+coeff*b;
    i = i+1;
end
u(1,:) = boundCondLeft(0:ht:tend);


% Time animation of the temperature :)
% for i=1:length(t)
%     titre = sprintf('Time t=%f',(i-1)*ht);
%     bar = u(:,i)*[1 1];
%     contourf(x,1:2,bar',0:0.01:1);title(titre);colorbar;caxis([0 1]);
%     F(i)=getframe;
% end
% movie(F);

u = flipud(u);
end

function U = boundCondLeft(t)
%Returns the boundary condition at the given time
U = (t<=1);
end
