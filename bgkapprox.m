% tsteps: number of time steps Nx: number of x points Nv: number of v
% points, here we are going to use only Nv=1 nu is used to find deltat a is
% the left end of the space interval b is the right end of the end interval
% v is the vector velocity, for this example v =1 Tfinal is the final time
% Nv is the number of grid points in the velocity

function [Fmatrix] = bgkapprox(Nx,Nv,nu,a,b,Tfinal)
 A = [0,0,0;0,0,0;0,1,0];
 w = [0,0.5,0.5];
 %% Don't use vpa **********************************************
 v = GaussLegendre(Nv);
 alpha=40; % alpha is a parameter used in the initial distribution f0
 stages=3; 
 deltax = (b-a)/(Nx-1);
 deltatini = nu*(deltax/max(v)); 
 NumTimesteps = ceil(Tfinal/deltatini);
 deltatfin = Tfinal/NumTimesteps; % I have this one in my notes but not sure
 %% Remove Brackets**********************************************
 x= a:deltax:b;
 
 
%% Initial Distribution
fo = exp(-alpha*(x-0.5).^2); % initial distribution

%% NEW TIME STEP
Fmatrix = zeros(Nv,Nx);
for k=1:Nv
    fold(k,:) = fo;
end
 disp('hello')
while NumTimesteps>0
    NumTimesteps=NumTimesteps-1;
    for k=1:Nv
    %[fnew]=bgkapproxtstep(fold,stages,Nx,deltax,v,deltatfin,A,w);
    [fnew]=bgkapproxtstepperiodic(fold(k,:),stages,Nx,deltax,v(k),deltatfin,A,w); 
    Fmatrix(k,:) = fnew;
    end
    disp('hello1')
    fold = Fmatrix;
end
  
 Fmatrix;
 %% Talk to stephanie about this Plotting routine**********************************************
 figure
 for k=1:Nv
     plot(x,fo,'g',x,Fmatrix(k,:),'r');
     hold on
 end
 
 figure()
 plot(x,fo,'g',x,Fmatrix)
%for k=1:Nx
%    e(k)=abs(fnew(k)-fo(k));
%end

%error=max(e)
 

end
 

  