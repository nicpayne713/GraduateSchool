function [fnew]=bgkapproxtstepperiodic(fold,stages,Nx,deltax,v,dt,A,w)

%% OLD TIME STEP
fold;

%% STAGES
% The stages are given by equation (22).
fstage=zeros(stages,Nx,Nv);
sigma=zeros(stages,Nx,Nv);
Flux=zeros(stages,Nx,Nv);

for i = 1:stages
    for l = i:stages
        if i == 1
            
        else
            % calculate flux
            for j = 2:Nx-1
                sigma(i,l,j) = (fstage(i,l,j+1)-fstage(i,l,j-1))/(2*deltax);
            end
            sigma(i,l,1)=(fstage(i,l,2)-fstage(i,l,Nx-1))/(2*deltax);
            sigma(i,l,Nx)= sigma(i,l,1);
            
            for j=2:Nx
                Flux(i,l,j) = max(v,0)*(fstage(i,l,j-1) +...
                    (deltax/2)*sigma(i,l,j-1))+ min(v,0)*(fstage(i,l,j) + (deltax/2)*sigma(i,l,j));
            end
            Flux(i,l,1) = max(v,0)*(fstage(i,l,Nx-1) + (deltax/2)*sigma(i,l,Nx-1))...
                + min(v,0)*(fstage(i,l,Nx) + (deltax/2)*sigma(i,l,Nx)); % Set Flux value at left endpoint
            Flux(i,l,Nx+1)=  max(v,0)*(fstage(i,l,Nx) + ...
                (deltax/2)*sigma(i,l,Nx))+ min(v,0)*(fstage(i,l,2) +...
                (deltax/2)*sigma(i,l,2)); % these lines also are part of the definition of the F.
            
            % moments
            for j = 1:Nx+1
                for k = 1:Nv
                    
                end
            end

%------STAGES 1---------------------------------------
i=1;

% Calculate stage 1. 
for j=1:Nx
    f1(j,k)=fold(j,k) + dt * (A(1,1)/tau1(j,K) * (;
end

%------STAGE 2---------------------------------------
% For this problem the collision term is zero so in the 2nd order method
% the 1st stage is equal to the previous time step (or the initial
% distribution if we are calculating the first time step) since A(2,1)=0.

i=2;
l=1;
  
 % Calculate sigma for stage 2. 
 for j = 2:Nx-1
        sigma(l,j) = (fstage(l,j+1)-fstage(l,j-1))/(2*deltax);
 end
%  Boundary condtions: sigma
% sigma(l,1)=0; % set sigma value at left endpoint
% sigma(l,Nx)=0; % set sigma value at right endpoint
 
%  % Boundary conditions: periodic
  sigma(l,1)=(fstage(l,2)-fstage(l,Nx-1))/(2*deltax);
  sigma(l,Nx)= sigma(l,1);
 
 
 % Calculate the Flux for stage 2. 
 for j=2:Nx
     Flux(l,j) = max(v,0)*(fstage(l,j-1) + (deltax/2)*sigma(l,j-1))+ min(v,0)*(fstage(l,j) + (deltax/2)*sigma(l,j)); 
 end
 Flux(l,1) = max(v,0)*(fstage(l,Nx-1) + (deltax/2)*sigma(l,Nx-1))+ min(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx)); % Set Flux value at left endpoint
 Flux(l,Nx+1)=  max(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx))+ min(v,0)*(fstage(l,2) + (deltax/2)*sigma(l,2)); % these lines also are part of the definition of the F.
 
 %Calculate stage 2. (Equal to the initial distribution in the 2nd order
 % method.)
 for j=1:Nx
     fstage(i,j)=fold(j)-(dt/deltax)*A(i,l)*(Flux(l,j+1)-Flux(l,j)); 
 end
 
 
 %------STAGE 3---------------------------------------
 i=3;
 l=2;
 
 
 % Calculate sigma for stage 3. 
 for j = 2:Nx-1
        sigma(l,j) = (fstage(l,j+1)-fstage(l,j-1))/(2*deltax);
 end
% Boundary condtions: sigma
% sigma(l,1)=0; % set sigma value at left endpoint
% sigma(l,Nx)=0; % set sigma value at right endpoint
 
% Boundary conditions: periodic
  sigma(l,1)=(fstage(l,2)-fstage(l,Nx-1))/(2*deltax);
  sigma(l,Nx)= sigma(l,1);
 
 
 % Calculate Flux for stage 3.

 
 for j=2:Nx
     Flux(l,j) = max(v,0)*(fstage(l,j-1) + (deltax/2)*sigma(l,j-1))+ min(v,0)*(fstage(l,j) + (deltax/2)*sigma(l,j)); 
 end
 Flux(l,1) = max(v,0)*(fstage(l,Nx-1) + (deltax/2)*sigma(l,Nx-1))+ min(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx)); % Set Flux value at left endpoint
 Flux(l,Nx+1)=  max(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx))+ min(v,0)*(fstage(l,2) + (deltax/2)*sigma(l,2)); % these lines also are part of the definition of the F.
 %Flux(l,1) = max(v,0)*(fstage(l,Nx-1) + (deltax/2)*sigma(l,Nx-1)); % Set Flux value at left endpoint
 %Flux(l,Nx+1)=  max(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx)); % these lines also are part of the definition of the F.
 
 
 %Calculate stage 3.
 for j=1:Nx
     fstage(i,j)=fold(j)-(dt/deltax)*A(i,l)*(Flux(l,j+1)-Flux(l,j)); 
 end
 
 
 
 
 %------New Time Step---------------------------------------

 l=3;
  
 % We need to calculate F^(3) before we can calculate the next time step.
 
 
 % Calculate sigma for l=3. 
 for j = 2:Nx-1
        sigma(l,j) = (fstage(l,j+1)-fstage(l,j-1))/(2*deltax);
 end
%  Boundary condtions: sigma
% sigma(l,1)=0; % set sigma value at left endpoint
%sigma(l,Nx)=0; % set sigma value at right endpoint
 
%  % Boundary conditions: periodic
  sigma(l,1)=(fstage(l,2)-fstage(l,Nx-1))/(2*deltax);
  sigma(l,Nx)= sigma(l,1);
 
 
 % Calculate Flux for l=3.
 for j=2:Nx
     Flux(l,j) = max(v,0)*(fstage(l,j-1) + (deltax/2)*sigma(l,j-1))+ min(v,0)*(fstage(l,j) + (deltax/2)*sigma(l,j)); 
 end
 Flux(l,1) = max(v,0)*(fstage(l,Nx-1) + (deltax/2)*sigma(l,Nx-1))+ min(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx)); % Set Flux value at left endpoint
 Flux(l,Nx+1)=  max(v,0)*(fstage(l,Nx) + (deltax/2)*sigma(l,Nx))+ min(v,0)*(fstage(l,2) + (deltax/2)*sigma(l,2)); % these lines also are part of the definition of the F.
 
 
% NEW TIME STEP
 
 for j=1:Nx
     fnew(j) = fold(j)-(dt/deltax)*(w(2)*(Flux(2,j+1)-Flux(2,j))+w(3)*(Flux(3,j+1)-Flux(3,j)));
 end