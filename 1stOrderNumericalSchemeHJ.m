clc
clear
% First order ENO Scheme using LLLF Hamiltonian in One-Dimension
% Test equation for convex case
    H = @(u) ((u+1).^2) ./ 2; % Burders' Equation with alpha = 1
    
% For exact solution, implement numerical scheme on super-fine mesh, then
% compare to sheme implemented on a coarse mesh

dx = 2*1e-2; dxx = 2*1e-1;
x = -1:dx:1;
x = x(1:length(x)-1); % fine grid for "exact" solution, not including x = 1
xx = -1:dxx:1;
xx = xx(1:length(xx)-1); % coarse grid for approximation
deltat = .1/(pi^2);
tfinal = 2/(pi^2);
t = 0: deltat :tfinal;


Phiexact = zeros(length(t),length(x)); % matrix for exact
Phihat = zeros(length(t),length(xx)); % matrix for approximation
% initial condition
    Phiexact(1,:) = (-1)*cos(pi.*x);
    Phihat(1,:) = (-1)*cos(pi.*xx);
   
tic
    
% first order scheme
    for k = 2:length(t) % tme step
        % left boundary
            deltaxphi1 = Phiexact(k-1,2)-Phiexact(k-1,1);
            u1 = deltaxphi1/dx;

            %U1 = (1/dx)*(phi1 - phi0);
            Phiexact(k,1) = Phiexact(k-1,1) + deltat*LLFHamil1D(u1,u1);

            hu1 = (Phihat(k-1,2)-Phihat(k-1,1))/dxx;
            %hU1 = (1/dxx)*(hphi1 - hphi0);
            Phihat(k,1) = Phihat(k-1,1) + deltat*LLFHamil1D(hu1,hu1);            
        
        % middle spacial points
        for i = 2:length(x)-1 % spatial step
            %phim1 = Phiexact(k-1,i-1); % 
            deltaxphiiplus = Phiexact(k-1,i+1) - Phiexact(k-1,i);
            %phi0 = Phiexact(k-1,i); % 
            deltaxphiiminus = -(Phiexact(k-1,i-1)-Phiexact(k-1,i));
            %phi1 = Phiexact(k-1,i+1);
            %U1 = (1/dx)*(phi1 - phi0);
            %U2 = (1/dx)*(phi0 - phim1);
            u1 = deltaxphiiplus/dx;
            u2 = deltaxphiiminus/dx;
            Phiexact(k,i) = Phiexact(k-1,i) + deltat*LLFHamil1D(u1,u2);
            if i < length(xx) 
                %hphim1 = Phihat(k-1,i-1);
                %hphi0 = Phihat(k-1,i);
                %hphi1 = Phihat(k-1,i+1);
                %hU1 = (1/dxx)*(hphi1 - hphi0);
                %hU2 = (1/dxx)*(hphi0 - hphim1);
                hdeltaxphiiplus = Phihat(k-1,i+1) - Phihat(k-1,i);
                hdeltaxphiiminus = -(Phihat(k-1,i-1)-Phihat(k-1,i));
                hu1 = deltaxphiiplus/dxx;
                hu2 = deltaxphiiminus/dxx;
                Phihat(k,i) = Phihat(k-1,i) + deltat*LLFHamil1D(hu1,hu2);
            end
        end
            
            
        % right boundary
            deltaxphiN = -(Phiexact(k-1,length(x)-1)-Phiexact(k-1,length(x)));
            u1 = deltaxphiN/dx;

            %U1 = (1/dx)*(phi1 - phi0);
            Phiexact(k,1) = Phiexact(k-1,1) + deltat*LLFHamil1D(u1,u1);

            hu1 = -(Phihat(k-1,length(xx)-1)-Phihat(k-1,length(xx)))/dxx;
            %hU1 = (1/dxx)*(hphi1 - hphi0);
            Phihat(k,length(xx)) = Phihat(k-1,length(xx)) + deltat*LLFHamil1D(hu1,hu1);
        
    end
    
toc

% 1 norm at t = .5/(pi^2), which is where the solution is smooth
[id] = find(t==.5/(pi^2));
Lin1 = norm(Phiexact(id,:),1) - norm(Phihat(id,:),1);

% 1 norm at t = 1.5/(pi^2), where the solution has a discontinuous
%derivative
[id2] = find(t==1.5/(pi^2));
Lin1d = norm(Phiexact(id2,:),1) - norm(Phihat(id2,:),1);

disp(Lin1)
disp(Lin1d)
figure()
plot(x, Phiexact(end,:),'r')
figure(2)
plot(xx, Phihat(end,:),'c')