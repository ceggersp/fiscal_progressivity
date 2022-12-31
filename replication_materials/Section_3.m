clear all; clc; close all;
tic

I= 1000;                % Assets (a) grid size. As Figure 1 represents a derivative on the a dimension, we use a finer grid.
J = 400;                % z grid size.

ga = 2;                 % CRRA utility with parameter gamma
alpha = 0.35;           % Production function F = K^alpha * L^(1-alpha)
delta = log(1.07);      % Capital depreciation
sig2vec = [0.10, 0.12]; % Diffusion parameters for z
MPC = zeros(I,J,2);     % MPCs will be stored here
PSav = zeros(I,J,2);    % Precautionary savings will be stored here
L_ave = 0.9404;

for sigit = 1:length(sig2vec)
    sig2 = sig2vec(sigit);
    rho = log(1.057);   %(Kaplan, Moll & Violante, 2017) discount rate
    eta = sig2vec(1)/2; % The process with sig2 = 0.1 will have a stationary variance of 1.
    Var = sig2/(2*eta);
    sdev = sqrt(Var);
    
    w = 1.4675;         % Equilibrium wage
    amin = -w/4;
    amax = 45*w;
    zmean = -Var/2;
    zmin = zmean - 2.5*sdev;
    zmax = zmean + 2.5*sdev;
    a = linspace(amin,amax,I)';
    z = linspace(zmin,zmax,J);   % productivity vector
    da = (amax-amin)/(I-1);
    
    dz = (zmax-zmin)/(J-1);
    dz2 = dz^2;
    aa = a*ones(1,J);
    zz = ones(I,1)*z;
    expzz = exp(zz);
    pdfzz = normpdf((zz-zmean*ones(I,J))/sdev);
    pdfzz = pdfzz./sum(dz*pdfzz(1,:));
    cdfzz = normcdf((zz-zmean*ones(I,J))/sdev);
    
    mu = eta*(zmean - z);
    s2 = sig2.*ones(1,J);
    
    maxit= 100;
    crit = 10^(-6);
    Delta = 1000;
    
    %Finite difference approximation of the partial derivatives
    Vaf = zeros(I,J);
    Vab = zeros(I,J);
    Vzf = zeros(I,J);
    Vzb = zeros(I,J);
    Vzz = zeros(I,J);
    c = zeros(I,J);
    
    L = L_ave;
    %INITIAL GUESS
    K = 9;
    %r = ((1-tau)*alpha*Aprod*K^(alpha-1))*L^(1-alpha) - delta; %interest rates %0.011693820316052
    w = (1-alpha)*(K^alpha)*L^(-alpha);
    wexpzz = w.*expzz;
    r = alpha*(K^(alpha-1))*(L^(1-alpha)) - delta;
    
    %CONSTRUCT MATRIX Aswitch SUMMARIZING EVOLUTION OF z
    % yy = - s2/dz2 - mu/dz;
    % chi =  s2/(2*dz2);
    % zeta = mu/dz + s2/(2*dz2);
    chi =  - min(mu,0)/dz + s2/(2*dz2);
    yy =  min(mu,0)/dz - max(mu,0)/dz - s2/dz2;
    zeta = max(mu,0)/dz + s2/(2*dz2);
    
    Ir = 35;
    crit_S = 0.5*10^(-5);
    v0 = ((L/L_ave)*wexpzz + r.*aa).^(1-ga)/(1-ga)/rho;
    r0 = 0.03;
    rmin = 0.01*rho;
    rmax = 0.99*rho;
    %%%%%%%%%%%%%%%%
    % STEADY STATE %
    %%%%%%%%%%%%%%%%
    for ir=1:Ir
        
        r_r(ir)=r;
        rmin_r(ir)=rmin;
        rmax_r(ir)=rmax;
        
        KD(ir) = (alpha/(r + delta))^(1/(1-alpha))*L;
        w = (1-alpha)*(KD(ir)^alpha)*L^(-alpha);
        wexpzz = w.*expzz;
              
        %This will be the upperdiagonal of the matrix Aswitch
        updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
        for j=1:J
            updiag=[updiag;repmat(zeta(j),I,1)];
        end
        %This will be the center diagonal of the matrix Aswitch
        centdiag=repmat(chi(1)+yy(1),I,1);
        for j=2:J-1
            centdiag=[centdiag;repmat(yy(j),I,1)];
        end
        centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];
        %This will be the lower diagonal of the matrix Aswitch
        lowdiag=repmat(chi(2),I,1);
        for j=3:J
            lowdiag=[lowdiag;repmat(chi(j),I,1)];
        end
        %Add up the upper, center, and lower diagonal into a sparse matrix
        Aswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);
        
        if ir>1
            v0 = V_r(:,:,ir-1);
        end
        
        v = v0;
        
        for n=1:maxit
            V = v;
            V_n(:,:,n)=V;
            % forward difference
            Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
            Vaf(I,:) = ((L/L_ave)*w*exp(z) + r.*amax).^(-ga);
            % backward difference
            Vab(2:I,:) = (V(2:I,:,:)-V(1:I-1,:,:))/da;
            Vab(1,:) = ((L/L_ave)*w.*exp(z) + r.*amin).^(-ga);  %state constraint boundary condition
            I_concave = Vab > Vaf;
            %consumption and savings with forward difference
            cf = Vaf.^(-1/ga);
            sf = (L/L_ave)*wexpzz + r.*aa - cf;
            %consumption and savings with backward difference
            cb = Vab.^(-1/ga);
            sb = (L/L_ave)*wexpzz + r.*aa - cb;
            %consumption and derivative of value function at steady state
            c0 = (L/L_ave)*wexpzz + r.*aa;
            Va0 = c0.^(-ga);
            
            % dV_upwind makes a choice of forward or backward differences based on
            % the sign of the drift
            If = sf > 0; %positive drift --> forward difference
            Ib = sb < 0; %negative drift --> backward difference
            I0 = (1-If-Ib); %at steady state
            %make sure backward difference is used at amax
            %Ib(I,:) = 1; If(I,:) = 0;
            %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
            %already taken care of automatically
            
            Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term
            c = Va_Upwind.^(-1/ga);
            u = c.^(1-ga)/(1-ga);
            
            %CONSTRUCT MATRIX A
            X0 = - min(sb,0)/da;
            Y0 = - max(sf,0)/da + min(sb,0)/da;
            Z0 = max(sf,0)/da;
            updiag0=[0]; %This is needed because of the peculiarity of spdiags.
            for j=1:J
                updiag0=[updiag0;Z0(1:I-1,j);0];
            end
            centdiag0=reshape(Y0,I*J,1);
            lowdiag0=X0(2:I,1);
            for j=2:J
                lowdiag0=[lowdiag0;0;X0(2:I,j)];
            end
            AA0=spdiags(centdiag0,0,I*J,I*J)+spdiags([updiag0;0],1,I*J,I*J)+spdiags([lowdiag0;0],-1,I*J,I*J);
            
            AA = AA0;
            A = AA + Aswitch;
            
            if max(abs(sum(A,2)))>10^(-9)
                disp('Improper Transition Matrix')
                break
            end
            
            B = (1/Delta + rho)*speye(I*J) - A;
            
            u_stacked = reshape(u,I*J,1);
            V_stacked = reshape(V,I*J,1);
            
            b = u_stacked + V_stacked/Delta;
            
            V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
            
            V = reshape(V_stacked,I,J);
            
            Vchange = V - v;
            v = V;
            
            dist(n) = max(max(abs(Vchange)));
            if dist(n)<crit
                %disp('Value Function Converged, Iteration = ')
                %disp(n)
                break
            end
        end
        toc;
        
        % FOKKER-PLANCK EQUATION %
        AT = A';
        b = zeros(I*J,1);
        
        %need to fix one value, otherwise matrix is singular
        i_fix = 1;
        b(i_fix)=.1;
        row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
        AT(i_fix,:) = row;
        
        %Solve linear system
        gg = AT\b;
        g_sum = gg'*ones(I*J,1)*da*dz;
        gg = gg./g_sum;
        
        g = zeros(I,J);
        g = reshape(gg,I,J);
        
        g_r(:,:,ir) = g;
        sav = (L/L_ave)*wexpzz + r.*aa - c;
        adot(:,:,ir) = sav;
        inc(:,:,ir) = (L/L_ave)*wexpzz + r.*aa;
        V_r(:,:,ir) = V;
        
        KS(ir) = sum(sum(g.*aa))*da*dz;
        S(ir) = KS(ir) - KD(ir);
        
        %UPDATE INTEREST RATE
        disp(ir);
        S_K(ir) = S(ir)/abs(KD(ir));
        disp(S_K(ir));
        if S_K(ir)<-crit_S
            disp('Excess Demand');
            rmin = r;
            r = 0.5*(r+rmax);
        elseif S_K(ir)>crit_S;
            disp('Excess Supply');
            rmax = r;
            r = 0.5*(r+rmin);
        elseif abs(S_K(ir))<crit_S;
            display('Equilibrium Found, Unemployment rate =')
            disp(100-100*L/L_ave);
            disp(KD(ir));
            disp(S_K(ir));
            break
        end
        
    end
    
    sum(g(1,:))*da*dz
    
    Izf = repmat(mu,I,1) > 0;
    Izb = repmat(mu,I,1) < 0;
    Iz0 = repmat(mu,I,1) == 0;
    
    G = zeros(I,J);
    Gaf = G;
    Gab = G;
    Ga0 = G;
    Gzf = G;
    Gzb = G;
    Gz0 = G;
    Gzz = G;
    
    steps = 2000;
    T = 1;
    dt = T/steps;
    
    MU = repmat(mu,I,1);
    
    for t = 1:steps
        
        if floor(t/100) == t/20
            display('MPC iteration progress:')
            strcat(num2str(1000*t/steps),'%')/10
        else
        end

        Gaf(1:I-1,:) = (G(2:I,:)-G(1:I-1,:))/da;
        Gaf(I,:) = Gaf(I-1,:);
        Gab(2:I,:) = (G(2:I,:)-G(1:I-1,:))/da;
        Gab(1,:) = Gab(2,:);
        Ga0(2:I-1,:) = (G(3:I,:)-G(1:I-2,:))/(2*da);
        Ga0(I,:) = Ga0(I-1,:);
        Ga0(1,:) = Ga0(2,:);
        Ga_Upwind = If.*Gaf + Ib.*Gab + I0.*Ga0;
        
        Gzf(1:I-1,:) = (G(2:I,:)-G(1:I-1,:))/dz;
        Gzf(I,:) = Gzf(I-1,:);
        Gzb(2:I,:) = (G(2:I,:)-G(1:I-1,:))/dz;
        Gzb(1,:) = Gzb(2,:);
        Gz0(2:I-1,:) = (G(3:I,:)-G(1:I-2,:))/(2*dz);
        Gz0(I,:) = Gz0(I-1,:);
        Gz0(1,:) = Gz0(2,:);
        Gz_Upwind = Izf.*Gzf + Izb.*Gzb + Iz0.*Gz0;
        
        Gzz(2:I-1,:) = (G(3:I,:)-2*G(2:I-1,:)+G(1:I-2,:))/(dz^2);
        Gzz(I,:) = Gzz(I-1,:);
        Gzz(1,:) = Gzz(2,:);
        
        G = G + dt*(c + sav.*Ga_Upwind + MU.*Gz_Upwind + (sig2/2)*Gzz);
        
    end
    
    MPC(1:I-1,:,sigit) = (G(2:I,:) - G(1:I-1,:))/da;
    MPC(I,:,sigit) = MPC(I-1,:,sigit);
 
    pval = aa + repmat(PVal(z,zmean,w/L_ave,eta,r,1000,0.02),I,1);
    cdet = (r-(r-rho)/ga)*pval;
    PSav(:,:,sigit) = (cdet - c)./pval;
    
end

%% Figure (1)
plot(100*da*dz*cumsum(sum(g')),MPC(:,0.25*J,1),'Color',[0.9 0.0 0.4]); hold on
plot(100*da*dz*cumsum(sum(g')),MPC(:,0.75*J,1),'Color',[0.9 0.0 0.4],'LineStyle','--'); hold on
plot(100*da*dz*cumsum(sum(g')),MPC(:,0.25*J,2),'Color',[0.9 0.6 0.4]); hold on
plot(100*da*dz*cumsum(sum(g')),MPC(:,0.75*J,2),'Color',[0.9 0.6 0.4],'LineStyle','--'); hold off
hl = legend('Percentile 10 of $z$ $(\sigma^2=0.10)$', 'Percentile 90 of $z$ $(\sigma^2=0.10)$', 'Percentile 10 of $z$ $(\sigma^2=0.12)$', 'Percentile 90 of $z$ $(\sigma^2=0.12)$')
set(hl, 'FontSize',11, 'Interpreter','latex','EdgeColor',[1 1 1])
xlabel('Wealth percentile','FontSize',14,'interpreter','latex')
ylabel('MPC','FontSize',14,'interpreter','latex')
xlim([0 100])
ylim([0 0.19])
saveas(gcf,'Graphs/Figure_1.eps')
saveas(gcf,'Graphs/Figure_1.png')

%% Figure (2)
plot(da*dz*cumsum(sum(g')'),PSav(:,0.25*J,1),'Color',[0.9 0.0 0.4],'LineWidth',1); hold on
plot(da*dz*cumsum(sum(g')'),PSav(:,0.50*J,1),'Color',[0.9 0.4 0.4],'LineWidth',1); hold on
plot(da*dz*cumsum(sum(g')'),PSav(:,0.75*J,1),'Color',[0.9 0.8 0.4],'LineWidth',1); hold on
h2 = legend('Percentile 10 of $z$', 'Percentile 50 of $z$', 'Percentile 90 of $z$')
set(h2, 'FontSize',11, 'Interpreter','latex','EdgeColor',[1 1 1])
ylabel('PS as fraction of $a+E(PV)$','FontSize',14,'interpreter','latex')
xlabel('Wealth percentile','FontSize',14,'interpreter','latex')
xlim([0 1]); hold off
saveas(gcf,'Graphs/Figure_2.eps')
saveas(gcf,'Graphs/Figure_2.png')

%% Figure (3)
plot(da*dz*cumsum(sum(g')'),PSav(:,0.50*J,1),'Color',[0.9 0.0 0.4],'LineWidth',1); hold on
plot(da*dz*cumsum(sum(g')'),PSav(:,0.50*J,2),'Color',[0.9 0.5 0.4],'LineWidth',1); hold on
h2 = legend('$\sigma^2=0.10$', '$\sigma^2=0.12$')
set(h2, 'FontSize',11, 'Interpreter','latex','EdgeColor',[1 1 1])
ylabel('PS as fraction of $a+E(PV)$','FontSize',14,'interpreter','latex')
xlabel('Wealth percentile','FontSize',14,'interpreter','latex')
xlim([0 1]); hold off
saveas(gcf,'Graphs/Figure_3.eps')
saveas(gcf,'Graphs/Figure_3.png')

%%
function [PV,disc,zcorr] = PVal(z,zmean,w,eta,r,T,dt)
    J = length(z);
    zdemeaned = z - zmean;
    N = T/dt;
    for t = 2:N
        if floor(t*dt) == t*dt
            display('PVal iteration progress:')
            strcat(num2str(round(1000*t/N)/10),'%')
        else
        end
        zdemeaned(t,:) = exp(-eta*dt*(t-1))*zdemeaned(1,:); 
    end
    disc = exp(-r*dt*[0:N-1]');
    PV = zeros(N,J);
    znew = zdemeaned + zmean;
    znew = w*exp(znew);
    for j = 1:J
        PV(:,j) = disc.*znew(:,j);
    end
    PV = dt*sum(PV); 
end