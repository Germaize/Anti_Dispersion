close all;
clear all;


%addpath(genpath('/Users/yiminggan/Desktop/Paper_AntiDispersion/'));

Coeff1=1; %% No tracer can escape if Coeff=1


D=1e-11;% Diffusivity for taylor dispersion
D=1e-10;% Diffusivity for anti dispersion


Nz=500;


%% Define the differential matrixs
[Dzed,zed]  = Differentiation_Matrix(Nz-1);
l_PVS=5e-3; %5000 mum
    Z           = (zed+1)/2*1; %Z [10 9 ... 1,0];
    Dz          = 2/1*Dzed; 


%% Parameters
mu=7e-4; % Dynamic viscosity
h=10e-6; % 10 um, width of the channel


k=1e-7; % permeability of the leaky wall
%k=0;

B=sqrt(3*k*mu/h^3)*l_PVS; 

epsilon=h/l_PVS; % aspect ratio


P0=10; % The actual pressure that drives such a flow velocity.



u_0=P0/(2*mu/epsilon/h); % flow velocity
t_list=linspace(0.1,0.9,10)/epsilon*3/2; % simulate 10 time points



P_left=1;
P_right=0;
%P_right=0;

Pe=0;

P_term1=((P_left-Pe)*exp(-1*B)-(P_right-Pe))*exp(B*Z)...
    /(exp(-B*(1))-exp(B*(1)));
P_term2=((P_right-Pe)-(P_left-Pe)*exp(1*B))*exp(-B*Z)...
    /(exp(-B*(1))-exp(B*(1)));
P=Pe+P_term1+P_term2;

dPdX=B*P_term1+(-B)*P_term2;

d2PdX2=B^2*(P-Pe);

d3PdX3=B^3*P_term1+(-B^3)*P_term2;

U_bar=-2/3*dPdX;

V_h=-2*d2PdX2*(1/6-1/2);
dV_hdX=-2*d3PdX3*(1/6-1/2);

DU_barDX=-d2PdX2;


if k==0 % for nonpermeable wall simulation

    P=P_left*(1-Z);
    dPdX=-1;
    d2PdX2=0;
    d3PdX3=0;
    U_bar=-2/3*dPdX;
    
    V_h=-2*d2PdX2*(1/6-1/2);
    dV_hdX=-2*d3PdX3*(1/6-1/2);
    
    DU_barDX=-d2PdX2;


end


%% Plot the flow profile
figure; plot(Z,P); %xlim([0,0.01])
xlabel('X (m)');
ylabel('P (Pa)');

%figure; plot(Z,U_bar);
%figure; plot(Z,du_bardx);
%figure; plot(Z,v_h);
figure; plot(Z*l_PVS,U_bar*u_0); xlabel('x (m)'); ylabel('axial velocity (m/s)');
yyaxis right; plot(Z*l_PVS,V_h*u_0*epsilon); ylabel('radial velocity (m/s)');

disp('Pe number is:');
Pe_number=u_0*h/D;

for kk=1:numel(t_list)
    t=t_list(kk)
    
         a=6; ns=30; nd=30;      % implicit parameters for laplace transform
        %radt=linspace(tini,tend,nnt); % time vector
    
        for n=1:ns+1+nd               % prepare necessary coefficients
           alfa(n)=a+(n-1)*pi*1i;
           beta(n)=-exp(a)*(-1)^n;
        end
        n=1:nd;
        bdif=fliplr(cumsum(gamma(nd+1)./gamma(nd+2-n)./gamma(n)))./2^nd;
        beta(ns+2:ns+1+nd)=beta(ns+2:ns+1+nd).*bdif;
        beta(1)=beta(1)/2;
    
           for ii=1:numel(t)
               tt=t(ii);
               s(ii,:)=alfa/tt;                 % complex frequency s
               bt(ii,:)=beta/tt;
           end
    
           s_all=s.';
           s_all_vec=s_all(:).';
    
    
         % fs_VectorInSpace1=
         % A bolous with a width of 0.01

          %C0=10*exp(-(zed+0.6).^2/(0.015));

          C0=1*exp(-(Z-0.2).^2/(0.015/4));
         %+10*exp(-(zed-0.6).^2/(0.1*0.1));

       %  C0=10*exp(-(zed+0.7).^2/(0.1*0.1))+10*exp(-(zed-0.3).^2/(0.1*0.1));
          RHS=0*Z+C0;
      %    figure; plot(Z,C0)
    
    parfor ii=1:numel(s_all) % solve all time points parallelly
           
         sigma=s_all(ii);
    
          %  RHS(end)=10/sigma;
    
         C=0*Z;
        
         Matrix=zeros(Nz,Nz);


  
            Matrix=sigma.*eye(Nz,Nz)+...
                epsilon*U_bar.*(1-epsilon*2/105*Pe_number.*DU_barDX-Coeff1*epsilon*Pe_number/15*V_h).*Dz+Coeff1*epsilon^2*3/35*V_h.*U_bar.*Pe_number.*Dz...
                +Coeff1*epsilon^2.*Pe_number.*(2/5*V_h.^2-1/15*U_bar.*dV_hdX).*eye(Nz,Nz)-...
             epsilon^2./Pe_number.*((1+2/105.*Pe_number.^2.*U_bar.^2))...
             .*Dz^2-Coeff1*epsilon*V_h.*eye(Nz,Nz);


% Periodic BC    
        % Matrix(1,:)=Dz(1,:); %Matrix(1,1)=1;
         Matrix(1,:)=0; Matrix(1,1)=1; 
         Matrix(end,:)=0; Matrix(end,end)=1;
        % 

        % 
          %% The Neuman BC
         Matrix(1,:)=Dz(1,:); Matrix(end,:)=Dz(end,:);

         fs_VectorInSpace1(ii,:)=Matrix\RHS;



         Y=linspace(h,1000*h,100)';

    end


       btF1=bt'.*fs_VectorInSpace1;          % functional value F(s)
       conc1(kk,:) = sum(real(btF1),1);



end

        disp(['Effective enhancement:' num2str(max(2/105.*Pe_number^2))]);


 Y=linspace(h,1000*h,100)';

figure; 
[ZZ,YY]=meshgrid(Z,Y);

%% Reconstruct the pertubation term

clear C_prime
    [ZZ,YY]=meshgrid(Z,Y);
    [ZZ1,YY1]=meshgrid(Z,linspace(-1,1,100));

for ll=1:size(conc1,1)
    C_prime(ll,:,:)=epsilon*(Pe_number.*U_bar.*Dz*conc1(ll,:)')'.*(YY1.^2/4-YY1.^4/8-7/120)...
        +V_h'.*(conc1(ll,:)).*epsilon.*Pe_number.*(YY1.^2/2-1/6);

end



%% 


figure;
plot(Z,C0,'--','LineWidth',2); hold on;
plot(Z,conc1');
title('average concentration, strong dispersion');

%for jj=1:numel(t_list)
legend('t=t0','t=t1','t=t2','t=t3','t=t4');


Z0=linspace(0,l_PVS,Nz)';


plott=1;
if plott==1
    figure; 

    for jj=1:numel(t_list)

    
    subplot(2,1,1); contourf(ZZ1,YY1,conc1(jj,:).*ones(size(ZZ1))...
        +1*squeeze(C_prime(jj,:,:)),'EdgeColor','none');
    colorbar;
    set(gca,'clim',[0,1.2]);
    title(['$\bar{C}+C_{prime}$'],'Interpreter','latex');

    subplot(2,1,2); contourf(ZZ1,YY1,squeeze(C_prime(jj,:,:)),'EdgeColor','none');
    colorbar;
    title(['$C_{prime}$'],'Interpreter','latex');
   % set(gca,'clim',[-1,1]);

    pause(0.5);
    
    end
end
%end



%30*epsilon*U_bar(1)*t_list(ii)
if k~=0
    save("Traer_3rd_ReducedBolus.mat");
else
    save("Traer_3rd_IncreasedBolus.mat");
end
%save("Traer_3rd_IncreasedBolus.mat");

