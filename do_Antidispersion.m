close all;
clear all;

Coeff1=0; %% No tracer can escape if Coeff=1


Coeff1=1; %% No tracer can escape if Coeff=1


D=1e-11;% Diffusivity for taylor dispersion
D=1e-10;% Diffusivity for taylor dispersion


D1=1e-9;


%D1=5e-13;

%D=5e-8;% Diffusivity for taylor dispersion

Nz=500;


[Dzed,zed]  = Differentiation_Matrix(Nz-1);
l_PVS=5e-3; %5000 mum

    Z           = (zed+1)/2*1; %Z [10 9 ... 1,0];
    Dz          = 2/1*Dzed; 

%k=1e-10;
mu=7e-4;

%mu=1;
h=10e-6; % 10 um


k=1e-7; % From comsol simulation
%k=0;
%k=1e-10; % From comsol simulation
%k=0;
k=0;

B=sqrt(3*k*mu/h^3)*l_PVS;
%P0=500/factor; % A larger pressure


% u_0=1e-4;
% 
% u_0=1e-3;

%u_0=1e-10;
epsilon=h/l_PVS;


P0=10; % The actual pressure that drives such a flow velocity.



%u_0 =1.4286e-04; % The 
%u_0=1.5e-04; % The value to be used
u_0=P0/(2*mu/epsilon/h);
%P0=u_0*(2*mu/epsilon/h);



t_list=linspace(0.1,0.9,10)/epsilon*3/2;
%t_list=linspace(0.01,0.6,30)/epsilon*3/2;



P_left=1;
P_right=0;
%P_right=0;




P=P_left*(exp(B*(Z-1))+exp(-B*(Z-1)))...
    /(exp(-B*(1))-exp(B*(1)));




% dPdX=P0*(B*exp(B*(Z-l_PVS))+(-B)*exp(-B*(Z-l_PVS)))...
%     /(exp(-B*(l_PVS))-exp(B*(l_PVS)));
% d2PdX2=P0*(B^2*exp(B*(Z-l_PVS))+(B^2)*exp(-B*(Z-l_PVS)))...
%     /(exp(-B*(l_PVS))-exp(B*(l_PVS)));

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

if k==0

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
Pe_number_D1=u_0*h/D1;


for kk=1:numel(t_list)
    t=t_list(kk)
    
         a=6; ns=30; nd=30;      % implicit parameters
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

          C0=10*exp(-(Z-0.2).^2/(0.015/4));
         %+10*exp(-(zed-0.6).^2/(0.1*0.1));

       %  C0=10*exp(-(zed+0.7).^2/(0.1*0.1))+10*exp(-(zed-0.3).^2/(0.1*0.1));
          RHS=0*Z+C0;
      %    figure; plot(Z,C0)
    
    parfor ii=1:numel(s_all)
           
         sigma=s_all(ii);
    
          %  RHS(end)=10/sigma;
    
         C=0*Z;
        
         Matrix=zeros(Nz,Nz);
         Matrix1=zeros(Nz,Nz);
         phi=0.4;
         % Matrix=sigma*eye(Nz,Nz)+(u_bar.*(1-2*h^2/105/D*du_bardx)).*Dz-...
         %     (D*(1+2*h^2/105/D^2.*u_bar.^2)+...
         %     -5*11.*u_bar*h/60/phi^2*h.*(v_h/2/D-sqrt(v_h.^2+4*sigma*D)/2/D))...
         %     .*Dz^2;

         % The 0.1 models the tissue that only allows fluid out

  
            Matrix=sigma.*eye(Nz,Nz)+...
                epsilon*U_bar.*(1-epsilon*2/105*Pe_number.*DU_barDX-Coeff1*epsilon*Pe_number/15*V_h).*Dz+Coeff1*epsilon^2*3/35*V_h.*U_bar.*Pe_number.*Dz...
                +Coeff1*epsilon^2.*Pe_number.*(2/5*V_h.^2-1/15*U_bar.*dV_hdX).*eye(Nz,Nz)-...
             epsilon^2./Pe_number.*((1+2/105.*Pe_number.^2.*U_bar.^2))...
             .*Dz^2-Coeff1*epsilon*V_h.*eye(Nz,Nz);



            Matrix_D1=sigma.*eye(Nz,Nz)+...
                epsilon*U_bar.*(1-epsilon*2/105*Pe_number_D1.*DU_barDX-Coeff1*epsilon*Pe_number_D1/15*V_h).*Dz+Coeff1*epsilon^2*3/35*V_h.*U_bar.*Pe_number_D1.*Dz...
                +Coeff1*epsilon^2*Pe_number_D1.*(2/5*V_h.^2-1/15*U_bar.*dV_hdX).*eye(Nz,Nz)...
             -epsilon^2./Pe_number_D1.*((1+2/105.*Pe_number_D1.^2*U_bar.^2))...
             .*Dz^2-Coeff1*epsilon*V_h.*eye(Nz,Nz);

% Periodic BC    
        % Matrix(1,:)=Dz(1,:); %Matrix(1,1)=1;
         Matrix(1,:)=0; Matrix(1,1)=1; 
         Matrix(end,:)=0; Matrix(end,end)=1;
        % 
        % 
     Matrix_D1(1,:)=0; Matrix_D1(1,1)=1;
     Matrix_D1(end,:)=0; Matrix_D1(end,end)=1;

        % 
          %% The Neuman BC
           Matrix(1,:)=Dz(1,:); Matrix(end,:)=Dz(end,:);
           Matrix_D1(1,:)=Dz(1,:); Matrix_D1(end,:)=Dz(end,:);


           

         fs_VectorInSpace1(ii,:)=Matrix\RHS;
         fs_VectorInSpace1_D1(ii,:)=Matrix_D1\RHS;



         Y=linspace(h,1000*h,100)';
         %v_Y=1/mu*B^2*P0*exp(-B*Z)'.*(Y.^3/6-h*Y.^2/2);
         [ZZ,YY]=meshgrid(Z,Y);
         % fs_parenchyma(ii,:,:)=fs_VectorInSpace1(ii,:)./exp(h*(v_h/2/D-sqrt(v_h.^2+4*sigma*D)/2/D))'...
         %     .*ones(100,Nz).*exp(YY.*(v_h'/2/D-sqrt(v_h'.^2+4*sigma*D)/2/D))...
         %     ; % The superficial concentration


         % For ECS
           % fs_parenchyma(ii,:,:)=fs_VectorInSpace1(ii,:)...
           %   .*ones(100,Nz).*exp((YY-h).*(v_h'/2/D-sqrt(v_h'.^2+4*sigma*D)/2/D))...
           %   ; % The superficial concentration
           % 


    end


       btF1=bt'.*fs_VectorInSpace1;          % functional value F(s)
       conc1(kk,:) = sum(real(btF1),1);


       btF1_D1=bt'.*fs_VectorInSpace1_D1;          % functional value F(s)
       conc1_D1(kk,:) = sum(real(btF1_D1),1);


        % btF2=bt'.*fs_parenchyma;
        % conc_parenchyma(kk,:,:) = sum(real(btF2),1);




end

        disp(['Effective enhancement:' num2str(max(2/105.*Pe_number^2))]);
        disp(['Effective enhancement:' num2str(max(2/105.*Pe_number_D1^2))]);


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

        C_prime_D1(ll,:,:)=epsilon*(Pe_number_D1.*U_bar.*Dz*conc1_D1(ll,:)')'.*(YY1.^2/4-YY1.^4/8-7/120)...
        +V_h'.*(conc1_D1(ll,:)).*epsilon.*Pe_number_D1.*(YY1.^2/2-1/6);
end



%% 
figure;
subplot(2,1,1); 
plot(Z,C0,'--','LineWidth',2); hold on;
plot(Z,conc1_D1');
title('average concentration, weak dispersion');
subplot(2,1,2); 
plot(Z,C0,'--','LineWidth',2); hold on;
plot(Z,conc1');
title('average concentration, strong dispersion');

%for jj=1:numel(t_list)
legend('t=t0','t=t1','t=t2','t=t3','t=t4');


Z0=linspace(0,l_PVS,Nz)';

% 
phi=0.4; % the porosity
plott=1;
if plott==1
    figure; 

    for jj=1:numel(t_list)

    
    subplot(4,1,1); contourf(ZZ1,YY1,conc1_D1(jj,:).*ones(size(ZZ1))...
        +squeeze(C_prime_D1(jj,:,:)),'EdgeColor','none');
    colorbar;
    set(gca,'clim',[0,12]);
    title(['D=' num2str(D1*1e12) '\mum^2/s, C_{mean}+C']);

    subplot(4,1,2); contourf(ZZ1,YY1,squeeze(C_prime_D1(jj,:,:)),'EdgeColor','none');
    colorbar;
    title(['C']);
  %  set(gca,'clim',[-1,1]);
    
    subplot(4,1,3); contourf(ZZ1,YY1,conc1(jj,:).*ones(size(ZZ1))...
        +1*squeeze(C_prime(jj,:,:)),'EdgeColor','none');
    colorbar;
    set(gca,'clim',[0,12]);
    title(['D=' num2str(D*1e12) '\mum^2/s, C_{mean}+C']);

    subplot(4,1,4); contourf(ZZ1,YY1,squeeze(C_prime(jj,:,:)),'EdgeColor','none');
    colorbar;
    title(['C']);
   % set(gca,'clim',[-1,1]);

    pause(0.5);
    
    end
end
%end



% 
figure; plot(trapz(Z0,conc1_D1,2)/l_PVS);
hold on;
plot(trapz(Z0,conc1,2)/l_PVS)
% (trapz(Z0,conc1,2)/l_PVS)
% 

%C0=10*exp(-(zed+0.6).^2/(0.15*0.1));



%%%   Compare the theoretical solution in infinte space to the numerical
%%%   solution
Pe_number_D11=u_0*h/(D1*(1+2*(u_0*h/D1)^2/105));

T0=0.015/4/(4*epsilon^2/Pe_number_D11);


C00=10*exp(-(Z-0.2).^2/(0.015/4));
close all;


figure; 
%subplot(2,1,1);
%t_list=[0,1,2,3,4,5]*0.1;

for ii=1:1:numel(t_list)
    plot(Z,10*sqrt(T0/(t_list(ii)+T0)).*exp(-Pe_number_D11/epsilon^2*(Z-0.2-epsilon*U_bar(1)*t_list(ii)).^2/4/(t_list(ii)+T0))...
        );
    c_theoretical(ii,:)=10*sqrt(T0/(t_list(ii)+T0)).*exp(-Pe_number_D11/epsilon^2*(Z-0.2-epsilon*U_bar(1)*t_list(ii)).^2/4/(t_list(ii)+T0));
    hold on;
end
    % plot(Z,10*sqrt(T0/(0+T0)).*exp(-Pe_number_D11/epsilon*(Z-0.2-epsilon*U_bar(1)*0).^2/4/(0+T0))...
    %     );
        hold on;
%subplot(2,1,2);
plot(Z,conc1_D1,'*');



figure;
plot(Z,(conc1_D1-c_theoretical)./10);


%30*epsilon*U_bar(1)*t_list(ii)
if k~=0
    save("Traer_3rd_ReducedBolus.mat");
else
    save("Traer_3rd_IncreasedBolus.mat");
end
%save("Traer_3rd_IncreasedBolus.mat");

