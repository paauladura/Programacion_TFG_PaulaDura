%% EXAMPLE for Slepian functions within a spherical cap             - 2021-10-20
clear all; clc

% checking available pathes
if exist('gruenbaum_evp','file') ~= 2,error('Please add path of ''SHBundle'''),   end; 
if exist('circ4sph','file') ~= 2,     error('Please add path of ''visBundle'''),  end; 
if exist('checkcoor','file') ~= 2,    error('Please add path of ''uberall'''),    end; 


%% Gruenbaum's eigenvalue problem
Theta0 = 30; lmax = 60;
[GB_vectors, GB_values, GB_orders, GB_INFO__]= gruenbaum_evp(Theta0*pi/180,lmax);

%% selection of Slepian functions
[MAX_chi,ind5]  = max(GB_values);       % worst concentration
[ZERO_chi,ind4]  = min(abs(GB_values)); % eigen value close to zero   
if ZERO_chi < 0, ind4 = ind+1; end
ind3 = min([fix(ind4/3+1),60]);         % a well concentrated basis function
ind2 = min([fix(ind4/6+1),15]);         % a well concentrated basis function
select = [1,ind2,ind3,ind4,ind5];       % index of visualized functions

%% zero-padding of the coefficients + rotation of center
alpha  = 30; beta   = 50; gamma  = 55;
[slepian3D,SL_INFO__] = getSlepianSHC(GB_vectors, GB_values, GB_orders,GB_INFO__,select,alpha*pi/180,beta*pi/180,gamma*pi/180);


%% spherical harmonic synthesis per Slepian basis function
theta  = linspace(0,pi,100);
lambda = linspace(-pi,pi,200);
F = gshs_grid(slepian3D, lambda, pi/2-theta, 6371e3,'quant','none','GM',1,'max_lm',lmax,'sub_WGS84',false);

%% visualization of selected Slepian basis functions 
% boundaries of circles on the sphere
[lambound,phibound] = circ4sph(alpha, 90-beta,Theta0);
thebound = 90-phibound;

%% visualization
figure 
colormap(colbrew(1,129))
lambdaG =  lambda*180/pi;thetaG = theta*180/pi;
for ii = 1: numel(select)
    f = F(:,:,ii);
    % ... subfigures
    subplot(2,3,ii);
    maxf = max(max(abs(f)));
    imagesc(lambdaG,thetaG,f);hold on; plot(alpha,beta,'ko');
    plot(lambound,thebound,'k','LineWidth',2);
    colorbar('Location', 'SouthOutside'); caxis(maxf*[-1,1])
    set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[0:30:180]); 
    axis equal; axis tight; 
    title(['Slepian function with eigen value \chi = ' num2str(SL_INFO__.eigenvalue(ii),'%6.2f') ])
end


% Gruenbaum eigenvalues
subplot(236)
plot(1:numel(GB_values), GB_values,'m','LineWidth',2); hold on; plot([1,numel(GB_values)],[0,0],'k-.')
plot(select,GB_values(select),'mo','MarkerSize',12)
title('Eigenvalues of Gruenbaum''s problem')
