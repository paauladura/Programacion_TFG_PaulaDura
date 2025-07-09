%% EXAMPLE for the use of WIGNER_ALL &  ROTATE_SHC - 2021-04-12
% (rotation of spherical harmonics by Wigner-d-functions)
clear all; clc


% checking available path
if exist('wigner_all','file') ~= 2,   error('Please add path of ''SHBundle'''),   end; 
if exist('colbrew','file') ~= 2,      error('Please add path of ''visBundle'''),  end; 
if exist('checkcoor','file') ~= 2,    error('Please add path of ''uberall'''),    end; 

% spherical harmonic function to be rotated
lmax   = 15;
m      = 11;


% rotation angles in radian
alpha  = 1.2;
beta   = 2.3;
gamma  = 3.4;


%% rotation of function by Wigner-d-functions
% all function up to degree 'lmax'
DLMK      = wigner_all(lmax,beta,'dlmk');
% ... selection of degree and order
dlmk      = DLMK{lmax+1};
dlmk      = dlmk(lmax+1+m,:);
dlmk_plus = dlmk(lmax+1:end);
dlmk_minus= [0, fliplr(dlmk(1:lmax))];

% Legendre functions 
thetaG = 0:1:180; theta  = thetaG*pi/180;
lambdaG = -180:180; lambda = lambdaG*pi/180;
[Pnm,tmp, tmp, DegreeOrder] = legendreP(lmax,theta);
% ... selection of degree and order
index = DegreeOrder(:,1) == lmax;
Pnm   = Pnm(index,:);


%% Summation over k
% for-loop: not the best method, but easy to understand
ykm_Wignerd = 0;
for kk = 0:lmax
    ykm_Wignerd = ykm_Wignerd + LeNorm(kk,'geocomplex') * dlmk_plus(kk+1) * (Pnm(kk+1,:))'*exp(1i*kk*(lambda-alpha));
    ykm_Wignerd = ykm_Wignerd + LeNorm(-kk,'geocomplex') * dlmk_minus(kk+1) * (Pnm(kk+1,:))'*exp(-1i*kk*(lambda-alpha));
end
ykm_Wignerd = ykm_Wignerd *exp(-1i*m*gamma);               % rotation around gamma
ykm_Wignerd = ykm_Wignerd/LeNorm(m,'geocomplex');          % back to real normalization
ykm_Wignerd = real(ykm_Wignerd);                           % ensure real values

%% rotation of the coordinate system
R             = R3(gamma,'3x3')*R2(beta,'3x3')*R3(alpha,'3x3');
[LAMBDA,PHI]  = meshgrid(lambda,pi/2-theta);
[xyz(1,:),xyz(2,:),xyz(3,:)] = sph2cart(LAMBDA(:),PHI(:),1);
xyz           = R*xyz;
[LAMBDAr,PHIr]= cart2sph(xyz(1,:),xyz(2,:),xyz(3,:));
Ynm_coord     = plm(lmax,m,pi/2-PHIr(:)) .* exp(1i*m*LAMBDAr(:));
Ynm_coord     = reshape(real(Ynm_coord), size(ykm_Wignerd )); % viusalization of the cosine part in geodetic normalization


%% rotation of the (real) coefficients by Wigner-d-function
shc_cs = zeros(lmax+1,lmax+1); 
shc_cs(lmax+1,1+m) = 1;
shc_cs   = rotate_shc(shc_cs,alpha,beta,gamma);
shs_rot  =  gshs_ptw(shc_cs, LAMBDA(:), PHI(:), 1, 1,'quant', 'none','GM',1,'max_lm',lmax,'sub_WGS84',false);
shs_rot  = reshape(real(shs_rot), size(ykm_Wignerd ));

%% visualisation
alphaG = alpha*180/pi;
betaG = beta*180/pi;
figure
colormap(colbrew(1,129))
subplot(231);
imagesc(lambdaG,thetaG,ykm_Wignerd); 
hold on; plot(alphaG,betaG,'ko');
colorbar('Location', 'SouthOutside'); 
caxis([-3,3])
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[0:30:180]); axis equal; axis tight
title({'Rotation of spherical harmonics','by Wigner-d-functions'})

subplot(232);
imagesc(lambdaG,thetaG,Ynm_coord); 
hold on; plot(alphaG,betaG,'ko'); 
colorbar('Location', 'SouthOutside'); 
caxis([-3,3])
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[0:30:180]); axis equal; axis tight
title('Rotation of the coordinate system')

subplot(233);
imagesc(lambdaG,thetaG,shs_rot); 
hold on; plot(alphaG,betaG,'ko');
colorbar('Location', 'SouthOutside'); 
caxis([-3,3]) 
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[0:30:180]); axis equal; axis tight
title({'Rotation of SH-coefficients','by Wigner-d-functions'})


subplot(234);
imagesc(lambdaG,thetaG,Ynm_coord-ykm_Wignerd); 
hold on; plot(alphaG,betaG,'ko');
colorbar('Location', 'SouthOutside');  caxis(3e-14*[-1,1])
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[0:30:180]); axis equal; axis tight
title('Check: ''rotated coordinates'' minus ''rotated spherical harmoncis''')


subplot(236);
imagesc(lambdaG,thetaG,Ynm_coord-shs_rot); 
hold on; plot(alphaG,betaG,'ko'); 
colorbar('Location', 'SouthOutside');  caxis(3e-14*[-1,1])
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[0:30:180]); axis equal; axis tight
title('Check: ''rotated coordinates'' minus ''rotated SH-coefficients''')
