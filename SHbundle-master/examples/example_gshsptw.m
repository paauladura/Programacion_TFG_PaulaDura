%% EXAMPLE for the use of GSHS_PTW                    - 2015-09-02
%  removed the demonstration of depricated functions  - 2021-04-12
%  added comparison with gshs_grid + GUI for selection- 2021-10-19
clear all;
clc;



% checking available path
if exist('clm2sc','file') ~= 2,       error('Please add path of ''SHBundle'''),   end; 
if exist('checkcoor','file') ~= 2,    error('Please add path of ''uberall'''),    end; 

%% coefficients and constants
constants; % load constants from UBERALL/constants
[gsm, lmax, lmin, info] = parse_icgem('./data/example.icgem', 'max_lm', 10); % parse ICGEM-file up to max_lm = 10
[field, lmax]           = clm2sc(gsm, 'max_lm', lmax); % convert l,m-indexed list into /S|C\-format


%% selection of field quantity by user
quantity = {'none','potential','geoid','tr','trr','water','smd'};
[select,ok] = listdlg('PromptString','quanitity:','SelectionMode','single','ListString',quantity);
if ok == false,  select = 2; end
quantity = quantity{select};


%% defining coordinates
row = 100; col = 200;
[lamRAD,phiRAD] = meshgrid(linspace(-pi, pi, col), linspace(-pi/2, pi/2, row)); % define coordinates
r      = ones(row,col) * 7e6;

%% spherical harmonic synthesis
% ... by the speed optimized function for point-wise locations...
tic; F_pwt = gshs_ptw(field, lamRAD(:), phiRAD(:), r(:), ae, 'GM', GM, 'max_lm', lmax, 'quant', quantity, 'waitbar', false); toc; % calculate potential for points (lamRAD, phiRAD)
F_pwt = reshape(F_pwt,row,col);
% ... by the speed optimized function for grid-wise locations...
lamRAD = lamRAD(1,:);  phiRAD =  phiRAD(:,1);
tic; f_grid   = gshs_grid(field,lamRAD, phiRAD, ae, 'GM', GM, 'height', r(1)-ae, 'max_lm', lmax, 'quant', quantity, 'sub_wgs84', true, 'curvature', false, 'waitbar', true, 'legendre', 'mex'); 


%% visualization
figure
subplot(221);
imagesc(lamRAD*180/pi, phiRAD*180/pi,F_pwt); shading flat; 
xlabel('longitude'); ylabel('latitude'); title(['SHS(' quantity ') by ''gshs\_ptw'' for point set and on' ]); 
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[-90:30:90]); 
axis equal; axis tight; colorbar

subplot(222);
imagesc(lamRAD*180/pi, phiRAD*180/pi,f_grid);shading flat; 
xlabel('longitude'); ylabel('latitude'); title(['SHS(' quantity ') by ''gshs\_grid'' on a grid and on a sphere']); 
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[-90:30:90]); 
axis equal; axis tight; colorbar


subplot(223); 
imagesc(lamRAD*180/pi, phiRAD*180/pi,F_pwt-f_grid); shading flat; 
xlabel('longitude'); ylabel('latitude'); title('difference: ''gshs\_ptw'' minus ''gshs\_grid'''); 
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[-90:30:90]); 
axis equal; axis tight; colorbar