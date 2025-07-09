%% TEST FOR ORDINARY DIFFERENTIAL EQUATION (ODE) - 2014-10-27

close all;
clear all;
clc;

addpath('..'); % path to 'SHbundle'

theta = 45 * pi / 180;

% Lmax = 100;
Lmax = 1000;

LEGENDRE = 'legendre_mex_plm';                 % plm
% LEGENDRE = 'Legendre_mex_plm';  % plm-behaviour of Legendre_mex
% LEGENDRE = 'Legendre_mex';      % Legendre_mex (all coeffs)

switch lower(LEGENDRE)
    case 'plm'
        R = nan(Lmax+1);
        for m = 0:Lmax
            [P, dP, ddP] = plm(m:Lmax, m, theta);   
            for l = m:Lmax
                R(m + 1, l + 1) = ddP(l - m + 1) + dP(l - m + 1) * cos(theta) / sin(theta) + P(l - m + 1) * (l^2 + l - (m / sin(theta))^2); % nico
            end
        end

    case 'legendre_mex_plm'
        R = nan(Lmax+1);
        for m = 0:Lmax
            [P, dP, ddP] = Legendre_mex(m:Lmax, m, theta);   
            for l = m:Lmax
                R(m + 1, l + 1) = ddP(l - m + 1) + dP(l - m + 1) * cos(theta) / sin(theta) + P(l - m + 1) * (l^2 + l - (m / sin(theta))^2);
            end
        end

    case 'legendre_mex'
        [P, dP, ddP] = Legendre_mex('speed_od', Lmax, theta); 
        i = 0;
        R = nan(Lmax+1);
        for m = 0:Lmax
            for l = m:Lmax
                i = i + 1;
                R(m + 1, l + 1) = ddP(i) + dP(i) * cos(theta) / sin(theta) + P(i) * (l^2 + l - (m / sin(theta))^2);
            end
        end

    otherwise
        error('not implemented');
end


figure;
pcolor(log10(abs(R)));
shading flat;
colorbar;
caxis([-17 -8]);
title(sprintf('Diff. Eq.: %s, max(v) = %1.2e', LEGENDRE, max(abs(R(:)))), 'Interpreter', 'none');
xlabel('degree l');
ylabel('order m');
axis equal;
axis tight;        
