close all;
clear all;
clc;

addpath('..');            % path to 'SHbundle'
addpath('../../uberall'); % path to 'uberall'

max_l = 2500;
% max_l = 100;

% LEGENDRE = 'plm';                      % plm
% LEGENDRE = 'legendre_mex_plm_noxnum';  % plm-behaviour of Legendre_mex
LEGENDRE = 'legendre_mex_plm_xnum';    % Legendre_mex with X-numbers

% NOT USABLE (memory overflow for higher degree/order): Legendre_mex, LegendreP, LegendreP_Xnum

switch lower(LEGENDRE)
    case 'plm' % PLM of Nico
        tic;
        v = nan(max_l + 1, max_l + 1);
        w = nan(max_l + 1, max_l + 1);
        x = nan(max_l + 1, max_l + 1);
        for l = 0:max_l
            [w((0:l) + 1, l + 1), x((0:l) + 1, l + 1)] = neumann(l + 1); % w ... Neumann weight elements
        end
        theta = acos(x);
        
        for m = 0:max_l
            fprintf('m = %d\n', m);
            for l = m:max_l
                 P = plm(l, m, theta((0:l) + 1, l + 1))';
            
                 if m == 0
                     v(m + 1, l + 1) = abs(P.^2 * w((0:l) + 1, l + 1) - 2);
                 else
                     v(m + 1, l + 1) = abs(P.^2 * w((0:l) + 1, l + 1) - 4);
                 end
            end
        end
        toc;
        
    case 'legendre_mex_plm_noxnum' % PLM-Klon of Matthias (disable X-numbers)
        tic;
        v = nan(max_l + 1, max_l + 1);
        w = nan(max_l + 1, max_l + 1);
        x = nan(max_l + 1, max_l + 1);
        for l = 0:max_l
            [w((0:l) + 1, l + 1), x((0:l) + 1, l + 1)] = neumann(l + 1); % w ... Neumann weight elements
        end
        theta = acos(x);
        
        for m = 0:max_l
            fprintf('m = %d\n', m);
            for l = m:max_l
                 P = Legendre_mex(l, m, theta((0:l) + 1, l + 1), 'plm', 'noxnum')'; 
                 if m == 0
                     v(m + 1, l + 1) = abs(P.^2 * w((0:l) + 1, l + 1) - 2);
                 else
                     v(m + 1, l + 1) = abs(P.^2 * w((0:l) + 1, l + 1) - 4);
                 end
            end
        end
        toc;
    case 'legendre_mex_plm_xnum' % PLM-Klon of Matthias (enable X-numbers)
        tic;
        v = nan(max_l + 1, max_l + 1);
        w = nan(max_l + 1, max_l + 1);
        x = nan(max_l + 1, max_l + 1);
        for l = 0:max_l
            [w((0:l) + 1, l + 1), x((0:l) + 1, l + 1)] = neumann(l + 1); % w ... Neumann weight elements
        end
        theta = acos(x);
        
        for m = 0:max_l
            fprintf('m = %d\n', m);
            for l = m:max_l
                 P = Legendre_mex(l, m, theta((0:l) + 1, l + 1), 'plm', 'xnum')'; 
                 if m == 0
                     v(m + 1, l + 1) = abs(P.^2 * w((0:l) + 1, l + 1) - 2);
                 else
                     v(m + 1, l + 1) = abs(P.^2 * w((0:l) + 1, l + 1) - 4);
                 end
            end
        end
        toc;

     otherwise
        error('not implemented');
end

%% output
figure;
pcolor(log10(v));
shading flat;
colorbar;
caxis([-17 -10]);
title(sprintf('Orthogonality: %s, max(v) = %1.2e', LEGENDRE, max(v(:))), 'Interpreter', 'none');
xlabel('degree l');
ylabel('order m');
axis equal;
axis tight;

