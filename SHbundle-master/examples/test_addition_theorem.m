%% TEST OF ADDITION THEOREM - 2014-10-27
close all;
clear all;
clc;

addpath('..'); % path to 'SHbundle'

theta_Q   = (0:1:90) / 180 * pi;

max_l = 30000; 
% max_l = 100;

% LEGENDRE = 'plm';                     % plm
LEGENDRE = 'Legendre_mex_plm_noxnum'; % standard plm-behaviour of Legendre_mex
% LEGENDRE = 'Legendre_mex_plm_xnum';   % plm-behaviour of Legendre_mex, X-numbers stabilized
% LEGENDRE = 'Legendre_mex';            % Legendre_mex (all coeffs)
% LEGENDRE = 'LegendreP';               % LegendreP (all coeffs but different order)
% LEGENDRE = 'LegendreP_Xnum';          % LegendreP with X-numbers (like LegendreP)

switch lower(LEGENDRE)
    case 'plm' % PLM of Nico
        P_b = zeros(length(theta_Q), max_l+1);
        tic; 
        for m = 0:max_l
            P_b = P_b + plm(0:max_l, m, theta_Q) .^ 2;
        end
        toc;
        A_b = ones(length(theta_Q), 1) * (2 * (0:max_l) + 1);
        R = (abs(A_b - P_b) ./ A_b);

    case 'legendre_mex_plm_noxnum' % PLM-Klon of Matthias (disable X-numbers)
        P_b = zeros(length(theta_Q), max_l+1);
        tic; 
        for m = 0:max_l
%             P_b = P_b + Legendre_mex(0:max_l, m, theta_Q) .^ 2; % if Lmax > 1500: enable X-numbers
            P_b = P_b + Legendre_mex(0:max_l, m, theta_Q, 'plm', 'noxnum') .^ 2;
        end
        toc;
        A_b = ones(length(theta_Q), 1) * (2 * (0:max_l) + 1);
        R = (abs(A_b - P_b) ./ A_b);
        
    case 'legendre_mex_plm_xnum' % PLM-Klon of Matthias (enable X-numbers)
        P_b = zeros(length(theta_Q), max_l+1);
        tic; 
        for m = 0:max_l
            P_b = P_b + Legendre_mex(0:max_l, m, theta_Q, 'plm', 'xnum') .^ 2;
        end
        toc;
        A_b = ones(length(theta_Q), 1) * (2 * (0:max_l) + 1);
        R = (abs(A_b - P_b) ./ A_b);
        
    case 'legendrep' % LegendreP of Hailong
        tic; P_theta_Q = legendreP(max_l, theta_Q); toc; 
        R(length(theta_Q), max_l + 1) = 0;
        for idx_q = 1:size(theta_Q, 2)
            for l = 0:max_l
                idx_s = sum(1:l) + 1; % for element (l, 0)
                idx_e = idx_s + l;    % for element (l, l)
                R(idx_q, l + 1) = abs(1 - (2 * l + 1)^(-1) * sum(P_theta_Q(idx_s:idx_e, idx_q) .^ 2));                
            end
        end    
        
    case 'legendrep_xnum' % LegendreP with X-numbers of Hailong
        tic; P_theta_Q  = legendreP_Xnum(max_l, theta_Q); toc;
        R(length(theta_Q), max_l + 1) = 0;
        for idx_q = 1:size(theta_Q, 2)
            for l = 0:max_l
                idx_s = sum(1:l) + 1; % for element (l, 0)
                idx_e = idx_s + l;    % for element (l, l)
                R(idx_q, l + 1) = abs(1 - (2 * l + 1)^(-1) * sum(P_theta_Q(idx_s:idx_e, idx_q) .^ 2));                
            end
        end    
        
    case 'legendre_mex' % fast mex-Version of Matthias (Lmax > 1500: uses X-numbers)
        tic; P_theta_Q_ = Legendre_mex(max_l, theta_Q, 'speed_od'); toc;
        % rearrange
        P_theta_Q(size(P_theta_Q_, 1), size(P_theta_Q_, 2)) = 0;%   = P_phi_Q_;
        idx_old(max_l * ((max_l + 1) / 2 + 1)) = 0; % initialize array
        m = 0:max_l;
        m_pre = m * (max_l + 1) - m .* (m - 1)/2 - m + 1;
        i = 1;
        for l = 0:max_l
            idx_old(i:(i + max_l)) = m_pre + l;
            i = i + l + 1;
        end
        idx_new = 1:(max_l * ((max_l + 1) / 2 + 1) + 1);
        P_theta_Q(idx_new, :) = P_theta_Q_(idx_old, :);
      
        R(length(theta_Q), max_l + 1) = 0;
        for idx_q = 1:size(theta_Q, 2)
            for l = 0:max_l
                idx_s = sum(1:l) + 1; % for element (l, 0)
                idx_e = idx_s + l;    % for element (l, l)
%                 R(idx_q, l + 1) = abs(1 - 1 / (2 * l + 1) * sum(P_theta_Q(idx_s:idx_e, idx_q).^2));
                R(idx_q, l + 1) = abs(1 - (2 * l + 1)^(-1) * sum(P_theta_Q(idx_s:idx_e, idx_q) .^ 2));                
            end
        end    
 
    otherwise
        error('not implemented');
end

figure;
pcolor(log10(R));
shading flat;
title(sprintf('Addition Theorem: %s, max(R) = %1.2e', LEGENDRE, max(R(:))), 'Interpreter', 'none');
xlabel('degree l');
ylabel('co-latitude \theta');
% axis equal;
axis tight;
colorbar;
caxis([-17 -10]);


