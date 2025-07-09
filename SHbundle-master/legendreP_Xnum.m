function [P, dP, ddP] = legendreP_Xnum(kmax, thetaRAD)

% Compute normalised legendre function using X number
% Refer to paper: Numerical computation of spherical harmonics of arbitrary
% degree and order by extending exponent of floating point numbers
% Any number can be writen in:  X = x*B^expo,  B is 2^960

% IN:
%    kmax ....... maximum degree of the Legendre functions          [1 x 1]
%    thetaRAD ... co-latitude in [rad]                              [N x 1]     
%
% OUTPUT
%    P .......... normalized Legendre functions                     [N x M]
%    dP ......... first order derivatives of Legendre functions     [N x M]   
%    ddP ........ second order derivatives of Legendre functions    [N x M]   
%
% REMARKS:
%    1. P, dP, dPP will be 2D matrices
%    2. all 2D matrix indices are [ N(length of thetaRAD) * M(kmax+1) ]
%    3. if dP, dPP are not specified as outputs, they won't be computed
%    4. if thetaRAD is not N x 1, it will be reshaped

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Hailong FU (HF), Uni Stuttgart
%    Matthias ROTH (MR), GIS, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-09-29: MR, remove unnecessary output parameters, update help text
%    2013-06-30: HF, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

%% check number of output parameters
if nargout >= 2
    test_dP = true;
else
    test_dP = false;
end

if nargout >= 3
    test_ddP = true;
else
    test_ddP = false;
end

%% 2D parallel recursion
dmax               = kmax + 1;
Wlm_1_pre(1, dmax) = 0; % faster than "Wlm_1_pre = zeros(1, dmax);"
Wlm_1(1, dmax)     = 0; % faster than "Wlm_1 = zeros(1, dmax);"
Wlm_2(1, dmax)     = 0; % dito
Wlm_1_pre(2)       = sqrt(3);
Wlm_1_pre(1)       = sqrt(3);

siz      = size(thetaRAD); % reshape col-latitude to be one column
lenLat   = siz(1) * siz(2);
thetaRAD = reshape(thetaRAD, lenLat, 1);
coslat   = cos(thetaRAD);
sinlat   = sin(thetaRAD); 
cotlat   = cot(thetaRAD);
 
P(lenLat, (1 + dmax) * dmax / 2) = 0.0; % faster than "P = zeros(lenLat, (1 + dmax) * dmax / 2);"
P(:, 1) = 1;
P(:, 2) = sqrt(3) * coslat;
P(:, 3) = sqrt(3) * sinlat;

% initialize x-numbers
IND   = 960;
BIG   = 2^IND;

P_expo = zeros(lenLat, (1 + dmax) * dmax / 2, 'int8'); % exponent of X number

if test_dP
    dP(lenLat, (1 + dmax) * dmax / 2) = 0.0; % faster than "dP = zeros(lenLat, (1 + dmax) * dmax / 2);"
    dP(:, 2) = -P(:, 3);
    dP(:, 3) =  P(:, 2);
    dP_expo  =  P_expo;
    if test_ddP
        ddP(lenLat, (1 + dmax) * dmax / 2) = 0.0; % faster than "ddP = zeros(lenLat, (1 + dmax) * dmax / 2);"
        ddP(:, 2) = -P(:, 2);
        ddP(:, 3) = -P(:, 3);
        ddP_expo  =  P_expo;
    end    
end

preCols = 3; % previously finished columns

for l = 2:kmax % if l >= 2
    lp1 = l + 1;  
    m = 0:(l - 1);  
    Wlm_1(m + 1)   = sqrt((4 * l^2 - 1) ./ (l^2 - m.^2));
    Wlm_1(lp1)     = sqrt(1 + 0.5 / l);  
    m_             = 1:(l - 1);  % 1 based index, for whom pre-pre P exists
    Wlm_2(m_)      = Wlm_1(m_) ./ Wlm_1_pre(m_); 
    accum_m_       = preCols + m_;    
    accum_m_pre    = accum_m_ - l;   
    accum_m_prepre = accum_m_pre - (l - 1);  
    
%% calculate P
    P_expo_dif_pre_prepre = P_expo(:, accum_m_pre) - P_expo(:, accum_m_prepre);
    idx1 = P_expo_dif_pre_prepre > 0;
    idx2 = P_expo_dif_pre_prepre < 0;
    
    C1 = ones(lenLat, l - 1); % coefficients for X number
    C2 = C1;

    C2(idx1) = BIG.^(-double(P_expo_dif_pre_prepre(idx1)));
    C1(idx2) = BIG.^(double(P_expo_dif_pre_prepre(idx2)));
    
    % prepare Exponent array
    P_expo_cur          = P_expo(:, accum_m_pre); 
    P_expo_prepre       = P_expo(:, accum_m_prepre);
    P_expo_cur(idx2)    = P_expo_prepre(idx2);
    P_expo(:, accum_m_) = P_expo_cur;

    P_expo(:, preCols + l)   = P_expo(:, preCols);
    P_expo(:, preCols + lp1) = P_expo(:, preCols);    
 
    % Legendre
    P(:, accum_m_)      = bsxfun(@times, bsxfun(@times, P(:, accum_m_pre) .* C1, coslat),  Wlm_1(m_)) ...
                          - bsxfun(@times, P(:, accum_m_prepre) .* C2, Wlm_2(m_));         
    P(:, preCols + l)   = bsxfun(@times, bsxfun(@times, P(:, preCols), coslat),  Wlm_1(l));                   
    P(:, preCols + lp1) = Wlm_1(lp1) * bsxfun(@times, P(:, preCols), sinlat);

    accum_m = preCols + (1:(l+1));
    [P(:, accum_m), P_expo(:, accum_m) ] = xnorm(P(:, accum_m), P_expo(:, accum_m));
    
%% calculate dP
    if test_dP
        accum_m_                = preCols + (1:l); 
        accum_m_pre             = accum_m_ - l; 
        
        P_expo_dif_cur_pre      = P_expo(:,accum_m_) - P_expo(:,accum_m_pre);
        idx1 = P_expo_dif_cur_pre>0;
        idx2 = P_expo_dif_cur_pre<0;
        
        C1 = ones(lenLat, l); % coefficients for X number
        C2 = C1;

        C2(idx1) = BIG.^(-double(P_expo_dif_cur_pre(idx1)));
        C1(idx2) = BIG.^(double(P_expo_dif_cur_pre(idx2)));
        
        % prepare Exponent array
        dP_expo_cur = P_expo(:,accum_m_); 
        P_expo_pre = P_expo(:,accum_m_pre);
        dP_expo_cur(idx2) = P_expo_pre(idx2);
        dP_expo(:,accum_m_) = dP_expo_cur;
        dP_expo(:,preCols+lp1) = P_expo(:,preCols+lp1); 

        e           = sqrt((2*l+1)/(2*l-1)*(l*l-(0:l-1).*(0:l-1))); 
        dP(:, accum_m_) = bsxfun(@times, P(:, accum_m_).*C1, l*cotlat) - (1./sinlat*e).* (P(:, accum_m_pre).*C2);
        dP(:, preCols+lp1) = bsxfun(@times, P(:, preCols+lp1), l*cotlat);
        
        [dP(:, accum_m), dP_expo(:, accum_m) ] = xnorm(dP(:, accum_m), dP_expo(:, accum_m));
        
%% calculate ddP
        if test_ddP
            dP_P_expo_dif      = dP_expo(:,accum_m) - P_expo(:,accum_m);
            idx1 = dP_P_expo_dif>0;
            idx2 = dP_P_expo_dif<0;
            
            C1 = ones(lenLat, l + 1); % coefficients for X number
            C2 = C1;
            
            C2(idx1) = BIG.^(-double(dP_P_expo_dif(idx1)));
            C1(idx2) = BIG.^(double(dP_P_expo_dif(idx2)));
            
            % prepare Exponent array
            ddP_expo_cur = dP_expo(:,accum_m); 
            P_expo_cur = P_expo(:,accum_m);
            ddP_expo_cur(idx2) = P_expo_cur(idx2);
            ddP_expo(:,accum_m) = ddP_expo_cur;
            
            ddP(:, accum_m) = bsxfun(@times, dP(:, accum_m).*C1, -cotlat) + (1./sinlat./sinlat*((0:l).*(0:l))-l*(l+1)).*( P(:, accum_m).*C2 );
            [ddP(:, accum_m), ddP_expo(:, accum_m) ] = xnorm(ddP(:, accum_m), ddP_expo(:, accum_m));
        end
    end        
    
    Wlm_1_pre = Wlm_1;
    preCols = preCols + lp1;
end

%% transpose that output is the same as Markus'
P      = x2f(P, P_expo)';
if test_dP
    dP = x2f(dP, dP_expo)';
    if test_ddP
        ddP = x2f(ddP, ddP_expo)';        
    end
end

end

%% internal function X2F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = x2f(x, ix)

IND   = 960;
BIG   = 2^IND;
BIGI  = 2^(-IND);

idx = ix < 0;
x(idx) = x(idx) * BIGI;

idx = ix > 0;
x(idx) = x(idx) * BIG;

end

%% internal function XNORM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, ix] = xnorm(x, ix)

IND   = 960;
BIG   = 2^IND;
BIGI  = 2^(-IND);
BIGS  = 2^(IND / 2);
BIGSI = 2^(-IND / 2);

w = abs(x);

idx = w >= BIGS;
 x(idx) = x(idx) * BIGI; ix(idx) = ix(idx) + 1;
 
idx = w < BIGSI;
 x(idx) = x(idx) * BIG;  ix(idx) = ix(idx) - 1;

end
