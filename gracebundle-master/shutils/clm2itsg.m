function sh_coeff = clm2itsg(clm)

% CLM2ITSG converts the spherical harmonic coefficients given in four
% column order-leading format to the ITSG format
%
% ITSG format           clm format
% [l    m   klm]        [l  m   Clm     Slm]
% [0    0   C00]        [0  0   C00     0]
% [1    0   C10]        [1  0   C10     0]
% [1    1   C11]        [2  0   C20     0]
% [1    -1  S11]        [1  1   C11     S11]
% [2    0   C20]        [2  1   C21     S21]
% [2    1   C21]        [2  2   C22     S22]
% [2    -1  S21]
% [2    2   C22]
% [2    -2  S22]
%
% INPUT
% clm       - Spherical harmonic coefficients in clm-format
%
% OUTPUT
% sh_coeff  - Rearranged coefficients in ITSG-format
%

narginchk(2, 2)

sh_coeff = itsg_sh_ordering(max(clm(:,1)), 0);
sh_coeff = [sh_coeff, zeros(length(sh_coeff), 1)];

for k = 1:length(sh_coeff)
    idx = find((clm(:,1) == sh_coeff(k,1)) & (clm(:,2) == abs(sh_coeff(k,2))))
    if sh_coeff(k, 2) < 0
        sh_coeff(k, 3) = clm(idx, 4);
    else
        sh_coeff(k, 3) = clm(idx, 3);
    end
end
