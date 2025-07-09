function gfc_ = removetemporalmean(gfc, mean_gfc)
    
% REMOVETEMPORALMEAN removes the temporal mean computed by the function
% GETGRACEMEAN. 
%
% gfc_ = removetemporalmean(gfc, mean_gfc)
%
%
% INPUT
% gfc       - GRACE spherical harmonic coefficients in the 10-column
%             cell-array format
% mean_gfc  - Temporal mean of the GRACE coefficients
%
% OUTPUT
% gfc_      - Anomalies of the monthly GRACE spherical harmonic 
%             coefficients given in the same format as 'gfc'
%
%-----------------------------------------------------------------------

% Balaji Devaraju, IIT Kanpur, 12/7/2019

gfc_ = gfc;
for k = 1:size(gfc, 1)
    gfc_{k, 9} = gfc{k, 9} - mean_gfc;
    gfc_{k, 10} = sqrt(2)*gfc{k,10};
end