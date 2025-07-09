function fil = prpfltr(fil,lmax)

% PRPFLTR prepares the filter to applied on the SH co-efficients, supplied
% with the syn structure
%
% fil = prpfltr(fil,lmax)
%
% I/P
% fil   -   Filter in CS, SC or look-up-table [l m Clm Slm] formats
% lmax  -   Maximum degree of the given SH development
%
% O/P
% fil   -   Prepared filter
%--------------------------------------------------------------------------

if ~isint(lmax) || ~isequal(max(size(lmax)),1)
    error('LMAX must be a scalar integer.')
end

if isequal(size(fil),[lmax+1,2*lmax+1])
    fil = sc2cs(fil);
elseif isequal(size(fil),[sum(1:lmax+1),4])
    fil = sc2cs(clm2sc(fil));
elseif ~isequal(size(fil),[lmax+1 lmax+1])
    fil = [];
    fprintf('WARNING: The filter provided does not conform to CS, SC or Colombo format and hence will not be applied \n')
end
