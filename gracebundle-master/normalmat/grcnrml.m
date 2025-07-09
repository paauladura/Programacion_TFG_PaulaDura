function [N,Qxx,eigval] = grcnrml(pos1,pos2,lmax)

% GRCNRML(POS,MAXDEG) - GRACE Normal matrix Simulation, simulates the
% variance-covariance matrices based on the position of the two-satellites
% and the maximum degree specified. This function is designed for desktop
% calculations and hence, this function might not work for higher orders
% of degrees. The simulation is done based on observation equations for the
% difference in disturbing potential observed by the two satellites.
% 
% [N,Qxx,eigval] = grcnrml(pos1,pos2)
% [N,Qxx,eigval] = grcnrml(pos1,pos2,lmax)
%
% INPUT
% pos(1,2) - position in spherical co-ordinates with time in Julian days,
% longitude, co-latitude [radians], and radius [m] for both satellites as 
% two individual nx4 position vectors. [time,lambda,theta,r].
% Radius has to be given with the height of satellite from the center of
% the earth.
% lmax - maximum degree upto which the errors have to be simulated.
%
% OUTPUT
% N - Normal matrix
% Qxx - inverse of the normal matrix
% eigval - eigen values of the normal matrix
%
%--------------------------------------------------------------------------
% USES  GRACEBundle/blddsgn, blockinv 
%       uberall/constants_grace
%--------------------------------------------------------------------------

% Author: Balaji Devaraju (BD)
% Created on: 4 February 2008, Stuttgart
%
% Revision:
% 2014-03-05    BD  Added const as an input
%--------------------------------------------------------------------------

if nargin < 2
    error('Insufficient input arguments. Minimum two arguments required for executing the function \n')
elseif nargin < 3
    lmax = 60;
end

[lpos1,cpos1] = size(pos1);
[lpos2,cpos2] = size(pos2);

if cpos1 < 4 || cpos2 < 4
    error('Insufficient orbit information')
end

[nanr1,nanc1] = find(isnan(pos1));
[nanr2,nanc2] = find(isnan(pos2));

pos1(unique([nanr1;nanr2]),:) = [];
pos2(unique([nanr1;nanr2]),:) = [];

if lpos1<lpos2
    indx = zeros(lpos2,1);
    for i = 1:lpos1
        indx = indx + (pos1(i,1)==pos2(:,1));
    end
    indx = logical(indx);
    pos2 = pos2((indx==1),:);
elseif lpos1>lpos2
    indx = zeros(lpos1,1);
    for i = 1:lpos2
        indx = indx + (pos2(i,1)==pos1(:,1));
    end
    indx = logical(indx);
    pos1 = pos1(indx,:);
end

if isequal(pos1(:,1),pos2(:,1))
    display(' ')
    display('------------------------------------------------------')
    display('Provided co-ordinates match each other temporally')
    display('------------------------------------------------------')
    display(' ')
else
    temp = find(pos1(:,1)==pos2(:,1));
    if ~isempty(temp) && (length(temp)>2*(lmax+1)^2)
        fprintf('\n Only %f out of %f co-ordinates match \n',size(pos1,1),lpos1)
        pos1 = pos1(temp,:);
        pos2 = pos2(temp,:);
    else
        keyboard
        error('Check co-ordinates')
    end
end

lpos1 = size(pos1,1);

clmcnt = sum(1:lmax+1);
slmcnt = clmcnt-lmax-1;

cc = zeros(clmcnt);
ss = zeros(slmcnt);
sc = zeros(slmcnt,clmcnt);
elmnts = (lmax+1)^2;

k = 1:elmnts:lpos1;
cnt = [k(2:end)-1, lpos1];

lpos2 = [];
clear lpos2

constants_grace
const = [GM ae];

hwb = waitbar(0,'Percentage of normal matrix calculated ...');
set(hwb,'NumberTitle','off','Name','GRACE covariance simulation')
for i = 1:length(k)
    [c1,s1] = blddsgn(pos1(k(i):cnt(i),2:4),lmax,'potential',const);
    [c2,s2] = blddsgn(pos2(k(i):cnt(i),2:4),lmax,'potential',const);
    c1 = c1-c2;
    s1 = s1-s2;
    cc = cc + c1'*c1;
    sc = sc + s1'*c1;
    ss = ss + s1'*s1;
    waitbar(i/length(k))
end
close(hwb)
c1 = []; c2 = []; s1 = []; s2 = [];

idx = [1:2, lmax+2];

if nargout >= 2
    tmp = 1 - cc(3:end,3:end)./cc(3:end,3:end)';
    fprintf('Symmetry of cc: range %g mean %g std %g \n', range(1-tmp(:)), mean(tmp(:)), std(tmp(:)))
    tmp = 1 - ss(2:end,2:end)./ss(2:end,2:end)';
    fprintf('Symmetry of ss: range %g mean %g std %g \n', range(1-tmp(:)), mean(tmp(:)), std(tmp(:)))
    c1 = cc(3:end,3:end);
    c1(lmax,:) = [];
    c1(:,lmax) = [];
    s1 = sc(2:end,3:end);
    s1(:,lmax) = [];
    
    Qxx = blockinv(c1,ss(2:end,2:end),s1,'struct');
    c1 = []; s1 = [];
    
    Qxx.NW = [zeros(2,clmcnt-3);Qxx.NW(1:lmax-1,:); ...
                zeros(1,clmcnt-3);Qxx.NW(lmax:end,:)];
    Qxx.NW = [zeros(clmcnt,2) Qxx.NW(:,1:lmax-1) ...
                zeros(clmcnt,1) Qxx.NW(:,lmax:end)];
            
    Qxx.SE = [zeros(1,slmcnt-1);Qxx.SE];
    Qxx.SE = [zeros(slmcnt,1) Qxx.SE];
    
    Qxx.SW = [zeros(1,clmcnt-3);Qxx.SW];
    Qxx.SW = [zeros(slmcnt,2) Qxx.SW(:,1:lmax-1) ...
                zeros(slmcnt,1) Qxx.SW(:,lmax:end)];
end
N = [cc sc'; sc ss];
if nargout == 3, eigval = eig(N); end
