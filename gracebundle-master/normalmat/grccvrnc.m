function grccvrnc(dvec,mnpth,lmax,cntr,rls)

% GRCCVRNC manages the GRACE L1B data to compute the simulated covariance
% matrices for each month of the time-variable gravity field solution.
%
% grccvrnc(ymdd,mnpth,lmax,cntr,rls)
%
% INPUT
% dvec      -   [Year Month StartDay EndDay]
% mnpth     -   Path where the L1B navigation position data are stored
% lmax      -   Maximum degree of spherical harmonic expansion
% cntr      -   Centre that provides the data [string]
% rls       -   Release number of the dataset [string]
%
% OUTPUT
% Simulated covariance matrices will be generated and will be stored in the
% same folder as the 'mnpth/year/' 
%--------------------------------------------------------------------------
% USES GRACE/grcnrml
%--------------------------------------------------------------------------

% Created on: 4 Spetember 2009, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

tic

if nargin==2
    lmax = 60;
    cntr = 'GRC';
    rls = 'xx';
elseif nargin==3
    cntr = 'GRC';
    rls = 'xx';
elseif nargin==4
    rls = 'xx';
elseif nargin<2
    error('Insufficient input arguments')
elseif lmax<0
    error('Degree of the spherical harmonics must always be positive')
end

[n,m] = size(dvec);

if m<4
    error('Start and end days of observation for each month missing')
end

for i = 1:n
    fname = [mnpth,num2str(dvec(i,1)),'/','GNV1B',num2str(dvec(i,1))];
    grcl1b = load(fname);
    fld = fieldnames(grcl1b);
    grcl1b = grcl1b.(fld{1});
    
    tmp = cell2mat(grcl1b(:,1));
    ind1 = find(tmp(:,4) == dvec(i,3));
    if isempty(ind1)
        itmp = find(tmp(:,1)==dvec(i,1) & tmp(:,2)==dvec(i,2));
        ind1 = itmp(1);
    end
    ind2 = find(tmp(:,4) == dvec(i,4));
    if isempty(ind2)
        itmp = find(tmp(:,1)==dvec(i,1) & tmp(:,2)==dvec(i,2));
        ind2 = itmp(end);
    end
    
    posA = cell2mat(grcl1b(ind1:ind2,2));
    posB = cell2mat(grcl1b(ind1:ind2,3));
    
    posA(:,3) = pi/2 - posA(:,3);
    posB(:,3) = pi/2 - posB(:,3);
    
    [N,Q,e] = grcnrml(posA,posB,lmax);
    
    fname = [mnpth,num2str(dvec(i,1)),'/', ...
        cntr,rls,num2str(dvec(i,1)),num2str(dvec(i,2)),'N',num2str(lmax)];
    save(fname,'N','e')
    fname = [mnpth,num2str(dvec(i,1)),'/', ...
        cntr,rls,num2str(dvec(i,1)),num2str(dvec(i,2)),'Q',num2str(lmax)];
    save(fname,'Q')
    fprintf('%s done \n',mat2str(dvec(i,:)))
end

fprintf('done \n')

toc
