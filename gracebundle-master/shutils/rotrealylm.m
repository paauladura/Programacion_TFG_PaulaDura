function [Ynmrot,Ynm] = rotrealylm(l,m,theta,lambda,alpha,beta,gamma)

% ROTREALYLM rotates real spherical harmonics using the Euler rotations
%
% [Ylmrot,Ylm] = rotrealylm(l,m)
% [Ylmrot,Ylm] = rotrealylm(l,m,theta,lambda)
% [Ylmrot,Ylm] = rotrealylm(l,m,theta,lambda,alpha)
% [Ylmrot,Ylm] = rotrealylm(l,m,theta,lambda,alpha,beta)
% [Ylmrot,Ylm] = rotrealylm(l,m,theta,lambda,alpha,beta,gamma)
%
% INPUT 
% l     -   spherical harmonic degree
% m     -   spherical harmonic order
% theta -   Co-latitude points of the desired grid! Must be a vector.
%                                                               [radians]
% lambda-   Longitude points of the desired grid. Must be a vector.
%                                                               [radians]
% alpha,beta,gamma -   Euler rotation angles                    [radians]
%
% OUTPUT
% Ynmrot -  Rotated spherical harmonics for cosine and sine components of
%           degree n and order m.
% Ynm    -  Unrotated spherical harmonics for cosine and sine components of
%           degree n and order m.
% -------------------------------------------------------------------------
%
% USES wignerD, ylm
%
% See also eulerYnm
%--------------------------------------------------------------------------

if nargin == 2
    Ynm = ylm(l,m);
    [r,c] = size(Ynm);
    Ynm= zeros(2*l+1,size(Ynm(:)',2));
    k = (l:-1:1);
    ylmc = ylm(l,0);
    Ynm(l+1,:) = ylmc(:)';
    for i = 1:l
        [ylmc,ylms] = ylm(l,k(i));
        Ynm(i,:) = ylms(:)';
        Ynm(2*l+2-i,:) = ylmc(:)';
    end
    
    D = wignerD(l,pi/6,pi/4,pi/3,'r');
    
    tmp = D*Ynm;
    Ynmrot.s = reshape(tmp(l-m+1,:),r,c);
    Ynmrot.c = reshape(tmp(l+m+1,:),r,c);
    
    tmp = Ynm;
    Ynm = [];
    Ynm.s = reshape(tmp(l-m+1,:),r,c);
    Ynm.c = reshape(tmp(l+m+1,:),r,c);
elseif nargin == 4
    Ynm = ylm(l,m,theta,lambda);
    [r,c] = size(Ynm);
    Ynm= zeros(2*l+1,size(Ynm(:)',2));
    k = (l:-1:1);
    ylmc = ylm(l,0,theta,lambda);
    Ynm(l+1,:) = ylmc(:)';
    for i = 1:l
        [ylmc,ylms] = ylm(l,k(i),theta,lambda);
        Ynm(i,:) = ylms(:)';
        Ynm(2*l+2-i,:) = ylmc(:)';
    end
    
    D = wignerD(l,pi/2,pi/2,pi/2,'r');
    
    tmp = D*Ynm;
    Ynmrot.s = reshape(tmp(l-m+1,:),r,c);
    Ynmrot.c = reshape(tmp(l+m+1,:),r,c);
    
    tmp = Ynm;
    Ynm = [];
    Ynm.s = reshape(tmp(l-m+1,:),r,c);
    Ynm.c = reshape(tmp(l+m+1,:),r,c);
elseif nargin == 5
    Ynm = ylm(l,m,theta,lambda);
    [r,c] = size(Ynm);
    Ynm= zeros(2*l+1,size(Ynm(:)',2));
    k = (l:-1:1);
    ylmc = ylm(l,0,theta,lambda);
    Ynm(l+1,:) = ylmc(:)';
    for i = 1:l
        [ylmc,ylms] = ylm(l,k(i),theta,lambda);
        Ynm(i,:) = ylms(:)';
        Ynm(2*l+2-i,:) = ylmc(:)';
    end
    
    D = wignerD(l,alpha,0,0,'r');
    
    tmp = D*Ynm;
    Ynmrot.s = reshape(tmp(l-m+1,:),r,c);
    Ynmrot.c = reshape(tmp(l+m+1,:),r,c);
    
    tmp = Ynm;
    Ynm = [];
    Ynm.s = reshape(tmp(l-m+1,:),r,c);
    Ynm.c = reshape(tmp(l+m+1,:),r,c);
elseif nargin == 6
    Ynm = ylm(l,m,theta,lambda);
    [r,c] = size(Ynm);
    Ynm= zeros(2*l+1,size(Ynm(:)',2));
    k = (l:-1:1);
    ylmc = ylm(l,0,theta,lambda);
    Ynm(l+1,:) = ylmc(:)';
    for i = 1:l
        [ylmc,ylms] = ylm(l,k(i),theta,lambda);
        Ynm(i,:) = ylms(:)';
        Ynm(2*l+2-i,:) = ylmc(:)';
    end
    
    D = wignerD(l,alpha,beta,0,'r');
    % D = D{l+1};
    
    tmp = D*Ynm;
    Ynmrot.s = reshape(tmp(l-m+1,:),r,c);
    Ynmrot.c = reshape(tmp(l+m+1,:),r,c);
    
    tmp = Ynm;
    Ynm = [];
    Ynm.s = reshape(tmp(l-m+1,:),r,c);
    Ynm.c = reshape(tmp(l+m+1,:),r,c);
elseif nargin == 7
    Ynm = ylm(l,m,theta,lambda);
    [r,c] = size(Ynm);
    Ynm= zeros(2*l+1,size(Ynm(:)',2));
    k = (l:-1:1);
    ylmc = ylm(l,0,theta,lambda);
    Ynm(l+1,:) = ylmc(:)';
    for i = 1:l
        [ylmc,ylms] = ylm(l,k(i),theta,lambda);
        Ynm(i,:) = ylms(:)';
        Ynm(2*l+2-i,:) = ylmc(:)';
    end
    
    D = wignerD(l,alpha,beta,gamma,'r');
    
    tmp = D*Ynm;
    Ynmrot.s = reshape(tmp(l-m+1,:),r,c);
    Ynmrot.c = reshape(tmp(l+m+1,:),r,c);
    
    tmp = Ynm;
    Ynm = [];
    Ynm.s = reshape(tmp(l-m+1,:),r,c);
    Ynm.c = reshape(tmp(l+m+1,:),r,c);
end


if nargout<1
    [Theta,Lambda] = meshgrid(theta,lambda);
    X = cos(Lambda).*sin(Theta);
    Y = sin(Lambda).*sin(Theta);
    Z = cos(Theta);

    figure
    axH(1) = subplot('position',[0.05 0.4 0.5 0.5]); surf(X,Y,Z,real(Ynm.c)')
    axis equal; shading interp; 
    title('original spherical harmonic'); axis off
    axH(2) = subplot('position',[0.55 0.4 0.5 0.5]); surf(X,Y,Z,real(Ynmrot.c)')
    axis equal; shading interp; axis off;
    title('rotated spherical harmonic') 
    
    MAX = max(max(Ynm.c));
    clim = [-MAX,MAX];
    
    caxis(clim);
    set(axH,'CLim',clim);
    h=colorbar;
    set(h,'position',[0.5 0.4 0.05 0.4])
end