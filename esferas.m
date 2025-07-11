degree=8;
order=8;
grid=40;
radius=5;

% Create the grid
delta = pi/grid;
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : 2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta)

% Calculate the harmonic
Ymn = legendre(degree,cos(theta(:,1)))
Ymn = Ymn(order+1,:)';
yy = Ymn;
for kk = 2: size(theta,1);
    yy = [yy Ymn]
end;
yy = yy.*cos(order*phi);

order = max(max(abs(yy)));
rho = radius + 2*yy/order;

% Apply spherical coordinate equations
r = rho.*sin(theta);
x = r.*cos(phi);  % spherical coordinate equations
y = r.*sin(phi);
z = rho.*cos(theta);

% Plot the surface
clf
surf(x,y,z)
light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.5)
colormap("jet")