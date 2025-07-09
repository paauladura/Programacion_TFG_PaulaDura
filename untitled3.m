% PASO 1 Crear rejilla esférica
[theta, phi] = meshgrid(linspace(0, pi, 100), linspace(0, 2*pi, 100));

% Calcular el valor del armónico esférico
n=0;
m=0;
Y = compute_spherical_harmonic(n, m, theta, phi);

% Convertir a coordenadas cartesianas
r = 1; % Radio constante para esfera perfecta
x = r * sin(theta) .* cos(phi);
y = r * sin(theta) .* sin(phi);
z = r * cos(theta);
figure(1)
% Graficar superficie
surf(x, y, z, Y, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet")

% PASO 2 Calcular el valor del armónico esférico
n=4;
m=0;
Y40 = compute_spherical_harmonic(n, m, theta, phi);
% Graficar superficie
figure(2)
surf(x, y, z, Y40, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet")
figure(3)
surf(x, y, z, Y+Y40, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet")

% PASO 3 Calcular el valor del armónico esférico
n=5;
m=5;
Y55 = compute_spherical_harmonic(n, m, theta, phi);
% Graficar superficie
figure(4)
surf(x, y, z, Y55, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet")
figure(5)
surf(x, y, z, Y+Y40+Y55, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet")

% PASO 4 Calcular el valor del armónico esférico
n=8;
m=4;
Y84 = compute_spherical_harmonic(n, m, theta, phi);
% Graficar superficie
figure(6)
surf(x, y, z, Y84, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet")
figure(7)
surf(x, y, z, Y+Y40+Y55+Y84, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
axis equal
colormap("jet") 