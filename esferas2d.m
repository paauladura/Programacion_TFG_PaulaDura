function spherical_harmonics_3D_interactive
    % Crear figura principal
    fig = figure('Position', [100, 100, 1200, 900], 'Name', 'Armónicos Esféricos 3D', 'NumberTitle', 'off');
    
    % Ejes 3D
    ax = axes('Position', [0.1, 0.25, 0.8, 0.7]);
    axis(ax, 'equal');
    axis(ax, 'off');
    hold(ax, 'on');
    view(ax, 3);
    
    % Controles deslizantes
    uicontrol('Style', 'text', 'Position', [50, 60, 100, 20], 'String', 'Grado (n):', 'FontSize', 11, 'FontWeight', 'bold');
    n_slider = uicontrol('Style', 'slider', 'Position', [160, 60, 350, 25],...
        'Min', 0, 'Max', 8, 'Value', 0, 'SliderStep', [1/8, 1/8],...
        'Callback', @update_plot);
    
    uicontrol('Style', 'text', 'Position', [50, 20, 100, 20], 'String', 'Orden (m):', 'FontSize', 11, 'FontWeight', 'bold');
    m_slider = uicontrol('Style', 'slider', 'Position', [160, 20, 350, 25],...
        'Min', 0, 'Max', 8, 'Value', 0, 'SliderStep', [1/8, 1/8],...
        'Callback', @update_plot);
    
    % Etiquetas de valores
    n_label = uicontrol('Style', 'text', 'Position', [520, 60, 80, 25], 'String', 'n = 0', 'FontSize', 11, 'BackgroundColor', [0.9 0.9 0.9]);
    m_label = uicontrol('Style', 'text', 'Position', [520, 20, 80, 25], 'String', 'm = 0', 'FontSize', 11, 'BackgroundColor', [0.9 0.9 0.9]);
    
    % Botón para alternar colormap
    uicontrol('Style', 'pushbutton', 'Position', [620, 60, 150, 25],...
        'String', 'Alternar Colormap', 'FontSize', 10, 'Callback', @toggle_colormap);
    
    % Información del patrón
    pattern_label = uicontrol('Style', 'text', 'Position', [620, 20, 350, 25],...
        'String', 'Tipo: Zonal | Polos: Norte=1, Sur=1 | Ceros θ: 0, Ceros λ: 0', ...
        'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]);
    
    % Variables de estado
    current_colormap = 1; % 1 = rojo/azul, 2 = rojo/negro
    pole_north = 1;
    pole_south = 1;
    zeros_theta = 0;
    zeros_lambda = 0;
    
    % Función para alternar colormap
    function toggle_colormap(~,~)
        current_colormap = mod(current_colormap, 2) + 1;
        update_plot();
    end

    % Función principal de actualización
    function update_plot(~,~)
        n = round(get(n_slider, 'Value'));
        m = round(get(m_slider, 'Value'));
        
        % Validar m <= n
        if m > n
            m = n;
            set(m_slider, 'Value', m);
        end
        
        set(n_label, 'String', sprintf('n = %d', n));
        set(m_label, 'String', sprintf('m = %d', m));
        
        % Determinar propiedades según n y m
        if m == 0
            pattern_type = 'Zonal';
            pole_north = 1;
            pole_south = (-1)^n;
            zeros_theta = n;
            zeros_lambda = 0;
        elseif m == n
            pattern_type = 'Sectorial';
            pole_north = 0;
            pole_south = 0;
            zeros_theta = 0;
            zeros_lambda = 2*m;
        else
            pattern_type = 'Teseral';
            pole_north = 0;
            pole_south = 0;
            zeros_theta = n - m;
            zeros_lambda = 2*m;
        end
        
        set(pattern_label, 'String', sprintf(...
            'Tipo: %s | Polos: Norte=%.1f, Sur=%.1f | Ceros θ: %d, Ceros λ: %d', ...
            pattern_type, pole_north, pole_south, zeros_theta, zeros_lambda));
        
        % Graficar el armónico esférico
        plot_spherical_harmonic(n, m);
    end

    % Función para graficar el armónico esférico en 3D
    function plot_spherical_harmonic(n, m)
        % Limpiar ejes
        cla(ax);
        
        % Crear rejilla esférica
        [theta, phi] = meshgrid(linspace(0, pi, 100), linspace(0, 2*pi, 100));
        
        % Calcular el valor del armónico esférico
        Y = compute_spherical_harmonic(n, m, theta, phi);
        
        % Convertir a coordenadas cartesianas
        r = 1; % Radio constante para esfera perfecta
        x = r * sin(theta) .* cos(phi);
        y = r * sin(theta) .* sin(phi);
        z = r * cos(theta);
        
        % Graficar superficie
        surf(ax, x, y, z, Y, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        
        % Configurar colormap
        if current_colormap == 1
            % Rojo (positivo) a azul (negativo)
            colormap(ax, [linspace(0, 0.2, 64)' linspace(0, 0.2, 64)' linspace(1, 0.8, 64)';  % Azules
                         linspace(1, 0.8, 64)' linspace(0.2, 0, 64)' linspace(0.2, 0, 64)']);  % Rojos
            caxis(ax, [-1 1]);
        else
            % Rojo (positivo) a negro (negativo)
            colormap(ax, [linspace(0, 0, 128)' linspace(0, 0, 128)' linspace(0, 0, 128)';  % Negros
                         linspace(1, 0.8, 128)' linspace(0, 0.2, 128)' linspace(0, 0.2, 128)']);  % Rojos
            caxis(ax, [-1 1]);
        end
        
        % Configuración de visualización
        axis(ax, 'equal');
        shading(ax, 'interp');
        material(ax, 'dull');
        
        % Iluminación
        light('Position', [1 1 1], 'Style', 'infinite');
        light('Position', [-1 -1 -1], 'Style', 'infinite');
        lighting(ax, 'gouraud');
        
        % Añadir líneas de referencia
        plot3(ax, [-1.2 1.2], [0 0], [0 0], 'k-', 'LineWidth', 0.5); % Eje X
        plot3(ax, [0 0], [-1.2 1.2], [0 0], 'k-', 'LineWidth', 0.5); % Eje Y
        plot3(ax, [0 0], [0 0], [-1.2 1.2], 'k-', 'LineWidth', 0.5); % Eje Z
        
        % Marcar polos
        plot3(ax, 0, 0, 1.05, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
        plot3(ax, 0, 0, -1.05, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
        text(ax, 0, 0, 1.15, 'N', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        text(ax, 0, 0, -1.15, 'S', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        
        % Título
        title(ax, sprintf('Armónico Esférico R_{%d,%d} = P_{%d,%d}·cos(%dλ)', n, m, n, m, m), ...
            'FontSize', 16, 'FontWeight', 'bold');
        
        % Añadir colorbar
        cb = colorbar(ax);
        cb.Position = [0.92, 0.25, 0.02, 0.5];
        cb.Label.String = 'Valor del Armónico';
        cb.Label.FontSize = 12;
    end

    % Función para calcular armónico esférico en coordenadas esféricas
    function Y = compute_spherical_harmonic(n, m, theta, phi)
        % Calcular polinomio asociado de Legendre
        P = legendre(n, cos(theta));
        
        if n == 0
            Plm = P;
        else
            Plm = squeeze(P(m+1, :, :));
        end
        
        % Factor de normalización
        norm_factor = sqrt((2*n+1)/(4*pi) * factorial(n-m)/factorial(n+m));
        
        % Aplicar la parte angular
        angular = cos(m * phi);
        
        % Combinar componentes
        Y = norm_factor * Plm .* angular;
        
        % Normalizar para mantener valores entre -1 y 1
        max_val = max(abs(Y(:)));
        if max_val > 0
            Y = Y / max_val;
        end
    end

    % Inicializar
    update_plot();
end