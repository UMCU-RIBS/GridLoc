function cmap = eb_visualizer_make_colormap(mapstr, lim, base_col)
% (c) Julia Berezutskaya, https://github.com/Immiora

if nargin < 2, lim = 64; end
if nargin < 3, base_col = [0.5 0.5 0.5]; end

switch mapstr
    case 'w2g'
        g = ones(lim, 1);
        r = linspace(0.7, 0, lim)'; % linspace(1, 0, lim)';
        b = linspace(0.7, 0, lim)';
        cmap = [r, g, b];
    case 'w2r'
        g = linspace(0.7, 0, lim)';
        r = ones(lim, 1);
        b = linspace(0.7, 0, lim)';
        cmap = [r, g, b];
    case 'w2b'
        g = linspace(0.7, 0, lim)';
        r = linspace(0.7, 0, lim)';
        b = ones(lim, 1);
        cmap = [r, g, b];
    case 'w2m'
        g = linspace(0.7, 0, lim)';
        r = ones(lim, 1);
        b = ones(lim, 1);
        cmap = [r, g, b];
    case 'w2c'
        g = ones(lim, 1);
        r = linspace(0.7, 0, lim)';
        b = ones(lim, 1);
        cmap = [r, g, b];
    case 'w2y'
        g = ones(lim, 1);
        r = ones(lim, 1);
        b = linspace(0.7, 0, lim)';
        cmap = [r, g, b];
    case 'custom_jet'
        m = colormap(jet(64));
        m(26:lim/2-1, 1) = linspace(0, base_col(1), 6); 
        m(26:lim/2-1, 2) = linspace(1, base_col(2), 6);
        m(26:lim/2-1, 3) = linspace(1, base_col(3), 6);

        m(lim/2, :)     = 0.5;
        m(lim/2 + 1, :) = 0.5;

        m(lim/2 + 2:39, 1) = linspace(base_col(1), 1, 6); 
        m(lim/2 + 2:39, 2) = linspace(base_col(2), 1, 6);
        m(lim/2 + 2:39, 3) = linspace(base_col(3), 0, 6);
        cmap = m;
    case 'custom_hot'
        m = colormap(hot(64));
        m(1, :)     = base_col;

        m(2:4, 1) = linspace(base_col(1), 0.2, 3); 
        m(2:4, 2) = linspace(base_col(2), 0, 3);
        m(2:4, 3) = linspace(base_col(3), 0, 3);
        cmap = m;

    case 'custom_red'
        m = zeros(64, 3);
        m(1:end, 1) = linspace(0.3, 0.5, 64);
        m(1:end, 2:3) = zeros(64, 2);
        cmap = m;
    case 'gr2g'
        g = linspace(0.7, 1, lim)';
        r = linspace(0.7, 0.3, lim)';
        b = linspace(0.7, 0.3, lim)';
        cmap = [r, g, b];
    case 'gr2r'
        g = linspace(0.7, 0.3, lim)';
        r = linspace(0.7, 1, lim)';
        b = linspace(0.7, 0.3, lim)';
        cmap = [r, g, b];
    case 'cold'
        cmap = colormap('hot');
        cmap = cmap(:, [3 2 1]);
    case 'blues'
        r = linspace(0, 0.7, lim)';
        g = linspace(0, 0.9, lim)';
        b = linspace(0.6, 1, lim)';
        cmap = [r, g, b];
    case 'pinks'
        r = linspace(0.6, 1, lim)';
        g = linspace(0, 0.8, lim)';
        b = linspace(0.4, 0.8, lim)';
        cmap = [r, g, b];
    case 'yels'
        r = linspace(0.8, 1, lim)';
        g = linspace(0.2, 1, lim)';
        b = linspace(0, 0.2, lim)';
        cmap = [r, g, b];
    case 'greens'
        r = linspace(0, 0.8, lim)';
        g = linspace(0.3, 1, lim)';
        b = linspace(0, 0.8, lim)';
        cmap = [r, g, b];
        
    case 'custom_cool'
        m = colormap(cool(64));
        m(lim/2, :)     = 0.5;
        m(lim/2 + 1, :) = 0.5;
        cmap = m;
        
    case 'custom_winter'
        m = colormap(winter(64));
        m(lim/2, :)     = 0.5;
        m(lim/2 + 1, :) = 0.5;
        cmap = m;
    case 'custom_red_blue'
        r = linspace(0, 1, lim)';
        g = linspace(0, 0, lim)';
        b = linspace(1, 0, lim)';
        cmap = [r, g, b];
        cmap(lim/2, :)     = 0.5;
        cmap(lim/2 + 1, :) = 0.5;
        cmap(12:lim/2-1, 1) = linspace(0, 0.5, 20); 
        cmap(12:lim/2-1, 2) = linspace(0, 0.5, 20);
        cmap(12:lim/2-1, 3) = linspace(0.85, 0.5, 20);

        cmap(lim/2 + 2:53, 1) = linspace(0.5, 0.85, 20); 
        cmap(lim/2 + 2:53, 2) = linspace(0.5, 0, 20);
        cmap(lim/2 + 2:53, 3) = linspace(0.5, 0, 20);        
end

