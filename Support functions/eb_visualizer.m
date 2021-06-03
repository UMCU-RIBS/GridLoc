function eb_visualizer(surface, surface_data, electrodes, electrode_data, opts)

% (c) Julia Berezutskaya, https://github.com/Immiora

if nargin < 1, error('No input'); end
if nargin < 2, surface_data = []; electrodes = []; electrode_data = []; end
if nargin < 3, electrodes = []; electrode_data = []; end
if nargin < 4, electrode_data = []; end
if nargin < 5, opts = struct; end

eb_visualizer_set_defauts(opts);

plot_surface    = 1;
plot_electrodes = 1;
plot_s_data     = 1;
plot_e_data     = 1;
e2s_data        = 0; 

if isempty(surface),        plot_surface = 0; end
if isempty(electrodes),     plot_electrodes = 0; end
if isempty(surface_data),   plot_s_data = 0; s_data_colorbar_show = 0; end
if isempty(electrode_data), plot_e_data = 0; e_data_colorbar_show = 0; end    

if ~isempty(surface)
    faces    = surface.tri;
    vertices = surface.vert;
else
    faces    = [];
    vertices = [];
end

els          = electrodes;


if plot_e_data == 1,
    if size(electrode_data, 2) == 1,   
        [e_data, e_data_norm] = eb_visualizer_x_data2triplets(electrode_data, e_data_cmap, e_data_span, e_data_thresh, 0.5);
    else
        e_data                = electrode_data;
    end
    e_data_nans = isnan(e_data_norm);

    switch e_data_plot_form
        case 'size'
            e_mkr_size = zeros(size(e_data_norm)); 
            e_mkr_size(~e_data_nans) = e_data_scale(e_data_norm(~e_data_nans));
            e_color = repmat(e_data_color_default, length(e_data), 1);

        case 'color'
            e_mkr_size = repmat(e_mkr_sise_default, 1, length(e_data));
            e_color = zeros(size(e_data)); 
            e_color(~e_data_nans, :) = e_data_cmap(e_data_norm(~e_data_nans), :);

        case 'size_color'
            e_mkr_size = zeros(size(e_data_norm));           
            e_mkr_size(~e_data_nans) = e_data_scale(e_data_norm(~e_data_nans));
            e_color = e_data; 
    end
end

if plot_s_data == 1,
    if size(surface_data, 1) == size(els, 1)
        s_data_span    = e_data_span;
        s_data_thresh  = e_data_thresh;
        s_data_cmap    = e_data_cmap;
        s_data         = eb_visualizer_e_data2s_data(vertices, els, surface_data, show_hemi_default);
        e2s_data = 1;       
    else
        s_data = surface_data;
        e2s_data = 0;
    end
    
    if size(s_data, 2) == 1,       
       e_data2s_data_cmap = eb_visualizer_make_colormap(['custom_', s_data_cmap], 64);
       s_data = eb_visualizer_x_data2triplets(s_data, e_data2s_data_cmap, ...
           s_data_span, s_data_thresh, surf_color_default(1)); % eb_visualizer_s_data2colors

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_surface == 1
P = patch('faces', faces, 'vertices', vertices,  ...
    'facecolor', surf_color_default, ...
    'edgecolor', 'none');
    if plot_s_data == 1
        set(P, 'facevertexCData', s_data);
        set(P, 'facecolor', 'interp');
    end
end

if plot_electrodes == 1  && plot_e_data ~= 1 && e2s_data ~= 1
    hold on
    e_data_color_default = [0 0 0];
    switch els_form
        case 'full'
            eM = '.';
            e_mkr_size = e_mkr_sise_default;
            plot3(els(:, 1) + els_bringup, els(:, 2), els(:, 3), 'marker', eM,'color', ...
                    e_data_color_default, 'linestyle', 'none', 'markersize', e_mkr_size + 5);
            plot3(els(:, 1) + els_bringup, els(:, 2), els(:, 3),'marker', eM, 'color', ...
                    e_data_color_default, 'linestyle', 'none', 'markersize', e_mkr_size);
        
        case 'empty'
            eM = 'o';
            e_mkr_size = e_mkr_sise_default * 0.2;
            plot3(els(:, 1) + els_bringup, els(:, 2), els(:, 3), 'marker', eM,'color', ...
                    e_data_color_default, 'linestyle', 'none', 'markersize', e_mkr_size, ...
                     'linewidth', e_data_line_width);
    end
end


view(show_view_pos)

if plot_electrodes == 1 && plot_e_data == 1   
    for i = 1:size(els, 1)
        hold on
        if any(e_data(i, :) ~= 0.5) % ~isnan(e_data(i)) && e_data(i) ~= 0 
            
            plot3(els(i, 1) + els_bringup, els(i, 2), els(i, 3),'.','color', ...
                [0 0 0],'markersize', e_mkr_size(i) + 5);
            plot3(els(i, 1) + els_bringup, els(i, 2), els(i, 3),'.','color', ...
                e_color(i, :),'markersize', e_mkr_size(i));
            
            if e_data_add_gloss == 1
                plot3(els(i, 1) + els_bringup, els(i, 2), els(i, 3), ...
                    'marker', '.', 'color', [0.6 0.8 1], ...
                    'linestyle', 'none', ...
                    'markersize', e_mkr_sise_default * 0.2);
            end
            
        elseif e_data_nodata_show == 1
             plot3(els(i, 1) + els_bringup, els(i, 2), els(i, 3), ...
                 'marker', e_data_nodata_style_char, ...
                 'color', [0.1 0.1 0.1], ...
                 'markersize', e_mkr_sise_default * 0.2, ...
                 'linewidth', e_data_line_width);
        end
    end
end


l = light;
lighting gouraud;
if show_view_pos(1)     > 200,
    set(l,'position',[-1 0 1]) % [-1 0 1] for left
elseif show_view_pos(1) < 200,
    set(l,'position',[1 0 1]) % [-1 0 1] for left
end
camlight('headlight','infinite')
material dull; 
axis off; axis equal; axis tight
hold on

if plot_e_data == 1 && e_data_colorbar_show == 1,
     colormap(e_data_cmap)
     % set(gca, 'clim', e_data_span)
     caxis(e_data_span)
     e_data_cbar = colorbar('east', 'fontsize', 10, 'axislocation', 'out'); 
     if ~isempty(e_data_colorbar_label)
        e_data_cbar.Label.String = e_data_colorbar_label;
        e_data_cbar.Label.FontSize = 16;
     end
end


if plot_s_data == 1 && s_data_colorbar_show == 1,
     s_data_cbar = newcolorbar('west', 'fontsize', 10, 'axislocation', 'in');     
     colormap(gca, s_data_cmap)

     caxis(s_data_span)
     if ~isempty(s_data_colorbar_label)
        s_data_cbar.Label.String = s_data_colorbar_label;
        s_data_cbar.Label.FontSize = 16;
     end
end








