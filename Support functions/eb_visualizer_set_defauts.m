function opts = eb_visualizer_set_defauts(opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% general %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts, 'fig_handle'),
    opts.fig_handle = figure('units', 'normalized', 'outerposition', [0 0 1 1]); end

if ~isfield(opts, 'show_hemi_default'), 
    opts.show_hemi_default      = 'l'; end

if ~isfield(opts, 'show_view_pos'), 
    if strcmp(opts.show_hemi_default, 'l')
        opts.show_view_pos      = [270, 0]; 
    elseif strcmp(opts.show_hemi_default, 'r')
        opts.show_view_pos      = [90, 0];
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts, 'surf_color_default'), 
    opts.surf_color_default    = [0.5 0.5 0.5]; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% surface data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts, 's_data_cmap'), 
    opts.s_data_cmap           = colormap('jet'); end

if ~isfield(opts, 's_data_span'), 
    opts.s_data_span           = [-10 10]; end 

if ~isfield(opts, 's_data_thresh'), 
    opts.s_data_thresh         = [-6 6]; end

if ~isfield(opts, 's_data_colorbar_show'), 
    opts.s_data_colorbar_show  = 0; end

if ~isfield(opts, 's_data_colorbar_label'), 
    opts.s_data_colorbar_label  = ''; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts, 'els_bringup'), 
    if strcmp(opts.show_hemi_default, 'l')
        opts.els_bringup       = -7; 
    elseif strcmp(opts.show_hemi_default, 'r')
        opts.els_bringup       = 7; 
    end
end

if ~isfield(opts, 'els_form'), 
    opts.els_form              = 'empty'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% electrode data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts, 'e_mkr_sise_default'), 
    opts.e_mkr_sise_default    = 40; end

if ~isfield(opts, 'e_data_span'), 
    opts.e_data_span           = [-1 1]; end

if ~isfield(opts, 'e_data_thresh'), 
    opts.e_data_thresh         = [0 0]; end

if ~isfield(opts, 'e_data_cmap_rng'),
    opts.e_data_cmap_rng       = 64; end

if ~isfield(opts, 'e_data_scale_step'),
    opts.e_data_scale_step     = 2; end

if ~isfield(opts, 'e_data_color_default'),
    opts.e_data_color_default  = [0 1 0]; end  

if ~isfield(opts, 'e_data_cmap'),
    opts.e_data_cmap           = colormap('hot'); end%eb_visualizer_make_colormap('w2g', opts.e_data_cmap_rng); end

if ~isfield(opts, 'e_data_scale'),
    opts.e_data_scale          = round(linspace(0.5, opts.e_data_scale_step, opts.e_data_cmap_rng) * opts.e_mkr_sise_default); end

if ~isfield(opts, 'e_data_plot_form'),
    opts.e_data_plot_form      = 'size_color'; end

if ~isfield(opts, 'e_data_line_width'), 
    opts.e_data_line_width     = 2; end

if ~isfield(opts, 'e_data_nodata_style'), 
    opts.e_data_nodata_style  = 'empty'; end

if ~isfield(opts, 'e_data_nodata_show'), 
    opts.e_data_nodata_show  = 1; end

if ~isfield(opts, 'e_data_nodata_style_char'), 
    if strcmp(opts.e_data_nodata_style, 'empty');
        opts.e_data_nodata_style_char  = 'o';
    elseif strcmp(opts.e_data_nodata_style, 'full');
        opts.e_data_nodata_style_char  = '.';
    end
end

if ~isfield(opts, 'e_data_add_gloss'), 
    opts.e_data_add_gloss = 0; end

if ~isfield(opts, 'e_data_colorbar_show'), 
    opts.e_data_colorbar_show  = 0; end

if ~isfield(opts, 'e_data_colorbar_label'), 
    opts.e_data_colorbar_label  = ''; end




fields = fieldnames(opts);
for o = 1:length(fields)
    assignin('caller', fields{o}, opts.(fields{o}))
end