function [x_data_triplets, x_data_norm] = eb_visualizer_x_data2triplets(x_data, cmap, rng, th, base_col)
% (c) Julia Berezutskaya, https://github.com/Immiora

if ischar(cmap), cmap = colormap(cmap); end

indices_x_colors= x_data < th(1) | x_data > th(2);
x_data_triplets = zeros(size(x_data, 1), 3);
d               = x_data(indices_x_colors);
s_data_range    = range(rng);
s_data_min      = rng(1);

x_data_norm     = round((d - s_data_min) ./ s_data_range * 63 + 1);
x_colors        = ones(size(d, 1), 3) * base_col;
x_data_norm(x_data_norm > 64) = 64;
x_data_norm(x_data_norm < 1) = 1;

for i =  1:size(cmap, 1)
    temp = size(x_colors(x_data_norm == i, :), 1);
    x_colors(x_data_norm == i, :) = repmat(cmap(i, :), temp, 1);
end

x_data_triplets(indices_x_colors, :) = x_colors;
x_data_triplets(~indices_x_colors, :) = base_col;

x_data_norm     = round((x_data - s_data_min) ./ s_data_range * 63 + 1);
x_data_norm(x_data_norm > 64) = 64;
x_data_norm(x_data_norm < 1) = 1;
