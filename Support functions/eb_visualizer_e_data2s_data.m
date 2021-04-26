function s_data = eb_visualizer_e_data2s_data(vertices, els, edata, hemi)
% (c) Julia Berezutskaya, https://github.com/Immiora

len     = size(vertices, 1);
s_data  = zeros(len, 1);
fade    = 10;

nonanind = find(~isnan(edata));

for el = nonanind';
    b_z = abs(vertices(:,3)-els(el,3));
    b_y = abs(vertices(:,2)-els(el,2));
    if strcmp(hemi, 'l'), 
        b_x = abs(vertices(:,1) - els(el,1));  
    elseif strcmp(hemi, 'r'), 
        b_x = abs(els(el,1) - vertices(:,1)); %put everything on left hemisphere
    end
    [~, im] = min(b_x.^2+b_z.^2+b_y.^2);

    b_z = abs(vertices(:,3)-vertices(im,3));
    b_y = abs(vertices(:,2)-vertices(im,2));
    if strcmp(hemi, 'l'), 
        b_x = abs(vertices(:,1) - vertices(im,1));  
    elseif strcmp(hemi, 'r'), 
        b_x = abs(vertices(im,1) - vertices(:,1)); %put everything on left hemisphere
    end
    d = edata(el) * exp((-(b_x.^2+b_z.^2+b_y.^2))/(fade)); %exponential fall off % electrode_weight is a measure: cor, power, lag or whatever else
    s_data = s_data + d;
end

