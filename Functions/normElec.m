function [normal] = normElec( M, coords_Closest, normway)
%NORMELEC       Calculate the normals from a given point on a surface model.
%
%   Calculate the average normals for every point in the struct coords_Closest 
%   in an area specified by normdist. Normway = 25.
%   
%   CALLING SEQUENCE:
%   [ subjstructs ] = projectElectrodes( M, subjstructs, normway, interstype, intersval)
%
%   [ subjstructs ] = normElec( M, subjstructs, normway )
%
%   INPUT:
%       M:              struct('vert', Vx3matrix, 'tri', Tx3matrix) - brain model 
%       subjstructs:    field of structures, for each subject: struct('electrodes', Nsubjx3matrix)
%       normway:        controls the way the normal vector to the surface is computed
%           possible values: [x, y, z] - every electrode is projected onto the surface using the normal vector line [electrode - [x, y, z]]
%                            normdist (double) - only those vertices vert whose distance(vert, electrode) < normdist are used to compute the average normal vector from the normals at these vertices
%   OUTPUT:
%       %The average normal of all the normals in an area around the ROI
%       point specified by normdist.
%
%   Example:
%       [ subjstructs ] = normElec( M, subjstructs, 25)
%           projection using the surrounding normal vectors; most commonly used for models
%           that contain low number of vertices (< 10000) - every brain
%           model can be reduced to about this number without affecting the
%           accuracy of the projection - and this is normally done (see
%           also demo.m)

%   Author: Jan Kubanek, edited for normal calculation by Michael Leibbrand
%   Institution: Wadsworth Center, NYSDOH, Albany, NY; UMC Utrecht,
%   Utrecht, Utrecht. 
%   Date: August 2005; August 2016


%THE LOOP FOR ALL ELECTRODES OF ALL SUBJETS-------------
%one subject per time
Ss = 1;
Vv = length(M.vert);
Tt = length(M.tri);

normal = zeros(length(coords_Closest),3);


for subj = 1 : Ss,
    Ee = size(coords_Closest, 1);
   
   
   for eg = 1 : Ee,
     
%compute the distance of all vertices from the electrode being processed
       dvect = zeros(1, Vv);
       vert = M.vert; %reallocate to speed up the loop
       electrode = coords_Closest(eg, :); %reallocate to speed up the loop
       for v = 1 : Vv,
           %compute the distance ||eg - v||^2 and store into a vector dvect
           delta = electrode - vert(v, :);
           dvect(v) = delta * delta';
       end


%NORMAL VECTOR PART------------------------------------
%normal vector computation (computation of line points p1 and p2):

       if length(normway) == 1, %the normdist is specified, compute average normals from near vertices
           normdist = normway;
           closevert = find(dvect < normdist^2);
           if isempty(closevert),
               error('The normdist distance from electrode %d of subject %d to compute the normal vector is too large', eg, subj);
           end
       
       %compute the normal vector
           %surf the brain so that the normals get computed:
           if verLessThan('matlab', '8.4'), %<= R2014a
               h = trisurf(M.tri, M.vert(:, 1), M.vert(:, 2), M.vert(:, 3), 'FaceVertexCData', 1, 'LineStyle', 'none');
               normals = get(h, 'VertexNormals');
               close(gcf);
           else %>= R2014b
               tr = triangulation(M.tri, M.vert);
               [~, MSGID] = lastwarn();
               warning('off', MSGID);
               normals = vertexNormal(tr);
           end
           normals2av = normals(closevert, :);
           [row, col] = find(isnan(normals2av)); %find NaN values
           normals2av(row, :) = []; %and remove them
           if length(normals2av) < 1,
               error('There is less than one normal vector used for the normal vector averaging, please increase the parameter normway to more than %d', normway);
               normal(Ss) = [NaN, NaN, NaN];
           else
               normal(eg,:) = mean(normals2av, 1);
           end
           p1 = coords_Closest(eg, :);
           p2 = p1 + normal(eg,:);
       elseif length(normway) == 3,
           pointorig = normway(:)';
           p1 = coords_Closest(eg, :);
           p2 = pointorig;
       else
           error('The attribute normway is specified incorrectly.');
       end


   end
end