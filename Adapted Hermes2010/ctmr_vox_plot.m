function c = ctmr_vox_plot(cortex,xyz,weights,ssize,v,cm,varargin)

% ctmr_vox_plot(xyz,weights,ssize)
% xyz and weights output from voxplot_func_gm
% ssize=2; % size of voxel squares on surface
% v='l'; % 'l' or 'r' left or right
% functions has a popup to load cortex
%
%     Copyright (C) 2009  D. Hermes & K.J. Miller, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% edits:
% dh 2009: ctmr_gauss_plot function adjusted for fMRI maps with square voxels

% adapted by MPBranco, Nov 2020:
% - added first input variable to function
% - added cm as input

%%%OLD CODE:
%[cname]=spm_select(1,'mat','select cortex.mat');
%%%
if isempty(varargin{1})
    facealpha = 1;
else
    facealpha = varargin{1};
end

brain=cortex.vert;

if length(weights)~=length(xyz(:,1))
    error('you sent a different number of weights than xyz (perhaps a whole matrix instead of vector)')
end

c=zeros(length(cortex.vert),1);
for k=1:length(xyz(:,1))
    b_z=abs(brain(:,3)-xyz(k,3));
    b_y=abs(brain(:,2)-xyz(k,2));
    b_x=abs(brain(:,1)-xyz(k,1));
    
    d=b_z<ssize & b_y<ssize & b_x<ssize;
    d=d*weights(k); % no smoothing
    
    c=max(c,d); %overlap is going to maximum
end

% a=tripatch(cortex, '', c);
% shading interp;
% a=get(gca);
% %%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
% %%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
% d=a.CLim;
% set(gca,'CLim',[0 max(weights)]);
% 
% l=light;
% 
% %%% NEW CODE
% %load in colormap
% if isempty(cm)
%     cm=colormap('hot');
%     cm=[.7 .7 .7;cm];
% else
%     cm=[.7 .7 .7;cm];
% end
% %%%%
% colormap(cm);
% lighting gouraud; %play with lighting...
% % material dull;
% material([.3 .8 .1 10 1]);
% axis off
% set(gcf,'Renderer', 'zbuffer')
% 
% if v=='l'
%     view(270, 0);
%     %view(-92,32);
%     set(l,'Position',[-1 0 1])
% elseif v=='r'
%     view(90, 0);
%     set(l,'Position',[1 0 1])
% end
end