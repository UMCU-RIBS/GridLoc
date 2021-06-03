function [ forward_ind, inverse_ind ] = indexFuncLegacy( subj_info )
%Get the correct indices for the grids based upon their orientation and
%shape.
% forward_ind - GRID --> GRIDLOC: indexes taht map the grid layout into
%                       the (simulated) electrodes coordinates order.
% inverse_ind - GRIDLOC --> GRID: indexes taht map the (simulated)
%                       electrode coordinates into the grid layout

global dims
global numElec
global inputGrid

side = subj_info.hemiVect.side;
hemi = subj_info.hemiVect.hemi;

%rehsape the grid layout input by user into a NxM matrix.
baseGrid = reshape(inputGrid,dims(2),dims(1));

%For each hemisphere and lead side there is a mapping
switch hemi
    
    case 'l'
        
        switch side
            
            case 'u'
                %forward
                forward_ind      = inputGrid;
                
                %inverse;
                baseGrid         = reshape(forward_ind,dims(2),dims(1));
                [~,inverse_ind]  = sort(baseGrid(:));
                
            case 'd'
                %forward
                forward_ind      = inputGrid(numElec:-1:1);
                
                %inverse
                baseGrid         = reshape(forward_ind,dims(2),dims(1));
                [~,inverse_ind]  = sort(baseGrid(:));
                
            case 'l'
                %forward
                indGrid          = rot90(baseGrid,1);
                forward_ind      = indGrid(:);
                
                %inverse
                indx             = flipud(rot90(reshape(inputGrid,dims(2),dims(1))',2));
                indx             = indx(:);
                baseGrid         = rot90(reshape(indx,dims(2),dims(1)),2);
                [~,inverse_ind]  = sort(baseGrid(:));
                
            case 'r'
                %forward
                indGrid          = rot90(baseGrid,-1);
                forward_ind      = indGrid(:);
                
                %inverse
                indx             = flipud(reshape(inputGrid,dims(2),dims(1))');
                indx             = indx(:);
                baseGrid         = rot90(reshape(indx,dims(2),dims(1)),2);
                [~,inverse_ind]  = sort(baseGrid(:));
        end
        
    case 'r'
        
        switch side
            
            case 'u' 
                %forward
                indGrid             = fliplr(baseGrid);
                forward_ind         = indGrid(:);
                
                %inverse
                baseGrid            = reshape(forward_ind,dims(2),dims(1));
                [~,inverse_ind]     = sort(baseGrid(:));
                
            case 'd'
                %forward
                indGrid             = fliplr(baseGrid);
                auxInd              = indGrid(:);
                forward_ind         = auxInd(numElec:-1:1);
                
                %inverse
                baseGrid            = reshape(forward_ind,dims(2),dims(1));
                [~,inverse_ind]     = sort(baseGrid(:));
                
            case 'r'
                %forward
                indGrid             = fliplr(baseGrid);
                auxIndGrid          = rot90(indGrid,1);
                forward_ind         = auxIndGrid(:);
                
                %inverse
                indx                = reshape(inputGrid,dims(2),dims(1))';
                indx                = indx(:);
                baseGrid            = reshape(indx,dims(2),dims(1));
                [~,inverse_ind]     = sort(baseGrid(:));
                
            case 'l'
                %forward
                indGrid             = fliplr(baseGrid);
                auxIndGrid          = rot90(indGrid,-1);
                forward_ind         = auxIndGrid(:);
                
                %inverse
                indx                = rot90(reshape(inputGrid,dims(2),dims(1))',2);
                indx                = indx(:);
                baseGrid            = reshape(indx,dims(2),dims(1));
                [~,inverse_ind]     = sort(baseGrid(:));
        end
end

end


