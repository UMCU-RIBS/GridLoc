function [ ind ] = indexFuncLegacy( subj_info )
%Get the correct indices for the grids based upon their orientation and
%shape
%   Detailed explanation goes here
global dims
global numElec
global inputGrid

side = subj_info.hemiVect.side;
hemi = subj_info.hemiVect.hemi;

baseGrid=reshape(inputGrid,dims(2),dims(1));
switch hemi
    case 'l'
        switch side
            case 'u'
                ind = inputGrid;
            case 'd'
                ind = inputGrid(numElec:-1:1);
            case 'l'
                indGrid=rot90(baseGrid,1);
                ind = indGrid(:);
            case 'r'
                indGrid=rot90(baseGrid,-1);
                ind = indGrid(:);
        end
    case 'r'
        switch side
            case 'u'
                indGrid=fliplr(baseGrid);
                ind = indGrid(:);
            case 'd'
                indGrid=fliplr(baseGrid);
                auxInd=indGrid(:);
                ind = auxInd(numElec:-1:1);
            case 'r'
                indGrid=fliplr(baseGrid);
                auxIndGrid=rot90(indGrid,1);
                ind = auxIndGrid(:);
            case 'l'
                indGrid=fliplr(baseGrid);
                auxIndGrid=rot90(indGrid,-1);
                ind = auxIndGrid(:);
        end
end
end



