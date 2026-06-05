function [ errorInput ] = checkInput( inputGrid )
global numElec
%check input is greater than the amount of electrodes present
for ee = 1:numElec
    inputGrid(ee);
    elecsPresent = [1:numElec];
    z1 = ismember(inputGrid(ee),elecsPresent);
    if z1 == 0
       errorInput = 1;
       fprintf(2,'ERROR. The input of Electrode %d is incorrect (%d) because it is greater than the maximal amount of electrodes present. Please check and correct this value in the inputGrid variable.  \n \n',ee,inputGrid(ee))
    else
       errorInput = 0;
    end
end
%check for duplicate values
[n, bin] = histc(inputGrid, unique(inputGrid));
multiple = find(n > 1);
index    = find(ismember(bin, multiple));
if isempty(index) == 0
        fprintf(2,'ERROR. The input of Electrode %d contains a duplicate value. Please check and correct the value in the inputGrid variable.  \n \n',index)
        errorInput = 1;
end
    
end

