function Createplots(pcase)


%--------------------------------------------------------------------------
% Directory settings
%--------------------------------------------------------------------------

filepath = strcat(pwd,'\Figures');
addpath(filepath)

switch pcase
    case 0

        CreateFig1C0
        CreateFig2C0

    case 1

        CreateFig1C1 
        CreateFig2C1 

    case 2

        CreateFig1C2 
        CreateFig2C2 

    case 3

        CreateFig1C3 
        CreateFig2C3

    case 4

        CreateFig1C4 
        CreateFig2C4 

    case 5

        CreateFig1C5 
        CreateFig2C5 

    case 6

        CreateFig1C6

    case 7

        CreateFig1C7  


end

close all

end