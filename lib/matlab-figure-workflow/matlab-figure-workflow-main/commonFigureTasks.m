% this should be a script

switch taskflag
%--------------------------------------------------------------------------
    case 'axes'
    try
        ha = gca; % get current axis handle
        ha.XAxis.Color = 'k'; % change the x axis color to black
        ha.YAxis.Color = 'k'; % change the y axis color to black
        ha.XAxis.FontSize = fonttick; % change x tick font size
        ha.YAxis.FontSize = fonttick; % change y tick font size
        ha.XAxis.Label.FontSize = fontlabel; % change x label font size
        ha.YAxis.Label.FontSize = fontlabel; % change y label font size
        ha.Layer = 'top'; % place the axes on top of the data
        ha.LineWidth = 1; % increase line width
    catch
        error("Needs (fonttick,fontlabel) variables defined.")
    end
%--------------------------------------------------------------------------
    case 'axes-z' % NEED
    try
        ha = gca; % get current axis handle
        ha.ZAxis.Color = 'k'; % change the x axis color to black
        ha.ZAxis.FontSize = fonttick; % change x tick font size
        ha.ZAxis.Label.FontSize = fontlabel; % change x label font size
    catch
        error("Needs (fonttick,fontlabel) variables defined.")
    end
%--------------------------------------------------------------------------
    case 'legend'
    try
        hl.FontSize = fontlegend; % change legend font size
        hl.EdgeColor = 'k'; % change the legend border to black
    catch
        error("Needs (fontlegend) variables defined.")
    end
%--------------------------------------------------------------------------
    case 'export'
    try
        pathpdf = mfoldername(mfilename('fullpath'),'pdf');
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)
    catch
        error("Needs (savename) variables defined.")
    end
%--------------------------------------------------------------------------
end