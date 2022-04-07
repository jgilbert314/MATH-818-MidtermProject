function [  ] = print2pdf(F, filename)
    
    set(F,'Units','Inches');
    pos = get(F,'Position');
    set(F,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print([filename, '.pdf'], '-dpdf', '-fillpage');
    
end