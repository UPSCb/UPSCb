function pltstmp(FunctionName)
%% PLTSTMP - plot stample with date and name of function in the bottom of figure
% Function Plots additional pannel on figure. Pannel contains information
% about subject, digit, function name and time 
% Original file source:
% https://se.mathworks.com/matlabcentral/fileexchange/4999-pltstmp
h=gcf;
fn = sprintf('Created by :: %s :: %s', FunctionName,datestr(now,'yyyy-mm-dd HH:MM:SS'));

in = 'tex';
% create axes for text string and add string to the bottom of the
% figure
a = axes('parent',h,...
    'posi',[-0 -0 1 1],...
    'ydir','norm','xdir','norm',...
    'visible','off',...
    'Handlevisibility','off',...
    'Tag','PLTSTMP-AXES');
text(.5,0.01,fn,'hor','center','vert','bot','fontsize',9,...
     'FontWeight', 'bold','interpre',in,'parent',a);
return