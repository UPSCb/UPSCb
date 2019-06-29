% calculatecurvature(x,y)
% 2018-06-21 
% Script loads result structure and calculates curvature for each day
% using calccurvature.m
% Kamil @Umeå University base on Matthew Zinkgraf R script.

function CurvatureFinalTable = calculatecurv_final(meta, result, nPoints)
%% calculatecurv_final - calculate curvature for all trees within one experiment
% INPUTS:
% meta struct with metadata:
% Contains meta information about each tree line, example:
%          FolderName             genotype       tree     data     treatment     hook     length
%     ________________________    _____________    ____    ______    _________    ______    ______
% 
%     'EXP5-T89_1_1-28_H2O'       'EXP5-T89'       '1'     '1-28'      'H2O'      69.229    121.24
%     'EXP5-T89_2_1-28_H2O'       'EXP5-T89'       '2'     '1-28'      'H2O'       68.21    113.24
%     'EXP5-T89_5_1-28_H2O'       'EXP5-T89'       '5'     '1-28'      'H2O'      71.679    120.11
%     'EXP5-T89_6_1-28_H2O'       'EXP5-T89'       '6'     '1-28'      'H2O'       69.51    118.03
%     'EXP5-T89_7_1-28_H2O'       'EXP5-T89'       '7'     '1-28'      'H2O'      61.191    112.93
%     'EXP5-T89_9_1-28_H2O'       'EXP5-T89'       '9'     '1-28'      'H2O'      79.294    125.16
%     'EXP5-etr1L6_2_1-28_H2O'    'EXP5-etr1L6'    '2'     '1-28'      'H2O'      64.329    109.69
%     'EXP5-etr1L6_3_1-28_H2O'    'EXP5-etr1L6'    '3'     '1-28'      'H2O'      63.171    122.59
%     'EXP5-etr1L6_4_1-28_H2O'    'EXP5-etr1L6'    '4'     '1-28'      'H2O'      70.974    124.06
%     'EXP5-etr1L6_5_1-28_H2O'    'EXP5-etr1L6'    '5'     '1-28'      'H2O'      69.494    129.24
%     'EXP5-etr1L6_6_1-28_H2O'    'EXP5-etr1L6'    '6'     '1-28'      'H2O'      71.888     144.4
%     'EXP5-etr1L6_7_1-28_H2O'    'EXP5-etr1L6'    '7'     '1-28'      'H2O'      67.168    137.93
%
% result Contains x,y, and angles data for all trees, example:
% result = 
% 
%   struct with fields:
% 
%        EXP5_T89_1_1_28_H2O: [1×1 struct]
%        EXP5_T89_2_1_28_H2O: [1×1 struct]
%        EXP5_T89_5_1_28_H2O: [1×1 struct]
%        EXP5_T89_6_1_28_H2O: [1×1 struct]
%        EXP5_T89_7_1_28_H2O: [1×1 struct]
%        EXP5_T89_9_1_28_H2O: [1×1 struct]
%     EXP5_etr1L6_2_1_28_H2O: [1×1 struct]
%     EXP5_etr1L6_3_1_28_H2O: [1×1 struct]
%     EXP5_etr1L6_4_1_28_H2O: [1×1 struct]
%     EXP5_etr1L6_5_1_28_H2O: [1×1 struct]
%     EXP5_etr1L6_6_1_28_H2O: [1×1 struct]
%     EXP5_etr1L6_7_1_28_H2O: [1×1 struct]
% Each tree has 3 fields, for example:
% >> result.EXP4_139L4_3_1_28_H2O
%   struct with fields: 
%     output: [28×4 table] - contains curvature and lift for each day
%        out: {14140×3 cell} - contains XY position for each measuring
%        point (several during one day)
%        vec: [14084×5 table] - contains difference between x,y and angle
%        at each measuring time points.
% nPoints - how many points from each day are used for calculation
% curvature, example: 40

FolderNames      = meta.FolderName;
ResultStructName = meta.FolderName;
% Adjust Struct Names - it can't have '-'
for i = 1:numel(FolderNames)
  ResultStructName{i} = strrep(FolderNames{i}, '-', '_');
end

CurvatureFinalTable = table(meta.FolderName, ResultStructName, meta.hook,cell(size(meta.hook)), cell(size(meta.hook)),...
    'VariableNames', {'MetaFolderName', 'ResultFieldName', 'Hook', 'Curvature', 'Lift'});
  
 %% Main Loop
for iFolder = 1:numel(ResultStructName)
  fig = figure('Position', [969 49 944 948]);
  cData = result.(CurvatureFinalTable.ResultFieldName{iFolder}).out;
  CurvLift = result.(CurvatureFinalTable.ResultFieldName{iFolder}).output;
  cData = cell2table(cData, 'VariableNames', {'Dir', 'X', 'Y'});
  G = findgroups(cData.Dir);
  
  Hook = CurvatureFinalTable.Hook(iFolder);
  CurvatureFinalTable.Lift{iFolder} = CurvLift.lift;
  
  for i = 1:28 % 28 days
    X = cData.X(G==i);
    Y = cData.Y(G==i);
    X = X-X(1); % Normalize to 0
    Y = Y-Y(1); % Normalize to 0
    Y = -Y; % reverse Y axis
    X = X(round(linspace(1,numel(X), nPoints))); % take NPOINTS samples
    Y = Y(round(linspace(1,numel(Y), nPoints))); % take NPOINTS samples
    curvature(i) = calccurvature(X,Y,Hook); %#ok<*AGROW>
    subplot(6,5,i), plot(X,Y, '.'), axis tight
    
    XL(:,i) = xlim;
    YL(:,i) = ylim;
    title({sprintf('Curv: %1.1f',curvature(i)),...
        sprintf('Lift: %1.1f', CurvLift.lift(i))})
  end
  
  %% Fix axes
  xl = [min(XL(1,:)) max(XL(2,:))];
  yl = [min(YL(1,:)) max(YL(2,:))];
  % 28 days
  for i = 1:28, subplot(6,5,i), xlim(xl), ylim(yl), end
  subplot(6,5,26),  xlabel('cm'),  ylabel('cm')
  fname = CurvatureFinalTable.ResultFieldName{iFolder};
  subplot(6,5,28), XLls = xlim;
  text(XLls(2)+20,0,fname, 'FontSize', 14, 'Interpreter', 'none')
 
  pltstmp('calculatecurv final ')
  saveas(fig,['Results\' fname '10sampl_fixedaxes.png'])
  
  CurvatureFinalTable.Curvature{iFolder} = curvature;
end
