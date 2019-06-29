% StemCurve.R
% Script rewritted from Matthew Zinkgraf R script:
% 2018-05 Kamil Antos, Umeå Uniersity

clear all

directory = 'DATA\EXP5 Full\';
d = dir([directory 'EXP5*']);
%% Create MetaData table:
disp('Creating meta table...')
for i = 1:numel(d)
fName = d(i).name;
qq = strfind(fName, '_');
FolderName{i} = fName;    %#ok<*SAGROW>
genotype{i}    = fName(1:qq(1)-1);
tree{i}        = fName(qq(1)+1:qq(2)-1);
data{i}        = fName(qq(2)+1:qq(3)-1);
treatment{i}   = fName(qq(3)+1:end);
end

meta = table(FolderName', genotype', tree', data',treatment', 'VariableNames',...
  {'FolderName', 'genotype', 'tree','data', 'treatment'});
meta.hook = nan(size(meta,1), 1);
meta.length = nan(size(meta,1), 1);
clear fName qq genotype tree data treatment FolderName i d 
% directory and meta stay at this point
%% Load Distance to hook data
% R:
% hook<-read.table("Data/distance.txt",sep="\t",header=T,stringsAsFactors = F)
% meta<-merge(meta,hook,by.x="dir",by.y="dir")
distanceTable = readtable([directory 'distance.txt']);
distanceTable.dir = categorical(distanceTable.dir);
% Find the same folder names in distance table and copy hook and length
for i = 1:size(meta,1)
  meta.hook(i)  = distanceTable.hook(  distanceTable.dir == meta.FolderName(i));
  meta.length(i) = distanceTable.length(distanceTable.dir == meta.FolderName(i));
end

disp(meta)
save('meta5.mat','meta')
disp('Meta table (meta5.mat) created')

%% Apply stemcurve function
disp('Applying stemcurve function..')
grep_key = '.txt';
for i = 1:size(meta,1)
  disp(meta.FolderName(i))
  fname(i) = strrep(meta.FolderName(i), '-', '_');
  temp = fname(i);
  result.(temp{:}) = stemcurve_fun([directory meta.FolderName{i,1}],grep_key,meta.hook(i), meta.length(i)) ;
end
disp(result)
save('result5.mat','result')
disp('Result saved at result5.mat')
%% Plot Bending Analysis
%  pdf(file="Data/Bending_analysis.pdf")
disp('Plotting Bending Analysis...')
fig = figure('Position', [9 49 944 948]);
for i = 1:size(meta,1)
  c = cell2mat(result.(fname{i}).out(:,2:3));
  subplot(6,2,i)
  CM = colormap;
  col = randi(size(CM,1));
  plot(c(:,1), c(:,2), 'ro')
  title(strrep(fname{i}, '_', ' '))
   set(gca, 'YDir', 'reverse');
end
pltstmp('StemCurve EXP5')
saveas(fig, 'Results\Bending_all5', 'png')
%% Create Bending Stats
% bending_stats<-data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("dir","curv","height","time"))), stringsAsFactors=F)
disp('Creating bending stats..')
bending_stats = [];
allFields = fields(result);
allDir = [];
allHook = [];
allLength = [];
allGenotype = [];
alltree = [];
allData = [];
allTreatment = [];
meta.FolderName = categorical(meta.FolderName);
for b = 1:size(allFields,1)
temp = result.(allFields{b}).output;
temp.dir = cellstr(temp.dir);
bending_stats = [bending_stats; temp]; %#ok<*AGROW>
 allHook      = [allHook; repmat(meta.hook(b),size(temp,1),1)];
 allGenotype  = [allGenotype; repmat(meta.genotype(b),size(temp,1),1)];
 alltree      = [alltree; repmat(meta.tree(b),size(temp,1),1)];
 allData      = [allData; repmat(meta.data(b),size(temp,1),1)];
 allTreatment = [allTreatment; repmat(meta.treatment(b),size(temp,1),1)];
end

bending_stats.genotype     = allGenotype;
bending_stats.tree         = alltree;
bending_stats.Data         = allData;
bending_stats.allTreatment = allTreatment;
bending_stats.hook         = allHook;

bending_stats.dir = categorical(cellstr(bending_stats.dir));
bending_stats.time = categorical(bending_stats.time);

% Add rate and cum_height to the bending_stats
bending_stats.rate = zeros(size(bending_stats,1),1);
bending_stats.cum_lift = zeros(size(bending_stats,1),1);

for i = 1:size(bending_stats,1)
  if bending_stats.time(i)=='0'
    zr = bending_stats.lift(i);
  else
    bending_stats.rate(i) = bending_stats.lift(i)-bending_stats.lift(i-1);
    bending_stats.cum_lift(i) = bending_stats.lift(i)-zr;
  end
end

%% Create box plot
disp('Plotting box plot..')
[G, Time, gen] = findgroups(bending_stats.time, bending_stats.genotype);
for i = 1:length(Time)-1
  gr(:,i) = bending_stats.cum_lift(bending_stats.time == (Time(i)));
end
fig = figure('Position', [9 57 1906 921]);
boxplot(gr), title('Cum Lift')
saveas(fig, 'Results\Cum_lift5', 'png')
% B = findgroups(bending_stats.time);
lift_mean = splitapply(@mean,bending_stats.lift,G);
hn = splitapply(@numel, bending_stats.lift,G);
lift_std = splitapply(@std,bending_stats.lift,G);
lift_se = lift_std./sqrt(hn);

curv_mean = splitapply(@mean,bending_stats.curvature,G);
curv_std = splitapply(@std,bending_stats.curvature,G);
cn = splitapply(@numel,bending_stats.curvature,G);
curv_se = curv_std./sqrt(cn);

rate_mean = splitapply(@mean,bending_stats.rate,G);
rate_std = splitapply(@std,bending_stats.rate,G);
rn = splitapply(@numel,bending_stats.rate,G);
rate_se = rate_std./sqrt(rn);

cum_mean = splitapply(@mean,bending_stats.cum_lift,G);
cum_std = splitapply(@std,bending_stats.cum_lift,G);
cumn = splitapply(@numel,bending_stats.cum_lift,G);
cum_se = cum_std./sqrt(cumn);

%% Create SU Excell file
disp('Creating SU Excell file..')
SU = table(Time, lift_mean,lift_se, hn,curv_mean,curv_se,cn,rate_mean,rate_se, cum_mean,cum_se,gen);
writetable(SU, 'Results\SU5.xlsx')
SU.gen = categorical(SU.gen);
%%  FINAL PLOTS:
% file="Stem_curve.pdf"
disp('Finall Plots')
allgen = unique(SU.gen);
fig = figure('Position', [680 49 1136 948]);
for i = 1:2 
    qq = SU.gen == allgen(i);
subplot(3,1,1)
errorbar(SU.Time(qq), SU.curv_mean(qq), SU.curv_se(qq), 'o'), hold on
 xlabel('Time (days)')
 ylabel('Stem Curvature (degree)')
title(sprintf('Curv mean %s', string(allgen(i))))

% file="Stem_rate.pdf"
subplot(3,1,2)
errorbar(SU.Time(qq),SU.rate_mean(qq), SU.rate_se(qq), 'o'), hold on
xlabel('Time (days)')
ylabel({'Rate of Bending', '(lift\day\stem length)'}, 'Interpreter', 'none')
title(sprintf('rate of bending %s', string(allgen(i))))
ylim([-5 11])
% #pdf(file="Stem_height.pdf",width=8,height=6)
subplot(3,1,3)
h(i) = errorbar(SU.Time(qq),SU.cum_mean(qq), SU.cum_se(qq), 'o'); hold on 
xlabel('Time (days)')
ylabel({'Normalized Stem Lift', '(lift / stem length)'})
title(sprintf('Cum mean %s', string(allgen(i))))
ylim([-10 60])
end
legend(h', string(allgen), 'Location', 'SouthEast')
pltstmp('StemCurve EXP5')
saveas(fig, 'Results\Stem_curve5', 'png')
