function result = stemcurve_fun(directory, grep_key, hook, len)
% STEMCURVE_FUN(DIRECTORY,GREP_KEY, HOOK,LEN)
% Processing X,Y coordinates from Poplar tension wood bending to extract curvature
% INPUTS
% directory - [STRING] 
% grep_key - [STRING]
% OUTPUTS: result struct with fields: 
%     output:  - contains curvature and lift for each day
%        out:  - contains XY position for each measuring point (several during one day)
%        vec:  - contains difference between x,y and angle at each measuring time points.
% 
  all_files = dir([directory, '\*.txt']);
  files     = all_files;
  n         = numel(all_files);

  out = []; % c("file","x","y")))
  vec = []; % c("file","xis","yis","angle")
  slashqq = strfind(directory, '\') ;
  folderName = directory(slashqq+1:end);
  output = table([0:n-1]', repmat(folderName,n,1), nan(n,1), nan(n,1),...
    'VariableNames', {'time', 'dir', 'curvature', 'lift'});
  
for j = 1:n
  outi = [];
  veci = [];
  XYdata = readtable([directory '\' all_files(j).name]);  
  cSum = 0;
  % Check distance from base of tree, up to hook distance. Save k for that
  % distance
  for k = 1:size(XYdata,1)-1
    X = XYdata{k:k+1, :};
    d(k) =  pdist(X,'euclidean');
    cSum = cSum +pdist(X,'euclidean');
    if(cSum>=hook), break, end
  end
  for kk = 1:k
    XYBase2Hook(kk,:) = {[directory files(j).name], XYdata{kk,1}, XYdata{kk,2}};
  end
  outi = XYBase2Hook;
  XYBase2Hook = cell2table(XYBase2Hook, 'VariableNames', {'file', 'x', 'y'});
% Normalize Base2hook by first values of X and Y
XYBase2Hook.x = XYBase2Hook.x-XYBase2Hook.x(1);
XYBase2Hook.y = XYBase2Hook.y-XYBase2Hook.y(1);

%%   # calculate vectors
xis = XYBase2Hook.x(2:end) - XYBase2Hook.x(1:end-1); % distance between each x (x1:x2)
yis = XYBase2Hook.y(2:end) - XYBase2Hook.y(1:end-1); % distance between each y
vectors = table([1:numel(xis)]', xis, yis, 'VariableNames', {'ID', 'xis', 'yis'});
 
%  #http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/index.htm
% Calculate Angle between every pair of pixels:
for i = 1:size(vectors,1)-1
    angle(i) = rad2deg(atan2(yis(i+1), xis(i+1))- atan2(yis(i), xis(i)));
end

vectors = vectors(1:numel(angle),:);
vectors.angle = angle';
vectors.file = repmat({[directory files(j).name]},size(vectors,1),1);
% curvature = sum of incremental deviations from a straight line
 output.curvature(j) =  sum(vectors.angle);
 output.lift(j) = XYdata{1,2}-XYdata{k,2}; %  [YHook] -[YBase]

out = [out;outi]; %#ok<*AGROW>
vec = [vec;vectors];
clear XYBase2Hook angle
end
% I assume that all are txt
result.output = output;
result.out = out;
result.vec = vec;
