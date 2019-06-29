% 2018-06-18 
% This file is using the testing curves from 
% Working through the Stem curve function with Example datasets####
% Ex1.flat line
% Ex2 line with unifrom slope
% Ex3 line with changing slope (ie curved)
% Ex4 line with changigne slope
% Ex5 negative slope line

% 1) The goal of the curvature calculation is determine how much a curved 
%   line deviates from a straight line. Where a straight line has a uniform 
% slope or constant change in degrees along a line. My previous statement
% about sum(degrees) will not tell you how curved a line is but instead that
% calculation would represent the incline angle of a plant and this is going 
% to be highly correlated with how much a plant lifted.

% 2) To calculate curvature, I calculated how much the angle (degrees) 
% changed from one segment to the next. 
% At the smallest level, curvature = Segment 1 angle - segment 0.
% At a plant level, the code loops across all the segments and then sums 
% all the deviations from a straight line. To illustrate this, 
% I created several example data sets ranging from flat lines with constant
% slopes to lines that curve up and down. In the case of a flat line (Ex1 & Ex2),
% the angle from one segment to the next does not change which results in a 
% curvature value of zero. It is important to note that in Ex2 the line still 
% lifts but it does not curve. For lines that curve upward the curvature 
% values are greater than zero because the the slope keep increasing from 
% one segment to the next and negative curvature represents lines that bend down. 

clear all
close all
clc

NumberOfExamples = 8;
for i = 1:NumberOfExamples
    filename = sprintf('Ex%i.txt', i);
    ex = loadexample(['Testing\' filename]);
    subplot(1,NumberOfExamples,i), plot(ex.X,ex.Y), hold on
    plot(0:10, 0:10,'--r') % Straight line
    axis tight
    xlabel(filename)
end
result = stemcurve_fun('Testing\', '.txt', 10, []);
result.output;
for i = 1:NumberOfExamples
    subplot(1,NumberOfExamples,i)
    title(result.output.curvature(i))
end

%% Function for importing data from the text file
function Ex2 = loadexample(filename)
%% Initialize variables.
delimiter = '\t';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Create output variable
Ex2 = table;
Ex2.X = cell2mat(raw(:, 1));
Ex2.Y = cell2mat(raw(:, 2));

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp;
end