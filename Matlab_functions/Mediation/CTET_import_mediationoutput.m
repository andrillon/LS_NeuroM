function mediation_output = CTET_import_mediationoutput(filename)

%% Initialize variables.
delimiter = ' ';
startRow = 2;
endRow = inf;

%% Format for each line of text:
%   column1: text (%q)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
mediation_output = table;
mediation_output.Channels = dataArray{:, 1};
mediation_output.ACME_control = dataArray{:, 2};
mediation_output.ACME_treated = dataArray{:, 3};
mediation_output.ADE_control = dataArray{:, 4};
mediation_output.ADE_treated = dataArray{:, 5};
mediation_output.Prop_Mediated_control = dataArray{:, 6};
mediation_output.Prop_Mediated_treated = dataArray{:, 7};
mediation_output.ACME_average = dataArray{:, 8};
mediation_output.ADE_average = dataArray{:, 9};
mediation_output.Prop_Mediated_average = dataArray{:, 10};

