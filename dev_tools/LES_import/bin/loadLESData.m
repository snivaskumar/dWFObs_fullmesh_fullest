function [filesInFolder,flowData,turbData,dataSource] = loadLESData(scriptOptions,rawTurbData)
% Sort and import files from data folder
disp('Sorting and importing files from source folder...')
rawFileList = dir(scriptOptions.sourcePath); % Gather all files from PALM/SOWFA folder
if length(rawFileList) < 3
    error('The specified directory does not exist/has no LES files.');
end
for j = 1:length(rawFileList)-2
    filesInFolder{j} = [rawFileList(j+2).folder '/' rawFileList(j+2).name];   % Remove '.' and '..'
end
filesInFolder = natsortfiles(filesInFolder); % Sort numerically

% Decide whether this is SOWFA or PALM data
if nnz(cell2mat(strfind(filesInFolder,'.nc'))) > 0
    disp('Found a .nc file. Assuming this is PALM data.')
    % Load all flow data and turbine data at once
    [ flowData,turbData ] = loadPALMdata(filesInFolder,rawTurbData.hubHeight);
    dataSource            = 'palm';
elseif nnz(cell2mat(strfind(filesInFolder,'.vtk'))) > 0
    disp('Found a .vtk file. Assuming this is SOWFA data.')
    % Load only flow data mesh and complete turbine data
    [ flowData,turbData ] = loadSOWFAdata(filesInFolder,scriptOptions);
    dataSource            = 'sowfa';
else
    error('Did not find any SOWFA/PALM data in the folder.');
end
end