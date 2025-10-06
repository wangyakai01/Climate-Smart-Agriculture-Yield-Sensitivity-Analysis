function loadAllVars(folderPath)
      % Ensure the folder path ends with a separator
    % Ensure folder path ends with a separator
    if folderPath(end) ~= filesep
        folderPath = [folderPath, filesep];
    end

    % Get list of .mat files
    matFiles = dir([folderPath, '*.mat']);

    % Load each file's variables into the base workspace
    for i = 1:length(matFiles)
        filePath = fullfile(matFiles(i).folder, matFiles(i).name);
        disp(['Loading: ', matFiles(i).name]);

        % Load into a structure first
        dataStruct = load(filePath);

        % Push each variable into the base workspace
        varNames = fieldnames(dataStruct);
        for j = 1:length(varNames)
            assignin('base', varNames{j}, dataStruct.(varNames{j}));
        end
    end
end
