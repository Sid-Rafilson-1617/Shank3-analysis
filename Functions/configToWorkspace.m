function configToWorkspace(config)
    % configToWorkspace assigns each field of the config structure
    % to a variable in the base workspace

    fields = fieldnames(config);
    for i = 1:length(fields)
        assignin('base', fields{i}, config.(fields{i}));
    end
end