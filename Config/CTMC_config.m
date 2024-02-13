function config = CTMC_config(main_config)

    % number of states
    config.NS = 7;
    
    %minimum and maximum sniff frequencies for binning
    config.min_frequency = 4;
    config.max_frequency = 12;
    
    % number of first guess' to avoid local minima in the gradient descent
    config.Ntrials = 50;
    
    config.empirical = cell(1,main_config.Nmice);
    config.stationary = cell(1,main_config.Nmice);
    config.Q = cell(1,main_config.Nmice);
end