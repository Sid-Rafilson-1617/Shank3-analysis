function config = main_config()

    %graphical visualization
    config.visualization = true;
    
    %number of trials and mice to run
    config.Nsignals = 7;
    config.Nmice = 2;
    
    %Mouse names and sessions
    config.mouse_names = {'Mouse 3001', 'Mouse 3005'};
    config.signal_names = {'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5', 'Session 6','Session 7'};

    %oscilliscope sampling rate
    config.f = 800;
end