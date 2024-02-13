function config = sniff_analysis_config()

    config.mouse_0 = 1;
    config.mouse_1 = 2;
    
    config.signal_0 = 1;
    config.signal_1 = 2;
    
    config.Nmice = config.mouse_1 - config.mouse_0 + 1;
    config.Nsignals = config.signal_1 - config.signal_0 + 1;

end