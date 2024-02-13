function [times, states] = simulate_CTMC(Q, max_time)
    % Using the Gillespie algorithm to simulate finite homogeneous CTMC
    
    %storing transition rate for each state in 1D array
    NS = length(Q);
    transition_rate = zeros(1,NS);
    for state = 1:NS 
        transition_rate(state) = -Q(state, state);
    end
    
    %finding the stationary distribution of the CTMC
    sigma = linsolve([Q; ones(1, length(Q))], [zeros(NS, 1); 1]);
    
    % draw initial state based on stationary distribution
    X_0 = randsample(NS, 1, true, sigma);

    times = [0];
    states = [X_0];
    
    while times(end) < max_time
        X = states(end);
    
        % wait time until next transition
        wait_time = exprnd(1/transition_rate(X));
        
        %draw new state
        X = randsample(NS, 1, true, transition_rate/wait_time);
        
        % update time and states
        times(end + 1) = times(end) + wait_time;
        states(end + 1) = X;
    end
end