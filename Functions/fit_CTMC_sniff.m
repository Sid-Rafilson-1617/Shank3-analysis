%% Fit CTMC to sniff data

function [best_fitted_Q, empirical_proportions, sigma] = fit_CTMC_sniff(NS, data, Ntrials, max_bin, load_Q)
    % Fits a continuous time markov chain. Input arguements: number of states(NS), and an instantaneous sniff frequency file
    % containing the structures "inhalations_times_seconds" and
    % "sniff_frequencies". The function outputs a transition rate matrix
    % (best_fitted_Q), the empirically calculated proportions of state
    % occupance (empirical_proportions), and the stationary distibution (sigma). 
   
    N = size(data,2);

    for session = 1:N
        single_series = cell2mat(data(1,session));
        data_times{session} = single_series(:,1);
        freqs{session} = single_series(:,2);
    end
    %binning sniff frequencies into NS number of bins
    fprintf('___Sniff Frequency State Space___\n')
    for sesh = 1:N
        for state = 1:NS-1
            bin_indices = freqs{sesh} > ((state-1)*max_bin/(NS-1)) & freqs{sesh} <= (state*max_bin/(NS-1));
            data_values{sesh}(bin_indices,1) = state;
            if sesh == 1
                fprintf("frequencies %u to %u becomes state %u\n", (state-1)*max_bin/(NS-1),(state*max_bin/(NS-1)),state)
            end
        end
        data_values{sesh}(freqs{sesh} > (NS-1)*max_bin/(NS-1), 1) = NS;
        if sesh == 1
            fprintf("frequencies %u to %u becomes state %u\n", (NS-1)*max_bin/(NS-1),Inf,NS)
        end
    end
    
    % fit null model to loaded time series
    fprintf('Fitting null model to loaded data... \n');
    
    % Initializing parameters for training
    if load_Q
        first_q = load("Q_matrix.mat").best_fitted_Q;
        count = 1;
        for row = 1:NS
            for col = 1:NS
                guess_q(count) = first_q(row, col);
                count = count + 1;
            end
        end
    else
        guess_q = 0.01*ones(1,NS*(NS-1));
    end

    guess_q = ones([NS^2-NS,1]);
    fit_options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'Display', 'iter');
    best_logL = -Inf;
    
    for r = 1:Ntrials
        randomized_guess_q = guess_q .* (10.^(1*2-rand(size(guess_q))));
        initial_logL = loglikelihood(Qfunction(randomized_guess_q), data_times, data_values);
        if(~isfinite(initial_logL))
            continue;
        end
    
        % now we maximize the loglikelihood
        fitted_q = fminunc(@(q) - loglikelihood(Qfunction(q), data_times, data_values), randomized_guess_q, fit_options);
    
        % calculate loglikelihood for maximized parameters
        logL = loglikelihood(Qfunction(fitted_q), data_times, data_values);
    
        if((r==1) || (logL > best_logL))
            best_logL = logL;
            best_fitted_q = fitted_q;
        end
    end
    
    
    fprintf('Maximized log likelihood: %.10g\n', best_logL)
    best_fitted_Q = Qfunction(best_fitted_q);
    fprintf('Estimated transition rate matrix: \n')
    best_fitted_Q
    
    
    
    % determine stationary distribution
    fprintf('Stationary distribution of fitted null model:\n')
    sigma = linsolve([best_fitted_Q; ones([1,NS])], [zeros([NS,1]);1])'
    
    % determine empirical proportions for all series'
    fprintf('Empirical state proportion in data: \n')
    empirical_proportions = zeros(NS,N);
    for series = 1:N
        for state = 1:NS
            empirical_proportions(state, series) = mean(data_values{series} == state);
        end
    end

    



 

    
    
    % mapping free parameters q to the transition matrix
    function Q = Qfunction(q)
        Q = zeros(NS);
        count = 1;
        for i = 1:NS
            for j = 1:NS
                if i == j
                    Q(i,j) = 0;
                else
                    Q(i,j)=q(count);
                    count = count + 1;
                end
            end
        end
        Q = Q - diag(sum(Q,1));
    end
    
    
    % computes log likelihood of data given the transition rate matrix Q
    % and cell arrays containing transition times and states for N series'
    function logL = loglikelihood(Q, times, values)
        if(any(Q(~eye(NS)) <0, 'all'))
            logL = -Inf;
            return;
        end
    
        %diagonolize transition rate matrix,
        [V,D,W] = eig(Q);
        W = W';
        Winv = inv(W);
        diagonalizable = (rank(W) == NS);
    
        logL = 0;
        % loop over all sniff signals
        for session = 1:N
            M = length(data_times{session});
            for m =2:M
                delta = (times{session}(m)-times{session}(m-1));
                if(diagonalizable)
                    Qexp = real(Winv * diag(exp(diag(delta*D)))*W);
                else
                    Qexp = expm(delta.*Q);
                end
                logL = logL + log(Qexp(values{session}(m), values{session}(m-1)));
            end
        end
    end
end



