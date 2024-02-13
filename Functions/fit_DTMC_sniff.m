function [best_fitted_P, empirical_proportions, sigma] = fit_DTMC_sniff(NS, min_frequency, max_frequency, data, Ntrials, max_bin, load_P)

    N = size(data,2)

    for session = 1:N
        single_series = cell2mat(data(1,session));
        data_times{session} = single_series(:,1);
        freqs{session} = single_series(:,2);
    end

    %binning sniff frequencies into NS number of bins from min to max
    %frequency
    freq_range = max_frequency - min_frequency;
    SS_bin_size = freq_range / (NS - 2);
    fprintf('___Sniff Frequency State Space___\n')
    for sesh = 1:N
        for state = min_frequency:SS_bin_size:max_frequency - SS_bin_size
            bin_indices = freqs{sesh} >= state & freqs{sesh} < state + SS_bin_size;
            data_values{sesh}(bin_indices,1) = 2 + (state - min_frequency) /SS_bin_size;
            if sesh == 1
                fprintf("frequencies %u to %u becomes state %u\n", state, state + SS_bin_size, 2 + (state - min_frequency) /SS_bin_size)
            end
        end
        data_values{sesh}(freqs{sesh} >= max_frequency) = NS;
        data_values{sesh}(freqs{sesh} < min_frequency) = 1;
        if sesh == 1
            fprintf("frequencies %u to %u becomes state %u\n", max_frequency,Inf,NS)
            fprintf("frequencies %u to %u becomes state %u\n", 0, min_frequency,1)
        end
    end


    % Initializing parameters for training
    if load_P
        first_p = load("P_matrix.mat").P;
        first_p = cell2mat(first_p);
        count = 1;
        for row = 1:NS
            for col = 1:NS
                if row ~= col
                    guess_q(count) = first_p(row, col);
                    count = count + 1;
                end
            end
        end
    else
        guess_q = 0.01*ones(1,NS*(NS-1));
    end

    fit_options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
    
	best_logL = -Inf; % keep track of best log-likelihood so far
    for r=1:Ntrials
		% randomly modify each parameter guess
		% within 10% and 1000% of its original value
		randomized_guess_q = guess_q .* (10.^(1-2*rand(size(guess_q))));
		
		% skip this trial if randomized_guess_q is problematic
		initial_logL = loglikelihood(randomized_guess_q);
		if(~isfinite(initial_logL))
            continue; end

		% maximize log-likelihood starting with the new parameter guess
		fitted_q = fminunc(@(q) -loglikelihood(q), randomized_guess_q, fit_options);
	
		% calculate maximized log-likelihood for this fit
		logL = loglikelihood(fitted_q);
        fprintf('Fit # %u, LogL = %e \n', r, logL)

		if((r==1) || (logL>best_logL))
			% found a better fit, so keep track of it
			best_logL     = logL;
			best_fitted_q = fitted_q;
		end
    end
    
    fprintf('Maximized log-likelihood: %.10g\n',best_logL);

	% calculate corresponding estimated transition matrix
	best_fitted_P = Pfunction(best_fitted_q);
	
	fprintf('Estimated transition matrix:\n');
	best_fitted_P

    % empirical proportions
    fprintf('Empirical state proportion in data: \n')
    empirical_proportions = zeros(NS,N);
    for session = 1:N
        for state = 1:NS
            empirical_proportions(state, session) = mean(data_values{session} == state);
        end
    end
    empirical_proportions'

    % determine stationary distribution
    fprintf('Stationary distribution of fitted null model:\n')
    [V,D] = eig(best_fitted_P);
    lambdas = diag(D);
    u = V(:,abs(lambdas-1)<1e-10);
    sigma = bsxfun(@rdivide, u, sum(u,1))







    %% Auxillary Functions

    % mapping free parameters q to the transition matrix 
    function P = Pfunction(q)
        P = zeros(NS);
        count = 1;
        for i = 1:NS
            for j = 1:NS
                if i == j
                    P(i,j) = 0;
                else
                    P(i,j)=q(count);
                    count = count + 1;
                end
            end
        end
        P = P + diag(1-sum(P,1));
    end



    % log likelihood function
    function logL = loglikelihood(q)
		P = Pfunction(q); % get transition matrix for the specific parameters q
		if(any(P<0,'all') || any(P>1,'all'))
			logL = -Inf;
			return;
		end
        logL = 0;
        for session = 1:N
		    % pre-calculate all required powers of P
		    % and store them in a 3D array
            M = length(data_values{session});
		    max_power = max(diff(data_times{session})); % maximum power to which we must raise P
		    Ppowers = zeros([NS NS max_power]);
		    Ppowers(:,:,1) = P;
		    for p=2:max_power
			    Ppowers(:,:,p) = P * Ppowers(:,:,p-1);
		    end
    
		    % compute log-likelihood using pre-computed matrix powers
            for m = 2:M
                p = data_times{session}(m) - data_times{session}(m-1);
	            logL = logL + log(Ppowers(data_values{session}(m), data_values{session}(m-1),p));
            end
        end
    end


end
