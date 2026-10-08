function LL = calc_DPH_logL(dwell_times,T_mat,pi_vec)
    % CALC_DPH_LOGL Calculates the exact log-likelihood of discrete dwell times.
    %
    % Inputs:
    %   pi_vec      : 1 x K initiation probability vector
    %   T_mat       : K x K transient transition matrix
    %   dwell_times : N x 1 array of strictly positive integers (discrete times)
    %
    % Output:
    %   LL          : Scalar log-likelihood value
    
    % Force dimensions
    pi_vec = pi_vec(:)'; 
    K = length(pi_vec);
    
    % Calculate exit probability vector
    T_mat = T_mat(1:K,1:K);
    t_vec = (eye(K) - T_mat) * ones(K, 1);
    
    if size(dwell_times,2)==1
        % Find the maximum discrete dwell time
        max_k = max(dwell_times);
        bin_values = 1:max_k;
        
        % Bin the data to get frequencies (counts) of each dwell time.
        % histcounts is highly optimized in MATLAB.
        counts = histcounts(dwell_times, 0.5:(max_k + 0.5));
    else
        bin_values = dwell_times(:,1)';
        counts = dwell_times(:,2)';
    end

    bin_values = bin_values(counts>0);
    counts = counts(counts>0);
    
    LL = 0;
    
    % Sequential evaluation up to the maximum observed time
    for k = 1:length(bin_values)
        % Compute P(X = k)
        prob = pi_vec * mpower(T_mat,(bin_values(k)-1)) * t_vec;
        
        % Guard against log(0) numerical underflow 
        if prob <= 0
            prob = realmin; % MATLAB's smallest positive normalized float
        end
        
        % Add weighted log-probability
        LL = LL + counts(k) * log(prob);
    end
end