function [K_trace, lambda_iter, w_iter] = iEMM(y, hypers, K0, niter, thin)
%iEMM - Infinite Exponential Mixture Model via Collapsed Gibbs Sampling
%   This function performs Bayesian nonparametric clustering on continuous, 
%   positive data using a Dirichlet Process Mixture Model (DPMM). It 
%   assumes an Exponential likelihood for the observations and a conjugate 
%   Gamma prior for the cluster rate parameters.
%
%   Syntax
%     [K_trace, lambda_iter, w_iter] = iEMM(y, hypers, K0, niter, thin)
%
%   Input Arguments
%     y - Array of observed data
%       column vector
%     hypers - model hyperparameters
%       structure
%         hypers.lambda_a - Shape parameter for the Gamma prior
%         hypers.lambda_b - Rate parameter for the Gamma prior
%         hypers.gamma - Concentration parameter for the Dirichlet Process
%     K0 - Initial number of clusters
%       integer
%     niter - Total number of MCMC iterations to run
%       integer
%     thin - Thinning interval for saving (save every 'thin' iterations)
%       integer
%
%   Output Arguments
%     K_trace - Number of active clusters at each saved iteration
%       column vector
%     lambda_iter - Rate parameters of active clusters at each saved iteration.
%       column cell vector
%     w_iter - Mixture weights of active clusters at each saved iteration.
%       column cell vector
%
%   Algorithm
%     Written by MCASH with the help of Gemini 3 (01/2026).
%     The function implements a Collapsed Gibbs Sampler utilizing the 
%     Chinese Restaurant Process (CRP). It analytically integrates out the 
%     mixture weights to sample individual data point assignments, applies 
%     the log-sum-exp trick for numerical stability during Softmax 
%     sampling, and concludes each loop with a global Gibbs update for the 
%     cluster rates. The first half of the iterations are automatically 
%     discarded as burn-in. Cluster parameters are sorted prior being saved
%     according to rate constants.
%
%   Reference
%     Hines, K. E., et al. (2015). "Analyzing Single-Molecule Time Series 
%     via Nonparametric Bayesian Inference." Biophysical Journal
%...


% Hyperparameters
lambda_a = hypers.lambda_a;
lambda_b = hypers.lambda_b;
gamma = hypers.gamma;

% Initialize randomizer
rng(1);

% Storing
burnin = round(niter/2); % burn the first half of iteration
nsave = floor((niter - burnin)/thin);
K_trace = zeros(nsave, 1);
lambda_iter = cell(nsave, 1);
w_iter = cell(nsave, 1);
save_idx = 0;

% Vectorize raw data and get sample size
y = y(:);
N = length(y);

% Pre-compute constant for new cluster
LOG_GAMMA = log(gamma);
LOG_PRED_CONST = log(lambda_a) + lambda_a*log(lambda_b); % log(Marginal Predictive)

% Initialization
z = randi(K0, [N, 1]); % random assignment of dwell times to clusters
K = K0;
Nz_k = zeros(1,K); % cluster population
lambda_k = zeros(1,K); % cluster lifetime
sum_y_k = zeros(1,K); % cluster time fraction
for k = 1:K
    Nz_k(k) = nnz(z==k);
    sum_y_k(k) = sum(y(z==k));
    lambda_k(k) = gamrnd(lambda_a + Nz_k(k), 1./(lambda_b + sum_y_k(k)));
end

% Collapsed Gibbs sampling loop
for it = 1:niter
    for i = 1:N
        % Remove data point i
        old_k = z(i);
        Nz_k(old_k) = Nz_k(old_k) - 1;
        sum_y_k(old_k) = sum_y_k(old_k) - y(i);

        % Delete cluster if empty after removal of point i
        if Nz_k(old_k) == 0
            lambda_k(old_k) = [];
            Nz_k(old_k) = [];
            sum_y_k(old_k) = [];
            z(z > old_k) = z(z > old_k) - 1;
            K = K - 1;
        end

        % Calculate re-assingment probabiltiies (CRP)
        % Log-Likelihood for existing clusters
        log_lik = log(Nz_k) + log(lambda_k) - (lambda_k .* y(i));

        % Log-Likelihood for a new cluster
        log_lik_new = LOG_GAMMA + LOG_PRED_CONST - ...
            (lambda_a+1) * log(lambda_b+y(i));

        % Sampling via Softmax (Log-Sum-Exp trick pour la stabilité)
        all_logs = [log_lik(:); log_lik_new];
        probs = exp(all_logs - max(all_logs));
        new_k = find(cumsum(probs/sum(probs)) >= rand(), 1);

        % 3. Mise à jour
        z(i) = new_k;
        if new_k > K
            % Création d'un nouveau cluster
            K = K + 1;
            Nz_k(new_k) = 1;
            sum_y_k(new_k) = y(i);
            lambda_k(new_k) = gamrnd(lambda_a + 1, 1./(lambda_b + y(i)));
        else
            Nz_k(new_k) = Nz_k(new_k) + 1;
            sum_y_k(new_k) = sum_y_k(new_k) + y(i);
        end
    end

    % Global update of lambda for each clusters (Gibbs Step)
    lambda_k = gamrnd(lambda_a + Nz_k, 1./(lambda_b + sum_y_k));

    % Save iteration
    if it > burnin && mod(it, thin) == 0
        save_idx = save_idx + 1;
        [~,ord] = sort(lambda_k);
        K_trace(save_idx) = K;
        lambda_iter{save_idx} = lambda_k(ord);
        w_iter{save_idx} = Nz_k(ord) / N;
    end
end

% % Calculate log-posterior
% LP_iter = calculate_log_posterior(q_iter, w_iter, y, A, B, alpha);
end


% function [logPost] = calculate_log_posterior(q_iter, w_iter, y, A, B, alpha)
%     % q_iter : cell array des vecteurs de taux [1 x K]
%     % w_iter : cell array des vecteurs de poids [1 x K]
%     % y      : données (dwell times)
%     % A, B   : hyperparamètres de la prior Gamma
%     % alpha  : paramètre de concentration du processus de Dirichlet
% 
%     num_samples = numel(q_iter);
%     logPost = zeros(num_samples, 1);
%     N = numel(y);
% 
%     for s = 1:num_samples
%         q = q_iter{s}(:); % Taux de l'itération s
%         w = w_iter{s}(:); % Poids (n_k / N)
%         K = numel(q);
% 
%         % 1. Log-Vraisemblance (Incomplète)
%         % sum_{i=1}^N log( sum_{k=1}^K w_k * q_k * exp(-q_k * y_i) )
%         % Utilisation du Log-Sum-Exp pour la stabilité numérique
%         log_comp = log(w) + log(q) - (q * y'); % Matrice [K x N]
%         logL = sum(max(log_comp, [], 1) + ...
%             log(sum(exp(log_comp - max(log_comp, [], 1)), 1)));
% 
%         % 2. Log-Prior des paramètres (Gamma distribution)
%         % sum_{k=1}^K [ (A-1)*log(q_k) - B*q_k ]
%         logPrior_params = sum((A - 1) * log(q) - B * q);
% 
%         % 3. Log-Prior de la structure (Dirichlet Process / Chinese Restaurant Process)
%         % log(alpha^K) + sum(log(factorial(n_k - 1))) - log(Pochhammer(alpha, N))
%         % n_k = w_k * N
%         n_k = w * N;
%         logPrior_structure = K * log(alpha) + sum(gammaln(n_k)) - ...
%             gammaln(alpha + N) + gammaln(alpha);
% 
%         % Log-Postériore Totale
%         logPost(s) = logL + logPrior_params + logPrior_structure;
%     end
% end


% function plotiter(ax,bincnt,binedg,K,iter)
% 
%     chld = ax(1).Children;
%     if isempty(chld)
%         line(ax(1),iter,K);
%     else
%         chld.XData = [chld.XData,iter];
%         chld.YData = [chld.YData,K];
%     end
% 
%     oldxlim = ax(2).XLim;
%     oldylim = ax(2).YLim;
%     histogram(ax(2),'bincounts',bincnt,'binedges',binedg);
%     if ax(2).XLim(1)>oldxlim(1)
%         ax(2).XLim(1) = oldxlim(1);
%     end
%     if ax(2).XLim(2)<oldxlim(2)
%         ax(2).XLim(2) = oldxlim(2);
%     end
%     if ax(2).YLim(2)<oldylim(2)
%         ax(2).YLim(2) = oldylim(2);
%     end
%     title(ax(2),sprintf('iteration %i, K=%i',iter,K));
%     ax(2).XTick = 1:floor(ax(2).XLim(2));
%     drawnow;
% end
