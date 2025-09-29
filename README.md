# AMS595
Project 1

% AMS 595 Programming, John Walsh, Fall 2025
%Monte Carlo Estimation of Pi, Project 1, Part 1

clear; clc; close all;

% Define number of points (N) to test (increasing)
N_values = round(logspace(2, 6, 10)); % from 1e2 to 1e6 (10 steps)
pi_estimates = zeros(size(N_values));
errors = zeros(size(N_values));
times = zeros(size(N_values));

true_pi = pi;

for k = 1:length(N_values)
    N = N_values(k);
    tic; % start timing
    
    inside = 0; % counter for points inside circle
    for i = 1:N
        x = rand;
        y = rand;
        if x^2 + y^2 <= 1
            inside = inside + 1;
        end
    end
    
    pi_est = 4 * inside / N;
    pi_estimates(k) = pi_est;
    errors(k) = abs(pi_est - true_pi);
    
    times(k) = toc; % end timing
end

%% Plot results
figure;

% Plot convergence of Pi estimate
subplot(2,1,1);
semilogx(N_values, pi_estimates, 'b-o', 'LineWidth', 2.5);
hold on;
yline(true_pi, 'r--', 'LineWidth', 2);
xlabel('Number of Points (N)');
ylabel('Estimated \pi');
title('Project 1: Monte Carlo Estimation of \pi');
legend('Estimated \pi','True \pi','Location','best');
grid on;

% Plot error vs computation time
subplot(2,1,2);
loglog(N_values, errors, 'k-o', 'LineWidth', 1.5);
hold on;
loglog(N_values, times, 'g-s', 'LineWidth', 1.5);
xlabel('Number of Points (N)');
ylabel('Error / Time (s)');
title('Error and Computation Time vs. N');
legend('Error','Execution Time','Location','best');
grid on;




% AMS 595 Programming, John Walsh, Fall 2025
%Monte Carlo Estimation of Pi, Project 1, Part 2
% Modification to Part 1, Monte Carlo Pi with while loop and precision by significant figures

clear; clc; close all;

% === Settings ===
target_sigs = [2 3 4];     % target significant figures to achieve
alpha = 0.05;              % confidence level for the precision check (95% CI)
batch = 1e4;               % draw points in batches for speed, but keep a while loop

% Storage
N_needed   = zeros(size(target_sigs));
pi_hat_final = zeros(size(target_sigs));
time_taken = zeros(size(target_sigs));

Z = norminv(1 - alpha/2);  % 1.96 for 95% CI  %Normalization

for t = 1:numel(target_sigs)
    s = target_sigs(t);

    % Reset counters per precision target
    N = 0;
    inside = 0;
    tic;

    % While loop until required precision (in significant figures) is met
    done = false;
    while ~done
        % --- Monte Carlo draws (in a batch to maximize efficiency) ---
        x = rand(batch,1);
        y = rand(batch,1);
        inside = inside + sum(x.^2 + y.^2 <= 1);
        N = N + batch;

        % --- Current estimate and CI-based precision check (NO true pi used) ---
        p_hat = inside / N;        % probability a point falls in the quarter circle
        pi_hat = 4 * p_hat;        % estimate of pi

        % Standard error for p, then delta for pi
        se_p  = sqrt(p_hat * (1 - p_hat) / N);
        halfwidth = Z * 4 * se_p;  % CI half-width for pi_hat

        % Translate "s significant figures" into an absolute tolerance
        % If x ~ 3.14..., the first significant digit is the ones place (10^0).
        % Absolute tolerance corresponding to s sig figs is 0.5 * 10^(d - s),
        % where d = floor(log10(|x|)) + 1 (number of digits to left of decimal).
        % Estimate's order of magnitude (not the true value).
        d = floor(log10(abs(pi_hat))) + 1;
        tol_sf = 0.5 * 10^(d - s);

        % Stop when the *confidence half-width* fits inside the sig-fig tolerance
        % (i.e., the estimate is statistically precise enough to fix s significant digits)
        if halfwidth <= tol_sf
            done = true;
        end
    end

    time_taken(t) = toc;
    N_needed(t)   = N;
    pi_hat_final(t) = pi_hat;
end

% === Create a Results table ===
T = table(target_sigs(:), N_needed(:), pi_hat_final(:), time_taken(:), ...
    'VariableNames', {'SigFigs','Iterations_N','Pi_Estimate','Time_sec'});
disp(T);

% === Plots ===
figure;

% Iterations required vs significant figures
subplot(2,1,1);
semilogy(target_sigs, N_needed, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on;
xlabel('Significant Figures Target');
ylabel('Iterations N (log scale)');
title('Iterations Required to Reach Target Significant Figures (95% CI)');


% AMS 595 Programming, John Walsh, Fall 2025
%Monte Carlo Estimation of Pi, Project 1, Part 3

function pi_hat = mcPiSigFigs(target_sigs)  % I am having difficulty implementing the call function.  I have not fully debugged.

% mcPiSigFigs  Monte Carlo estimation of pi to a user-specified number of significant figures.
%   pi_hat = mcPiSigFigs(target_sigs)
%
% Features from Part 3 modification:
%   • Uses a while loop and a CI-based precision criterion (no true pi used).
%   • Plots points live: inside quarter circle vs outside (different colors).
%   • Prints the final estimate to the command window and onto the plot, rounded
%     to the requested number of significant figures.
%   • Returns the computed pi estimate.
%
% Inputs:
%   target_sigs (positive integer): number of significant figures to achieve.
%
% Notes:use Confidence Interval and batch method from Part 2
%   - Uses a 95% CI by default (modifiable via 'alpha' below).
%   - Draws points in batches for efficiency, while still updating the plot.

    arguments
        target_sigs (1,1) {mustBeInteger, mustBePositive}
    end

    % ---- Tunable parameters ----
    alpha    = 0.05;   % 95% confidence interval for the stopping rule
    batch    = 5e3;    % number of points per iteration (plot update granularity)
    showEveryBatch = 1;     % update plot every 'showEveryBatch' batches

    % ---- Init counters ----
    N = 0;                % total samples
    inside = 0;           % hits inside quarter circle
    Z = norminv(1 - alpha/2);   % ~1.96 for 95% CI

    % ---- Figure setup ----
    figure('Name','Monte Carlo \pi Estimation','Color','w');
    ax = axes; hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
    xlim([0 1]); ylim([0 1]);
    xlabel('x'); ylabel('y');
    title(sprintf('Monte Carlo \\pi Estimation (Target: %d s.f.)', target_sigs));

    % Quarter circle outline for reference
    th = linspace(0, pi/2, 400);
    plot(cos(th), sin(th), 'k-', 'LineWidth', 1.5);

    % Pre-create two plot handles for fast incremental updates
    insidePlot  = plot(nan, nan, '.', 'MarkerSize', 8, 'Color', [0 0.45 0.74], 'DisplayName','Inside');
    outsidePlot = plot(nan, nan, '.', 'MarkerSize', 8, 'Color', [0.85 0.33 0.1], 'DisplayName','Outside');

    legend('Quarter circle','Inside','Outside','Location','southoutside','Orientation','horizontal');

    % For live plotting, keep running buffers
    Xin = []; Yin = [];
    Xout = []; Yout = [];

    % Text handle for status updates
    statusText = text(0.02, 0.98, '', 'Units','normalized', 'VerticalAlignment','top', ...
                      'FontName','Consolas', 'FontSize',10, 'Color',[0.25 0.25 0.25]);

    % ---- While loop until precision target met ----
    done = false;
    iter = 0;
    while ~done
        iter = iter + 1;

        % --- Draw a batch of points ---
        x = rand(batch,1);
        y = rand(batch,1);
        r2 = x.^2 + y.^2;
        in = r2 <= 1;

        % Update counters
        inside = inside + sum(in);
        N = N + batch;

        % --- Update plot buffers ---
        Xin  = [Xin;  x(in)];
        Yin  = [Yin;  y(in)];
        Xout = [Xout; x(~in)];
        Yout = [Yout; y(~in)];

        % Update plot every 'showEveryBatch' batches for responsiveness
        if mod(iter, showEveryBatch) == 0
            set(insidePlot,  'XData', Xin,  'YData', Yin);
            set(outsidePlot, 'XData', Xout, 'YData', Yout);
            drawnow limitrate;
        end

        % --- Estimate and CI-based stopping rule (NO true pi used) ---
        p_hat  = inside / N;     % probability of landing inside quarter circle
        pi_hat = 4 * p_hat;

        % Binomial standard error for p, mapped to pi by factor 4
        se_p = sqrt(max(p_hat*(1 - p_hat), eps) / N); % guard for numerical safety
        halfwidth = Z * 4 * se_p;  % CI half-width for pi_hat

        % Translate "target_sigs significant figures" into absolute tolerance
        d = floor(log10(abs(pi_hat))) + 1;  % digits to left of decimal of current estimate
        tol_sf = 0.5 * 10^(d - target_sigs);

        % Update status text
        set(statusText, 'String', sprintf( ...
            'N = %d\n\\pî = %.6f\nCI half-width = %.3g\nSig-fig tol = %.3g', ...
            N, pi_hat, halfwidth, tol_sf));

        % Stop when CI half-width fits inside the sig-fig tolerance
        if halfwidth <= tol_sf
            done = true;
        end
    end

    % ---- Final reporting at requested precision ----
    pi_str = sprintf(['%.' num2str(target_sigs) 'g'], pi_hat);

    % Command window output
    fprintf('Final estimate of pi (to %d significant figures): %s\n', target_sigs, pi_str);
    fprintf('Total samples used: %d\n', N);

    % Print on the plot
    finalText = sprintf('\\pi \\approx %s  (N = %d)', pi_str, N);
    text(0.5, -0.08, finalText, 'Units','normalized', 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'FontSize', 12, 'Color', [0.1 0.1 0.1]);

end

