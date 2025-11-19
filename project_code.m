%% Enable Parallel Pool
p = gcp('nocreate');   % If no pool, do NOT auto-start yet
if isempty(p)
    p = parpool;       % Start pool with default workers
end
numWorkers = p.NumWorkers;
fprintf('Using %d parallel workers.\n', numWorkers);

%% Normal (Non-parallel) Execution
fprintf('\n---- Running SERIAL version ----\n');
tic;

noise_dbm=-80;
noise=10^((noise_dbm-30)/10);

P_ts_d = linspace(-30, 10, 30);
P_ts = 10.^((P_ts_d-30) / 10);

N = 4;
iterations = 1e6;
omega = 1; 
delta_n_square = noise;
threshold_dB = 3;
threshold=10^(threshold_dB/10);
P_dBm=-60;
P=10.^((P_dBm-30) / 10);

r=2;
theta=1/2;
beta_u=0.5;
d_1=3;
tau=-3;
h_si=gamrnd(r, theta, 1, iterations);

shape_hs = [3,4,5];
mean_outage_serial = zeros(length(shape_hs), length(P_ts));

z = (threshold*P) ./(d_1^(4*tau).* beta_u^2 .* P_ts);

% SERIAL loop
for k=1:length(shape_hs)
    shape_h = shape_hs(k);
    scale = omega/shape_h;

    hs = sqrt(gamrnd(shape_h, scale, N, iterations));
    gs = sqrt(gamrnd(shape_h, scale, N, iterations));
    ns = sqrt(gamrnd(shape_h, scale, N, iterations));
    ms = sqrt(gamrnd(shape_h, scale, N, iterations));

    H = sum(hs .* gs .* ms .* ns).^2;

    for i = 1:length(P_ts)
        SNR = d_1^(4*tau).* beta_u^2 .* P_ts(i) .* H ./ (P.*h_si + noise);
        count = sum(SNR < threshold);
        mean_outage_serial(k,i) = count / iterations;
    end
end

serialTime = toc;
fprintf('Serial Execution Time: %.4f sec\n', serialTime);



%% PARALLEL Execution
fprintf('\n---- Running PARALLEL version ----\n');

tic;

mean_outage_parallel = zeros(length(shape_hs), length(P_ts));  

parfor k = 1:length(shape_hs)
    shape_h = shape_hs(k);
    scale = omega/shape_h;

    % Each worker generates its own random matrices
    hs = sqrt(gamrnd(shape_h, scale, N, iterations));
    gs = sqrt(gamrnd(shape_h, scale, N, iterations));
    ns = sqrt(gamrnd(shape_h, scale, N, iterations));
    ms = sqrt(gamrnd(shape_h, scale, N, iterations));

    H = sum(hs .* gs .* ms .* ns).^2;

    tmp_outage = zeros(1, length(P_ts));

    for i = 1:length(P_ts)
        SNR = d_1^(4*tau).* beta_u^2 .* P_ts(i) .* H ./ (P.*h_si + noise);
        tmp_outage(i) = sum(SNR < threshold) / iterations;
    end

    mean_outage_parallel(k,:) = tmp_outage;
end

parallelTime = toc;
fprintf('Parallel Execution Time: %.4f sec\n', parallelTime);



%% Compute Speedup
speedup = serialTime / parallelTime;
fprintf('\n=====================================\n');
fprintf('   Workers Used   : %d\n', numWorkers);
fprintf('   Serial Time    : %.4f sec\n', serialTime);
fprintf('   Parallel Time  : %.4f sec\n', parallelTime);
fprintf('   Speedup        : %.2f × faster\n', speedup);
fprintf('=====================================\n');
