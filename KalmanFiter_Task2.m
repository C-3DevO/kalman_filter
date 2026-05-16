%% Lab 2 - Kalman Filter
clear;clc;close all

N = 100;
p = 2;
A = [0.99 0.01;0.01 0.999];
Q = [0.001 0 ; 0 0.001];
sigma_w= 0.1;
%% Fig 13.16 (channel realization)
Hn = zeros(p, N+2);
hn = [1; 0.9];    % h[-1]
Hn(:,1) = hn;
for n=1:N+1
    hn = A*hn+(mvnrnd(zeros(p,1),Q))'; %dynamic model
    Hn(:, n+1) = hn;
end
figure;
subplot(2,1,1)
plot(0:N,Hn(1,1:N+1));
ylim([0 2])
ylabel('tap weight,hn[0]]')
xlabel('Sample number, n')

subplot(2,1,2)
plot(0:N,Hn(2,1:N+1));
ylim([0 2]);
ylabel('tap weight,hn[1]]')
xlabel('Sample number, n')

%% Signal generation (Fig 13.17)
%%%% input v[n] %%%%
T=10;
vn = zeros(N+2,1);
for n=T/2+1:T+1
    vn(n:T:N+1) = 1;
end
figure;
subplot(3,1,1)
plot(0:N,vn(1:N+1)); %Fig 13.17
ylim([-1 3])
ylabel('channel input,v[n]')
xlabel('Sample number, n')
%%%% noiseless channel output  y[n] %%%%

yn = zeros(N+1,1);


yn(1) = Hn(1,2)*vn(1);

for n = 1:N
    yn(n+1) = Hn(1,n+2)*vn(n+1) + Hn(2,n+2)*vn(n);
end

subplot(3,1,2)
plot(0:N, yn(1:N+1) ) % plot y[n] (Fig 13.17b) 
ylim([-1 3]);
ylabel('noiseless channel output, y[n]')
xlabel('Sample number, n')
%%%% channel output x[n] %%%%
wn = sigma_w*randn(N+1,1);
xn = yn + wn; 


subplot(3,1,3)
plot(0:N, xn(1:N+1)) % plot received signal with noise (Fig 13.17b) 
ylim([-1 3]);
ylabel('channel output, x[n]')
xlabel('Sample number, n')
%% KALMAL FILTER %%%%%%%%%%%%%%%%%%%

sigma2 = 0.1;        % observation noise variance

h_hat = zeros(p,1);  % h[-1|-1] initial estimate (2x1)
M = 100*eye(p);      % M[-1|-1] initial error covariance (2x2)

% storage for plots (like Kay Fig 13.18/13.19/13.20)
Hhat = zeros(p, N+1);    % estimated taps over time: h[n|n]
Khist = zeros(p, N+1);   % Kalman gains K[n]
mse_trace = zeros(1, N+1); % trace of M[n|n]
Mdiag = zeros(p, N+1);      % diag of M[n|n] for each n (2 entries)
for n = 0:N
    % Build measurement vector v[n]^T = [v[n]  v[n-1]]
    if n == 0
        vvec = [vn(1); 0];      % v[-1]=0
    else
        vvec = [vn(n+1); vn(n)];
    end

    % ---- Prediction step (time update) ----
    h_pred = A*h_hat;
    M_pred = A*M*A' + Q;

    % ---- Kalman gain ----
    denom = sigma2 + vvec.'*M_pred*vvec;   % scalar
    K = (M_pred*vvec) / denom;             % (2x1)

    % ---- Measurement update ----
    innov = xn(n+1) - vvec.'*h_pred;       % innovation: x - predicted x
    h_hat = h_pred + K*innov;
    M = (eye(p) - K*vvec.')*M_pred;

    % store
    Hhat(:, n+1) = h_hat;
    Khist(:, n+1) = K;
    mse_trace(n+1) = trace(M);
    Mdiag(:, n+1) = diag(M);
end

figure;

subplot(4,1,1)
plot(0:N, Hn(1,2:N+2), 'k', 0:N, Hhat(1,:), 'r--');
ylim([0 2])
xlabel('Sample number, n')
ylabel('tap 0')
legend('True h_n[0]', 'Estimated \hat{h}_n[0|n]')

subplot(4,1,2)
plot(0:N, Hn(2,2:N+2), 'k', 0:N, Hhat(2,:), 'r--');
ylim([0 2])
xlabel('Sample number, n')
ylabel('tap 1')
legend('True h_n[1]', 'Estimated \hat{h}_n[1|n]')
%% Plotting Kalmar Gain
% === Fig 13.19: Kalman filter gains ===
figure;

subplot(5,1,1)
plot(0:N, Khist(1,:))
xlabel('Sample number, n')
ylabel('K_1[n]')

subplot(5,1,2)
plot(0:N, Khist(2,:))
xlabel('Sample number, n')
ylabel('K_2[n]')
%% Plotting MSE
% === Fig 13.20: Minimum MSE (diagonal of M) ===
figure;
subplot(6,1,1)
plot(0:N, Mdiag(1,:))
ylim([0 0.2])
xlabel('Sample number, n')
ylabel('M_{11}[n|n]')

subplot(6,1,2)
plot(0:N, Mdiag(2,:))
ylim([0 0.2])
xlabel('Sample number, n')
ylabel('M_{22}[n|n]')


