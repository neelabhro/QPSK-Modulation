%% PA 4: PSD with different pulses

%%% Periodogram

% Prr_raised = periodogram(r_raised);
Prr_rect = periodogram(r);
% Prr_trig = periodogram(r);

% 
% N = length(r);
% R = fft(r);
% Pr = 1/N*abs(R).^2;
% 

%%
figure(1)
plot(1:length(Prr_rect), Prr_trig, 'LineWidth',3)
hold on
title('Periodogram with re')
plot(1:length(Prr_rect), Prr_rect) 
legend('Raised cosine pulse', 'Rectangular pulse');
hold off


% 
% % figure(2)
% % plot(Pr)
% % title('Periodogram')
% 
% figsaver(1,'PA4')

% %%% Spectrogram:
% r_raised = r;
% % s = spectrogram(r);
% figure(1)
% spectrogram(r_rect,'yaxis')
% 
% figure(2)
% spectrogram(r_raised,'yaxis')


% plot(abs(s))

%% Simulation
% Skeleton code for simulation chain

% History:
%   2000-06-28  written /Stefan Parkvall
%   2001-10-22  modified /George Jongren

clc;

% Initialization

%%% Variable
EbN0_db = -20:10;                     % Eb/N0 values to simulate (in dB)

%%% Static
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 10;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 1000;                % Size of each data sequence (in nr bits)
nr_training_bits = 100;             % Size of training sequence (in nr bits)
nr_blocks = 50;                     % The number of blocks to simulate
Q = 8;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.

%%% Variable
% Step function
pulse_shape = ones(1, Q);

% Raised cosine
% pulse_shape = root_raised_cosine(Q);

% Triangular?
% pulse_shape = bartlett(Q);

figure(66)
plot(pulse_shape)


% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.
nr_errors = zeros(1, length(EbN0_db));   % Error counter

%%% ADDITIONS: list of t_samp to measure accuracy of sync
tsamp_list = zeros(1,length(EbN0_db));
nr_errors_pre = zeros(  1, length(EbN0_db));   % Pre-phase est error counter
nr_errors_psamp = zeros(1, length(EbN0_db)); % Perfect sampling error counter
nr_errors_psamp_pre = zeros(1,length(EbN0_db)); % Perfect pre-phase est error counter

for snr_point = 1:length(EbN0_db)
  
  %%% tsamp sum
  tsamp_sum = 0;
  
  % Loop over several blocks to get sufficient statistics.
  for blk = 1:nr_blocks

    %%%
    %%% Transmitter
    %%%

    % Generate training sequence.
    b_train = training_sequence(nr_training_bits);
    
    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);

    % Generate guard sequence.
    b_guard = random_data(nr_guard_bits);
 
    % Multiplex training and data into one sequence.
    b = [b_guard b_train b_data b_guard];
    
    % Map bits into complex-valued QPSK symbols.
    d = qpsk(b);
    %d = differential_qpsk(b);

    % Upsample the signal, apply pulse shaping.
    tx = upfirdn(d, pulse_shape, Q, 1);
    
    %%%%%%%%%%%%%% Try to use this?
%     tx = txfilter(d);

    %%%
    %%% AWGN Channel
    %%%
    
    % Compute variance of complex noise according to report.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);

    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(tx))+j*randn(size(tx)));

    % Received signal.
    rx = tx + n;
    
    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start=1+Q*nr_guard_bits/2;
    t_end=t_start+50;
    t_samp = sync(mf, b_train, Q, t_start, t_end);
%     t_samp = 48;
    %%% PA1 Perfect t_samp:
    t_psamp = 48;            % 10 / 2 * 8 = 40 guard bits, 40-48 is first 
                            % training bit
    
    
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r_pre = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);

    r_psamp_pre = mf(t_psamp:Q:t_psamp+Q*(nr_training_bits+nr_data_bits)/2-1);
    
    % Phase estimation and correction.
    phihat = phase_estimation(r_pre, b_train);
    
    r = r_pre * exp(-1i*phihat);
    
    r_psamp = r_psamp_pre * exp(-1i*phihat);      % with or without phase corrrection
    
    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat = detect(r);
    %bhat = differential_detect(r);
    bhat_pre = detect(r_pre);           % pre-phase correction bhat
    bhat_psamp = detect(r_psamp);
    bhat_psamp_pre = detect(r_psamp_pre);
    
    % Reference points:   
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
    temp=bhat(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors(snr_point) = nr_errors(snr_point) + sum(temp);
    
    %%% Pre-phase correction
    temp_pre = bhat_pre(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors_pre(snr_point) = nr_errors_pre(snr_point) + sum(temp_pre);
    
    %%% Perfect sampling
    temp_psamp = bhat_psamp(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors_psamp(snr_point) = nr_errors_psamp(snr_point) + sum(temp_psamp);
    
    %%% Perfect sampling pre-phase correction
    temp_psamp_pre = bhat_psamp_pre(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors_psamp_pre(snr_point) = nr_errors_psamp_pre(snr_point) + sum(temp_psamp_pre);
    
    %%% tsamp averaging
    tsamp_sum = tsamp_sum + t_samp;

    % Next block.
  end
  
  %%% tsamp
  tsamp_list(snr_point) = tsamp_sum / nr_blocks;    % contains the 
                                                    % different tsamp
  
  
%   two_scatter(0, r, EbN0_db(snr_point), 0)
  % Next Eb/No value.
end

% Compute the BER. 
BER = nr_errors / nr_data_bits / nr_blocks;
BER_pre = nr_errors_pre / nr_data_bits / nr_blocks;
BER_psamp = nr_errors_psamp / nr_data_bits / nr_blocks;
BER_psamp_pre = nr_errors_psamp_pre / nr_data_bits / nr_blocks;

disp('Done!')


%%% BER 
figure(2)
plot(EbN0_db, BER)
hold on

%%% Perfect BER: 
EbN0 = 10.^(EbN0_db/10);        % db conversion.. made a simple mistake.
BER_0 = qfunc(sqrt(2*EbN0));    % BER rate of QPSK, M page 128
BerT = 0.5 * erfc( sqrt(10 .^ (EbN0_db / 10)) ); % now equal to BER_0
plot(EbN0_db, BER_0);
hold off

% title('BER')
axis([EbN0_db(1) 20 10^-5 max(max(BER_0,BER))])
legend('Simulation', 'Theoretical value')
xlabel('SNR (dB)')
ylabel('BER (Pr)')


% figsaver(1:2,'PA1')
