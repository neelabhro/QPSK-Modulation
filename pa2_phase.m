%% PA 2: Phase and timing sensitivity to noise

%%% Performance: Measure phase est error / timing est error
    % What is correct phase est? 
    % What is error in phase est?
        % COMPARE BER BEFORE AND AFTER PHASE EST! 
    
%%% What is correct timing?
    % 48 samples in with current setup.    
        % Q: Can other samplepoints be better when noise is introduced?
    % tsamp_list contains different sampling points chosen

    
% Plot for different SNR (start with very large)
% Which is most sensitive?
% Why?
% Improved?

%%% TODO: 
% Fix SNR, plot different sample points and different phase errors.

% Gradually change synch and plot BER over time



% 
%%% BER 
figure(99)
for i = 1:length(EbN0_db)
    plot(error, BER(:,i))
    hold on
end
hold off
%%% Perfect BER: 
% EbN0 = 10.^(EbN0_db/10);        % db conversion.. made a simple mistake.
% BER_0 = qfunc(sqrt(2*EbN0));    % BER rate of QPSK, M page 128
% BerT = 0.5 * erfc( sqrt(10 .^ (EbN0_db / 10)) ); % now equal to BER_0
% plot(EbN0_db, BER_0);
% hold off

% title('BER')
% axis([EbN0_db(1) 20 10^-5 max(max(BER_0,BER))])
legend('Simulation with SNR 0' ,'Simulation with SNR 2','Simulation with SNR 4', ...
'Simulation with SNR 6', 'Simulation with SNR 8', 'Simulation with SNR 10')
xlabel('Phase offset rad')
ylabel('BER (Pr)')


%%% Other plot if you want it
figure(100)
for i = 1:length(EbN0_db)
    plot(error, BER_psamp_pre(:,i))
    hold on
end
hold off
%%% Perfect BER: 
% EbN0 = 10.^(EbN0_db/10);        % db conversion.. made a simple mistake.
% BER_0 = qfunc(sqrt(2*EbN0));    % BER rate of QPSK, M page 128
% BerT = 0.5 * erfc( sqrt(10 .^ (EbN0_db / 10)) ); % now equal to BER_0
% plot(EbN0_db, BER_0);
% hold off

% title('BER')
% axis([EbN0_db(1) 20 10^-5 max(max(BER_0,BER))])
legend('Simulation with SNR 5' ,'Simulation with SNR 6','Simulation with SNR 7', ...
'Simulation with SNR 8', 'Simulation with SNR 9', 'Simulation with SNR 10')
xlabel('Phase offset (rad)')
ylabel('BER (Pr)')



%% Simulation 

%%%%%
%%%%% MODIFIED TO LOOP OVER PHASE OR SYNC ERRORS
%%%%%

clc;
close all;

% Initialization

%%% Sync error NOT USED
% sync_error = -10:10;
% error = sync_error;

%%% Phase error 
phase_error = [-12:12] .* (pi/10);
error = phase_error; 

%%% NEEL
% sq_phase_error = zeros(1, length(EbN0_db));
% 
% sq_phase_error(snr_point) = sq_phase_error(snr_point)+abs(phihat)^2;
% 
% sq_phase_error = sqrt(sq_phase_error/nr_blocks);


%%% Variable
EbN0_db = 0:2:10;                     % Eb/N0 values to simulate (in dB)

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



% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);



BER = zeros(length(error),length(EbN0_db));
BER_pre = zeros(length(error),length(EbN0_db));
BER_psamp = zeros(length(error),length(EbN0_db));
BER_psamp_pre = zeros(length(error),length(EbN0_db));

for snr_point = 1:length(EbN0_db)
    % Loop over different values of Eb/No.
    nr_errors = zeros(1, length(error));   % Error counter

    %%% ADDITIONS: list of t_samp to measure accuracy of sync
    tsamp_list = zeros(1,length(error));
    nr_errors_pre = zeros(  1, length(error));   % Pre-phase est error counter
    nr_errors_psamp = zeros(1, length(error)); % Perfect sampling error counter
    nr_errors_psamp_pre = zeros(1,length(error)); % Perfect pre-phase est error counter

    
    % Neel phase code
    sq_phase_error = zeros(1, length(EbN0_db));
    
    for err_point = 1:length(error)

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
    
        %%% PA1 Perfect t_samp:
        t_psamp = 48;            % 10 / 2 * 8 = 40 guard bits, 40-48 is first 
                                % training bit
        %%% PA4 Sync error
%         t_samp = t_psamp + sync_error(err_point);


        % Down sampling. t_samp is the first sample, the remaining samples are all
        % separated by a factor of Q. Only training+data samples are kept.
        r_pre = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);

        r_psamp_pre = mf(t_psamp:Q:t_psamp+Q*(nr_training_bits+nr_data_bits)/2-1);

        % Phase estimation and correction.
%         phihat = phase_estimation(r_pre, b_train);
%         r = r_pre * exp(-1i*phihat);
        r = r_pre * exp(-1i*phase_error(err_point));

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
        nr_errors(err_point) = nr_errors(err_point) + sum(temp);

        %%% Pre-phase correction
        temp_pre = bhat_pre(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
        nr_errors_pre(err_point) = nr_errors_pre(err_point) + sum(temp_pre);

        %%% Perfect sampling
        temp_psamp = bhat_psamp(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
        nr_errors_psamp(err_point) = nr_errors_psamp(err_point) + sum(temp_psamp);

        %%% Perfect sampling pre-phase correction
        temp_psamp_pre = bhat_psamp_pre(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
        nr_errors_psamp_pre(err_point) = nr_errors_psamp_pre(err_point) + sum(temp_psamp_pre);

        %%% tsamp averaging
        tsamp_sum = tsamp_sum + t_samp;

        % Next block.
      end

      %%% tsamp
      tsamp_list(err_point) = tsamp_sum / nr_blocks;    % contains the 
                                                        % different tsamp


    %   two_scatter(0, r, EbN0_db(snr_point), 0)
      % Next Eb/No value.
    end


    % Compute the BER. 
    BER(:,snr_point) = nr_errors / nr_data_bits / nr_blocks;
    BER_pre(:,snr_point) = nr_errors_pre / nr_data_bits / nr_blocks;
    BER_psamp(:,snr_point) = nr_errors_psamp / nr_data_bits / nr_blocks;
    BER_psamp_pre(:,snr_point) = nr_errors_psamp_pre / nr_data_bits / nr_blocks;
    
    BER
    disp('Done!')

end

