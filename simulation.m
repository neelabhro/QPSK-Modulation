% Skeleton code for simulation chain

% History:
%   2000-06-28  written /Stefan Parkvall
%   2001-10-22  modified /George Jongren

clear;
clc;
close all;

% Initialization
EbN0_db = -5:20;                     % Eb/N0 values to simulate (in dB)
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
pulse_shape = ones(1, Q);
%pulse_shape = root_raised_cosine(Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.
nr_errors = zeros(1, length(EbN0_db));   % Error counter
for snr_point = 1:length(EbN0_db)
  
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
%     t_samp = 48;            % 10 / 2 * 8 = 40 guard bits, 40-48 is first 
                            % training bit
    
    
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r_pre = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);

    % Phase estimation and correction.
    phihat = phase_estimation(r_pre, b_train);
    r = r_pre * exp(-1i*phihat);
    
    %%% PA1 Perfect phase: calculate phase from noise?
    
    
    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat = detect(r);
    %bhat = differential_detect(r);
    
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
    temp=bhat(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors(snr_point) = nr_errors(snr_point) + sum(temp);

    % Next block.
  end
  
  %%% PA3: Plot signal constellation
%   two_scatter(r,r_pre, snr_point)
  
  % Next Eb/No value.
end

% Compute the BER. 
BER = nr_errors / nr_data_bits / nr_blocks;

%% PA 1: BER
%%% BER 
figure(2)
plot(EbN0_db, BER)
hold on

%%% Perfect BER: 
EbN0 = 10.^(EbN0_db/20);        % db conversion
BER_0 = qfunc(sqrt(2*EbN0));    % BER rate of QPSK, M page 128
BerT = 0.5 * erfc( sqrt(10 .^ (EbN0_db / 10)) );
plot(EbN0_db, BerT);
hold off

title('BER')
axis([-5 20 10^-5 1])
legend('Simulation', 'Theoretical value')
xlabel('SNR (dB)')
ylabel('BER (Pr)')
%%% Exact Error probability of QPSK: why?

% P = 2*qfunc(sqrt(2*EbN0_db))-(qfunc(sqrt(EbN0_db))).^2;     % M page 118
% plot(EbN0_db, P)
% 
% hold off

% figsaver(1:2,'PA1')

%% PA 2: Phase and timing sensitivity to noise

%%% Measure phase est error / timing est error
    % What is correct phase est? 
    
peek = 4;
constallation_scatter(tx(peek),3)
hold on
constallation_scatter(rx(peek),3)
hold off
    
%%% What is correct timing?
    % 48 samples in
    
% Plot for different SNR (start with very large)
% Why?
% Improved?

%% PA 3: Signal constellations and noise


%% PA 4: PSD with different pulses

%%% Periodogram
% 
% function pa4(r)
%     Prr = periodogram(r);
%     figure(1)
%     plot(Prr)
%     title('Periodogram')
% end


%% Save figures

% enter figure numbers and folder to where to save them
% figsaver(1:3, 'testfig')

%% Functions

function one_scatter(data, fig)
% one_scatter(data,figure)
% 
% Plots real and imaginary part of data using scatter plot with figure 
% number according to input.
%
% Input: 
%  data         = the data you want to plot
%  fig       = the number of the figure plotted
%
% Output: 
%   None

    figure(fig)
    scatter(real(data), imag(data))
    
    title('Constallation plot')

end


function two_scatter(tx, rx, snr_point)
    fignum = snr_point + 5;
    figure(fignum)
    scatter(real(tx), imag(tx))
    hold on
    scatter(real(rx), imag(rx))
    hold off
    title('Recieved signal constallation with SNR = ',  snr_point)
%     axis([-3 3 -3 3]);

end
