
%% Clear all workspace variables except for sdr objects

% List of variables to keep
varsToKeep = {'sdr0', 'sdr1', 't_hat_no_bfm', 't_hat_bfm_IC', 't_hat_bfm_NIC', 'input_bits', 'nFFT',  'scMin', 'scMax'};

% Get the list of all variables in the workspace
allVars = whos;

% Construct the command to clear all variables except those in varsToKeep
clearCommand = 'clear ';
for k = 1:length(allVars)
    if ~ismember(allVars(k).name, varsToKeep)
        clearCommand = [clearCommand allVars(k).name ' '];
    end
end

% Execute the clear command
eval(clearCommand);

clear allVars clearCommand k

%% Choose the stream

t_hat_plot = squeeze(t_hat_bfm_NIC(:, 1));

%% BER

rx_symbols =  t_hat_plot(nFFT/2 + 1 + scMin:nFFT/2 + scMax, 1);
[rx_bits, received_bits] = qpsk_demodulation(rx_symbols);
errors = sum(received_bits ~= input_bits.');
BER = errors / size(input_bits,1)