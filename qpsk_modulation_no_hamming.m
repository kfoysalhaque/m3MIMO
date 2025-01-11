function tx_symbols = qpsk_mod_no_hamming(bits)
    % QPSK Modulation function without Hamming encoding
    % Input: bits - Input bits to be modulated
    % Output: tx_symbols - Modulated QPSK symbols

    % Ensure the input bits length is even for QPSK
    if mod(length(bits), 2) ~= 0
        error('Input bits length must be even for QPSK modulation.');
    end

    % Map bits to QPSK symbols
    symbols_I = 2 * bits(1:2:end) - 1; % Map bits to I-channel: 0 -> -1, 1 -> 1
    symbols_Q = 2 * bits(2:2:end) - 1; % Map bits to Q-channel: 0 -> -1, 1 -> 1

    % Combine I and Q channel to form complex symbols
    tx_symbols = symbols_I + 1j * symbols_Q;
end
