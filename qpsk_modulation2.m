function [tx_symbols, encoded_bits] = qpsk_modulation2(bits)
    % QPSK Modulation function with Hamming encoding
    % Input: bits - Input bits to be modulated
    % Output: tx_symbols - Modulated QPSK symbols
    %         encoded_bits - Encoded bits with Hamming code

    % Ensure the input bits length is a multiple of 4 for Hamming(7,4) encoding
    n = length(bits);
    if mod(n, 4) ~= 0
        error('Input bits length must be a multiple of 4');
    end

    % Encode input bits with Hamming code
    encoded_bits = encode(bits, 7, 4, 'hamming/binary');

    % Generate QPSK symbols from encoded bits
    symbols_I = 2*encoded_bits(1:2:end) - 1; % Map bits to I-channel: 0 -> -1, 1 -> 1
    symbols_Q = 2*encoded_bits(2:2:end) - 1; % Map bits to Q-channel: 0 -> -1, 1 -> 1

    % Combine I and Q channel to form complex symbols
    tx_symbols = symbols_I + 1j * symbols_Q;
end