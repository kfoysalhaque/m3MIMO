function [rx_bits, corrected_bits] = qpsk_demodulation2(rx_symbols)
    % QPSK Demodulation function with Hamming decoding
    % Input: rx_symbols - Received QPSK symbols
    % Output: rx_bits - Demodulated bits
    %         corrected_bits - Corrected bits after Hamming decoding

    % Demodulate QPSK symbols
    rx_bits_I = real(rx_symbols) > 0; % Demodulate I-channel: 1 -> 1, 0 -> 0
    rx_bits_Q = imag(rx_symbols) > 0; % Demodulate Q-channel: 1 -> 1, 0 -> 0

    % Combine demodulated bits
    rx_bits = zeros(1, 2 * length(rx_symbols));
    rx_bits(1:2:end) = rx_bits_I;
    rx_bits(2:2:end) = rx_bits_Q;

    % Ensure the received bits length is a multiple of 7 for Hamming(7,4) decoding
    n = length(rx_bits);
    if mod(n, 7) ~= 0
        % Pad with zeros if the length is not a multiple of 7
        padding_length = 7 - mod(n, 7);
        rx_bits = [rx_bits, zeros(1, padding_length)];
    end

    % Decode received bits using Hamming code
    corrected_bits = decode(rx_bits, 7, 4, 'hamming/binary');
end
