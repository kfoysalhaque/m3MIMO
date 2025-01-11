function rx_bits = qpsk_demod_no_hamming(rx_symbols)
    % QPSK Demodulation function without Hamming decoding
    % Input: rx_symbols - Received QPSK symbols
    % Output: rx_bits - Demodulated bits

    % Demodulate QPSK symbols
    rx_bits_I = real(rx_symbols) > 0; % Demodulate I-channel: 1 -> 1, 0 -> 0
    rx_bits_Q = imag(rx_symbols) > 0; % Demodulate Q-channel: 1 -> 1, 0 -> 0

    % Combine demodulated bits
    rx_bits = zeros(1, 2 * length(rx_symbols));
    rx_bits(1:2:end) = rx_bits_I;
    rx_bits(2:2:end) = rx_bits_Q;
end
