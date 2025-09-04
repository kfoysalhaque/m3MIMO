function rx_bits = qpsk_demod_no_hamming(rx_symbols)
    % QPSK Demodulation function without Hamming decoding
    % Input: rx_symbols - Received QPSK symbols
    % Output: rx_bits - Demodulated bits

    % Preallocate for speed
    rx_bits = zeros(1, 2 * length(rx_symbols));
    
    % Iterate over each received symbol
    for i = 1:length(rx_symbols)
        real_part = real(rx_symbols(i));
        imag_part = imag(rx_symbols(i));
        
        % Determine the bits based on the quadrant
        if real_part > 0 && imag_part > 0
            % First quadrant
            rx_bits(2*i-1:2*i) = [0, 0];
        elseif real_part > 0 && imag_part < 0
            % Fourth quadrant
            rx_bits(2*i-1:2*i) = [0, 1];
        elseif real_part < 0 && imag_part > 0
            % Second quadrant
            rx_bits(2*i-1:2*i) = [1, 0];
        else
            % Third quadrant
            rx_bits(2*i-1:2*i) = [1, 1];
        end
    end
end
