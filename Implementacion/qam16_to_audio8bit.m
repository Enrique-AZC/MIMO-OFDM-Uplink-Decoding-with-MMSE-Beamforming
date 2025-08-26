function xsint = qam16_to_audio8bit(qam_data)
    if mod(length(qam_data), 2) ~= 0
        qam_data = qam_data(1:end-1);
    end

    % LSB primero, MSB segundo
    LSB = qam_data(1:2:end);
    MSB = qam_data(2:2:end);

    byte = 16 * MSB + LSB;
    xsint = (double(byte) - 128) / 128;
end
