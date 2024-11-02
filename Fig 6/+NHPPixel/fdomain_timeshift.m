function X = fdomain_timeshift(x, shifts)
    % based on https://github.com/int-brain-lab/ibl-neuropixel/blob/a56a4953155bb4036a9d7f33e373d837684f56e6/src/neurodsp/fourier.py#L206   

    % x is input signal channels x time x ...
    % shift in samples, positive shifts forward

    arguments
        x % channels x time x ...
        shifts (:, 1)
    end

    sz = size(x);
    ns = sz(2);

    % create a vector that contains a 1 sample shift on the axis
    dephas = zeros(1, ns);
    dephas(1, 2) = 1;
    dephas = fft(dephas, [], 2);

    % fft the data along the axis
    X = fft(x, [], 2);
    
    % apply the shift (s) to the fft angle to get the phase shift and broadcast
    X = X .* exp(1j .* angle(dephas) .* shifts);
    
    X = real(ifft(X, ns, 2));
end