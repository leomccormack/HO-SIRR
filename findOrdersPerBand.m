function orderPerBin = findOrdersPerBand(centerFreqs, freqLimits, procOrder)

    nBands = size(centerFreqs,1);
    if length(procOrder)==1
        orderPerBin = ones(nBands,1)*procOrder;
    else       
        nLims = length(freqLimits);
        if length(procOrder)~=nLims+1, error('The band limits for the specified orders should be length(orders)-1'); end
        
        orderPerBin = zeros(nBands,1);
        % first band from DC to lim1
        selector = (centerFreqs<=freqLimits(1));
        orderPerBin(selector==1) = procOrder(1);
        % intermediate bands from lim1 to limN
        for band = 1:length(freqLimits)-1
            selector = ((centerFreqs>freqLimits(band)).*(centerFreqs<=freqLimits(band+1)));
            orderPerBin(selector==1) = procOrder(band + 1);
        end
        % last band from limN to Nyquist
        selector = (centerFreqs>freqLimits(end));
        orderPerBin(selector==1) = procOrder(end);
    end
end