function c = pb_constants()
c.R = 8.3143;     % universal gas constant    [J/(mol*K)]
c.T = 310.15;     % absolute temperature      [K]  (37 degC)
c.F = 96.4867;    % Faraday constant          [C/mmol]
c.RTF = c.R*c.T/c.F;   % ~26.73 mV at 37 degC
end
