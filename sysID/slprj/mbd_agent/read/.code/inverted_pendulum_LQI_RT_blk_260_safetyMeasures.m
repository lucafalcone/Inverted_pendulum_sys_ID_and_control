function u = safety_measures(u, x, theta, borderProtection, inputLimitation, coneLimitation, inputDisable)

% border protection parameters
forceFieldStrenght = 2;
leftborder = 1;
rightborder = -1;
safety_margin = 0.30;

% input limitation parameters
inputLimit = 2.5;

% cone limitation parameters
coneWidth = 35; % deg

% execution
if inputLimitation
    if u > inputLimit
        u = inputLimit;
    elseif u < -inputLimit
        u = -inputLimit;
    end
end

if coneLimitation
    if abs(theta) > coneWidth/180*pi
        u = 0;
    end
end

if borderProtection 
    if(x < rightborder + safety_margin)
        u = forceFieldStrenght;
    end
    if(x > leftborder - safety_margin)
        u = -forceFieldStrenght;
    end
end

if inputDisable
    u = 0;
end

end