function dz = eom_force(z, F, M, m, b, c, l)
    g = 9.81;
    v   = z(2);
    theta   = z(3);
    omega = z(4);

    Delta = M + m*sin(theta)^2;
    
    num1 = F - c*v + m*l*omega^2*sin(theta) - m*g*sin(theta)*cos(theta) ...
        + b * cos(theta) * omega / l;
    num2 = (M+m)*(m*g*l*sin(theta) - b*omega) ...
        - m*l*cos(theta)*(F - c*v + m*l*omega^2*sin(theta));

    dz = zeros(4,1);
    dz(1) = v;
    dz(2) = num1 / Delta;
    dz(3) = omega;
    dz(4) = num2 / (m*l^2 * Delta);
end