function dYdt = exampleODE(t,Y,k1,km1,k2)
    dsdt = Y(3) * km1 - Y(1) * Y(2) * k1;
    dpdt = k2 * Y(3);
    dcdt = Y(1) * Y(2) * k1 - Y(3) * (km1 + k2);
    dedt = Y(3) * (k2 + km1) - Y(1) * Y(2) * k1;

    dYdt = [dsdt dedt dcdt dpdt]';
end
