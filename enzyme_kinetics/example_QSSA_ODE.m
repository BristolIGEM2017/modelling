function dCdt = example_QSSA_ODE(t,C,e0,k1,km1,k2)
    p = C(2); s = C(1);
    Km = (km1+k2)/k1;
    
    dpdt = k2*e0*s/(Km+s);
    dsdt = -dpdt;
    
    dCdt = [dsdt; dpdt];
end