function dCdt = QSSA(t,C,k1,km1,k2)
    p = C(1); s = C(2);
    Km = (km1+k2)/k1;
   
    dpdt = k2*e0*s/(Km+s);
    dsdt = -dpdtl
end
