function dCdt = gene_expression_ODE(t,C,k_DNA_Nap,k_DNA_Nrf,d_mRNA,...
    k_mRNA,m_Nap,m_Nrf,d_Nap,d_Nrf)

    DNA = C(1); mRNA_Nap = C(2); mRNA_Nrf = C(3); Nasc_Nap = C(4);
    Nasc_Nrf = C(5); Nap = C(6); Nrf = C(7);
    
    dDNA = 0;
    dmRNA_Nap = k_DNA_Nap*DNA - d_mRNA*mRNA_Nap;
    dmRNA_Nrf = k_DNA_Nrf*DNA - d_mRNA*mRNA_Nrf;
    dNasc_Nap = k_mRNA*mRNA_Nap; % IS TRANSLATION RATE SPECIFIC TO ENZYME?
    dNasc_Nrf = k_mRNA*mRNA_Nrf;
    dNap = m_Nap*Nasc_Nap - d_Nap*Nap;
    dNrf = m_Nrf*Nasc_Nrf - d_Nrf*Nrf;
    
    dCdt = [dDNA, dmRNA_Nap, dmRNA_Nrf, dNasc_Nap, dNasc_Nrf, dNap, dNrf]';
end