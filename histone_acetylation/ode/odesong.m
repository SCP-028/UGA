function dy=odesong(t,y)
dy= zeros(11,1)
 dy(1) = ((69 * 13.2815118 * y(2) * y(3)) /((64.8 + y(2)) * (33 + y(3))))
         + ((1.9 * 0.7508832 * y(4) * (15 -y(1) )) / ((73 +  y(4)) * (11 + (15 - y(1) ))))
         - ((167 * 5.3787853 * y(6) * y(1)) / ((5.9 + y(6)) * (5 + y(1) ))) 
        % citrate = CS - ACO2 - ACLY
        % - ((1.48 * 11.3360529 * y(1)) / (47 + y(1))) 
 dy(7) = ((167 * 5.3787853 * y(6) * y(1)) / ((5.9 + y(6)) * (5 + y(1))))
          +((1.1 * 12.5295931 * y(8)) / (120 + y(8)) )
          - ((5.3 * 12.5295931 * y(7)) / (480 + y(7) ))
          - ((2.2* 5.9057259 * y(7)) / (78 + y(7)))
        % isocitrate = ACO2 - IDH2
 dy(8) = ((5.3 * 12.5295931 * y(7))/(480 + y(7)) )
        -((1.1 * 12.5295931 * y(8))/(120 + y(8)) )
        - ((30 * 92.8854856 * y(8) * y(3)) / ((2400 + y(8)) * (80 + y(3))))
        % alpha-KG = IDH2 - OGDH
 dy(9) = ((30 * 92.8854856 * y(8) * y(3)) / ((2400 + y(8)) * (80 + y(3))))
          - ((30 * 13.7951864 * y(9) * y(3)) / ((4000 + y(9)) * (80 + y(3))))
        % OXa = OGDH + PC - CS
 dy(6) = ((30 * 13.7951864 * y(9) * y(3)) / ((4000 + y(9)) * (80 + y(3)))) 
        + ((60 * 48.5013029 * y(2) * y(10)) / ((220 + y(2)) * (3000 + y(10))))
        - ((167 * 5.3787853 * y(6) * y(1)) / ((5.9 + y(6)) * (5 + y(1))))
        % Ac-CoA2 = ACLY + ACSS2 - ACACA - HMGCS1 - KAT2A - ACOT12
 dy(11) = ((2.2* 5.9057259 * y(7)) / (78 + y(7)))
        + ((1.9 * 14.7159301 * y(4) * (15 - y(11))) / ((73 + y(4)) * (11 + (15 - y(11)))))
        - ((2 * 10.1 * y(11) * 1.3937241 * y(10)) / ((34 + y(11)) * (2100 + y(10))))
        - ((0.041 * 26.2653794 * y(11) * y(11)) / ((14 + y(11)) * (14 + y(11))))
        - ((0.011 * 3.4616622 *y(11))/(6.7+y(11))) 
        - ((7.26 * 21.4535639 * y(11)) / (500 + y(11)))
        % acetate = SLC16A3 + HDAC1 + HDAC2 + HDAC3 + ACOT13 + ACOT12 - ACSS1 - ACSS2
 dy(4) = 88.07 + (2.8 * 9.6725321) + (2 * 1.3847293) + (1.5 * 5.3305112) 
        + ((7.26 * 21.4535639 * y(11)) / (500 + y(11))) 
        - ((1.9 * 0.7508832 * y(4) * (15 - y(1))) / ((73 + y(4)) * (11 + (15 - y(1)))))
        - ((1.9 * 14.7159301 * y(4) * (15 - y(11))) / ((73 + y(4)) * (11 + (15 - y(11)))))
        % NAD = v_NAD - PDHA1 - IDH2 - OGDH
        % + ((1.48 * 11.3360529 * y(1)) / (47 + y(1)))
 dy(3) = ((0.20844 - (0.0012 * y(3))) / 85) 
        - ((69 * 13.2815118 * y(2) * y(3)) / ((64.8 + y(2)) * (33 + y(3))))
        - ((30 * 92.8854856 * y(8) * y(3)) / ((2400 + y(8)) * (80 + y(3)))) 
        - ((30 * 13.7951864 * y(9) * y(3)) / ((4000 + y(9)) * (80 + y(3))))
       
              
%%  annotate
%{
y(1) : AC-COA_1
y(2) : Pyruvate
y(3) : NAD
y(4) : acetate
y(5) : coa
y(6) : oxa
y(7) : citrate
y(8) : isocitrate
y(9) : a-kg
y(10) :HCO3
y(11) :ac-coa_2

equations
dy(1) = ((Kcat_PHDA1 * PDHA1 * y(2) * y(3)) /((Km_PDHA1_Pyruvate + y(2)) * (Km_PDHA1_NAD + y(3))))
         + ((Kcat_ACSS1 * ACSS1 * y(4) * (15 -y(1) )) / ((Km_ACSS1_Acetate +  y(4)) * (Km_ACSS1_CoA + (15 - y(1) ))))
         - ((Kcat_ACOT13 * ACOT13 * y(1)) / (Km_ACOT13_AcCoA1 + y(1))) 
         - ((Kcat_CS * CS * y(6) * y(1)) / ((Km_CS_OXa + y(6)) * (Km_CS_AcCoA1 + y(1) )))
        % citrate = CS - ACO2 - ACLY
 dy(7) = ((Kcat_CS * CS * y(6) * y(1)) / ((Km_CS_OXa + y(6)) * (Km_CS_AcCoA1 + y(1))))
          +((Kcat_ACO2_2 * ACO2 * y(8)) / (Km_ACO2_Isocitrate + y(8)) )
          - ((Kcat_ACO2_1 * ACO2 * y(7)) / (Km_ACO2_Citrate + y(7) ))
          - ((Kcat_ACLY * ACLY * y(7)) / (Km_ACLY_Citrate + y(7)))
        % isocitrate = ACO2 - IDH2
 dy(8) = ((Kcat_ACO2_1 * ACO2 * y(7))/(Km_ACO2_Citrate + y(7)) )
        -((Kcat_ACO2_2 * ACO2 * y(8))/(Km_ACO2_Isocitrate + y(8)) )
        - ((Kcat_IDH2 * IDH2 * y(8) * y(3)) / ((Km_IDH2_Isocitrate + y(8)) * (Km_IDH2_NAD + y(3))))
        % alpha-KG = IDH2 - OGDH
 dy(9) = ((Kcat_IDH2 * IDH2 * y(8) * y(3)) / ((Km_IDH2_Isocitrate + y(8)) * (Km_IDH2_NAD + y(3))))
          - ((Kcat_OGDH * OGDH * y(9) * y(3)) / ((Km_OGDH_AlphaKG + y(9)) * (Km_OGDH_NAD + y(3))))
        % OXa = OGDH + PC - CS
 dy(6) = ((Kcat_OGDH * OGDH * y(9) * y(3)) / ((Km_OGDH_AlphaKG + y(9)) * (Km_OGDH_NAD + y(3)))) 
        + ((Kcat_PC * PC * y(2) * y(10)) / ((Km_PC_Pyruvate + y(2)) * (Km_PC_HCO3 + y(10))))
        - ((Kcat_CS * CS * y(6) * y(1)) / ((Km_CS_OXa + y(6)) * (Km_CS_AcCoA1 + y(1))))
        % Ac-CoA2 = ACLY + ACSS2 - ACACA - HMGCS1 - KAT2A - ACOT12
 dy(11) = ((Kcat_ACLY * ACLY * y(7)) / (Km_ACLY_Citrate + y(7)))
        + ((Kcat_ACSS2 * ACSS2 * y(4) * (15 - y(11))) / ((Km_ACSS2_Acetate + y(4)) * (Km_ACSS2_CoA + (15 - y(11)))))
        - ((2 * Kcat_ACACA * y(11) * ACACA * y(10)) / ((Km_ACACA_AcCoA2 + y(11)) * (Km_ACACA_HCO3 + y(10))))
        - ((Kcat_HMGCS1 * HMGCS1 * y(11) * y(11)) / ((Km_HMGCS1_AcCoA2 + y(11)) * (Km_HMGCS1_AcCoA2 + y(11))))
        - ((Kcat_KAT2A * KAT2A*y(11))/(Km_KAT2A_AcCoA2+y(11))) 
        - ((Kcat_ACOT12 * ACOT12 * y(11)) / (Km_ACOT12_AcCoA2 + y(11)))
        % acetate = SLC16A3 + HDAC1 + HDAC2 + HDAC3 + ACOT13 + ACOT12 - ACSS1 - ACSS2
 dy(4) = 88.07 + (Kcat_HDAC1 * HDAC1) + (Kcat_HDAC2 * HDAC2) + (Kcat_HDAC3 * HDAC3)
        + ((Kcat_ACOT12 * ACOT12 * y(11)) / (Km_ACOT12_AcCoA2 + y(11))) 
        - ((Kcat_ACSS1 * ACSS1 * y(4) * (15 - y(1))) / ((Km_ACSS1_Acetate + y(4)) * (Km_ACSS1_CoA + (15 - y(1)))))
        - ((Kcat_ACSS2 * ACSS2 * y(4) * (15 - y(11))) / ((Km_ACSS2_Acetate + y(4)) * (Km_ACSS2_CoA + (15 - y(11)))))
        % NAD = v_NAD - PDHA1 - IDH2 - OGDH  
        %  + ((Kcat_ACOT13 * ACOT13 * y(1)) / (Km_ACOT13_AcCoA1 + y(1))) 
 dy(3) = ((0.20844 - (0.0012 * y(3))) / 85) 
        - ((Kcat_PDHA1 * PDHA1 * y(2) * y(3)) / ((Km_PDHA1_Pyruvate + y(2)) * (Km_PDHA1_NAD + y(3))))
        - ((Kcat_IDH2 * IDH2 * y(8) * y(3)) / ((Km_IDH2_Isocitrate + y(8)) * (Km_IDH2_NAD + y(3)))) 
        - ((Kcat_OGDH * OGDH * y(9) * y(3)) / ((Km_OGDH_AlphaKG + y(9)) * (Km_OGDH_NAD + y(3))))
       
%}

        
        