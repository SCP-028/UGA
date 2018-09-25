function dy=odesong_201711(t,y)
dy= zeros(8,1)
%%  parameters
                Km_PDHA1_Pyruvate = 64.8; Km_PDHA1_NAD = 33; Km_PDHA1_CoA1 = 4; Kcat_PDHA1 = 69;
                Km_ACSS1_Acetate = 73; Km_ACSS1_CoA1 = 11; Kcat_ACSS1 = 1.9;
                Km_ACOT12_AcCoA2 = 47; Kcat_ACOT12 = 1.48;
                Km_CS_OXa = 5.9; Km_CS_AcCoA1 = 5; Kcat_CS = 167;
                Km_ACO2_Citrate = 480; Km_ACO2_Isocitrate = 120; Kcat_ACO2_1 = 5.3; Kcat_ACO2_2 = 1.1;
                Km_ACLY_Citrate = 78; Km_ACLY_CoA2 = 14; Kcat_ACLY = 2.2;
                Km_IDH2_Isocitrate = 320; Km_IDH2_NAD = 80; Kcat_IDH2 = 30;
                Km_OGDH_AlphaKG = 117; Km_OGDH_NAD = 80; Kcat_OGDH = 30;
                Km_PC_Pyruvate = 220; Km_PC_HCO3 = 3000; Kcat_PC = 60;
                Km_ACSS2_Acetate = 73; Km_ACSS2_CoA2 = 11; Kcat_ACSS2 = 1.9;
                Km_ACACA_AcCoA2 = 34; Km_ACACA_HCO3 = 2100; Kcat_ACACA = 10.1;
                Km_FASN_AcCoA2 = 7; Km_FASN_NADPH = 5; Km_FASN_HCO3 = 7; Kcat_FASN = 2.7;
                Kcat_HDAC1 = 2.8;
                Kcat_HDAC2 = 2;
                Kcat_HDAC3 = 1.5;
                Km_ACAT2_AcCoA2 = 29;
                Km_HMGCS1_AcCoA2 = 14; Kcat_HMGCS1 = 0.041;
                Km_KAT2A_AcCoA2 = 6.7; Kcat_KAT2A = 0.028;
                % Km_ACOT12_AcCoA2 = 47; Kcat_ACOT12 = 1.48;
                % CONSTANTS
                c_Pyruvate = 77;
                c_HCO3 = 11200;
                c_CoA_total_C = 15;
                c_CoA_total_M = 15;
                c_NAD_total = 46.3;  % NADH 22; NAD+ 24.3; FAD 0.078
                v_max_SLC16A3 = 0.14195;
                c_Acetate_blood = 125;
                k_T_Acetate = 0.157;
                n_ATP_NAD = 2.5;
                v_ATP_ss = 0.2;
                c_NADH_ss = 22;
                c_NADPH = 51;
                % % constraints
                % c_CoA1 = max(0; c_CoA_total - c_AcCoA1);
                % c_CoA2 = max(0; c_CoA_total - c_AcCoA2);
                % c_NADH = max(0; c_NAD_total - c_NAD);
                % expression from LIHC
                c_ACACA = 1.3937241;
                c_ACLY = 5.9057259;
                c_ACO2 = 12.5295931;
                c_ACOT12 = 21.4535639;
                %c_ACOT13 = 11.3360529;
                c_ACSS1 = 0.7508832;
                c_ACSS2 = 14.7159301;
                c_CS = 5.3787853;
                c_HDAC1 = 9.6725321;
                c_HDAC2 = 1.3847293;
                c_HDAC3 = 5.3305112;
                c_FASN = 26.7645829;
                c_HMGCS1 = 26.2653794;
                c_IDH2 = 92.8854856;
                c_KAT2A = 3.4616622;
                c_OGDH = 13.7951864;
                c_PC = 48.5013029;
                c_PDHA1 = 13.2815118;
                %c_SLC16A3 = 0.5819519
%% function
%{
y(1) : AcCoA1
y(2) : Citrate
y(3) : Isocitrate
y(4) : AlphaKG
y(5) : OXa
y(6) : AcCoA2
y(7) : Acetate
y(8) : NAD
%}

        dy(1) = ( (Kcat_PDHA1 * c_PDHA1 * c_Pyruvate * y(8) * max(0,c_CoA_total_M -y(1) ) ) / ( (Km_PDHA1_Pyruvate + c_Pyruvate) * (Km_PDHA1_NAD + y(8)) * (Km_PDHA1_CoA1 + max(0, c_CoA_total_M - y(1))) ) ) ...
                + ( (Kcat_ACSS1 * c_ACSS1 * y(7) * max(0,c_CoA_total_M -y(1))) / ((Km_ACSS1_Acetate + y(7)) * (Km_ACSS1_CoA1 + max(0,c_CoA_total_M - y(1))) ) ) ...
                - ( (Kcat_CS * c_CS * y(5) * y(1)) / ( (Km_CS_OXa + y(5)) * (Km_CS_AcCoA1 + y(1)) ) ) ;
        % dcitrate = v_CS - v_ACO2 - v_ACLY
        dy(2) = ( (Kcat_CS * c_CS * y(5) * y(1)) / ( (Km_CS_OXa + y(5)) * (Km_CS_AcCoA1 + y(1)) ) ) ...
                - ( ( (Kcat_ACO2_1 * c_ACO2 * y(2)) /  (Km_ACO2_Citrate + y(2)) ) - ( (Kcat_ACO2_2 * c_ACO2 * y(3)) / (Km_ACO2_Isocitrate + y(3)) ) ) ...
                - ( (Kcat_ACLY * c_ACLY * y(2) * max(0, c_CoA_total_M - y(6))) / ( (Km_ACLY_Citrate + y(2)) * (Km_ACLY_CoA2 + max(0, c_CoA_total_M - y(6))) ) ) ;
        % disocitrate = v_ACO2 - v_IDH2
        dy(3) = ( ( (Kcat_ACO2_1 * c_ACO2 * y(2)) / (Km_ACO2_Citrate + y(2)) ) - ( (Kcat_ACO2_2 * c_ACO2 * y(3)) / (Km_ACO2_Isocitrate + y(3)) ) ) ...
                - ( (Kcat_IDH2 * c_IDH2 * y(3) * y(8)) / ( (Km_IDH2_Isocitrate + y(3)) * (Km_IDH2_NAD + y(8)) ) ) ;
        % dalpha-KG = v_IDH2 - v_OGDH
        dy(4) = ( (Kcat_IDH2 * c_IDH2 * y(3) * y(8)) / ( (Km_IDH2_Isocitrate + y(3)) * (Km_IDH2_NAD + y(8)) ) ) ...
                - ( (Kcat_OGDH * c_OGDH * y(4) * y(8)) / ( (Km_OGDH_AlphaKG + y(4)) * (Km_OGDH_NAD + y(8)) ) ) ;
        % dOXa = v_OGDH + v_PC - v_CS
        dy(5) =  ( (Kcat_OGDH * c_OGDH * y(4) * y(8)) / ( (Km_OGDH_AlphaKG + y(4)) * (Km_OGDH_NAD + y(8)) ) ) ...
                + ( (Kcat_PC * c_PC * c_Pyruvate * c_HCO3) / ( (Km_PC_Pyruvate + c_Pyruvate) * (Km_PC_HCO3 + c_HCO3) ) ) ...
                - ( (Kcat_CS * c_CS * y(5) * y(1)) / ( (Km_CS_OXa + y(5)) * (Km_CS_AcCoA1 + y(1)) ) ) ;
        % dAc-CoA2 = v_ACLY + v_ACSS2 - v_ACACA - v_FASN - v_HMGCS1 - v_KAT2A - v_ACOT12
        dy(6) =  ( (Kcat_ACLY * c_ACLY * y(2) * max(0, c_CoA_total_C - y(6))) / ( (Km_ACLY_Citrate +  y(2)) * (Km_ACLY_CoA2 + max(0, c_CoA_total_C - y(6))) ) ) ...
                + ( (Kcat_ACSS2 * c_ACSS2 * y(7) *  max(0, c_CoA_total_C - y(6))) / ( (Km_ACSS2_Acetate + y(7)) * (Km_ACSS2_CoA2 + max(0, c_CoA_total_C - y(6))) ) ) ...
                - ( (Kcat_ACACA * c_ACACA * y(6) * c_HCO3) / ( (Km_ACACA_AcCoA2 + y(6)) * (Km_ACACA_HCO3 + c_HCO3) ) ) ...
                - ( (2 * Kcat_FASN * c_FASN * y(6) ^ 2 * c_HCO3 *  c_NADPH ^ 2) / ( (Km_FASN_AcCoA2 + y(6)) ^ 2 * (Km_FASN_HCO3 + c_HCO3) * (Km_FASN_NADPH + c_NADPH) ^ 2 ) ) ...
                - ( (3 * Kcat_HMGCS1 * c_HMGCS1 * y(6) ^ 3) / ( (Km_ACAT2_AcCoA2 + y(6)) ^ 2 * (Km_HMGCS1_AcCoA2 + y(6)) ) ) ...
                - ( (Kcat_KAT2A *  c_KAT2A * y(6)) /  (Km_KAT2A_AcCoA2 + y(6)) ) ...
                - ( (Kcat_ACOT12 * c_ACOT12 *  y(6)) / (Km_ACOT12_AcCoA2 + y(6)) ) ;
        % dacetate = v_SLC16A3 + v_HDAC1 + v_HDAC2 + v_HDAC3 + v_ACOT12 - v_ACSS1 - v_ACSS2
        dy(7) = ( v_max_SLC16A3 * ( ( c_Acetate_blood / (c_Acetate_blood + k_T_Acetate) ) - ( y(7) / (y(7) + k_T_Acetate) ) ) ) ...
                + (Kcat_HDAC1 * c_HDAC1) ...
                +  (Kcat_HDAC2 * c_HDAC2) ...
                + (Kcat_HDAC3 * c_HDAC3) ...
                + ( (Kcat_ACOT12 * c_ACOT12 * y(6)) / (Km_ACOT12_AcCoA2 + y(6)) ) ...
                - ( (Kcat_ACSS1 * c_ACSS1 * y(7) * max(0, c_CoA_total_C - y(1))) / ( (Km_ACSS1_Acetate + y(7)) * (Km_ACSS1_CoA1 + max(0, c_CoA_total_C - y(1))) ) ) ...
                - ( (Kcat_ACSS2 *  c_ACSS2 * y(7) * max(0, c_CoA_total_C - y(6))) / ( (Km_ACSS2_Acetate + y(7)) * (Km_ACSS2_CoA2 +  max(0, c_CoA_total_C - y(6))) ) ) ;
        % dNAD = v_NAD - v_PDHA1 - v_IDH2 - v_OGDH
        dy(8) =  ( (v_ATP_ss / n_ATP_NAD) * (max(0, c_NAD_total - y(8)) / c_NADH_ss) ) ...
                - ( (Kcat_PDHA1 * c_PDHA1 *  c_Pyruvate * y(8)) / ( (Km_PDHA1_Pyruvate + c_Pyruvate) * (Km_PDHA1_NAD + y(8))  ) ) ...
                - ( (Kcat_IDH2 * c_IDH2 * y(3) * y(8)) / ( (Km_IDH2_Isocitrate + y(3)) * (Km_IDH2_NAD + y(8)) ) ) ...
                - ( (Kcat_OGDH *  c_OGDH * y(4) * y(8)) / ( (Km_OGDH_AlphaKG + y(4)) * (Km_OGDH_NAD + y(8)) ) ) ;
              