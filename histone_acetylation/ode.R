library(deSolve)

parameters <- c(Vmax_SLC16A3=0.000156, Km_SLC16A3=5.4,
                Km_ACOT1=35.8, Vmax_ACOT1=0.024,
                Km_HDAC1_AcProtein=29.43, Vmax_HDAC1=18.3,
                Km_ACSS2_acetate=350, Vmax_ACSS2=199,
                #Km_ACSS2_CoA
                Km_ACLY=73.8, Vmax_ACLY=2,
                Km_FASN_AcCoA=12.48, Vmax_FASN=0.191,
                #Km_FASN_NADPH
                Km_HMGCS1_AcCoA=1408, Vmax_HMGCS1=0.0447,
                Km_NAT1_AcCoA=268, Vmax_NAT1=1.435
                )
state <- c(Acetate=1,
           AcCoA=1,
           AcProtein=1,
           CoA=1,
           Citrate=1,
           NADPH=1,
           Km_ACSS2_CoA=1,
           Km_FASN_NADPH=1
           )

michaelis.menten <- function(t, y, parameters) {
    with(as.list(c(y, parameters)), {
        dAcetate <- ((Vmax_SLC16A3 * Acetate) / (Km_SLC16A3 + Acetate)) +
                    ((Vmax_ACOT1 * AcCoA) / (Km_ACOT1 + AcCoA)) +
                    ((Vmax_HDAC1 * AcProtein) / (Km_HDAC1_AcProtein + AcProtein)) -
                    ((Vmax_ACSS2 * Acetate * CoA) / ((Km_ACSS2_acetate + Acetate) * (Km_ACSS2_CoA + CoA)))

        dAcCoA <- ((Vmax_ACSS2 * Acetate * CoA) / ((Km_ACSS2_acetate + Acetate) * (Km_ACSS2_CoA + CoA))) +
                  ((Vmax_ACLY * Citrate) / (Km_ACLY + Citrate)) -
                  ((Vmax_FASN * AcCoA * NADPH) / ((Km_FASN_AcCoA + AcCoA) * (Km_FASN_NADPH + NADPH))) -
                  ((Vmax_HMGCS1 * AcCoA) / (Km_HMGCS1_AcCoA + AcCoA)) -
                  ((Vmax_NAT1 * AcProtein) / (Km_NAT1_AcCoA + AcProtein))
        list(c(dAcetate, dAcCoA))
        }
    )
}

times <- seq(0, 100, by = 0.01)
out <- ode(y = state, times = times, func = michaelis.menten, parms = parameters)