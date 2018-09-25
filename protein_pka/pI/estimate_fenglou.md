# Estimating cellular pH with protein pK values
For each acidic amino acid (ASP and GLU), we have:

$$HAA \rightleftharpoons AA^- + H^+$$
$$\frac{[H^+][AA^-]}{[HAA]} = Ka $$

We also know that $[AA^-] + [HAA] = c_{AA}$. Similarly, for each alkaline amino acid (ARG, LYS and HIS), we have:

$$AAOH \rightleftharpoons AA^+ + OH^-$$
$$\frac{[OH^-][AA^+]}{[AAOH]} = Kb $$

In an aqueous solution, we know that $[H^+][OH^-] = 1 \times 10^{-14}$. Combined, we can get the following:

$$\begin{cases}
\frac{[H^+][AA^-]}{[HAA]} = Ka \\
\frac{[OH^-][AA^+]}{[AAOH]} = Kb
\end{cases}$$

---

We should take the pH-buffers into consideration. We start at the normal cellular $pH = 7.4$. 
