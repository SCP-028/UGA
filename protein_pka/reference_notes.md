# Background from the Literature
## Intracellular H+ buffering power and its dependency on intracellular pH, Saleh, 1991
There are three types of intracellular $H^+$ buffering processes: `physicochemical buffering`, `biochemical reactions`, and transport of acid-base equivalents across organellar membranes (that is, `organellar buffering`).

Cell $H^+$ buffering power increases with decreasing intracellular pH.

Relationship between intracellular proton buffering capacity and intracellular pH mentioned a caveat in this paper.


## pH sensing and regulation in cancer
To overcome low $pH_i$ due to elevated rates of glycolysis, cancer cells employ a large redundancy of mechanisms to remove acids in order to maintain physiological $pH_i$. In response to glycolytic acidosis, $pH_i$ can be be maintained via lactate and $H^+$ efflux by monocarboxylate transporters and $Na^+$-driven proton extrusion, respectively (Gillies, 2002; Gallagher et al., 2008). As a consequence, the pH of the extracellular space of tumors becomes acidic; forming a reversed pH gradient ($pH_e$ < $pH_i$) in comparison to normal physiological conditions ($pH_e$ > $pH_i$).

The relatively acidic $pH_e$ induces migration and invasion (Bradley et al., 2011; Hanahan and Weinberg, 2011). Increased $pH_i$ by itself also has effects on cancer cell function such as increased proliferation (Moolenaar et al., 1986), promoting cell survival by limiting apoptosis, which is associated with intracellular acidification (Matsuyama and Reed, 2000), and selective advantages of growth-factor independent proliferation.

# Discuss with Sha
1. Transformation from mRNA to protein abundancy - [陈炜](http://www.sustc.edu.cn/biology_04/f/chenwei)
   - Protein array data?

$$P_i = f(E_i) \approx aE_i^2 + bE_i + c$$

where $E_i$ is the expression level (TPM) of gene $i$, and $P_i$ is the concentration of the protein product of gene $i$.

2. Assuming the volume of the cell is $100 \mu m^3$ and the pH buffer coefficient is $2.0 \times 10^5$, the number of protons needed to change the intracellular pH ($pH_i$) from $pH_0$ to $pH_1$ is:

$$(10^{-pH_1} - 10^{-pH_0}) \times 100 \times 2 \times 10^5 \times 10^{-15} \times 6.02 \times 10^{23} \cong 1.1 \times 10^9$$

where $1L = 10^{15}\mu m^3$ and $6.02 \times 10^{23}$ is the Avogadro constant. 

3. The new "pH level" $pH_1$ could be represented as:

$$
\begin{aligned}
(10^{-pH_1} - 10^{-pH_0}) \cdot C &= n\left(H^+\right) \\
10^{-pH_1} &= \frac{n(H^+)}{C} + 10^{-pH_0}\\
pH_1 &= -\log_{10} \left(10^{-pH_0} + \frac{1}{C}\left(\sum_{basic}\frac{P_i}{1 + 10^{pH_0 - pKa}} - \sum_{acidic}\frac{P_i}{1 + 10^{pKa - pH_0}}\right)\right)
\end{aligned}
$$

where $C = 100 \times 2 \times 10^5 \times 10^{-15} \times 6.02 \times 10^{23} \approx 1.2 \times 10^{16}$.

4. The problem now is that this $pH_1$ is not the actual $pH_i$. It is merely considering the proteins as a pH buffer system and calculating their effects on the intracellular pH. However, protein pKa often acts as a pH sensor, and regulate the intracellular pH through various proton pumps and transporters.

$$pH_i = \alpha \cdot pH_1 + \beta \Longrightarrow \text{neural network} \Rightarrow \text{predicted pH}$$
