# Normalize by Initial Perturbation

to calculate delta(IRF): "normalize by initial perturbation":

After adding a pulse and re-running FaIR, compute the concentration change as:

$$
\Delta C(t) = C^{\text{pert}}(t) - C^{\text{base}}(t)
$$

Then define the impulse response function (IRF) as:

$$
IRF(H) = \frac{\Delta C(t)}{\Delta C(t_0)} =  \frac{\Delta C(t_0 + H)}{\Delta C(t_0)}
$$


where: 
- t / H denotes the time after the emission pulse at which the concentration response is evaluated, 
- t₀ is the reference time immediately after the pulse is realized (i.e., the first post-pulse concentration point used to normalize the response so that IRF(t₀)=1).

------------------------------------------------------------------------

## What This Does Mathematically

Suppose at pulse year $t_0$:

-   we inject **+1 unit emission**
-   FaIR gives a concentration jump:

$$
\Delta C(t_0) = 0.45 \text{ ppm}
$$

At later years:

$$
\Delta C(t_0 + 10) = 0.30
$$

$$
\Delta C(t_0 + 50) = 0.18
$$

Then:

$$
IRF(0) = \frac{0.45}{0.45} = 1
$$

$$
IRF(10) = \frac{0.30}{0.45} = 0.667
$$

$$
IRF(50) = \frac{0.18}{0.45} = 0.4
$$

------------------------------------------------------------------------

## Interpretation

The IRF becomes a **dimensionless decay curve** that:

-   Starts at **1** at $H = 0$
-   Declines over time
-   Represents the fraction of the initial perturbation remaining in the
    atmosphere

This normalization removes dependence on the absolute concentration jump and isolates the atmospheric decay behavior.
