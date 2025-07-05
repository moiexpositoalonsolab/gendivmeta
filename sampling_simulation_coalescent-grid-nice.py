import msprime
import numpy as np
import pandas as pd
from itertools import product
import math

# --------------------------------------------------------------------------------
# Nice 1–2–5 ladder
def nice_log_series(low, high):
    k_min = int(math.floor(math.log10(low)))
    k_max = int(math.ceil (math.log10(high)))
    series = [m * 10**k
              for k in range(k_min, k_max + 1)
              for m in (1, 2, 5)
              if low <= m * 10**k <= high]
    return series

# --------------------------------------------------------------------------------
# Simulation parameters
recomb_rate = 1e-8
mu         = 1e-5          # high μ keeps the run-times modest
Ne_grid    = [100, 1_000, 10_000, 100_000]
n_grid     = nice_log_series(1, 10_000)       # *single* grid for n0 & n1
L_grid     = nice_log_series(1, 1_000_000)
reps       = 10

results = []

# Cartesian product over all factors except replicate
for rep, Ne, n, L in product(range(reps), Ne_grid, n_grid, L_grid):

    # Two contemporaneous sample sets of size n each
    samples = [
        msprime.SampleSet(n, time=0, population=0, ploidy=1),
        msprime.SampleSet(n, time=0, population=0, ploidy=1)
    ]

    ts  = msprime.sim_ancestry(
            samples=samples,
            sequence_length=L,
            population_size=Ne,
            recombination_rate=recomb_rate)

    mts = msprime.sim_mutations(ts, rate=mu)

    # First n haplotypes → “early”; next n → “recent”
    early_idx  = range(0, n)
    recent_idx = range(n, 2 * n)

    pi0 = mts.diversity(early_idx)
    pi1 = mts.diversity(recent_idx)
    rel_change = 100 * (pi1 - pi0) / pi0   # percentage

    results.append({
        "rep": rep,
        "Ne": Ne,
        "n0": n,
        "n1": n,
        "L": L,
        "pi0": pi0,
        "pi1": pi1,
        "relative_change": rel_change,
        "num_sites": mts.num_sites
    })

df = pd.DataFrame(results)
df.to_csv("pi_grid_simulation_results.csv", index=False)
print(df.head())