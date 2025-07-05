import msprime
import numpy as np
import pandas as pd
import random

# Parameters
recomb = 1e-8
mu = 1e-5
Ne_grid = [100, 1000, 10000, 100000]
reps = 500

# Load empirical data
# df_meta = pd.read_csv("tmp/dfpicisubset.csv")
df_meta = pd.read_csv("tmp/dfpici.csv")

# Drop rows with missing values in any required column
df_meta = df_meta.dropna(subset=["Early.Sample.N", "Recent.Sample.N", "N.Loci"])

# Convert to integers
df_meta["Early.Sample.N"] = df_meta["Early.Sample.N"].astype(int)
df_meta["Recent.Sample.N"] = df_meta["Recent.Sample.N"].astype(int)
df_meta["N.Loci"] = df_meta["N.Loci"].astype(int)

results = []

for r in range(reps):
    for Ne in Ne_grid:
        for i, row in df_meta.iterrows():
            n0 = row["Early.Sample.N"]
            n1 = row["Recent.Sample.N"]
            L = row["N.Loci"]

            if n0 < 2 or n1 < 2 or L < 1:
                continue  # skip invalid values

            samples = [
                msprime.SampleSet(n0, time=0, population=0, ploidy=1),
                msprime.SampleSet(n1, time=0, population=0, ploidy=1)
            ]

            ts = msprime.sim_ancestry(samples=samples,
                                      sequence_length=L,
                                      population_size=Ne,
                                      recombination_rate=recomb)

            mts = msprime.sim_mutations(ts, rate=mu)

            sample_indices_time0 = list(range(0, n0))
            sample_indices_time1 = list(range(n0, n0 + n1))

            pi0 = mts.diversity(sample_indices_time0)
            pi1 = mts.diversity(sample_indices_time1)
            relative_change = (pi1 - pi0) / pi0 if pi0 > 0 else float("nan")

            results.append({
                "rep": r,
                "Ne": Ne,
                "n0": n0,
                "n1": n1,
                "L": L,
                "pi0": pi0,
                "pi1": pi1,
                "relative_change": relative_change * 100,
                "num_sites": mts.num_sites
            })

df = pd.DataFrame(results)
df.to_csv("pi_grid_simulation_results_realistic_loci.csv", index=False)
print(df)
