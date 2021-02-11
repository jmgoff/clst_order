# Paralell example for Pt/Pd surface alloy system

At one Pd surface fraction (0.08), 100 parallel Monte Carlo simulations were ran
as described in PRB ... The probability of a Pd-Pd-Pd-Pd "channel" cluster
occurring in the top 4-layers is quantified as an esemble-averaged quantity from
all simulations. From these, the cluster order parameter can be calculated as (1
- P(Pd-Pd-Pd-Pd)/P<sub>random</sub>) where P<sub>random</sub> = c<sub>Pd</sub>
c<sub>Pd</sub> c<sub>Pd</sub> c<sub>Pd</sub>. See Goff et al. 2021 for more
details. The measured cluster probability, P(Pd-Pd-Pd-Pd), is calculated in
<pre><code>run_sro.py</code></pre>.

## Running the example:

Run the script that does the parallelized execution of clst_prob. The trajectory
files of the simulations, .dc files, and the primitive cell structure,
'atoms.cif', are needed for execution. The resulting probabilities for the 100
parallel simulations are stored in .json files with an 'example' prefix.

<pre><code>python run_sro.py</code></pre>

This may take a while depending on your parallel environment. To speed it up,
run it for the first 10 trajectories instead of 100. After the cluster
probabilities have been calculated, the cluster order parameters can be
measured. 

<pre><code>python parallel_analysis.py example</code></pre>

where the name of the prefix for the .json probability files are given as the
first argument.
