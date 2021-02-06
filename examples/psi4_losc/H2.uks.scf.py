from psi4_losc import psi4_losc
import psi4

# Memory for Psi4 in GB
psi4.set_memory('2 GB')
psi4.core.set_output_file('H2.uks.scf.out', False)

# Memory for numpy in GB
numpy_memory = 2

mol = psi4.geometry("""
    0 1
    H
    H 1 1.2
symmetry c1
""")

psi4.set_options({'guess':      'core',
                  'basis':      '3-21g',
                  'scf_type':   'df',
                  'e_convergence': 1e-8,
                  'reference':  'uhf',
                  'maxiter': 50})
dfa_name = "b3lyp"

# psi4 SCF-B3LYP calculation.
SCF_E_psi, dfa_wfn = psi4.energy(dfa_name, return_wfn=True)

# psi4_losc SCF-LOSC-B3LYP calculation.
wfn = psi4_losc.scf_losc(dfa_wfn)

# Compare to Psi4. LOSC should have no corrections to this small molecule
psi4.compare_values(SCF_E_psi, wfn.energy(), 6, 'SCF Energy')
