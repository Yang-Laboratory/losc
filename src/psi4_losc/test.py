import psi4_losc
import psi4

# Memory for Psi4 in GB
psi4.set_memory('2 GB')
psi4.core.set_output_file('output.dat', False)

# Memory for numpy in GB
numpy_memory = 2

# Triplet O2
mol = psi4.geometry("""
    0 3
    O
    O 1 1.2
symmetry c1
""")

psi4.set_options({'guess':      'core',
                  'basis':      '3-21g',
                  'scf_type':   'df',
                  'e_convergence': 1e-8,
                  'reference':  'uhf',
                  'maxiter': 50})
dfa_name = "hf"
dfa_name = "blyp"
dfa_name = "b3lyp"


SCF_E_psi, dfa_wfn = psi4.energy(dfa_name, return_wfn=True)
# Compare to Psi4
SCF_E = psi4_losc.scf_losc(dfa_wfn).energy()
psi4.compare_values(SCF_E_psi, SCF_E, 6, 'SCF Energy')
