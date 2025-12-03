import numpy as np
from scipy.linalg import eigh


def get_hcore(mol):
    return mol.intor("int1e_kin_sph") + mol.intor("int1e_nuc_sph")


def get_ovlp(mol):
    return mol.intor("int1e_ovlp_sph")


def get_eri(mol):
    return mol.intor("int2e_sph")


def get_core_guess(mol, hcore=None, ovlp=None):
    if ovlp is None:
        ovlp = get_ovlp(mol)
    if hcore is None:
        hcore = get_hcore(mol)
    mo_energy, mo_coeff = np.linalg.eigh(hcore)
    mo_occ = get_occ(mol, mo_energy)
    return make_rdm1(mo_coeff, mo_occ)


def make_rdm1(mo_coeff, mo_occ):
    mocc = mo_coeff[:, mo_occ > 0]
    # dm = (mocc * mo_occ[mo_occ > 0]).dot(mocc.T.conj())
    dm = np.einsum(
        "pi,i,iq->pq", mocc, mo_occ[mo_occ > 0], mocc.T.conj(), optimize=True
    )
    return dm


def get_fock(hcore, eri, dm):
    veff = get_veff(eri, dm)
    return hcore + veff


def get_veff(eri, dm):
    vj, vk = get_jk(eri, dm)
    return vj - 0.5 * vk


def get_jk(eri, dm):
    J = np.einsum("pqrs,rs->pq", eri, dm, optimize=True)
    K = np.einsum("prqs,rs->pq", eri, dm, optimize=True)
    return J, K


def get_occ(mol, mo_energy):
    e_idx = np.argsort(mo_energy)
    e_sort = mo_energy[e_idx]
    nmo = len(mo_energy)
    mo_occ = np.zeros(nmo)
    nocc = mol.nelectron // 2
    mo_occ[e_idx[:nocc]] = 2
    print(f" HOMO {e_sort[nocc - 1]:.15g}, LUMO {e_sort[nocc]:.15g}")
    return mo_occ


def energy_tot(mol, hcore, fock, dm):
    e_elec = energy_elec(hcore, fock, dm)
    e_nuc = mol.energy_nuc()
    return e_elec + e_nuc


def energy_elec(hcore, fock, dm):
    e_elec = 0.5 * np.einsum("pq,pq->", hcore + fock, dm, optimize=True)
    return e_elec


def scf(mol, max_iter=100, tol=1e-6):
    H = get_hcore(mol)
    S = get_ovlp(mol)
    I = get_eri(mol)
    D = get_core_guess(mol, H, S)

    hf_energy_old = 0.0
    for iter in range(max_iter):
        F = get_fock(H, I, D)
        hf_energy = energy_tot(mol, H, F, D)
        mo_energy, mo_coeff = eigh(F, S)
        mo_occ = get_occ(mol, mo_energy)
        D_new = make_rdm1(mo_coeff, mo_occ)

        E_diff = np.abs(hf_energy - hf_energy_old)
        D_diff = np.linalg.norm(D_new - D)

        print(
            f"Iter {iter:3d}: E = {hf_energy:.10f}, "
            f"dE = {E_diff:.3e}, dD = {D_diff:.3e}"
        )

        if E_diff < tol and D_diff < tol:
            print(f"SCF converged after {iter} iterations")
            print(f"E = {hf_energy:.10f}")
            return hf_energy

        hf_energy_old = hf_energy
        D = D_new

    raise RuntimeWarning(f"SCF did not converge after {iter} iterations")


if __name__ == "__main__":
    from pyscf import gto

    mol = gto.M(
        atom="""
        O 0.0 0.0 0.0
        H 0.0 0.0 1.0
        H 0.0 1.0 0.0
        """,
        basis="cc-pvtz",
        spin=0,
        charge=0,
        symmetry=False,
    )
    my_hf = scf(mol)
    # from pyscf import scf

    # pyscf_hf = scf.rhf.RHF(mol).kernel()
    # assert np.allclose(my_hf, pyscf_hf)
