"""Generator-facing fixed challenge model. No changes to abaqus_ufl itself."""
import numpy as np
import abaqus_ufl as au
from abaqus_ufl.core.tensor import log


class ChemoFatigue(au.Material):
    props = dict(E=10.0, nu=0.28, beta=0.035, c0=0.35, kappa=1e-5,
                 D=0.14, chi=0.25, ell=0.12, G_c=3.5e-6,
                 alpha_T=1e-4, k_fg=0.5)
    state_vars = dict(hist=0.0, acc=0.0, wprev=0.0)

    def _make_state(self, **kwargs):
        # The stress-free concentration is c0, not the framework default 0.
        defaults = dict(c=self.c0, damage=0.1, intact=0.9, p=0.01,
                        grad_c=np.array([0.1, -0.2, 0.0]),
                        grad_intact=np.array([-0.2, -0.1, 0.0]),
                        grad_p=np.array([0.03, -0.02, 0.0]))
        defaults.update(kwargs)
        return super()._make_state(**defaults)

    def _psi_plus(self, F, c):
        mu_lame = self.E / (2.0 * (1.0 + self.nu))
        lam = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        bulk = lam + 2.0 * mu_lame / 3.0
        exx = F[0, 0] - 1.0 - self.beta * (c - self.c0)
        eyy = F[1, 1] - 1.0 - self.beta * (c - self.c0)
        exy = 0.5 * (F[0, 1] + F[1, 0])
        trpos = exx + eyy
        if trpos.real < 0.0:
            trpos = 0.0
        return 0.5 * bulk * trpos * trpos + mu_lame * (0.5 * (exx - eyy)**2 + 2.0 * exy**2)

    def _trial_history(self, F, c, intact, hist_old, acc_old, wprev_old):
        psi_plus = self._psi_plus(F, c)
        degr = (1.0 - self.kappa) * intact**2 + self.kappa
        wnow = degr * psi_plus
        anew = acc_old
        if wnow.real > wprev_old.real:
            anew = acc_old + (wnow - wprev_old)
        fatigue = 1.0
        if anew.real > self.alpha_T:
            inner = 1.0 - self.k_fg * log(anew / self.alpha_T) / log(10.0)
            if inner.real < 0.0:
                inner = 0.0
            fatigue = inner * inner
            if fatigue.real < 1e-15:
                fatigue = 1e-15
        hnew = 2.0 * self.ell * psi_plus / (self.G_c * fatigue)
        if hnew.real < hist_old.real:
            hnew = hist_old
        return np.array([hnew, anew, wnow])

    def stress_PK1(self, F, c, intact, hist_old, acc_old, wprev_old, dt):
        mu_lame = self.E / (2.0 * (1.0 + self.nu))
        lam = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        exx = F[0, 0] - 1.0 - self.beta * (c - self.c0)
        eyy = F[1, 1] - 1.0 - self.beta * (c - self.c0)
        exy = 0.5 * (F[0, 1] + F[1, 0])
        degr = (1.0 - self.kappa) * intact**2 + self.kappa
        sxx = degr * (2.0 * mu_lame * exx + lam * (exx + eyy))
        syy = degr * (2.0 * mu_lame * eyy + lam * (exx + eyy))
        sxy = degr * 2.0 * mu_lame * exy
        stress = np.array([[sxx, sxy, 0.0], [sxy, syy, 0.0], [0.0, 0.0, 0.0]])
        hv = self._trial_history(F, c, intact, hist_old, acc_old, wprev_old)
        return stress, {"hist": hv[0], "acc": hv[1], "wprev": hv[2]}

    def pressure_resid(self, F, c, intact, p):
        mu_lame = self.E / (2.0 * (1.0 + self.nu))
        lam = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        etrace = F[0, 0] + F[1, 1] - 2.0 - 2.0 * self.beta * (c - self.c0)
        degr = (1.0 - self.kappa) * intact**2 + self.kappa
        return p - degr * 2.0 * (mu_lame + lam) * etrace / 3.0

    def solvent_storage(self, c, c_old, dt):
        return (c - c_old) / dt

    def solvent_flux(self, c, intact, grad_c, grad_p):
        degr = (1.0 - self.kappa) * intact**2 + self.kappa
        return -degr * self.D * grad_c + self.chi * c * grad_p

    def phase_storage(self, F, c, intact, hist_old, acc_old, wprev_old):
        hv = self._trial_history(F, c, intact, hist_old, acc_old, wprev_old)
        # Negative of the original damage equation, with s=1-d evaluated
        # independently to preserve tiny intactness at exhausted fatigue.
        return (1.0 + (1.0 - self.kappa) * hv[0]) * intact - 1.0

    def phase_flux(self, grad_intact):
        return -self.ell**2 * grad_intact

    def species_storage(self, damage, intact):
        return damage + intact - 1.0

    def species_flux(self, grad_intact):
        return 0.0 * grad_intact


class ChemoFatigueQuad(au.WeakForm):
    material = ChemoFatigue
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)
        self.c = au.ScalarField("c", degree=1, test="z")
        self.damage = au.ScalarField("damage", degree=1, test="eta")
        self.p = au.ScalarField("p", degree=1, test="q")
        self.intact = au.ScalarField("intact", degree=1, test="zint")

    def momentum_equation(self, v, F, c, intact, hist_old, acc_old, wprev_old, dt):
        return self.material.stress_PK1(F, c, intact, hist_old, acc_old, wprev_old, dt)

    def transport_equation(self, z, c, intact, grad_c, grad_p, c_old, dt):
        return (self.material.solvent_storage(c, c_old, dt),
                self.material.solvent_flux(c, intact, grad_c, grad_p))

    def phase_equation(self, eta, F, c, intact, grad_intact, hist_old, acc_old, wprev_old):
        return (self.material.phase_storage(F, c, intact, hist_old, acc_old, wprev_old),
                self.material.phase_flux(grad_intact))

    def pressure_equation(self, q, F, c, intact, p):
        return self.material.pressure_resid(F, c, intact, p)

    def species_transport_equation(self, zint, damage, intact, grad_intact):
        return (self.material.species_storage(damage, intact),
                self.material.species_flux(grad_intact))


class LumpedCoreMaterial(ChemoFatigue):
    """Generated bulk and history; the nodal reaction is in accepted_damage.f90."""

    def phase_storage(self):
        return 0.0


class LumpedCoreQuad(ChemoFatigueQuad):
    material=LumpedCoreMaterial

    def phase_equation(self,eta,grad_intact):
        return self.material.phase_storage(),self.material.phase_flux(grad_intact)


def verification_states(problem):
    """Smooth states inside each active regime; kink tests are separate."""
    F = np.array([[1.004, 0.002, 0.0], [-0.001, 1.001, 0.0], [0.0, 0.0, 1.0]])
    base = problem._mat._make_state(F=F, c=0.39, damage=0.3, intact=0.7, c_old=0.37,
                                    hist_old=0.0, acc_old=0.0, wprev_old=0.0,
                                    dt=0.1)
    unloading = dict(base, F=np.eye(3), hist_old=30.0,
                     acc_old=5e-4, wprev_old=0.001)
    fatigue = dict(base, acc_old=2e-4)
    exhausted = dict(base, acc_old=0.02)
    return [base, unloading, fatigue, exhausted]
