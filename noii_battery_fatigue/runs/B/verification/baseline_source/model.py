"""Q1 mixed chemo-mechanics/fatigue problem for the challenge's exact 2D laws.

The extra field q is an L2 projection of the degraded hydrostatic stress.
All history arguments are the committed start-of-increment values. Branches
are selected on real parts so that complex-step gives a semismooth tangent.
"""
import numpy as np
import abaqus_ufl as au
from abaqus_ufl.core.tensor import eye, log


class ChemoFatigue(au.Material):
    props = dict(E=10.0, nu=0.28, beta=0.035, c0=0.35, kappa=1e-5,
                 D=0.14, chi=0.25, ell=0.12, G_c=3.5e-6,
                 alpha_T=1e-4, k_fg=0.5)
    state_vars = dict(hhist=0.0, afat=0.0, wprev=0.0)

    def _make_state(self, **kwargs):
        # The generic verification reference assumes every scalar is zero.
        # The physical stress-free chemical reference here is c=c0.
        state = super()._make_state()
        state["c"] = self.c0
        state.update(kwargs)
        return state

    def _g(self, d):
        return (1.0 - self.kappa) * (1.0 - d) ** 2 + self.kappa

    def _psi(self, F, c):
        mu_e = self.E / (2.0 * (1.0 + self.nu))
        lam_e = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        bulk = lam_e + 2.0 * mu_e / 3.0
        exx = F[0, 0] - 1.0 - self.beta * (c - self.c0)
        eyy = F[1, 1] - 1.0 - self.beta * (c - self.c0)
        exy = 0.5 * (F[0, 1] + F[1, 0])
        tr_e = exx + eyy
        tr_plus = 0.0 * tr_e
        if tr_e.real > 0.0:
            tr_plus = tr_e
        return 0.5 * bulk * tr_plus ** 2 + mu_e * (0.5 * (exx - eyy) ** 2 + 2.0 * exy ** 2)

    def _alpha(self, F, c, d, afat_old, wprev_old):
        w_new = self._g(d) * self._psi(F, c)
        afat_new = afat_old
        if w_new.real > wprev_old.real:
            afat_new = afat_old + w_new - wprev_old
        return afat_new

    def _fatigue(self, alpha):
        fatigue = 1.0 + 0.0 * alpha
        if alpha.real > self.alpha_T:
            z = 1.0 - self.k_fg * log(alpha / self.alpha_T) / log(10.0)
            if z.real < 0.0:
                z = 0.0 * z
            fatigue = z ** 2
            if fatigue.real < 1e-15:
                fatigue = 1e-15 + 0.0 * alpha
        return fatigue

    def _H(self, F, c, d, hhist_old, afat_old, wprev_old):
        alpha = self._alpha(F, c, d, afat_old, wprev_old)
        fatigue = self._fatigue(alpha)
        candidate = 2.0 * self.ell * self._psi(F, c) / (self.G_c * fatigue)
        hhist_new = hhist_old
        if candidate.real > hhist_old.real:
            hhist_new = candidate
        return hhist_new

    def stress_PK1(self, F, c, d, hhist_old, afat_old, wprev_old, dt):
        # F=I+grad(u) is used only to carry the small displacement gradient.
        mu_e = self.E / (2.0 * (1.0 + self.nu))
        lam_e = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        exx = F[0, 0] - 1.0 - self.beta * (c - self.c0)
        eyy = F[1, 1] - 1.0 - self.beta * (c - self.c0)
        tr_e = exx + eyy
        degradation = self._g(d)
        # The generator carries 3x3 tensors; only the upper 2x2 block enters
        # 2D assembly. Its unused zz entry is not a model stress component.
        ident = eye(3)
        strain_carrier = 0.5 * (F + F.T) - ident - self.beta * (c - self.c0) * ident
        stress = degradation * (2.0 * mu_e * strain_carrier + lam_e * tr_e * ident)
        h_new = self._H(F, c, d, hhist_old, afat_old, wprev_old)
        a_new = self._alpha(F, c, d, afat_old, wprev_old)
        w_new = degradation * self._psi(F, c)
        return stress, {"hhist": h_new, "afat": a_new, "wprev": w_new}

    def pressure_resid(self, F, c, d, q):
        mu_e = self.E / (2.0 * (1.0 + self.nu))
        lam_e = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        tr_e = F[0, 0] + F[1, 1] - 2.0 - 2.0 * self.beta * (c - self.c0)
        return q - self._g(d) * 2.0 * (mu_e + lam_e) * tr_e / 3.0

    def species_storage(self, c, c_old, dt):
        return (c - c_old) / dt

    def species_flux(self, c, d, grad_c, grad_q):
        return -self._g(d) * self.D * grad_c + self.chi * c * grad_q

    def phase_storage(self, F, c, d, hhist_old, afat_old, wprev_old):
        h_trial = self._H(F, c, d, hhist_old, afat_old, wprev_old)
        return d - (1.0 - self.kappa) * (1.0 - d) * h_trial

    def phase_flux(self, grad_d):
        return -self.ell ** 2 * grad_d


class ChemoProblem(au.WeakForm):
    material = ChemoFatigue
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)
        self.c = au.ScalarField("c", degree=1, test="zeta")
        self.d = au.ScalarField("d", degree=1, test="eta_d")
        self.q = au.ScalarField("q", degree=1, test="eta_q")

    def momentum_equation(self, v, F, c, d, hhist_old, afat_old, wprev_old, dt):
        return self.material.stress_PK1(F, c, d, hhist_old, afat_old, wprev_old, dt)

    def pressure_equation(self, eta_q, F, c, d, q):
        return self.material.pressure_resid(F, c, d, q)

    def species_transport_equation(self, zeta, c, c_old, grad_c, d, grad_q, dt):
        return (self.material.species_storage(c, c_old, dt),
                self.material.species_flux(c, d, grad_c, grad_q))

    def phase_equation(self, eta_d, F, c, d, grad_d, hhist_old, afat_old, wprev_old):
        return (self.material.phase_storage(F, c, d, hhist_old, afat_old, wprev_old),
                self.material.phase_flux(grad_d))


def verification_states():
    """Nonzero trace, shear, concentration gradients, and both history branches."""
    base = dict(F=np.array([[1.021, .013, 0.], [.004, 1.006, 0.], [0., 0., 1.]]),
                c=.53, d=.27, q=.015, grad_c=np.array([.3, -.2, 0.]),
                grad_d=np.array([-.1, .2, 0.]), grad_q=np.array([.02, -.03, 0.]),
                c_old=.49, dt=.15, hhist_old=0., afat_old=0., wprev_old=0.)
    yield "active_undiminished", base
    yield "active_fatigue", {**base, "afat_old": .002, "wprev_old": 1e-7}
    yield "unloading_frozen", {**base, "hhist_old": 1e5, "wprev_old": .02, "afat_old": .002}
    yield "fatigue_floor", {**base, "afat_old": 2., "wprev_old": 1e-7}
    yield "compression", {**base, "F": np.array([[.98, .01, 0.], [.005, .99, 0.], [0., 0., 1.]])}

