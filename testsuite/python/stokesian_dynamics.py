import espressomd
from espressomd import constraints
import numpy as np
import unittest as ut
import unittest_decorators as utx
from tests_common import abspath
import scipy.optimize

s = espressomd.System(box_l=[1.0, 1.0, 1.0])


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDynamicsSetupTest(ut.TestCase):
    system = s

    def setUp(self):
        self.system.box_l = [10] * 3

        self.system.time_step = 1.0
        self.system.cell_system.skin = 0.4
        

        from espressomd.thermostat import flags
        self.system.thermostat.set_sd(viscosity=1.0,
                                      device="cpu",
                                      radii={0: 1.0},
                                      flags=flags.SELF_MOBILITY | flags.PAIR_MOBILITY | flags.FTS)

        # unset SD integrator so we can test whether set_sd fails
        # set_nvt() is the only way to ensure that integ_switch is 
        # set to a different value than INTEG_METHOD_SD
        self.system.integrator.set_nvt()

    def test_pbc_checks(self):

        self.system.periodicity = [0, 0, 1]
        with (self.assertRaises(Exception)): 
            self.system.integrator.set_sd()

        
        self.system.periodicity = [0, 0, 0]
        self.system.integrator.set_sd()
        with (self.assertRaises(Exception)):
            self.system.periodicity = [0, 1, 0]


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDynamicsTest(ut.TestCase):
    system = s

    # Digitized reference data of Figure 5b from
    #
    # Durlofsky et al., J. Fluid Mech. 180, 21 (1987)
    # https://doi.org/10.1017/S002211208700171X
    data = np.loadtxt(abspath('data/dancing.txt'))

    def setUp(self):
        self.system.box_l = [10] * 3
        self.system.periodicity = [0, 0, 0]

        self.system.time_step = 1.0
        self.system.cell_system.skin = 0.4

        self.system.part.add(pos=[-5, 0, 0], rotation=[1, 1, 1])
        self.system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])
        self.system.part.add(pos=[7, 0, 0], rotation=[1, 1, 1])

        from espressomd.thermostat import flags
        self.system.thermostat.set_sd(viscosity=1.0,
                                      device="cpu",
                                      radii={0: 1.0},
                                      flags=flags.SELF_MOBILITY | flags.PAIR_MOBILITY | flags.FTS)
        self.system.integrator.set_sd()

        gravity = constraints.Gravity(g=[0, -1, 0])
        self.system.constraints.add(gravity)

    def test(self):
        intsteps = int(8000 / self.system.time_step)
        pos = np.empty([intsteps + 1, 3 * len(self.system.part)])

        pos[0, :] = self.system.part[:].pos.flatten()
        for i in range(intsteps):
            self.system.integrator.run(1)
            for n, p in enumerate(self.system.part):
                pos[i + 1, 3 * n:3 * n + 3] = p.pos

        for i in range(self.data.shape[0]):
            for n in range(3):
                x_ref = self.data[i, 2 * n]
                y_ref = self.data[i, 2 * n + 1]
                x = pos[:, 3 * n]
                y = pos[:, 3 * n + 1]

                if y_ref < -555:
                    continue

                idx = np.abs(y - y_ref).argmin()
                dist = np.sqrt((x_ref - x[idx])**2 + (y_ref - y[idx])**2)
                self.assertLess(dist, 0.5)

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDiffusionTest(ut.TestCase):
    system = s

    kT = 1e-4
    R = 1.5
    eta = 1.0

    def setUp(self):
        #self.system.actors.clear()
        #self.system.part.clear()
        self.system.box_l = [10] * 3
        self.system.periodicity = [0, 0, 0]

        self.system.time_step = 1.0
        self.system.cell_system.skin = 0.4

        self.system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])

        from espressomd.thermostat import flags
        self.system.thermostat.set_sd(viscosity=self.eta,
                                      device="cpu",
                                      radii={0: self.R},
                                      kT=self.kT,
                                      seed=0,
                                      flags=flags.SELF_MOBILITY | flags.PAIR_MOBILITY | flags.FTS)
        self.system.integrator.set_sd()

    def test(self):
        intsteps = int(100000 / self.system.time_step)
        pos = np.empty([intsteps + 1, 3])
        orientation = np.empty((intsteps + 1, 3))

        pos[0, :] = self.system.part[0].pos
        orientation[0, :] = self.system.part[0].director
        for i in range(intsteps):
            self.system.integrator.run(1)
            pos[i + 1, :] = self.system.part[0].pos
            orientation[i + 1, :] = self.system.part[0].director

        t = np.arange(0, intsteps + 1)
        msd = np.linalg.norm(pos - pos[0, :], axis=1)**2
        costheta = np.dot(orientation[:, :], orientation[0, :])


        #translational diffusion coefficient
        D_expected = self.kT / (6 * np.pi * self.eta * self.R)

        # NOTE on steps per slice
        # The shorter these trajectories are, the more trajectories we get
        # and the more accurate the result is. 
        # However since we want to test that diffusion works for 
        # "long" trajectories we can't go too low.
        n_steps_per_slice = 200
        n_slices = int(intsteps / n_steps_per_slice)

        squared_displacement_per_slice = np.empty(n_slices)
        for i in range(n_slices):
            squared_displacement_per_slice[i] = np.linalg.norm(pos[i * n_steps_per_slice] - pos[(i+1) * n_steps_per_slice])**2
        D_measured = np.mean(squared_displacement_per_slice) / (6 * n_steps_per_slice * self.system.time_step)
        self.assertAlmostEqual(D_expected, D_measured, delta=D_expected * 0.2)


        #rotational diffusion coefficient
        Dr_expected = self.kT / (8 * np.pi * self.eta * self.R**3)
        tr_expected = int(1 / (2 * Dr_expected))
        fit = scipy.optimize.curve_fit(lambda t, Dr: np.exp(-2 * Dr * t),
                                       t[:2 * tr_expected],
                                       costheta[:2 * tr_expected],
                                       p0 = 1e-5)
        Dr_measured = fit[0][0]
        self.assertAlmostEqual(
            Dr_expected,
            Dr_measured,
            delta=Dr_expected * 1)

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()




if __name__ == '__main__':
    ut.main()
