#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np

# TODO: Feature guard
#@ut.skipIf(not espressomd.has_features("..."),"Generic potentials not enabled")
class GenericTest(ut.TestCase):
    def setUp(self):
        self.system = espressomd.System(box_l=[10,10,10])
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.1

        self.pos0 = self.system.box_l/2 - [.1,0,0]
        self.pos1 = self.system.box_l/2 + [.1,0,0]
        self.params = {
            'epsilon': 0.1,
            'sigma': 0.1,
            'cutoff': 1.0,
            'shift': 0.0
        }

        self.system.part.add(pos=self.pos0,type=0)
        self.system.part.add(pos=self.pos1,type=0)

    def test_non_bonded(self):
        # Run with LJ
        self.system.non_bonded_inter[0,0].lennard_jones.set_params(**self.params)

        self.system.integrator.run(10)

        reference = np.copy(self.system.part[:].pos)

        # Reset particles
        self.system.part[0].pos = self.pos0
        self.system.part[1].pos = self.pos1
        self.system.part[:].v = [0,0,0]
        self.system.part[:].f = [0,0,0]

        # Run with LJ from expression
        self.system.non_bonded_inter[0,0].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=self.params["cutoff"], shift=0.0)
    
        self.system.non_bonded_inter[0,0].generic.set_params(
            cutoff=1.0,
            energy="4*{epsilon}*(({sigma}/r)**12 - ({sigma}/r)**6)".format(**self.params),
            force="-48*{epsilon}*(({sigma}/r)**12 - .5*({sigma}/r)**6)/r".format(**self.params))

        self.system.integrator.run(10)

        np.testing.assert_allclose(np.copy(self.system.part[:].pos), reference)

if __name__ == "__main__":
    ut.main()

