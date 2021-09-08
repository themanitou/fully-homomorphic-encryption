##################################################################
#
#   Fully-Homomorphic Cryptography library,
#   based on Gentry-Halevi ideal lattice scheme.
#
#   Author: Quan Nguyen (https://github.com/themanitou)
#
#   This library is open-source software distributed under the
#   terms of the GNU Lesser General Public License (LGPL) version
#   2.1 or later.  See the file doc/copying.txt for complete
#   details on the licensing of this library.
#
###################################################################


#!/bin/sh
NP=`grep -c "^processor" /proc/cpuinfo`-1 && mpirun --np ${NP} --use-hwthread-cpus ../bin/fullhomcrypt
