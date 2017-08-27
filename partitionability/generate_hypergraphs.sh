#!/bin/bash

# Copyright (c) 2014-2017 Mitchel T. Keller and Stephen J. Young.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#PBS -lnodes=1:ppn=1,mem=200mb,walltime=24:00:00
/home/kellermt/nauty25r9/genbgL -z -d0:$k -D$e:$k $n $e $res/$mod /home/kellermt/partitionability/$n-hypergraphs/$n-$k-$e-hypergraphs-$res.g6
#gzip -9 /home/kellermt/sdepth/$n-hypergraphs/$n-$k-$e-hypergraphs-$res.g6
