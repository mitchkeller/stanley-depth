# stanley-depth
Code and data related to the Stanley depth of squarefree monomial
ideals as discussed in ["Combinatorial Reductions for the Stanley Depth
of I and S/I" by Mitchel T. Keller and Stephen J. Young](https://arxiv.org/abs/1702.00781).

Files in the `combcrit-split` directory relate to checking the
combinatorial criterion and splitting properties as discussed in
Section 7 of the paper. Files in the `partitionability` directory
relate to the computational proof of Theorem 14.

Data files are stored in Brendan McKay's .g6 format and generated
using [nauty](http://users.cecs.anu.edu.au/~bdm/nauty/). A version of
`genbg` which supports larger graphs was required to complete the
computational generation for n = 7.

The source code in this repository is free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
