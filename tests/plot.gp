#!/usr/bin/gnuplot --persist

p 'QChem_VG_small_basis/ref.specden.txt' u 1:3 w l lw 2, \
  'G09_VG_small_basis/ref.specden.txt' u 1:3 w l lw 2, \
  'G09_VG_big_basis/ref.specden.txt' u 1:3 w l lw 2, \
  'G09_AS_big_basis/ref.specden.txt' u 1:3 w l lw 2
