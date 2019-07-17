set term pngcairo
set size square
set pm3d map
load 'moreland.pal'
#load 'rdbu.pal'
#set palette model RGB
#set palette defined
set key off

unset xtics
unset ytics

unset ylabel

splot filename using ($1)