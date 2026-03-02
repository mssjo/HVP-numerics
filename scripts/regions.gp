# This GNUplot script generates the curves subdividing the monodromy regions.
# join_contours.py must be called to post-process the result for PGFplotting.

unset surface
set contour base
set cntrparam levels discrete 0.0

set xrange [-3:12]; set yrange [-3:3]

# The curve along which Im(the argument of the 4th root is 0 (plus at y=0)
im0(x,y) = -y**2 + 3*x**2 - 18*x + 3

# The curves along which Im(jsun) is 0 (plus at y=0)
im01(x,y) = 9*x - 6*x**2 + x**3 - 4*y**2 + x*y**2
im02(x,y) = -9 + 30*x**2 - 24*x**3 + 3*x**4 + 30*y**2 - 24*x*y**2 + 2*x**2*y**2 - y**4
im03(x,y) = 81 - 486*x + 1485*x**2 - 3456*x**3 + 5058*x**4 - 5868*x**5 + 4050*x**6 - 1680*x**7 + 333*x**8 - 30*x**9 + x**10 + 513*y**2 - 2808*x*y**2 + 4572*x**2*y**2 - 6768*x**3*y**2 + 5286*x**4*y**2 - 2616*x**5*y**2 + 444*x**6*y**2 - 32*x**7*y**2 + x**8*y**2 - 486*y**4 - 900*x*y**4 + 582*x**2*y**4 - 480*x**3*y**4 - 362*x**4*y**4 + 84*x**5*y**4 - 6*x**6*y**4 - 654*y**6 + 456*x*y**6 - 724*x**2*y**6 + 144*x**3*y**6 - 14*x**4*y**6 - 251*y**8 + 58*x*y**8 - 11*x**2*y**8 - 3*y**10

# The curves along which Re(jsun) is 1 (although re11 is strictly positive so it doesn't contribute)
re11(x,y) = 81 - 162*x + 135*x**2 - 60*x**3 + 87*x**4 - 18*x**5 + x**6 + 27*y**2 - 36*x*y**2 + 162*x**2*y**2 - 36*x**3*y**2 + 3*x**4*y**2 + 75*y**4 - 18*x*y**4 + 3*x**2*y**4 + y**6
re12(x,y) = 6561 - 39366*x + 216513*x**2 - 734832*x**3 - 761076*x**4 + 7540776*x**5 - 8373780*x**6 - 6955632*x**7 + 26726598*x**8 - 32699700*x**9 + 23810598*x**10 - 11618640*x**11 + 3944844*x**12 - 940536*x**13 + 156204*x**14 - 17616*x**15 + 1281*x**16 - 54*x**17 + x**18 - 98415*y**2 + 314928*x*y**2 - 7085880*x**2*y**2 + 26418960*x**3*y**2 - 44821836*x**4*y**2 + 25812432*x**5*y**2 + 37105128*x**6*y**2 - 84941136*x**7*y**2 + 79679214*x**8*y**2 - 45903024*x**9*y**2 + 17817624*x**10*y**2 - 4791312*x**11*y**2 + 891972*x**12*y**2 - 112464*x**13*y**2 + 9144*x**14*y**2 - 432*x**15*y**2 + 9*x**16*y**2 - 2335716*y**4 + 11739816*x*y**4 - 32819580*x**2*y**4 + 41698800*x**3*y**4 + 8669268*x**4*y**4 - 80740152*x**5*y**4 + 101393532*x**6*y**4 - 70805664*x**7*y**4 + 32218020*x**8*y**4 - 9994248*x**9*y**4 + 2126556*x**10*y**4 - 304848*x**11*y**4 + 28140*x**12*y**4 - 1512*x**13*y**4 + 36*x**14*y**4 - 10368324*y**6 + 16768944*x*y**6 - 2187000*x**2*y**6 - 34718544*x**3*y**6 + 59918940*x**4*y**6 - 52894944*x**5*y**6 + 29369232*x**6*y**6 - 10874592*x**7*y**6 + 2720340*x**8*y**6 - 453840*x**9*y**6 + 48552*x**10*y**6 - 3024*x**11*y**6 + 84*x**12*y**6 + 1201878*y**8 - 6344244*x*y**8 + 15507774*x**2*y**8 - 18891792*x**3*y**8 + 13803156*x**4*y**8 - 6463368*x**5*y**8 + 1988100*x**6*y**8 - 399600*x**7*y**8 + 51030*x**8*y**8 - 3780*x**9*y**8 + 126*x**10*y**8 + 1113750*y**10 - 2518128*x*y**10 + 2951640*x**2*y**10 - 1966608*x**3*y**10 + 808524*x**4*y**10 - 207216*x**5*y**10 + 33096*x**6*y**10 - 3024*x**7*y**10 + 126*x**8*y**10 + 172476*y**12 - 234360*x*y**12 + 159732*x**2*y**12 - 58224*x**3*y**12 + 12684*x**4*y**12 - 1512*x**5*y**12 + 84*x**6*y**12 + 9756*y**14 - 6768*x*y**14 + 2520*x**2*y**14 - 432*x**3*y**14 + 36*x**4*y**14 + 177*y**16 - 54*x*y**16 + 9*x**2*y**16 + y**18

set samples 1024
set isosamples 512,512
set table '../plot_data/im0.dat'; splot im0(x,y)
set table '../plot_data/im01.dat'; splot im01(x,y)
set table '../plot_data/im02.dat'; splot im02(x,y)
set table '../plot_data/im03.dat'; splot im03(x,y)


set isosamples 1024,1024
# set table '../plot_data/re11.dat'; splot re11(x,y)
set table '../plot_data/re12.dat'; splot re12(x,y)

unset table
plot 0 w l, '../plot_data/im01.dat' w l, '../plot_data/im02.dat' w l, '../plot_data/im03.dat' w l,  '../plot_data/re12.dat' w l

pause mouse close
