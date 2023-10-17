XLIBS = -L/usr/X11R6/lib -lX11
XINCLUDE = -I/usr/X11R6/include/

all: active directories plot

active: active.c random_numrec.c 
	gcc active.c random_numrec.c -o active -lm -O3 -mavx -ffast-math -freciprocal-math
directories:
	mkdir -p movies/gradient movies/opposite movies/perpendicular movies/seeds/gradient movies/seeds/opposite movies/seeds/perpendicular
	mkdir -p stat stat/seeds
	mkdir -p logs
	mkdir -p clusters/com_in_time/seeds/gradient clusters/com_in_time/seeds/opposite clusters/com_in_time/seeds/perpendicular clusters/com_in_time/gradient clusters/com_in_time/opposite clusters/com_in_time/perpendicular
	mkdir -p clusters/gradient clusters/opposite clusters/perpendicular clusters/seeds/gradient clusters/seeds/opposite clusters/seeds/perpendicular
	mkdir -p console_outputs povray_anim/seeds
	mkdir -p particles/gradient particles/opposite particles/perpendicular particles/seeds/gradient particles/seeds/opposite particles/seeds/perpendicular
	mkdir -p pinningsites/gradient pinningsites/opposite pinningsites/perpendicular pinningsites/seeds/gradient pinningsites/seeds/opposite pinningsites/seeds/perpendicular
	mkdir -p density_distribution/gradient density_distribution/opposite density_distribution/perpendicular density_distribution/seeds/gradient density_distribution/seeds/opposite density_distribution/seeds/perpendicular
	mkdir -p density_distribution/avg_gradient density_distribution/avg_opposite density_distribution/avg_perpendicular density_distribution/avg_seeds/gradient density_distribution/avg_seeds/opposite density_distribution/avg_seeds/perpendicular
plot: plot.c plot_grid.c
	gcc plot.c -o plot -lm  $(XLIBS) $(XINCLUDE)
	gcc plot_grid.c -o plot_grid -lm  $(XLIBS) $(XINCLUDE)
