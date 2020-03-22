import matplotlib
#matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from matplotlib import patches
import matplotlib.lines as mlines
#rc("text", usetex=True)
rc("font", family="serif")
rc("axes", grid=True)
rc("grid", linestyle="--")
rc("xtick", direction="in")
rc("ytick", direction="in")
rc("savefig", format="pdf", bbox="tight")
import numpy as np
import os

time = []
ninfected = []
vulnerable = []
recovered = []
dead = []

for fnum in range(200):
	print("fnum = %d" % fnum)

	# load discrete
	try:
		data = np.loadtxt("data/discrete_%05d.dat" % fnum)
	except:
		break
	time.append(data[0])
	population = data[1]
	contag_t0 = data[2]
	contag_t1 = data[3]
	ninfected.append(data[4])
	vulnerable.append(data[5])
	recovered.append(data[6])
	dead.append(data[7])

	# load infected
	data = np.transpose(np.loadtxt("data/infected_%05d.dat" % fnum))
	t = data[0]
	n = data[1]

	fig = plt.figure(figsize=(6, 3))
	gs = gridspec.GridSpec(1, 1)

	ax = fig.add_subplot(gs[0,0])
	ax.plot(t, n, "C1.-")
	ax.axvline(contag_t0, color="black")
	ax.axvline(contag_t1, color="black")
	ax.set_xlabel("days since infection")
	ax.set_ylabel("infected number per day")
	ax.set_title("t = %.3f days" % (time[-1]));

	gs.tight_layout(fig)
	#plt.show()
	plt.savefig("img/infected_%05d.png" % fnum)
	plt.close()

ninfected = np.array(ninfected)
vulnerable = np.array(vulnerable)
recovered = np.array(recovered)
dead = np.array(dead)

fig = plt.figure(figsize=(12, 15))
gs = gridspec.GridSpec(5, 2)

for i in range(2):
	ax = fig.add_subplot(gs[0,i])
	ax.plot(time, ninfected/population, "C1.-")
	ax.set_ylim([0,1])
	if i == 1:
		ax.set_ylim([1e-5,1])
		ax.set_yscale("log")
	ax.set_xlabel("day")
	ax.set_ylabel("ninfected fraction")

	ax = fig.add_subplot(gs[1,i])
	ax.plot(time, vulnerable/population, "C0.-")
	ax.set_ylim([0,1])
	if i == 1:
		ax.set_ylim([1e-5,1])
		ax.set_yscale("log")
	ax.set_xlabel("day")
	ax.set_ylabel("vulnerable fraction")

	ax = fig.add_subplot(gs[2,i])
	ax.plot(time, recovered/population, "C2.-")
	ax.set_ylim([0,1])
	if i == 1:
		ax.set_ylim([1e-5,1])
		ax.set_yscale("log")
	ax.set_xlabel("day")
	ax.set_ylabel("recovered fraction")

	ax = fig.add_subplot(gs[3,i])
	ax.plot(time, dead, "C3.-")
	if i == 1:
		ax.set_yscale("log")
	ax.set_xlabel("day")
	ax.set_ylabel("number dead")

	if i == 0:
		ax = fig.add_subplot(gs[4,i])
		ax.plot(time, (ninfected + vulnerable + recovered + dead) / population, "C4.-")
		ax.set_ylim([0.999,1.001])
		ax.set_xlabel("day")
		ax.set_ylabel("total frac (check conservation)")

gs.tight_layout(fig)
#plt.show()
plt.savefig("img/discrete_timeseries.png")
plt.close()
