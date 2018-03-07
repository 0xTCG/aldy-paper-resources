# 786

# Plot supplementary plots from Aldy

#%%
%matplotlib inline
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import re, os, collections
from matplotlib import rc

#%% Load data
profiles = collections.defaultdict(list)
for sample in 'NA11892 NA12878 NA12155'.split():
	with open('{}.bam.profile'.format(sample)) as f:
		pos = 0
		for l in f:
			l = l.strip().split()
			profiles[sample].append(float(l[0]))
			pos += 1

def load_profile(profile_path):
	result = []
	with open(profile_path) as f: # now it is NA19686
		for l in f:
			l = l.strip().split()
			if l[0] != 'chr22': continue
			l = map(int, l[1:])
			l[0] -= 42518969
			if l[0] < 0: continue
			#print l, len(result)
			result += [0] * (l[0] - len(result))
			result.append(l[1])
	return result
pgrnref = load_profile('pgrnseq.profile')
len(pgrnref)
xprofiles = {}
for pn, profile in profiles.iteritems():
	r = (28495, 29281)
	rc = 2 * float(sum(profile[slice(*r)])) / sum(pgrnref[slice(*r)])
	scale = 2.0 / rc
	xprofiles[pn] = [i for i in profile]
	print pn, scale
	profiles[pn] = [scale * i for i in profile]
CYP_REGIONS = { # intervals are inclusive on both sides
	'6.9e': ( 3605,  3783), '6.9i': ( 3300,  3605),  '6.0i': ( 7822,  8000),
	'6.8i': ( 3784,  3881),  '6.8e': ( 3882,  4023),
	'6.7i': ( 4024,  4477),  '6.7e': ( 4478,  4665),
	'6.6i': ( 4666,  4872),  '6.6e': ( 4873,  5014),
	'6.5i': ( 5015,  5204),  '6.5e': ( 5205,  5381),
	'6.4i': ( 5382,  5814),  '6.4e': ( 5815,  5975),
	'6.3i': ( 5976,  6063),  '6.3e': ( 6064,  6216),
	'6.2i': ( 6217,  6768),  '6.2e': ( 6769,  6940),
	'6.1i': ( 6941,  7642),  '6.1e': ( 7643,  7822),

	 '7.9i': (17000, 17316),  '7.9e': (17317, 17495),
	'7.8i': (17496, 17593),  '7.8e': (17594, 17735),
	'7.7i': (17736, 18189),  '7.7e': (18190, 18377),
	'7.6i': (18378, 18571),  '7.6e': (18572, 18713),
	'7.5i': (18714, 18905),  '7.5e': (18906, 19082),
	'7.4i': (19083, 19507),  '7.4e': (19508, 19668),
	'7.3i': (19669, 19756),  '7.3e': (19757, 19909),
	'7.2i': (19910, 20438),  '7.2e': (20439, 20610),
	'7.1i': (20611, 21312),  '7.1e': (21313, 21493),  '7.0i': (21493, 22000)
}

#%% Main function: enable usetex for publication-ready plots
#   TeX might need: tlmgr install type1cm dvipng
import matplotlib
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
matplotlib.rcParams['figure.autolayout'] = False

def xplot(title, mima=None, **kwargs):
	matplotlib.rcParams['text.usetex'] = True
	matplotlib.rcParams['font.size'] = 12

	if len(kwargs) > 1:
		matplotlib.rcParams['figure.figsize'] = (18, 5)
		x = kwargs.values()[1][0][:]
		y = [0] * len(x)
		for rn, r in CYP_REGIONS.iteritems():
			w = sum(x[r[0]:r[1]]) / float(sum(pgrnref[r[0]:r[1]]))
			for i in xrange(r[0], r[1]+1):
				y[i] = w
		f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
	else:
		matplotlib.rcParams['figure.figsize'] = (18, 2)
		f, (ax1, ax2, ax3) = plt.subplots(1, 3)
		ax4, ax5, y = None, None, None

	for axi, ax in enumerate([ax1, ax2, ax3, ax4, ax5]):
		if ax is None: continue
		for (a, b) in [(5483,5995),(3885,4667),(19335,19759),(17598,18380)]:
			ax.axvspan(a, b, color='orange', alpha=0.5, lw=0)
	#    ax.fill_between(range(len(mi)), av-2*st, av+2*st, facecolor='red', lw=0)
		if mima is not None:
			mi, ma = mima
			ax.fill_between(range(len(mi)), mi, ma, facecolor='red', lw=0)
		# ax.plot([1] * len(pgrnref), color='black')

		if axi < 3:
			ax.plot(pgrnref, '--', color='black')
			if len(kwargs) > 1:
				ax.plot(kwargs.values()[1][0], label='', **kwargs.values()[1][1])
			ax.set_ylim(0, 3000)
		elif len(kwargs) > 1:
			ax.plot([1] * len(pgrnref), color='black')
			ax.plot([.5] * len(pgrnref), **kwargs.values()[0][1])
			ax.plot([1.5] * len(pgrnref), **kwargs.values()[0][1])
			ax.plot([2.5] * len(pgrnref), **kwargs.values()[0][1])
			ax.plot(y, label='', **kwargs.values()[1][1])
			ax.set_ylim(0, 2.5)
			ax.set_yticklabels(['0', '1', '2', '3', '4', '5'])

		ax.xaxis.tick_top()

		cr = sorted([(a[1], v) for v, a in CYP_REGIONS.iteritems() if not v[2:] in ['9i', '0i', '6i', '8i', '7i', '5i', '3i']])
		ax.set_xticks([x for x, _ in cr])
		ax.set_xticklabels([('exon ' if x[-1] == 'e' else 'intron ') + x[-2] for _, x in cr])
		for tick in ax.get_xticklabels():
			tick.set_fontsize(8)
			tick.set_rotation(45)

	ax1.set_xlim(3300, 8000)
	if ax4 is not None:
		ax4.set_xlim(3300, 8000)
		ax5.set_xlim(17000, 22000)
		ax6.axis('off')
	ax2.set_xlim(17000, 22000)
	ax3.set_xlim(26900, 32000)

	# hide the spines between ax and ax2
	ax1.spines['right'].set_visible(False)
	#ax2.spines['right'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax3.spines['left'].set_visible(False)

	if ax4 is not None:
		ax4.spines['right'].set_visible(False)
		ax5.spines['left'].set_visible(False)
		ax5.tick_params(labeltop='off')
		ax4.tick_params(labeltop='off')
		ax5.yaxis.tick_right()


	ax1.yaxis.tick_left()
	ax1.tick_params(labelright='off')
	#ax2.tick_params(labelright='off')
	ax2.tick_params(labelleft='off')
	ax3.yaxis.tick_right()

	ax1.set_xlabel('CYP2D6')
	ax2.set_xlabel('CYP2D7')
	ax3.set_xlabel('CYP2D8')

	f.subplots_adjust(hspace=.7)

	d = .02
	kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
	ax1.plot((1-d,1+d), (-d,+d), **kwargs)
	ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)

	kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
	ax2.plot((-d,+d), (1-d,1+d), **kwargs)
	ax2.plot((-d,+d), (-d,+d), **kwargs)

	kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
	ax3.plot((-d,+d), (1-d,1+d), **kwargs)
	ax3.plot((-d,+d), (-d,+d), **kwargs)

	if ax4 is not None:
		kwargs = dict(transform=ax4.transAxes, color='k', clip_on=False)
		ax4.plot((1-d,1+d), (-d,+d), **kwargs)
		ax4.plot((1-d,1+d),(1-d,1+d), **kwargs)

		kwargs.update(transform=ax5.transAxes)  # switch to the bottom axes
		ax5.plot((-d,+d), (1-d,1+d), **kwargs)
		ax5.plot((-d,+d), (-d,+d), **kwargs)

		ax1.set_ylabel('Rescaled\nPGRN-Seq', rotation=90, size='12', x=1.1)
		ax4.set_ylabel('Normalized\nPGRN-Seq', rotation=90, size='12', x=1.2)
		# st=f.suptitle(title, fontsize=14, y=1.30)
	#else:
		# st=f.suptitle(title, fontsize=14, y=1.55)

	st = None
	return (ax1, ax2, ax3, ax4, ax5, y,f)


#%% Plot figure (1)
single = ([float(x)/2 for x in pgrnref], dict(color='black', linestyle='--', lw=.5))
p = xprofiles['NA11892']
ax1, ax2, ax3, ax4, ax5, y, fig = xplot('', None, insertion=(p, dict(color='red')), single=single)

ax1.annotate('$B_g$', xy=(155,148+55), xycoords='figure points', size='14')
ax1.annotate('$C_g$', xy=(155,190+55), xycoords='figure points', color='red', size='14')

w=sum(p[i] for i in range(28495, 29281))/sum(pgrnref[i] for i in range(28495, 29281))
print w
ax3.axvspan(28495, 29281, color='purple', alpha=0.5, lw=0)
ax3.annotate('$q$', xy=(870+12,225+55), xycoords='figure points', color='white', size='22')
ax3.annotate('$\\eta=\\frac{C_g(q)}{B_g(q)}=1.48$', xy=(905+12,220+55), xycoords='figure points', color='purple')

for ax in [ax4, ax5]:
	ax.cla()
	for (a, b) in [(5483,5995),(3885,4667),(19335,19759),(17598,18380)]:
		ax.axvspan(a, b, color='orange', alpha=0.5, lw=0)
	ax.plot(xprofiles['NA11892'], color='red')
	ax.plot([i*w for i in pgrnref], '--', color='purple', lw=2)
	ax.set_ylim(0, 3000)
ax4.set_xlim(3300, 8000)
ax5.set_xlim(17000, 22000)
ax4.set_ylabel('Rescaled\nPGRN-Seq', rotation=90, size='12', x=1.2)

ax4.annotate('$C_g$', xy=(140+12,35), xycoords='figure points', color='red', size='14')
ax5.annotate('$R_g = \\eta\\times B_g$', xy=(128+12,80), xycoords='figure points', color='purple', size='14')
plt.savefig('fig_sup_1.pdf', dpi=fig.dpi, bbox_inches='tight')


#%% Plot figure (2)/(ii)
p = profiles['NA12155']
ax1, _, ax3, ax4, ax5, y, st = xplot('(ii) \\textit{CYP2D6}: *1/*5 (whole gene deletion)', None, deletion=(p, dict(color='red')), single=single)
ax1.fill_between(range(len(pgrnref)), p[:len(pgrnref)], pgrnref, facecolor='red', lw=0, alpha=0.5, hatch='x')
ax4.fill_between(range(10000), y[:10000], [1] * 10000, facecolor='red', lw=0, alpha=0.5, hatch='x')

ax1.annotate('$R_g$', xy=(175,235), xycoords='figure points', size='13')
ax1.annotate('$C_g$', xy=(175,204), xycoords='figure points', color='red', size='13')

ax4.annotate('$\\frac{C_g}{R_g}$', xy=(300,28), xycoords='figure points', color='red', size='13')

ax3.cla()
ax3.axis('off')
plt.savefig('fig_sup_2b.pdf',  dpi=fig.dpi, bbox_inches='tight')


#%% Plot figure (2)/(iii)
p = profiles['NA12878']
ax1, ax2, ax3, ax4, ax5, y, st = xplot('(iv) \\textit{CYP2D6}: *10/*68+*4 (gene fusion in intron 1)', None, insertion=(p, dict(color='red')), single=single)
r = (6941, 8000)
ax1.fill_between(range(r[0], r[1]), p[r[0]:r[1]], pgrnref[r[0]:r[1]], facecolor='green', lw=0, alpha=0.5, hatch='x')
ax4.fill_between(range(r[0], r[1]), y[r[0]:r[1]], [1] * (r[1]-r[0]), facecolor='green', lw=0, alpha=0.5, hatch='x')
r = (17246, 21310)
ax2.fill_between(range(r[0], r[1]), p[r[0]:r[1]], pgrnref[r[0]:r[1]], facecolor='green', lw=0, alpha=0.5, hatch='x')
ax5.fill_between(range(r[0], r[1]), y[r[0]:r[1]], [1] * (r[1]-r[0]), facecolor='green', lw=0, alpha=0.5, hatch='x')

ax3.cla()
ax3.axis('off')

ax1.annotate('$R_g$', xy=(175,235), xycoords='figure points', size='13')
ax1.annotate('$C_g$', xy=(175,204), xycoords='figure points', color='red', size='13')

ax4.annotate('$\\frac{C_g}{R_g}$', xy=(300,43), xycoords='figure points', color='green', size='13')
plt.savefig('fig_sup_2c.pdf',  dpi=fig.dpi, bbox_inches='tight')
