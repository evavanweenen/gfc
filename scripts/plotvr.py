import gaia_fc as g
from matplotlib import pyplot as plt
import numpy as np

def RV_distribution(v_r, p_v_r, saveto=None, name=None, real_vr=None, entropy_text=True, labels=None, **kwargs):
    try:
        plt.figure()
        if len(p_v_r.shape) == 1: # only one distribution
            plt.plot(v_r, p_v_r, c='k', lw=2, label="$P(v_r)$")
        else:
            if labels is None:
                labels = ["${0}$".format(k) for k, p in enumerate(p_v_r)]
            for p, l in zip(p_v_r, labels):
                plt.plot(v_r, p, lw=2, label=l)
        if real_vr is not None:
            plt.axvline(real_vr, c='k', lw=2, ls="--", label="Observed")
        plt.xlabel(r"$v_r$")
        plt.ylabel(r"$P (v_r)$")
        plt.axis("tight")
        if name is not None:
            plt.text(0.05, 0.95, str(name), transform=plt.gca().transAxes, ha="left")
        if entropy_text:
            if len(p_v_r.shape) == 1:
                entropies = ["$H_1$: {0:.3f}".format(entropy(v_r, p_v_r))]
            else:
                entropies = ["$H_{0}$: {1:.3f}\n".format(k+1, entropy(v_r, p)) for k, p in enumerate(p_v_r)]
            entropiestr = reduce(add, entropies)
            plt.text(0.7, 0.95, entropiestr, transform=plt.gca().transAxes, va="top")
        plt.legend(loc="center right")
        if kwargs:
            plt.setp(plt.gca(), **kwargs)
        plt.tight_layout()
        show_or_save(saveto)
    finally:
        plt.close()
        
id0 = 3180473447906714368
vr = np.load("tgas/vr/xaxis.npy")
pvr = np.load("tgas/vr/{0}.npy".format(id0))
H_c = 3.15022324151
H_m = 4.72408456142
HRV =  41.684
eHRV =  0.879
lower_c = 10.9
upper_c = 64.5

name = "Gaia\n{0}".format(id0)

plt.figure(figsize=(10,10))

plt.plot(vr, pvr[0], lw=2, c='r', label = '$p(v_r|\mathbf{v}_t)$')
plt.plot(vr, pvr[1], lw=2, c='b', label = '$p(v_r)$')

plt.axvline(HRV, label = '$v_r$ (RAVE)', c='k', lw=2, ls='--')
plt.axvspan(HRV-eHRV, HRV+eHRV, color='0.5', alpha=0.5)
plt.axvline(lower_c, c='r', ls='--', lw=2)
plt.axvline(upper_c, c='r', ls='--', lw=2)

plt.axis("tight")
plt.ylim(0, pvr.max()*1.05)
plt.xlim(vr[0]-1, vr[-1]+1)
plt.xlabel("$v_r$ (km s$^{-1}$)")
plt.ylabel("$p(v_r)$")

entropiestr = "$H_C = {hc:.3f}$\n$H_M = {hm:.3f}$".format(hc = H_c, hm = H_m)
plt.text(0.9, 0.9, entropiestr, transform=plt.gca().transAxes, va="top", ha="right")
plt.text(0.05, 0.9, name, transform=plt.gca().transAxes, ha="left")

plt.show()
