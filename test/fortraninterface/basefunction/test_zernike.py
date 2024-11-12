from vacuum.fortraninterface import get_zernike
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True


rarr = np.linspace(0, 1, 64)
valuearr = list()
for _i, r in enumerate(rarr):
    valuearr.append(get_zernike(r, 4, 2))
valuearr = np.array(valuearr)
    
lradarr = [0, 1, 2, 2, 3, 4]
mpolarr = [0, 1, 0, 2, 1, 0]
    
fig, ax = plt.subplots()
for index, lrad in enumerate(lradarr):
    mpol = mpolarr[index]
    ax.plot(
        rarr,
        valuearr[:, lrad, mpol, 0],
        label = r'$\hat{\mathcal{R}}_'+str(lrad)+'^'+str(mpol)+'$'
    )
fig.legend()
fig.savefig('zernike.png', dpi=1000)