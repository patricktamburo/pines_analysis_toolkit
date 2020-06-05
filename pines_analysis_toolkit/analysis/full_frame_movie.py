import matplotlib.pyplot as plt
import glob
import numpy as np
import pdb
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

reduced_path = 'P:\\Photometry\\PINES\\202001\\2MASS 0916+1057\\domeflat_reduced\\'
reduced_files = np.sort(glob.glob(reduced_path+'*.fits'))

fig,ax = plt.subplots(1,1)
plt.ion()
for i in range(len(reduced_files)):
    image = fits.open(reduced_files[i])[0].data[200+50+60:350-10,612+150+55:912-50-20]
    avg,med,std = sigma_clipped_stats(image)
    ax.imshow(image,origin='lower',vmin=med+std,vmax=med+10*std)
    plt.pause(0.1)
    plt.cla()
pdb.set_trace()