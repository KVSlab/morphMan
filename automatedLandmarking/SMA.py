from scipy import signal
import numpy as np
xs = np.arange(0, 3*np.pi, 0.05)
data = np.sin(xs)
import matplotlib.pyplot as plt
peakind = signal.find_peaks_cwt(data, np.arange(1,10))
peakind, xs[peakind], data[peakind]
plt.plot(data)
plt.show()
print data
print np.arange(1,10)
print peakind
