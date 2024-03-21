import skimage.io
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
##put in files here
d = skimage.io.imread('*png')
img = skimage.io.imread('*.tiff')
fig, ax = plt.subplots(1,2, figsize=(20,10))
ax[0].text(50, 100, 'original image', fontsize=16, bbox={'facecolor': 'white', 'pad': 6})
ax[0].imshow(img)

ax[1].text(50, 100, 'depth map', fontsize=16, bbox={'facecolor': 'white', 'pad': 6})
ax[1].imshow(d)
d = np.flipud(d)
img = np.flipud(img)
fig = plt.figure(figsize=(15,10))
ax = plt.axes(projection='3d')

STEP = 5
for x in range(0, img.shape[0], STEP):
    for y in range(0, img.shape[1], STEP):
        ax.scatter(
            d[x,y], y, x,
            c=[tuple(img[x, y, :3]/255)], s=3)      
    ax.view_init(15, 165)
    
