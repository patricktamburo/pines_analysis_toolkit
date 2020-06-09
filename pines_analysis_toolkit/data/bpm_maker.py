import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import glob
import pickle

def bpm_maker():
    #THIS NEEDS TO BE RUN MULTIPLE TIMES, COMMENTING OUT APPROPRIATE BLOCKS EACH TIME.

    bias_path = 'P:\\Photometry\\PINES\\2MASS 0344+0111\\bias\\master_bias.fits'
    bias_image = fits.open(bias_path)[0].data
    bad_pixel_mask = np.ones(np.shape(bias_image))
    #Starts in bottom left, moves right across the row, then up. 

    # Read in a Kokopelli mask
    bad_pixel_x = []
    bad_pixel_y = []

    # #FIRST PASS: Just do the Kokopelli mask. 
    # crack_path = 'C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\crack_mask.fits'
    # crack_mask = fits.open(crack_path)[0].data[0:1024,:]
    # for y in range(1024):
    #     for x in range(1024):
    #         if crack_mask[y,x] == 0:
    #             bad_pixel_mask[y,x] = 0
    #             bad_pixel_x.append(x)
    #             bad_pixel_y.append(y)
    #     print('Done row ',y+1, ' of 1024, first pass.')

    # pdb.set_trace()
    # hdu = fits.PrimaryHDU(bad_pixel_mask)
    # hdu.writeto('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\'+'bpm1.fits')

    # #SECOND PASS:
    # threshold = 7
    # print('Restoring previous pass-1 bpm.')
    # bad_pixel_mask = fits.open('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\bpm1.fits')[0].data
    # #It's necessary to loop over the image again, this time ignoring previously-flagged pixels. 
    # bad_pixel_mask_2 = np.ones(np.shape(bad_pixel_mask))
    # bad_pixel_x = []
    # bad_pixel_y = []

    # for y in range(1024):
    #     for x in range(1024):
    #         if (3 < x < 1000) and (3 < y < 1000): #Don't worry about edge pixels for now.
    #             pixel_val = bias_image[y,x]
    #             adjacent_pixel_vals = []
    #             x_pos = np.linspace(x-4,x+4,9)
    #             y_pos = np.linspace(y-4,y+4,9)
    #             x_grid,y_grid = np.meshgrid(x_pos,y_pos)

    #             #Loop over the 80 adjacent pixels
    #             for i in range(len(x_grid)):
    #                 for k in range(len(x_grid[0])):
    #                     if (i == 4) and (k == 4):
    #                         continue
    #                     else:
    #                         #The second loop change: exclude previously flagged bad pixels from the adjacent pixels set. 
    #                         if bad_pixel_mask[int(y_grid[i][k]),int(x_grid[i][k])] == 1:
    #                             adjacent_pixel_vals.append(bias_image[int(y_grid[i][k]),int(x_grid[i][k])])
                
    #             median_adjacent_pixels = np.median(adjacent_pixel_vals)
    #             stdev_adjacent_pixels = np.std(adjacent_pixel_vals)
                
    #             #If the pixel's value deviates by more than the user-chosen threshold from its 80 neighbors, or
    #             #   if the pixel has already been flagged in crack_mask.fits, add it to the bad pixel mask!
    #             if ((abs(pixel_val-median_adjacent_pixels) > threshold * stdev_adjacent_pixels)) or (bad_pixel_mask[y,x] == 0):
    #                 bad_pixel_mask_2[y,x] = 0
    #                 bad_pixel_x.append(x)
    #                 bad_pixel_y.append(y)
    #     print('Done row ',y+1, ' of 1024, second pass.')

    # hdu = fits.PrimaryHDU(bad_pixel_mask_2)
    # hdu.writeto('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\'+'bpm2.fits')


    # #THIRD PASS:
    # threshold = 7
    # print('Restoring previous pass-2 bpm.')
    # bad_pixel_mask_2 = fits.open('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\bpm2.fits')[0].data
    # #Do another loop, ignoring those pixels already flagged in bpm2.fits.
    # bad_pixel_mask_3 = np.ones(np.shape(bad_pixel_mask_2))
    # bad_pixel_x = []
    # bad_pixel_y = []

    # for y in range(1024):
    #     for x in range(1024):
    #         if (3 < x < 1000) and (3 < y < 1000): #Don't worry about edge pixels for now.
    #             pixel_val = bias_image[y,x]
    #             adjacent_pixel_vals = []
    #             x_pos = np.linspace(x-4,x+4,9)
    #             y_pos = np.linspace(y-4,y+4,9)
    #             x_grid,y_grid = np.meshgrid(x_pos,y_pos)

    #             #Loop over the 80 adjacent pixels
    #             for i in range(len(x_grid)):
    #                 for k in range(len(x_grid[0])):
    #                     if (i == 4) and (k == 4):
    #                         continue
    #                     else:
    #                         if bad_pixel_mask_2[int(y_grid[i][k]),int(x_grid[i][k])] == 1:
    #                             adjacent_pixel_vals.append(bias_image[int(y_grid[i][k]),int(x_grid[i][k])])
                
    #             median_adjacent_pixels = np.median(adjacent_pixel_vals)
    #             stdev_adjacent_pixels = np.std(adjacent_pixel_vals)
                
    #             #If the pixel's value deviates by more than the user-chosen threshold from its 80 neighbors, or
    #             #   if the pixel has already been flagged in crack_mask.fits, add it to the bad pixel mask!
    #             if ((abs(pixel_val-median_adjacent_pixels) > threshold * stdev_adjacent_pixels)) or (bad_pixel_mask_2[y,x] == 0):
    #                 bad_pixel_mask_3[y,x] = 0
    #                 bad_pixel_x.append(x)
    #                 bad_pixel_y.append(y)
    #     print('Done row ',y+1, ' of 1024, third pass.')

    # #FOURTH PASS:
    # threshold = 5
    # print('Restoring previous pass-3 bpm.')
    # bad_pixel_mask_3 = fits.open('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\bpm3.fits')[0].data
    # #Do another loop, ignoring those pixels already flagged in bpm2.fits.
    # bad_pixel_mask_4 = np.ones(np.shape(bad_pixel_mask_3))
    # bad_pixel_x = []
    # bad_pixel_y = []

    # for y in range(1024):
    #     for x in range(1024):
    #         if (3 < x < 1000) and (3 < y < 1000): #Don't worry about edge pixels for now.
    #             pixel_val = bias_image[y,x]
    #             adjacent_pixel_vals = []
    #             x_pos = np.linspace(x-4,x+4,9)
    #             y_pos = np.linspace(y-4,y+4,9)
    #             x_grid,y_grid = np.meshgrid(x_pos,y_pos)

    #             #Loop over the 80 adjacent pixels
    #             for i in range(len(x_grid)):
    #                 for k in range(len(x_grid[0])):
    #                     if (i == 4) and (k == 4):
    #                         continue
    #                     else:
    #                         if bad_pixel_mask_3[int(y_grid[i][k]),int(x_grid[i][k])] == 1:
    #                             adjacent_pixel_vals.append(bias_image[int(y_grid[i][k]),int(x_grid[i][k])])
                
    #             median_adjacent_pixels = np.median(adjacent_pixel_vals)
    #             stdev_adjacent_pixels = np.std(adjacent_pixel_vals)
                
    #             #If the pixel's value deviates by more than the user-chosen threshold from its 80 neighbors, or
    #             #   if the pixel has already been flagged in crack_mask.fits, add it to the bad pixel mask!
    #             if ((abs(pixel_val-median_adjacent_pixels) > threshold * stdev_adjacent_pixels)) or (bad_pixel_mask_3[y,x] == 0):
    #                 bad_pixel_mask_4[y,x] = 0
    #                 bad_pixel_x.append(x)
    #                 bad_pixel_y.append(y)
    #     print('Done row ',y+1, ' of 1024, third pass.')

    # hdu = fits.PrimaryHDU(bad_pixel_mask_4)
    # hdu.writeto('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\'+'bpm4.fits')

    # #FIFTH PASS:
    # threshold = 5
    # print('Restoring previous pass-4 bpm.')
    # bad_pixel_mask_4 = fits.open('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\bpm4.fits')[0].data
    # #Do another loop, ignoring those pixels already flagged in bpm2.fits.
    # bad_pixel_mask_5 = np.ones(np.shape(bad_pixel_mask_4))
    # bad_pixel_x = []
    # bad_pixel_y = []

    # for y in range(1024):
    #     for x in range(1024):
    #         if (3 < x < 1000) and (3 < y < 1000): #Don't worry about edge pixels for now.
    #             pixel_val = bias_image[y,x]
    #             adjacent_pixel_vals = []
    #             x_pos = np.linspace(x-4,x+4,9)
    #             y_pos = np.linspace(y-4,y+4,9)
    #             x_grid,y_grid = np.meshgrid(x_pos,y_pos)

    #             #Loop over the 80 adjacent pixels
    #             for i in range(len(x_grid)):
    #                 for k in range(len(x_grid[0])):
    #                     if (i == 4) and (k == 4):
    #                         continue
    #                     else:
    #                         if bad_pixel_mask_4[int(y_grid[i][k]),int(x_grid[i][k])] == 1:
    #                             adjacent_pixel_vals.append(bias_image[int(y_grid[i][k]),int(x_grid[i][k])])
                
    #             median_adjacent_pixels = np.median(adjacent_pixel_vals)
    #             stdev_adjacent_pixels = np.std(adjacent_pixel_vals)
                
    #             #If the pixel's value deviates by more than the user-chosen threshold from its 80 neighbors, or
    #             #   if the pixel has already been flagged in crack_mask.fits, add it to the bad pixel mask!
    #             if ((abs(pixel_val-median_adjacent_pixels) > threshold * stdev_adjacent_pixels)) or (bad_pixel_mask_4[y,x] == 0):
    #                 bad_pixel_mask_5[y,x] = 0
    #                 bad_pixel_x.append(x)
    #                 bad_pixel_y.append(y)
    #     print('Done row ',y+1, ' of 1024, third pass.')

    # hdu = fits.PrimaryHDU(bad_pixel_mask_5)
    # hdu.writeto('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\'+'bpm5.fits')

    #SIXTH PASS:
    threshold = 5
    print('Restoring previous pass-5 bpm.')
    bad_pixel_mask_5 = fits.open('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\bpm5.fits')[0].data
    #Do another loop, ignoring those pixels already flagged in bpm2.fits.
    bad_pixel_mask_6 = np.ones(np.shape(bad_pixel_mask_5))
    bad_pixel_x = []
    bad_pixel_y = []

    for y in range(1024):
        for x in range(1024):
            if (3 < x < 1000) and (3 < y < 1000): #Don't worry about edge pixels for now.
                pixel_val = bias_image[y,x]
                adjacent_pixel_vals = []
                x_pos = np.linspace(x-4,x+4,9)
                y_pos = np.linspace(y-4,y+4,9)
                x_grid,y_grid = np.meshgrid(x_pos,y_pos)

                #Loop over the 80 adjacent pixels
                for i in range(len(x_grid)):
                    for k in range(len(x_grid[0])):
                        if (i == 4) and (k == 4):
                            continue
                        else:
                            if bad_pixel_mask_5[int(y_grid[i][k]),int(x_grid[i][k])] == 1:
                                adjacent_pixel_vals.append(bias_image[int(y_grid[i][k]),int(x_grid[i][k])])
                
                median_adjacent_pixels = np.median(adjacent_pixel_vals)
                stdev_adjacent_pixels = np.std(adjacent_pixel_vals)
                
                #If the pixel's value deviates by more than the user-chosen threshold from its 80 neighbors, or
                #   if the pixel has already been flagged in crack_mask.fits, add it to the bad pixel mask!
                if ((abs(pixel_val-median_adjacent_pixels) > threshold * stdev_adjacent_pixels)) or (bad_pixel_mask_5[y,x] == 0):
                    bad_pixel_mask_6[y,x] = 0
                    bad_pixel_x.append(x)
                    bad_pixel_y.append(y)
        print('Done row ',y+1, ' of 1024, third pass.')

    hdu = fits.PrimaryHDU(bad_pixel_mask_6)
    hdu.writeto('C:\\Users\\tambu\\Documents\\python_programs\\PINES\\Reduction\\Misc\\'+'bpm6.fits')


    plt.ion()
    plt.imshow(bias_image,origin='lower',vmin=-50,vmax=50)
    plt.plot(bad_pixel_x,bad_pixel_y,'rx',alpha=0.5)
    pdb.set_trace()