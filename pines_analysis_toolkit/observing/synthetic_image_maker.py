import numpy as np 

def synthetic_image_maker(x_centroids,y_centroids,fwhm=float):
    """Construct a synthetic image from centroid data

    :param x_centroids: array of x positions
    :type x_centroids: numpy array 
    :param y_centroids: array of y positions
    :type y_centroids: numpy array
    :param fwhm: fwhm of synthetic sources, defaults to 8
    :type fwhm: float, optional
    """
    #Construct synthetic images from centroid/flux data.
    synthetic_image = np.zeros((1024,1024))
    sigma = fwhm/2.355
    for i in range(len(x_centroids)):
        #Cut out little boxes around each source and add in Gaussian representations. This saves time. 
        int_centroid_x = int(np.round(x_centroids[i]))
        int_centroid_y = int(np.round(y_centroids[i]))
        y_cut, x_cut = np.mgrid[int_centroid_y-10:int_centroid_y+10,int_centroid_x-10:int_centroid_x+10]
        dist = np.sqrt((x_cut-x_centroids[i])**2+(y_cut-y_centroids[i])**2)
        synthetic_image[y_cut,x_cut] += np.exp(-((dist)**2/(2*sigma**2)+((dist)**2/(2*sigma**2))))
    return(synthetic_image)