from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

def bg_2d(data, box_size=50):
    """Removes large-scale changes in the background of an image. See https://photutils.readthedocs.io/en/stable/background.html.

    :param data: 2D array of pixel values
    :type data: numpy array
    :param box_size: Size of boxes used to do background estimation, defaults to 50
    :type box_size: int, optional
    :return: 2D array of pixel values with the background model subtracted
    :rtype: numpy array
    """
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (box_size, box_size), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    return data - bkg.background
