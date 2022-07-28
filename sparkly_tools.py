from astropy.io import fits
from astropy.visualization import PercentileInterval, ImageNormalize
from astropy.visualization import LinearStretch, SqrtStretch, LogStretch
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.wcs import WCS, FITSFixedWarning
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.table import Table

from photutils.centroids import centroid_sources, centroid_com
from photutils.segmentation import detect_threshold, detect_sources
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from photutils.aperture import aperture_photometry

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_zeropoint(filename):
    """
    Computer AB zeropoint from the PHOTFNU keyword in an image header.
    
    Parameters:
    -----------
    
        filename (string): 
            input FITS filename of the large image from the the sub-image is extracted
            
    Returns:
    --------
        zp (float):
            zeropoint (or None if it can't be computed from the header)
            
    """
    try:
        zp = -2.5*np.log10(fits.open(filename)[0].header['PHOTFNU']/3631)
    except:
        zp = None
    return zp
                           

def save_cutout(filename, position, size, postage_filename):
    """
    Extract a postage stamp image from a larger FITS image.

    Parameters:
    -----------
    
        filename (string): 
            input FITS filename of the large image from the the sub-image is extracted
            
        position (list):
            (ra,dec) of the center of object
        
        size (integer)
            size of the box in pixels
        
        postage_filename
            output FITS filename of the small postage stamp image.
            
     """

    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    hdu.writeto(postage_filename, overwrite=True)


def display_data(data, title=None, aperture=None, annulus=None, show_colorbar=False):
    """
    Generic display of an astronomical image that should look decent.

    Parameters:
    -----------
    
        data (Numpy array or string): 
            A two-dimensional array. If a string, try to read in a FITS file 
            using the string as a filename. The data must be in FITS extension 0.
            
        title (string):
            title for the plot (default = None)
        
        aperture (photutils CircularAperture)
            aperture object to display  (default = None)
        
        annulus (photutils CircularAnnulus)
            annulus object to display (default = None)

        show_colorbar (bool)
            display colorbar (default = False)
            
    """
   
    if isinstance(data, str):
        data = fits.getdata(data, ext=0)
 
    plt.figure()
    if title:
        plt.title(title)
    interval = PercentileInterval(95)
    interval.get_limits(data)
    stretch = SqrtStretch()
    norm = ImageNormalize(data, interval=interval, stretch=stretch)
    plt.imshow(data, cmap='gray', norm=norm)
    
    if aperture:
        aperture.plot(color='red', lw=0.5, label='Photometry aperture')
    if annulus:
        annulus.plot(color='yellow', lw=0.2, label='Background annulus')
    if show_colorbar:
    	plt.colorbar()



def add_at(big, pos, small):
    """
    Adds a small array into a big array at a specified position.

    Parameters:
    -----------
    
        big (Numpy array): 
            two-dimensional array that the small array is inserted into.
            
        pos (list):
            position to insert the small array
            
        small (Numpy array): 
            two-dimensional array to be inserted into the big array.

    Returns:
    --------
    
        output_image (Numpy array): 
            big array with small array embedded.
            
    """ 
    x1 = int(pos[0])
    y1 = int(pos[1])
    x2 = x1 + small.shape[0]
    y2 = y1 + small.shape[1]

    assert x2  <= big.shape[0], "the position will make the small matrix exceed the boundaries at x"
    assert y2  <= big.shape[1], "the position will make the small matrix exceed the boundaries at y"

    big[y1:y2,x1:x2] = big[y1:y2,x1:x2] + small

    return big


def create_starfield(original_image, psf_image, position_file, mag, zp, debug=False):
    """
    Return an array of PSF stars at specified positions.

    Parameters:
    -----------
    
        original_image (Numpy array): 
            two-dimensional array used to set the dimensions of the starfield
            
        psf_image (Numpy array): 
            two-dimensional array used as the point-spread functions
          
        position_file (int):
            text file with positions. This should be in the format returned by
            ds9 when the 'c' key is pressed on an image.

        mag (float): 
            magnitude of the star (mag)

        zp (float): 
            photometric zeropoint (mag/ADU)
            
        verbose (bool):
            print diagnostic information to the screen? Default = False.

    Returns:
    --------
    
        output_image (Numpy array): 
            two-dimensional array that can be added back to the original image
            to simulate stars.
            
    """
    
    # Read in the original image and zero it out.
    starfield = fits.getdata(original_image, ext=0)
    starfield_header = fits.getheader(original_image, ext=0)
    starfield_nx = starfield_header['NAXIS1']
    starfield_ny = starfield_header['NAXIS2']
    starfield = starfield * 0.0
    
    # Read in the PSF that we will add to the image
    psf = fits.getdata(psf_image, ext=0)
    psf_header = fits.getheader(psf_image, ext=0)
    psf_nx = psf_header['NAXIS1']
    psf_ny = psf_header['NAXIS2']
    
    # Normalize the PSF to unity
    normalized_psf = psf/(psf.sum())
    
    # Rescale using the zeropoint
    scale_factor = 10.**(-0.4*(mag - zp))
    scaled_psf =  scale_factor * normalized_psf
    
    # Load the list with the positions where the stars should go.
    df = pd.read_csv(position_file, sep=" ", 
                     header=None, usecols=[1,2,4,5], 
                     names=["ra","dec","x","y"])
    x_star = (df["x"]-1).to_numpy()
    y_star = (df["y"]-1).to_numpy()
    
    # Add the PSF to the image at the specified positions.
    for xxi, yyi in zip(x_star, y_star):        
        pos = ( xxi - int(psf_nx/2), yyi - int(psf_ny/2) ) 
        starfield = add_at(starfield, pos, scaled_psf)
        
    if debug:
        hdu = fits.PrimaryHDU(starfield)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto('debug_create_starfield.fits', overwrite=True)
        print("Saving test image: debug_create_starfield.fits")

    return starfield



def analyze_stars(data, zp, position_file, 
                  aperture_radius=5, sky_annulus=[7,10], 
                  centroid_box_size=5,
                  display = True,
                  title = None):
    
    # Load the guesstimated positions 
    df = pd.read_csv(position_file, sep=" ", 
                     header=None,
                     usecols=[1,2,4,5],
                     names=["ra","dec","x","y"])

    # Centroid the positions
    x_guess = (df["x"]-1).to_numpy()
    y_guess = (df["y"]-1).to_numpy()
    x, y = centroid_sources(data, x_guess, y_guess, box_size=centroid_box_size,
                            centroid_func=centroid_com)
    positions = []
    for new_x,new_y in zip(x,y):
        positions.append([new_x,new_y])

    # Define the photometric apertures
    aperture = CircularAperture(positions, r=aperture_radius)
    annulus_aperture = CircularAnnulus(positions, r_in=sky_annulus[0], r_out=sky_annulus[1])
    
    # Optionally display the data
    if display:
        display_data(data, title=title, aperture=aperture, annulus=annulus_aperture)
    
    # Photometer the data
    phot_table = aperture_photometry(data, aperture)

    # Compute the local background
    aperstats = ApertureStats(data, annulus_aperture)
    bkg_median = aperstats.median
    aperture_area = aperture.area_overlap(data)
    total_bkg = bkg_median * aperture_area

    # The final answer!
    mags = zp - 2.5*np.log10(phot_table['aperture_sum'] - total_bkg)
    phot_table.add_column(mags,name='mag')
    return(phot_table)
