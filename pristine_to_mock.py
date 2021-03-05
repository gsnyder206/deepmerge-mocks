#code for converting pristine simulated images to realistic mocks
import astropy
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import numpy as np
import astropy.units as u
import glob
import sys
import os
import photutils
from astropy.stats import gaussian_fwhm_to_sigma
import scipy.ndimage
import scipy as sp


#modify filename and metadata from input image, if desired
def output_pristine_fits_image(image_file,out_file):

    #erg/s/cm^2
    #probably a lambda*flambda
    #input images from Ray Sharma are in log this unit
    fo=fits.open(image_file,'readonly')

    outfo=fits.HDUList()
    outfo.append(fo[0])

        
    im=fo[0].data
    header=fo[0].header

    #hand-coded for now!!!!
    lambda_eff_microns=0.814*u.micron
    nu_eff_hz = astropy.constants.c.to(u.micron/u.s)/lambda_eff_microns
    
    #original image flux values
    #F_nu = lambda*flambda/nu
    image_flux = ((10.0**im)*u.erg/(u.s*u.cm*u.cm))/(nu_eff_hz)
    image_njy = image_flux.to('nJy')

    image_hdu=fits.ImageHDU(image_njy.value,header=header)

    image_hdu.header['EXTNAME']='SimulatedImage'
    image_hdu.header['BUNIT']='nanojanskies'



    #HAND CODED FOR NOW
    redshift=0.49253732071
    image_hdu.header['REDSHIFT']=redshift
    
    
    #hand-coded for now!!!
    fov_kpc=200.0
    npix=image_njy.shape[0]
    
    image_hdu.header['PIX_KPC']=(fov_kpc/npix,'kpc')
    
    
    outfo.append(image_hdu)
    outfo.writeto(out_file,overwrite=True)
    
    return

#Accepts PSF FWHM in arcseconds -- this FWHM can be estimated as the observatory diffraction limit at the associated wavelength
def convolve_with_fwhm(in_image, fwhm_arcsec=0.10):

    #open pristine image fits file
    in_fo=fits.open(in_image,'append')

    #load image data and metadata
    image_in=in_fo['SimulatedImage'].data
    header_in=in_fo['SimulatedImage'].header

    #Choi et al. 2018
    this_cosmology=astropy.cosmology.FlatLambdaCDM(0.72*100.0,0.26,Ob0=0.044)
    kpc_per_arcsec = this_cosmology.kpc_proper_per_arcmin(header_in['REDSHIFT']).value/60.0
    
    #calculate PSF width in pixel units
    pixel_size_arcsec=header_in['PIX_KPC']/kpc_per_arcsec
    sigma_arcsec=fwhm_arcsec*gaussian_fwhm_to_sigma
    sigma_pixels=sigma_arcsec/pixel_size_arcsec
    
    image_out=sp.ndimage.filters.gaussian_filter(image_in,sigma_pixels,mode='nearest')
    hdu_out = fits.ImageHDU(image_out,header=header_in)
    hdu_out.header['FWHMPIX']=(sigma_pixels/gaussian_fwhm_to_sigma,'pixels')
    hdu_out.header['FWHM']=(fwhm_arcsec,'arcsec')
    hdu_out.header['SIGMAPIX']=(sigma_pixels,'pixels')
    hdu_out.header['SIGMA']=(sigma_arcsec,'arcsec')
    hdu_out.header['EXTNAME']='MockImage_Noiseless'
    hdu_out.header['PIXSIZE']=(pixel_size_arcsec,'arcsec')
    
    in_fo.append(hdu_out)
    in_fo.flush()
    
    return

#accepts sb_maglim value which corresponds to magnitudes per square arcsecond
def add_simple_noise(in_image, sb_maglim=25.0, sb_label='25'):

    #algorithm from Snyder et al. (2019)
    #            sigma_nJy = (2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*maglim))*analysis_object.pixsize_arcsec[i]*(3.0*analysis_object.psf_fwhm_arcsec[i])

    in_fo=fits.open(in_image,'append')

    image_in=in_fo['MockImage_Noiseless'].data
    header_in=in_fo['MockImage_Noiseless'].header

    sigma_njy=(2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*sb_maglim))*header_in['PIXSIZE']*(3.0*header_in['FWHM'])
    
    npix=image_in.shape[0]
    
    noise_image = sigma_njy*np.random.randn(npix,npix)

    image_out=image_in + noise_image

    hdu_out = fits.ImageHDU(image_out,header=header_in)
    hdu_out.header['EXTNAME']='MockImage_SB'+sb_label
    hdu_out.header['SBLIM']=(sb_maglim,'mag/arcsec^2')
    hdu_out.header['RMSNOISE']=(sigma_njy,'nanojanskies')

    in_fo.append(hdu_out)
    in_fo.flush()
    
    
    return


