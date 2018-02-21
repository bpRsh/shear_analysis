# This is the config file for jedisim.py and associated programs.
#
# Convension: There should be no empty lines in this file.
# key=value  or key="value"    There should be no spaces around equal sign.
#
#------------------------------ bulge disk factor ------------------------------
# Ref: https://www.lsst.org/about/camera/features
# Ref: http://www.stsci.edu/hst/acs/documents/handbooks/current/c05_imaging2.html
LSST_r_blue_wavelength=5520 # Angstrom (bd_factors_scaled.py and bd_weights_psf.py)
LSST_r_red_wavelength=6910 # Angstrom (bd_factors_scaled.py and bd_weights_psf.py)
HST_ACS_WFC_central_wavelength=8333 # Angstrom for Broad I WFC camera.
HST_ACS_WFC_width=2511 # Angstrom
z_cutout=0.2
#
#
#------------------------------ jedicatalog --------------------------------------
# jedicatalog will read MAG,MAG0,PIXSCALE, and RADIUS from the input fitsfiles and
# will create three catalogs:  catalog.txt,convolved.txt,distortedlist.txt
# To create these files jedicatalog needs some physics settngs quantities
# which are given below.
# These setting will also be other programs.
#----------------------------- physics settings --------------------------------
nx=12288            # x pixels            (arg for jedidistort)
ny=12288            # y pixels            (arg for jedidistort)
pix_scale=0.06      # arseconds per pixel (arg for jedidistort)
lens_z=0.3          # lens redshift       (arg for jedidistort)
single_redshift=1   # 0 = not-fixed, 1=fixed
fixed_redshift=0.7  # the single source galaxy redshift to use (bdfactors, bdweights)
final_pix_scale=0.2 # LSST pixscale (arcsecords per pixel)
exp_time=6000       # exposure time in seconds
noise_mean=10       # mean for poisson noise
x_border=301        # must be large enough so that no image can overflow_
y_border=301        # must be large enough so that no image can overflow
x_trim=480          # larger than x_border to ensure no edge effects
y_trim=480          # larger than y_border to ensure no edge effects
num_galaxies=12420  # number of galaxies to simulate 138,000 default
min_mag=22          # minimum magnitude galaxy to simulate (inclusive)
max_mag=28          # maximum magnitude galaxy to simulate (inclusive)
power=0.33          # power for the power law galaxy distribution
# For the f814 filter for our 201 galaxy images
#minmag   =  19.4715 # f814w_gal_19.fits
#maxmag   =  25.9455 # f814w_gal_214.fits
# For the HST galaxies cut outs we have magnitude range from 19-26.
# Note that we do not have the file radius database 19.txt
# inside simdatabase/radius_db/ so jedicatalog will fail if we take min_mag 19.
# However, we are not interested in simulation of mag 19.
# We are interested in fainter LSST galaxies with range 22-28 right now.
#---------------------------- Lens for jedidistort  ----------------------------
# lens.txt has a single line with 5 parameters
# 6144 6144 1 1000.000000 4.000000
#  x    y  type p1       p2
#  x - x center of lens (in pixels)
#  y - y center of lens (in pixels)
#  type - type of mass profile
#         1. Singular isothermal sphere
#         2. Navarro-Frenk-White profile
#	        3. NFW constant distortion profile for grid simulations
#  p1 - first profile parameter
#         1. sigma_v [km/s]
#         2. M200 parameter [10^14 solar masses]
#		      3. Distance to center in px. M200 fixed at 20 default,
#            which can be modified in case 3
#  p2 - second profile parameter
#         1. not applicable, can take any numerical
#         2. c parameter [unitless]
#         3. c parameter [unitless]
lens_x=6144
lens_y=6144
lens_type=1
lens_p1=1000.0
lens_p2=4.0
lens_file="physics_settings/lens.txt"	  # arg for jedidistort
#----------------------- Output Folders ----------------------------------------
# Note that for scaled_bulge, output folder is jedisim_out/out0/scaled_bulge
output_folder="jedisim_out/out0/scaled_bulge/"         # jedicatalog etc.
90_output_folder="jedisim_out/out90/scaled_bulge/"  # jedicatalog etc.
prefix="trial1_" # used by jedicatalog etc.
#----------------------- HST images from jediconvolve --------------------------
# Note that HST and LSST will be updated to trial1_HST.fits and
# trial1_LSST_averaged.fits and so on.
HST_image="HST.fits"                      # jedipaste, jediconvolve
HST_convolved_image="HST_convolved.fits"  # jedipaste2, jedirescale
#----------------------- database folders --------------------------------------
# There are 10 radius database files 20.dat to 29.dat.
# which contains min and max radius to be used by jedicatalog.
# e.g. the file simdatabase/radius_db/20.dat has two lines: 36.72 3.51
# This must contain files "n.txt" for n= min_mag to max_mag
radius_db_folder="simdatabase/radius_db/"
# There are 15+2 redshift database files 19.dat to 33.dat with +- 99.dat.
# which contains min and max redshift to be used by jedicatalog.
# e.g. the file simdatabase/red_db/19.dat has two lines: 0.301000 0.138000
# For f814 filter images:
# minmag   =  19.4715 # f814w_gal_19.fits
# maxmag   =  25.9455 # f814w_gal_214.fits
# minrad   =  1.31    # f814w_gal_43.fits
# maxrad   =  21.895  # f814w_gal_9.fits
# For f606 filter images:
# minmag =  19.3882 # f606w_gal120.fits
# maxmag =  25.5268 # f606w_gal215.fits
# minrad =  1.157   # f606w_gal93.fits
# maxrad =  22.055  # f606w_gal110.fits
red_db_folder="simdatabase/red_db/"
#------------------- bulge to disk weights -------------------------------------
# This file will be updated by jedisim_ofiles.py to:
# bd_weights="physics_settings/bd_weights_z1.5.txt" or so on.
bd_weights="physics_settings/bd_weights_z0.7.txt" # jedisim_config, jedisimulate
bd_factors="physics_settings/bd_factors.txt" # jedisim_config, jedisimulate
bd_flux_rat="physics_settings/bd_flux_rat_z0.7.txt" # fr, frb, frd written by scaled_bd_flux_rat.py
#------------------- rescaled outfiles for jedirescale -------------------------
rescaled_outfileb="lsst_bulge.fits" # jedisimulate.py  (gcsb galaxy convolved rescaled)
rescaled_outfiled="lsst_disk.fits" # jedisimulate.py
rescaled_outfilem="lsst_bulge_disk.fits" # jedisimulate.py
#--------------------- Final outputs of jedisim -------------------------------
lsst_unnoised="lsst_unnoised.fits" # jedisimulate.py
lsst="lsst.fits" # jedisimulate.py
lsst_mono="lsst_mono.fits" # jedisimulate.py
#------------------- psf files for jediconvolve --------------------------------
psfb="psf/psfb.fits"
psfd="psf/psfd.fits"
psfm="psf/psfm.fits"
#------------------- catalog files for jedicatalog -----------------------------
# jedicatalog will read four headers MAG, MAG0, PIXSCALE, RADIUS from
# the images given in the end of this config.sh file.
#
# jedicatalog will write the following quantity needed by jeditransform in
# the file jedisim_out/out0/trial1_catalog.txt:
# name,x,y,angle,redshift,pixscale,old_mag,old_rad,new_mag,new_rad,stamp_name,dis_name
#
# jedicatalog also creates the catalog list for the program jediconvolve:
# jedisim_out/out0/trial1_convolvedlist.txt
#
# This text files have 6 rows:
# jedisim_out/out1/convolved/convolved_band_0.fits
# jedisim_out/out1/convolved/convolved_band_5.fits
#
# jedicatalog also creates the catalog list for the program jedidistort:
# jedisim_out/out0/trial1_distortedlist.txt
#
# This text file 12420 rows.
# There are 13 folder names and each folder name has 1000 filenames except last.
# The contents of distortedlist.txt is like this:
# line     1: jedisim_out/out0/distorted_0/distorted_0.fits
# line 12420: jedisim_out/out0/distorted_12/distorted_12419.fits
#
# The code snippets used in jedicatalog is:
# else if(strcmp(buffer3,"catalog_file")==0) sscanf(buffer2, "catalog_file=\"%[^\"]", temp_catalog_file);
#
# Three text files created by jedicatalog are:
# jedisim_out/out0/scaled_bulge/trail1_catalog.txt
# jedisim_out/out0/scaled_bulge/trail1_convolvedlist.txt
# jedisim_out/out0/scaled_bulge/trail1_distortedlist.txt
catalog_file="catalog.txt"
convolvedlist_file="convolvedlist.txt"
distortedlist_file="distortedlist.txt"
#----------------- catalog files for jedidistort -------------------------------
# jeditransform creates jedisim_out/out0/scaled_bulge/trial1_dislist.txt along with 12420 .gz stamps
# jedidistort will use these files.
# Input galaxy parameter file for jedidistort: x y nx ny zs file
#           x - x coord. of lower left pixel where galaxy should be embedded
#           y - y coord. of lower left pixel where galaxy should be embedded
#           nx - width of the galaxy in pixels
#           ny - height of the galaxy in pixels
#           zs - redshift of this galaxy
#           infile - filepath to the FITS file for this galaxy, 1024 chars max
#           outfile - filepath for the output FITS file for this galaxy, 1024 chars max
#
#  e.g. jedisim_out/out1/dislist.txt looks like this: (was created by jeditransform)
# 6813 888 10 23 1.500000 jedisim_out/out1/stamp_0/stamp_.fits.gz out1/distorted_0/distorted_0.fits
# x    y   nx ny zs       infile                      outfile
dislist_file="dislist.txt"  # arg for jedidistort
convlist_file="toconvolvelist.txt"
#-----------------------  source images for jedicatalog  -----------------------
num_source_images=201
# Jedicatalog will read this number of galaxies.
# It also needs to read the four header keys MAG, MAG0, PIXSCALE, and RADIUS.
#
# The code snippets in jedicatalog is: if(strcmp(buffer3,"image")==0 && nimage < num_source_images)
# so it will read these 201 fitsfiles from this config.sh file.
#
# After reading these fitsfiles jedicatalog will create 3 catalog files.
# jedisim_out/out0/trial1_catalog.txt
# jedisim_out/out0/trial1_convolvedlist.txt
# jedisim_out/out0/distortedlist.txt
# Note that jedicatalog do not create any fitsfiles, it will only read
# four header keys from them, and creates three text files.
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge0.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge1.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge2.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge3.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge4.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge5.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge6.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge7.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge8.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge9.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge10.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge11.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge12.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge13.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge14.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge15.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge16.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge17.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge18.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge19.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge20.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge21.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge22.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge23.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge24.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge25.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge26.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge27.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge28.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge29.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge30.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge31.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge32.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge33.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge34.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge35.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge36.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge37.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge38.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge39.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge40.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge41.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge42.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge43.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge44.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge45.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge46.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge47.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge48.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge49.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge50.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge51.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge52.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge53.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge54.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge55.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge56.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge57.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge58.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge59.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge60.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge61.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge62.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge63.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge64.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge65.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge66.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge67.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge68.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge69.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge70.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge71.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge72.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge73.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge74.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge75.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge76.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge77.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge78.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge79.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge80.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge81.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge82.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge83.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge84.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge85.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge86.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge87.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge88.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge89.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge90.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge91.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge92.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge93.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge94.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge95.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge96.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge97.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge98.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge99.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge100.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge101.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge102.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge103.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge104.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge105.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge106.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge107.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge108.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge109.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge110.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge111.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge112.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge113.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge114.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge115.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge116.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge117.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge118.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge119.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge120.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge121.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge122.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge123.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge124.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge125.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge126.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge127.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge128.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge129.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge130.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge131.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge132.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge133.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge134.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge135.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge136.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge137.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge138.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge139.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge140.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge141.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge142.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge143.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge144.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge145.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge146.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge147.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge148.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge149.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge150.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge151.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge152.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge153.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge154.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge155.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge156.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge157.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge158.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge159.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge160.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge161.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge162.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge163.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge164.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge165.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge166.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge167.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge168.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge169.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge170.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge171.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge172.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge173.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge174.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge175.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge176.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge177.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge178.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge179.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge180.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge181.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge182.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge183.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge184.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge185.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge186.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge187.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge188.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge189.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge190.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge191.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge192.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge193.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge194.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge195.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge196.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge197.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge198.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge199.fits"
image="simdatabase/scaled_bulge_f8/f814w_scaled_bulge200.fits"
