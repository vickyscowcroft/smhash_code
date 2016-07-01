# smhash_code
Python calibration scripts for SMHASH

See calibration_procedure.pdf for full documentation

- calculate_aperture_correction.py 
  - Calculate the aperture correction for a single epoch. Requires a .raw file with large aperture photometry first, then PSF photometry
- apply_aperture_correction.py
  - Run this after you've calculated the aperture correction with the previous script. Also applys the zero point correction and corrects to the Spitzer standard aperture.
- calculate_each_epoch_offset.py
  - Calculates the offset between epoch 1 and all other epochs.
  -Needs the .raw file and a .mtr file for each epoch.
- location_correction.py
  - Calculates and applies the location correction to each photometry file.
  - Must be run last.
  - Need to run once for each image
  - Give it the phot file name and correction image name as input. 


Other Useful scripts

- remove_bad_stars.py
  - Removes bad stars from a .als file (i.e. those with *** in the records, or the wrong number of columns)
  - WARNING: This will overwrite the original .als file
  - Need to update this so that it takes an optional second argument with a new filename for the output
  
- reddening_laws.py
  - Optical reddening laws from Cardelli, Clayton & Mathis 
  - IR reddening laws from Indebetouw

- reddening_plt_smc.py
  - Example of fitting the reddening laws to multiwavelength Cepheid data
  - This uses data from the SMC.
  - You'll need to change paths, wavelengths etc to make it fit your data.
  - email me for help