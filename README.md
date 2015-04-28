# smhash_code
Python calibration scripts for SMHASH

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
