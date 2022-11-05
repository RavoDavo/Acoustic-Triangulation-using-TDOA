# Acoustic-Triangulation-using-TDOA
This repository holds the code and report to the EEE4022S final year project for Electrical and Computer Engineering at the University of Cape Town.

In the repository is a Python script called Record.py, this script sends a record command to the 4 Raspberry Pi's. 
The script then plays the calibration sound, and follows by collecting the recordings from the devices.

These 4 recordings are then placed into TDOA_Triangulation.m. This Matlab script triangulates the signal source using TDOA.
