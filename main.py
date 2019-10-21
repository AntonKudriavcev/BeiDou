import initialize

# ./init.m

# Clean up the environment first =========================================
# clear
# close_('all')
# clc
# format('compact')
# format('long','g')
# --- Include folders with functions ---------------------------------------
# addpath('include')
# addpath('geoFunctions')
# Print startup ==========================================================
print ('Welcome to:  softGNSS\n')
# Initialize settings class=========================================
settings = initialize.Settings()

# Generate plot of raw data and ask if ready to start processing =========
try:
    print ('Probing data "%s"...' % settings.fileName)
    settings.probeData('TEST1.DAT')
finally:
    pass

print ('  Raw IF data plotted ')
print (' (run setSettings or change settings in "initialize.py" to reconfigure)')
print (' ')
gnssStart = True
# gnssStart = int(raw_input('Enter "1" to initiate GNSS processing or "0" to exit : ').strip())

if gnssStart:
    print (' ')
    settings.postProcessing()
