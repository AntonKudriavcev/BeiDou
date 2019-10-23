# ./initSettings.m

# Functions initializes and saves settings. Settings can be edited inside of
# the function, updated from the command line or updated using a dedicated
# GUI - "setSettings".

# All settings are described inside function code.

# settings = initSettings()

#   Inputs: none

#   Outputs:
#       settings     - Receiver settings (a structure).
import datetime

import numpy as np


class Result(object):
    def __init__(self, settings):
        self._settings = settings
        self._results = None
        self._channels = None

    @property
    def settings(self):
        return self._settings

    @property
    def channels(self):
        assert isinstance(self._channels, np.recarray)
        return self._channels

    @property
    def results(self):
        assert isinstance(self._results, np.recarray)
        return self._results

    @results.setter
    def results(self, records):
        assert isinstance(records, np.recarray)
        self._results = records

    def plot(self):
        pass


class True_Position(object):
    def __init__(self):
        self._E = None
        self._N = None
        self._U = None

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, e):
        self._E = e

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, n):
        self._N = n

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, u):
        self._U = u


class Settings(object):
    def __init__(self):
        # Processing settings ====================================================
        # Number of milliseconds to be processed used 36000 + any transients (see
        # below - in Nav parameters) to ensure nav subframes are provided
        self.ms_To_Process = 37000.0

        # Number of channels to be used for signal processing
        self.number_of_Channels = 8

        # Move the starting point of processing. Can be used to start the signal
        # processing at any point in the data record (e.g. for long records). fseek
        # function is used to move the file read point, therefore advance is byte
        # based only.
        self.skip_Number_of_Bytes = 0

        # Raw signal file name and other parameter ===============================
        # This is a "default" name of the data file (signal record) to be used in
        # the post-processing mode
        self.file_Name = 'TEST1.DAT'

        # Data type used to store one sample
        self.data_Type = 'int8'

        # Intermediate, sampling and code frequencies
        self.IF = 1250000.0

        self.sampling_Freq = 5000000.0

        self.code_Freq_Basis = 1023000.0 * 2 ## BeiDou B1I

        # Define number of chips in a code period
        self.code_Length = 1023 * 2 ## BeiDou B1I

        # Acquisition settings ===================================================
        # Skips acquisition in the script postProcessing.m if set to 1
        self.skip_Acquisition = False

        # List of satellites to look for. Some satellites can be excluded to speed
        # up acquisition
        self.acq_Satellite_List = range(1, 38) ## 37 BeiDou satellites

        # Band around IF to search for satellite signal. Depends on max Doppler
        self.acq_Search_Band = 14.0

        # Threshold for the signal presence decision rule
        self.acq_Threshold = 2.5

        # Tracking loops settings ================================================
        # Code tracking loop parameters
        self.dll_Damping_Ratio = 0.7

        self.dll_Noise_Bandwidth = 2.0

        self.dll_Correlator_Spacing = 0.5

        # Carrier tracking loop parameters
        self.pll_Damping_Ratio = 0.7

        self.pll_Noise_Bandwidth = 25.0

        # Navigation solution settings ===========================================

        # Period for calculating pseudoranges and position
        self.nav_Sol_Period = 500.0

        # Elevation mask to exclude signals from satellites at low elevation
        self.elevation_Mask = 10.0

        # Enable/dissable use of tropospheric correction
        self.use_Trop_Corr = True

        # 1 - On

        # True position of the antenna in UTM system (if known). Otherwise enter
        # all NaN's and mean position will be used as a reference .
        self.true_Position = True_Position()
        #         self.truePosition.E = np.nan

        #         self.truePosition.N = np.nan

        #         self.truePosition.U = np.nan

        # Plot settings ==========================================================
        # Enable/disable plotting of the tracking results for each channel
        self.plot_Tracking = True

        # 1 - On

        # Constants ==============================================================

        self._c = 299792458.0

        self._start_Offset = 68.802

    @property
    def c(self):
        return self._c

    @property
    def start_Offset(self):
        return self._start_Offset

    @property
    def samples_Per_Code(self):
        return np.long(np.round(self.sampling_Freq / (self.code_Freq_Basis / self.code_Length)))

    # makeCaTable.m
    def make_Ranging_Code_Table(self):

        # caCodesTable = makeCaTable(settings)

        #   Inputs:
        #       settings        - receiver settings
        #   Outputs:
        #       caCodesTable    - an array of arrays (matrix) containing C/A codes
        #                       for all satellite PRN-s

        # --- Find number of samples per spreading code ----------------------------
        samples_Per_Code = self.samples_Per_Code

        # --- Prepare the output matrix to speed up function -----------------------
        Ranging_Code_Table = np.zeros((37, samples_Per_Code))

        # --- Find time constants --------------------------------------------------
        ts = 1.0 / self.sampling_Freq

        tc = 1.0 / self.code_Freq_Basis

        # === For all satellite PRN-s ...
        for PRN in range(37):
            # --- Generate Ranging code for given PRN -----------------------------------
            Ranging_Code = self.generate_Ranging_Code(PRN)

            # --- Make index array to read Ranging code values -------------------------

            code_Value_Index = np.ceil(ts * np.arange(1, samples_Per_Code + 1) / tc) - 1
            code_Value_Index = np.longlong(code_Value_Index)

            code_Value_Index[-1] = 2045  ## is equal to 2045

            # The "upsampled" code is made by selecting values form the CA code
            # chip array (caCode) for the time instances of each sample.
            Ranging_Code_Table[PRN] = Ranging_Code[code_Value_Index]
        return Ranging_Code_Table

    # generateCAcode.m
    def generate_Ranging_Code(self, prn):

        assert prn in range(0, 37)
    ## list for phase assignment of G2 sequence         
        g2s = np.array([(1, 3), (1, 4), (1, 5), (1, 6), (1, 8), (1, 9), (1, 10), (1, 11), 
                (2, 7), 
                (3, 4),  (3, 5),  (3, 6),  (3, 8),  (3, 9),  (3, 10), (3, 11),
                (4, 5),  (4, 6),  (4, 8),  (4, 9),  (4, 10), (4, 11), 
                (5, 6),  (5, 8),  (5, 9),  (5, 10), (5, 11), 
                (6, 8),  (6, 9),  (6, 10), (6, 11),
                (8, 9),  (8, 10), (8, 11),
                (9, 10), (9, 11),
                (10, 11)])

        g2s = g2s - 1 ## reduction of phase coefficients to a convenient form

        # --- Pick right shift for the given PRN number ----------------------------
        g2_shift = g2s[prn]

        # --- Generate G1 code -----------------------------------------------------

        # --- Initialize g1 output to speed up the function ---
        g1 = np.zeros(2046)

        # --- Load shift register ---
        reg = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1]) ## init g1 register as 01010101010

        # --- Generate all G1 signal chips based on the G1 feedback polynomial -----
        for i in range(2046):
            g1[i] = reg[-1]

            save_Bit = reg[0] * reg[6] * reg[7] * reg[8] * reg[9] * reg[10]

            reg[1:] = reg[:-1]

            reg[0] = save_Bit

        # --- Generate G2 code -----------------------------------------------------

        # --- Initialize g2 output to speed up the function ---
        g2 = np.zeros(2046)

        # --- Load shift register ---
        reg = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1]) ## init g2 register as 01010101010

        # --- Generate all G2 signal chips based on the G2 feedback polynomial -----
        for i in range(2046):
            g2[i] = reg[g2_shift[0]] * reg[g2_shift[1]] ## XOR for different phase coefficients

            save_Bit = reg[0] * reg[1] * reg[2] * reg[3] * reg[4] * reg[7] * reg[8] * reg[10]

            reg[1:] = reg[:-1]

            reg[0] = save_Bit


        # --- Form single sample Ranging code by multiplying G1 and G2 -----------------
        Ranging_Code = -g1 * g2
        return Ranging_Code

    @staticmethod
    # calcLoopCoef.m
    def calc_Loop_Coef(LBW, zeta, k):
        # Function finds loop coefficients. The coefficients are used then in PLL-s
        # and DLL-s.

        # [tau1, tau2] = calcLoopCoef(LBW, zeta, k)

        #   Inputs:
        #       LBW           - Loop noise bandwidth
        #       zeta          - Damping ratio
        #       k             - Loop gain

        #   Outputs:
        #       tau1, tau2   - Loop filter coefficients

        # Solve natural frequency
        Wn = LBW * 8.0 * zeta / (4.0 * zeta ** 2 + 1)

        # solve for t1 & t2
        tau1 = k / (Wn * Wn)

        tau2 = 2.0 * zeta / Wn

        return tau1, tau2

    def probe_Data(self, file_Name_Str = None):

        import matplotlib.pyplot as plt
        from scipy.signal import welch
        from scipy.signal.windows.windows import hamming

        # Function plots raw data information: time domain plot, a frequency domain
        # plot and a histogram.

        # The function can be called in two ways:
        #   probeData(settings)
        # or
        #   probeData(fileName, settings)

        #   Inputs:
        #       fileName        - name of the data file. File name is read from
        #                       settings if parameter fileName is not provided.

        #       settings        - receiver settings. Type of data file, sampling
        #                       frequency and the default filename are specified
        #                       here.

        # Check the number of arguments ==========================================
        if file_Name_Str is None:
            file_Name_Str = self.file_Name
        if not isinstance(file_Name_Str, str):
            raise TypeError('File name must be a string')
        settings = self
        # Generate plot of raw data ==============================================

        try:
            with open(file_Name_Str, 'rb') as fid:
                # Move the starting point of processing. Can be used to start the
                # signal processing at any point in the data record (e.g. for long
                # records).
                fid.seek(settings.skip_Number_of_Bytes, 0)
                samples_Per_Code = settings.samples_Per_Code

                try:
                    data = np.fromfile(fid,
                                       settings.data_Type,
                                       10 * samples_Per_Code)

                except IOError:
                    # The file is too short
                    print ('Could not read enough data from the data file.')
                # --- Initialization ---------------------------------------------------
                plt.figure(100)
                plt.clf()
                time_Scale = np.arange(0, 0.005, 1 / settings.sampling_Freq)

                plt.subplot(2, 1, 1)
                plt.plot(1000 * time_Scale[1:int(samples_Per_Code / 1)],
                         data[1:int(samples_Per_Code / 1)])
                plt.axis('tight')
                plt.grid()
                plt.title('Time domain plot')
                plt.xlabel('Time (ms)')
                plt.ylabel('Amplitude')
                plt.subplot(2, 1, 2)
                f, Pxx = welch(data - np.mean(data),
                               settings.sampling_Freq / 1000000.0,
                               hamming(16384, False),
                               16384,
                               1024,
                               16384)
                plt.semilogy(f, Pxx)
                plt.axis('tight')
                plt.grid()
                plt.title('Frequency domain plot')
                plt.xlabel('Frequency (MHz)')
                plt.ylabel('Magnitude')
                plt.subplots_adjust(wspace = 0, hspace = 0.5)
                plt.show()
                plt.subplot(1, 1, 1)
                plt.hist(data, np.arange(- 128, 128))
                dmax = np.max(np.abs(data)) + 1

                plt.axis('tight')
                adata = plt.axis()

                plt.axis([-dmax, dmax, adata[2], adata[3]])
                plt.grid('on')
                plt.title('Histogram')
                plt.xlabel('Bin')
                plt.ylabel('Number in bin')
                plt.show()
            # === Error while opening the data file ================================
        except IOError as e:
            print ('Unable to read file "%s": %s' % (file_Name_Str, e))

    # ./postProcessing.m

    # Script postProcessing.m processes the raw signal from the specified data
    # file (in settings) operating on blocks of 37 seconds of data.

    # First it runs acquisition code identifying the satellites in the file,
    # then the code and carrier for each of the satellites are tracked, storing
    # the 1m sec accumulations.  After processing all satellites in the 37 sec
    # data block, then postNavigation is called. It calculates pseudoranges
    # and attempts a position solutions. At the end plots are made for that
    # block of data.

    #                         THE SCRIPT "RECIPE"

    # The purpose of this script is to combine all parts of the software
    # receiver.

    # 1.1) Open the data file for the processing and seek to desired point.

    # 2.1) Acquire satellites

    # 3.1) Initialize channels (preRun.m).
    # 3.2) Pass the channel structure and the file identifier to the tracking
    # function. It will read and process the data. The tracking results are
    # stored in the trackResults structure. The results can be accessed this
    # way (the results are stored each millisecond):
    # trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
    # XXX is a field name of the result (e.g. I_P, codePhase etc.)

    # 4) Pass tracking results to the navigation solution function. It will
    # decode navigation messages, find satellite positions, measure
    # pseudoranges and find receiver position.

    # 5) Plot the results.

    def post_Processing(self, file_Name_Str = None):
        # Initialization =========================================================
        import acquisition
        # import postNavigation
        # import tracking
        print ('Starting processing...')
        settings = self
        if not file_Name_Str:
            file_Name_Str = settings.file_Name
        if not isinstance(file_Name_Str, str):
            raise TypeError('File name must be a string')
        try:
            with open(file_Name_Str, 'rb') as fid:

                # If success, then process the data
                # Move the starting point of processing. Can be used to start the
                # signal processing at any point in the data record (e.g. good for long
                # records or for signal processing in blocks).
                fid.seek(settings.skip_Number_of_Bytes, 0)
                # Acquisition ============================================================
                # Do acquisition if it is not disabled in settings or if the variable
                # acqResults does not exist.
                if not settings.skip_Acquisition:  # or 'acqResults' not in globals():
                    # Find number of samples per spreading code
                    samples_Per_Code = settings.samples_Per_Code

                    # frequency estimation
                    data = np.fromfile(fid, settings.data_Type, 11 * samples_Per_Code)

                    print ('   Acquiring satellites...')
                    acq_Results = acquisition.Acquisition_Result(settings)
                    acq_Results.acquire(data)
                    # acqResults.plot()
                # Initialize channels and prepare for the run ============================
                # Start further processing only if a GNSS signal was acquired (the
                # field FREQUENCY will be set to 0 for all not acquired signals)
                if np.any(acq_Results.carr_Freq):
                    acq_Results.preRun()
                    acq_Results.show_Channel_Status()
                else:
                    # No satellites to track, exit
                    print ('No GNSS signals detected, signal processing finished.')
                    track_Results = None

                # # Track the signal =======================================================
                # startTime = datetime.datetime.now()

                # print ('   Tracking started at %s' % startTime.strftime('%X'))
                # trackResults = tracking.TrackingResult(acqResults)
                # try:
                #     trackResults.results = np.load('trackingResults_python.npy')
                # except IOError:
                #     trackResults.track(fid)
                #     np.save('trackingResults_python', trackResults.results)

                # print ('   Tracking is over (elapsed time %s s)' % (datetime.datetime.now() - startTime).total_seconds())
                # # Auto save the acquisition & tracking results to save time.
                # print ('   Saving Acquisition & Tracking results to storage')
                # # Calculate navigation solutions =========================================
                # print ('   Calculating navigation solutions...')
                # navResults = postNavigation.NavigationResult(trackResults)
                # navResults.postNavigate()

                # print ('   Processing is complete for this data block')
                # # Plot all results ===================================================
                # print ('   Plotting results...')
                # # TODO turn off tracking plots for now
                # if not settings.plotTracking:
                #     trackResults.plot()
                # navResults.plot()
                print ('Post processing of the signal is over.')
        except IOError as e:
            # Error while opening the data file.
            print ('Unable to read file "%s": %s.' % (settings.file_Name, e))
