import numpy as np

from initialize import Result


class Acquisition_Result(Result):
    def __init__(self, settings):
        self._settings = settings
        self._results = None
        self._channels = None

    @property
    def peak_Metric(self):
        assert isinstance(self._results, np.recarray)
        return self._results.peak_Metric

    @property
    def carr_Freq(self):
        assert isinstance(self._results, np.recarray)
        return self._results.carr_Freq

    @property
    def code_Phase(self):
        assert isinstance(self._results, np.recarray)
        return self._results.code_Phase

    def acquire(self, long_Signal):
        # ./acquisition.m
        # Function performs cold start acquisition on the collected "data". It
        # searches for GPS signals of all satellites, which are listed in field
        # "acqSatelliteList" in the settings structure. Function saves code phase
        # and frequency of the detected signals in the "acqResults" structure.

        # acqResults = acquisition(longSignal, settings)

        #   Inputs:
        #       longSignal    - 11 ms of raw signal from the front-end
        #       settings      - Receiver settings. Provides information about
        #                       sampling and intermediate frequencies and other
        #                       parameters including the list of the satellites to
        #                       be acquired.
        #   Outputs:
        #       acqResults    - Function saves code phases and frequencies of the
        #                       detected signals in the "acqResults" structure. The
        #                       field "carrFreq" is set to 0 if the signal is not
        #                       detected for the given PRN number.

        # Initialization =========================================================
        settings = self._settings

        # Find number of samples per spreading code
        samples_Per_Code = settings.samples_Per_Code

        # Create two 1m sec vectors of data to correlate with and one with zero DC
        signal1 = long_Signal[0:samples_Per_Code]

        signal2 = long_Signal[samples_Per_Code:2 * samples_Per_Code]

        signal_0DC = long_Signal - long_Signal.mean()

        # Find sampling period
        ts = 1.0 / settings.sampling_Freq

        # Find phase points of the local carrier wave
        phase_Points = np.arange(samples_Per_Code) * 2 * np.pi * ts

        # Number of the frequency bins for the given acquisition band (500Hz steps)
        number_of_Frq_Bins = np.int(np.round(settings.acq_Search_Band * 2) + 1)
        print(number_of_Frq_Bins)

        # Generate all C/A codes and sample them according to the sampling freq.
        Ranging_Code_Table = settings.make_Ranging_Code_Table()

        # --- Initialize arrays to speed up the code -------------------------------
        # Search results of all frequency bins and code shifts (for one satellite)
        results = np.zeros((number_of_Frq_Bins, samples_Per_Code))

        # Carrier frequencies of the frequency bins
        frq_Bins = np.zeros(number_of_Frq_Bins)

        # --- Initialize acqResults ------------------------------------------------
        # Carrier frequencies of detected signals
        carr_Freq = np.zeros(37)

        # C/A code phases of detected signals
        code_Phase_ = np.zeros(37)

        # Correlation peak ratios of the detected signals
        peak_Metric = np.zeros(37)

        print ('(')
        # Perform search for all listed PRN numbers ...
        for PRN in range(len(settings.acq_Satellite_List)):
            # Correlate signals ======================================================
            # --- Perform DFT of C/A code ------------------------------------------
            Ranging_Code_Freq_Dom = np.fft.fft(Ranging_Code_Table[PRN, :]).conj()

            for frq_Bin_Index in range(number_of_Frq_Bins):
                # --- Generate carrier wave frequency grid (0.5kHz step) -----------
                frq_Bins[frq_Bin_Index] = settings.IF - \
                                       settings.acq_Search_Band / 2 * 1000 + \
                                       500.0 * frq_Bin_Index

                sin_Carr = np.sin(frq_Bins[frq_Bin_Index] * phase_Points)

                cos_Carr = np.cos(frq_Bins[frq_Bin_Index] * phase_Points)

                I1 = sin_Carr * signal1

                Q1 = cos_Carr * signal1

                I2 = sin_Carr * signal2

                Q2 = cos_Carr * signal2

                IQ_freq_Dom1 = np.fft.fft(I1 + 1j * Q1)

                IQ_freq_Dom2 = np.fft.fft(I2 + 1j * Q2)

                # domain)
                conv_Code_IQ1 = IQ_freq_Dom1 * Ranging_Code_Freq_Dom

                conv_Code_IQ2 = IQ_freq_Dom2 * Ranging_Code_Freq_Dom

                acq_Res1 = abs(np.fft.ifft(conv_Code_IQ1)) ** 2

                acq_Res2 = abs(np.fft.ifft(conv_Code_IQ2)) ** 2

                # "blend" 1st and 2nd msec but will correct data bit issues
                if acq_Res1.max() > acq_Res2.max():
                    results[frq_Bin_Index, :] = acq_Res1

                else:
                    results[frq_Bin_Index, :] = acq_Res2

            # Look for correlation peaks in the results ==============================
            # Find the highest peak and compare it to the second highest peak
            # The second peak is chosen not closer than 1 chip to the highest peak
            # --- Find the correlation peak and the carrier frequency --------------
            peak_Size = results.max(1).max()
            frequency_Bin_Index = results.max(1).argmax()

            peak_Size = results.max(0).max()
            code_Phase = results.max(0).argmax()

            samples_Per_Code_Chip = np.long(round(settings.sampling_Freq / settings.code_Freq_Basis))

            exclude_Range_Index1 = code_Phase - samples_Per_Code_Chip

            exclude_Range_Index2 = code_Phase + samples_Per_Code_Chip

            # boundaries
            if exclude_Range_Index1 <= 0:
                code_Phase_Range = np.r_[exclude_Range_Index2 : samples_Per_Code + exclude_Range_Index1 + 1]

            elif exclude_Range_Index2 >= samples_Per_Code - 1:
                code_Phase_Range = np.r_[exclude_Range_Index2 - samples_Per_Code : exclude_Range_Index1]

            else:
                code_Phase_Range = np.r_[0 : exclude_Range_Index1 + 1, exclude_Range_Index2 : samples_Per_Code]

            # --- Find the second highest correlation peak in the same freq. bin ---
            second_Peak_Size = results[frequency_Bin_Index, code_Phase_Range].max()

            peak_Metric[PRN] = peak_Size / second_Peak_Size

            if (peak_Size / second_Peak_Size) > settings.acq_Threshold:
                # Fine resolution frequency search =======================================
                # --- Indicate PRN number of the detected signal -------------------
                print ('%02d ' % (PRN + 1))
                Ranging_Code = settings.generate_Ranging_Code(PRN)

                code_Value_Index = np.floor(ts * np.arange(1, 10 * samples_Per_Code + 1) / (1.0 / settings.code_Freq_Basis))

                long_Ranging_Code = Ranging_Code[np.longlong(code_Value_Index % 2046)]

                # (Using detected C/A code phase)
                xCarrier = signal_0DC[code_Phase : code_Phase + 10 * samples_Per_Code] * long_Ranging_Code

                fft_Num_Pts = 8 * 2 ** (np.ceil(np.log2(len(xCarrier))))

                # associated carrier frequency
                fftxc = np.abs(np.fft.fft(xCarrier, np.long(fft_Num_Pts)))

                uniq_Fft_Pts = np.long(np.ceil((fft_Num_Pts + 1) / 2.0))

                fft_Max = fftxc[4:uniq_Fft_Pts - 5].max()
                fft_Max_Index = fftxc[4:uniq_Fft_Pts - 5].argmax()

                fft_Freq_Bins = np.arange(uniq_Fft_Pts) * settings.sampling_Freq / fft_Num_Pts

                carr_Freq[PRN] = fft_Freq_Bins[fft_Max_Index]

                code_Phase_[PRN] = code_Phase

            else:
                # --- No signal with this PRN --------------------------------------
                print ('. ')

        # === Acquisition is over ==================================================
        print (')\n')
        acq_Results = np.core.records.fromarrays([carr_Freq, code_Phase_, peak_Metric],
                                                names='carr_Freq,code_Phase,peak_Metric')
        self._results = acq_Results
        return

    def plot(self):
        assert isinstance(self._results, np.recarray)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        # from scipy.io.matlab import loadmat

        # %% configure matplotlib
        mpl.rcdefaults()
        # mpl.rcParams['font.sans-serif']
        # mpl.rcParams['font.family'] = 'serif'
        mpl.rc('savefig',  bbox      = 'tight', transparent     = False, format    = 'png'          )
        mpl.rc('axes',     grid      = True,    linewidth       = 1.5,   axisbelow = True           )
        mpl.rc('lines',    linewidth = 1.5,     solid_joinstyle = 'bevel'                           )
        mpl.rc('figure',   figsize   = [8, 6],  autolayout      = False, dpi = 120                  )
        mpl.rc('text',     usetex    = True                                                         )
        mpl.rc('font',     family    = 'serif', serif           = 'Computer Modern Roman', size = 16)
        mpl.rc('mathtext', fontset   = 'cm'                                                         )

        # mpl.rc('font', size=16)
        # mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

        # ./plotAcquisition.m
        # Functions plots bar plot of acquisition results (acquisition metrics). No
        # bars are shown for the satellites not included in the acquisition list (in
        # structure SETTINGS).

        # plotAcquisition(acqResults)

        #   Inputs:
        #       acqResults    - Acquisition results from function acquisition.

        # Plot all results =======================================================
        f, hAxes = plt.subplots()

        plt.bar(range(1, 38), self.peak_Metric)
        plt.title('Acquisition results')
        plt.xlabel('PRN number (no bar - SV is not in the acquisition list)')
        plt.ylabel('Acquisition Metric ($1^{st}$ to $2^{nd}$ Correlation Peaks Ratio')
        old_Axis = plt.axis()

        plt.axis([0, 38, 0, old_Axis[-1]])
        plt.xticks(range(1, 38), size = 12)
        # plt.minorticks_on()
        hAxes.xaxis.grid()
        # Mark acquired signals ==================================================

        acquired_Signals = self.peak_Metric * (self.carr_Freq > 0)

        plt.bar(range(1, 38), acquired_Signals, FaceColor = (0, 0.8, 0))
        plt.legend(['Not acquired signals', 'Acquired signals'])
        plt.show()

    # preRun.m
    def pre_Run(self):
        assert isinstance(self._results, np.recarray)
        # Function initializes tracking channels from acquisition data. The acquired
        # signals are sorted according to the signal strength. This function can be
        # modified to use other satellite selection algorithms or to introduce
        # acquired signal properties offsets for testing purposes.

        # [channel] = preRun(acqResults, settings)

        #   Inputs:
        #       acqResults  - results from acquisition.
        #       settings    - receiver settings

        #   Outputs:
        #       channel     - structure contains information for each channel (like
        #                   properties of the tracked signal, channel status etc.).

        settings = self._settings
        # Initialize all channels ================================================
        PRN = np.zeros(settings.number_of_Channels, dtype = 'int64')
        acquired_Freq = np.zeros(settings.number_of_Channels)
        code_Phase = np.zeros(settings.number_of_Channels)
        status = ['-' for _ in range(settings.number_of_Channels)]

        # --- Copy initial data to all channels ------------------------------------

        # Copy acquisition results ===============================================

        # --- Sort peaks to find strongest signals, keep the peak index information
        PRN_indexes = sorted(enumerate(self.peak_Metric),
                            key = lambda x: x[-1], reverse = True)

        # --- Load information about each satellite --------------------------------
        # Maximum number of initialized channels is number of detected signals, but
        # not more as the number of channels specified in the settings.
        for ii in range(min(settings.number_of_Channels, sum(self.carr_Freq > 0))):
            PRN[ii] = PRN_indexes[ii][0] + 1

            acquired_Freq[ii] = self.carr_Freq[PRN_indexes[ii][0]]

            code_Phase[ii] = self.code_Phase[PRN_indexes[ii][0]]

            status[ii] = 'T'

        channel = np.core.records.fromarrays([PRN, acquired_Freq, code_Phase, status],
                                             names='PRN,acquired_Freq,code_Phase,status')
        self._channels = channel
        return

    def show_Channel_Status(self):
        # Prints the status of all channels in a table.

        # showChannelStatus(channel, settings)

        #   Inputs:
        #       channel     - data for each channel. It is used to initialize and
        #                   at the processing of the signal (tracking part).
        #       settings    - receiver settings

        channel  = self._channels
        settings = self._settings
        assert isinstance(channel, np.recarray)
        print ('\n*=========*=====*===============*===========*=============*========*')
        print ('| Channel | PRN |   Frequency   |  Doppler  | Code Offset | Status |')
        print ('*=========*=====*===============*===========*=============*========*')
        for channel_Nr in range(settings.number_of_Channels):
            if channel[channel_Nr].status != '-':
                print ('|      %2d | %3d |  %2.5e |   %5.0f   |    %6d   |     %1s  |' % (
                                    channel_Nr,
                                    channel[channel_Nr].PRN,
                                    channel[channel_Nr].acquired_Freq,
                                    channel[channel_Nr].acquired_Freq - settings.IF,
                                    channel[channel_Nr].code_Phase,
                                    channel[channel_Nr].status))
            else:
                print ('|      %2d | --- |  ------------ |   -----   |    ------   |   Off  |' % channel_Nr)

        print ('*=========*=====*===============*===========*=============*========*\n')


if __name__ == '__main__':
    pass
