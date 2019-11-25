
import numpy as np
import matplotlib.pyplot as plt

from initialize import Result
from graph_builder import Builder

class Acquisition_Result(Result):
    def __init__(self, settings):
        self._settings = settings
        self._results  = None
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

        # Initialization =========================================================
        settings = self._settings

        

        # Find number of samples per spreading code
        samples_Per_Code = settings.samples_Per_Code
        interleaved      = settings.interleaved

        # Create two 1m sec vectors of data to correlate with and one with zero DC
        signal1 = long_Signal[0 : (interleaved + 1) * samples_Per_Code]

        signal2 = long_Signal[(interleaved + 1) * samples_Per_Code : 2 * (interleaved + 1) * samples_Per_Code]

        # Find sampling period
        ts = 1.0 / settings.sampling_Freq

        # Find phase points of the local carrier wave
        phase_Points = np.arange(samples_Per_Code) * 2 * np.pi * ts

        # Number of the frequency bins for the given acquisition band (500Hz steps)
        number_of_Frq_Bins = np.int(np.floor(settings.acq_Search_Band / settings.acq_dopp_step) + 1)
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

        ## create exemplar plot graph


        print ('(')
        # Perform search for all listed PRN numbers ...

        for PRN in range(len(settings.acq_Satellite_List)):

            plot_builder = Builder(num_of_satellites = PRN + 1,
                                   sampl_freq        = settings.sampling_Freq, 
                                   samples_per_code  = settings.samples_Per_Code) 

            # Correlate signals ======================================================
            # --- Perform DFT of Ranging Code ----------------------------------------
            Ranging_Code_Freq_Dom = np.fft.fft(Ranging_Code_Table[PRN, :]).conj()

            for frq_Bin_Index in range(number_of_Frq_Bins):
                # --- Generate carrier wave frequency grid (0.5kHz step) -----------
                frq_Bins[frq_Bin_Index] = settings.IF - \
                                       settings.acq_Search_Band / 2 + \
                                       settings.acq_dopp_step * frq_Bin_Index

                IQ1 = signal1[0::2] + 1j * signal1[1::2]
                IQ2 = signal2[0::2] + 1j * signal2[1::2]

                IQ_freq_Dom1 = np.fft.fft(IQ1 * np.exp(-1j * phase_Points * (frq_Bins[frq_Bin_Index])))
                IQ_freq_Dom2 = np.fft.fft(IQ2 * np.exp(-1j * phase_Points * (frq_Bins[frq_Bin_Index])))

                # domain)
                conv_Code_IQ1 = IQ_freq_Dom1 * Ranging_Code_Freq_Dom
                conv_Code_IQ2 = IQ_freq_Dom2 * Ranging_Code_Freq_Dom

                acq_Res1 = abs(np.fft.ifft(conv_Code_IQ1))
                acq_Res2 = abs(np.fft.ifft(conv_Code_IQ2))

                # "blend" 1st and 2nd msec but will correct data bit issues
                if acq_Res1.max() > acq_Res2.max():
                    results[frq_Bin_Index, :] = acq_Res1

                else:
                    results[frq_Bin_Index, :] = acq_Res2

                plot_builder.add_to_plot(data           = results[frq_Bin_Index, :], 
                                         freq_deviation = frq_Bins[frq_Bin_Index] - settings.IF)

            plot_builder.show_plot()

            # Look for correlation peaks in the results ==============================
            # Find the highest peak and compare it to the second highest peak
            # The second peak is chosen not closer than 1 chip to the highest peak
            # --- Find the correlation peak and the carrier frequency --------------
            peak_Size           = results.max(1).max()
            frequency_Bin_Index = results.max(1).argmax()

            peak_Size  = results.max(0).max()
            code_Phase = results.max(0).argmax()

            # samples_Per_Code_Chip = np.long(round(settings.sampling_Freq / settings.code_Freq_Basis))

            # exclude_Range_Index1  = code_Phase - samples_Per_Code_Chip

            # exclude_Range_Index2  = code_Phase + samples_Per_Code_Chip

            # # boundaries
            # if exclude_Range_Index1 <= 0:
            #     code_Phase_Range = np.r_[exclude_Range_Index2 : samples_Per_Code + exclude_Range_Index1 + 1]

            # elif exclude_Range_Index2 >= samples_Per_Code - 1:
            #     code_Phase_Range = np.r_[exclude_Range_Index2 - samples_Per_Code : exclude_Range_Index1]

            # else:
            #     code_Phase_Range = np.r_[0 : exclude_Range_Index1 + 1, exclude_Range_Index2 : samples_Per_Code]

            # # --- Find the second highest correlation peak in the same freq. bin ---
            # second_Peak_Size = results[frequency_Bin_Index, code_Phase_Range].max()

            # peak_Metric[PRN] = peak_Size / second_Peak_Size

            if (peak_Size) > settings.acq_Threshold:
                # Fine resolution frequency search =======================================
                # --- Indicate PRN number of the detected signal -------------------
                print ('%02d ' % (PRN + 1))
                # Ranging_Code = settings.generate_Ranging_Code(PRN)

                # code_Value_Index = np.floor(ts * np.arange(1, 10 * samples_Per_Code + 1) / (1.0 / settings.code_Freq_Basis))

                # long_Ranging_Code = Ranging_Code[np.longlong(code_Value_Index % 2046)]

                # # (Using detected C/A code phase)
                # xCarrier = signal_0DC[code_Phase : code_Phase + 10 * samples_Per_Code] * long_Ranging_Code

                # fft_Num_Pts = 8 * 2 ** (np.ceil(np.log2(len(xCarrier))))

                # # associated carrier frequency
                # fftxc = np.abs(np.fft.fft(xCarrier, np.long(fft_Num_Pts)))

                # uniq_Fft_Pts = np.long(np.ceil((fft_Num_Pts + 1) / 2.0))

                # fft_Max = fftxc[4:uniq_Fft_Pts - 5].max()
                # fft_Max_Index = fftxc[4:uniq_Fft_Pts - 5].argmax()

                # fft_Freq_Bins = np.arange(uniq_Fft_Pts) * settings.sampling_Freq / fft_Num_Pts

                # carr_Freq[PRN] = fft_Freq_Bins[fft_Max_Index]

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
