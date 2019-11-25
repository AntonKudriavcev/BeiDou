

import numpy as np
import matplotlib.pyplot as plt

class Result(object):
    def __init__(self, settings):
        self._settings = settings
        self._results  = None
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


class Settings(object):
    def __init__(self):

        self.number_of_Channels = 8

        self.code_Freq_Basis = 2.046e6 ## BeiDou B1I
        self.code_Length     = 2046 ## BeiDou B1I 

        self.sampling_Freq = 10e6
        self.IF            = 1.098e6
        self.interleaved   = 1
        self.data_Type     = 'int8'
        self.file_Name     = 'D://study//5_year//Course_project//SoftGNSS_python//data.dat'

        self.ms_To_Process        = 7000.0
        self.skip_Number_of_Bytes = 4660


        # Acquisition settings ===================================================

        self.skip_Acquisition   = False
        self.acq_Satellite_List = range(1, 38) ## 37 BeiDou satellites
        self.acq_Search_Band    = 10e3 ## Hz
        self.acq_dopp_step      = 600
        self.acq_Threshold      = 500

    @property
    def samples_Per_Code(self):
        return np.long(np.round(self.sampling_Freq / (self.code_Freq_Basis / self.code_Length)))

    @property
    def code_period(self):
        return  self.code_Length/self.code_Freq_Basis

    def make_Ranging_Code_Table(self):

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

            code_Value_Index = np.floor(ts * np.arange(0, samples_Per_Code) / tc) 
            code_Value_Index = np.longlong(code_Value_Index)

            code_Value_Index[-1] = 2045  ## is equal to 2045

            Ranging_Code_Table[PRN] = Ranging_Code[code_Value_Index]

        return Ranging_Code_Table

    def generate_Ranging_Code(self, prn):

        assert prn in range(0, 37)
         
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

            g2[i]    = reg[g2_shift[0]] * reg[g2_shift[1]] ## XOR for different phase coefficients
            save_Bit = reg[0] * reg[1] * reg[2] * reg[3] * reg[4] * reg[7] * reg[8] * reg[10]
            reg[1:]  = reg[:-1]
            reg[0]   = save_Bit

        # --- Form single sample Ranging code by multiplying G1 and G2 -----------------
        Ranging_Code = g1 * g2
        return Ranging_Code


    def probe_Data(self, file_Name_Str = None):
        pass


    def post_Processing(self, file_Name_Str = None):
        import acquisition

        print ('Starting processing...')
        settings = self
        if not file_Name_Str:
            file_Name_Str = settings.file_Name
        if not isinstance(file_Name_Str, str):
            raise TypeError('File name must be a string')
        try:
            with open(file_Name_Str, 'rb') as fid:
                # Acquisition ============================================================
                    # Find number of samples per spreading code
                samples_Per_Code     = settings.samples_Per_Code
                skip_Number_of_Bytes = settings.skip_Number_of_Bytes
                interleaved          = settings.interleaved

                # frequency estimation
                data = np.fromfile(fid, 
                                   settings.data_Type, 
                                   count = (interleaved + 1) * (11 * samples_Per_Code + skip_Number_of_Bytes))
                data = data[2 * skip_Number_of_Bytes:]

                print ('   Acquiring satellites...')
                acq_Results = acquisition.Acquisition_Result(settings)
                acq_Results.acquire(data)
                # acqResults.plot()
                # Initialize channels and prepare for the run ============================

                if np.any(acq_Results.carr_Freq):
                    acq_Results.preRun()
                    acq_Results.show_Channel_Status()
                else:
                    # No satellites to track, exit
                    print ('No GNSS signals detected, signal processing finished.')
                    track_Results = None

                print ('Post processing of the signal is over.')
        except IOError as e:
            # Error while opening the data file.
            print ('Unable to read file "%s": %s.' % (settings.file_Name, e))
