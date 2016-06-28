class IdData:
    '''A class defined to store data for system identification

    This class stores the values of signals u (the inputs) and y (the outputs)
    as slices of data
    '''

    def __init__(self, u, y, Ts=1):
        '''Constructor of the iddata object

        The default constructor is iddata(u, y, Ts), u is the data for system
        inputs (N x nu) and y is the data for system outputs (N x ny) and Ts is
        the sampling time.

        N is the number of samples, nu is the number of inputs and ny is the
        number of outputs
        '''

        self.data_slice = [(u, y, Ts), ]; #slices stored in a list
        assert u.shape[0] == y.shape[0]

        self.input_count = u.shape[1]
        self.output_count = y.shape[1]
        self.Ts = Ts
        if Ts <= 0:
            raise ValueError('Value of sampling time Ts should be > 0')


    def merge(self, *args):
        '''Merge multiple IdData sets into the current data set

        merge(dat1, dat2, dat3, ... datN)

        The following conditions must be met:
            1) The number of inputs in all of datk should be the same
            2) The number of outputs in all of datk should be the same
            3) The sampling times of all datk should be the same

        Returns
        -------
        The updated IdData
        '''

        ndat = len(args)
        datout = self

        for idat in range(0,ndat):
            if args[idat].input_count == datout.input_count and args[idat].output_count == args[0].output_count and args[idat].Ts == datout.Ts:
                   datout.data_slice.append(*args[idat].data_slice)
            else:
                raise IdException('Input/Output count mismatch during merge')

        return datout

class IdException(Exception):
   def __init__(self, value):
      self.value = value
   def __str__(self):
     return repr(self.value)
