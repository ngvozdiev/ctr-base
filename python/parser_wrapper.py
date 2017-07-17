from ctypes import *
from sys import platform as _platform
import metrics_pb2
import matplotlib.pylab as plt
import numpy as np
import sys

c_double_p = POINTER(c_double)
c_ulonglong_p = POINTER(c_ulonglong)
c_char_p = POINTER(c_char)
c_byte_p = POINTER(c_byte)

def DictExtract(d, key, default=None):
    if key in d:
        return d.pop(key)
    return default

def ProcessDistributionBlob(metric_string, field_string, byte_buffer):
    output = metrics_pb2.PBDistribution()
    metrics_pb2.PBDistribution.ParseFromString(output, byte_buffer)
    return output

def Bin(x, y, bin_size=1):
    # Bins the x and y values. Each bin will have bin_size
    # values. Each y bin will be the average of the original y values,
    # each x bin will be the start of the bin.
    assert(len(x) == len(y))
    if bin_size == 1:
        return (x, y)

    if len(x) % bin_size != 0:
        padding = bin_size - (len(x) % bin_size)
        # Will pad with 0s
        x = np.pad(x, (0, padding), 'constant')
        y = np.pad(y, (0, padding), 'constant')
    
    new_y = np.mean(y.reshape(-1, bin_size), axis=1)
    new_x = x.reshape(-1, bin_size)[:,0]
    return new_x, new_y

class MetricsParser(object):
    PLOT_LINE = 1
    PLOT_SCATTER = 2
    PLOT_STEP = 3

    def __init__(self, metrics_file, sofile):
        self.lib = cdll.LoadLibrary(sofile)

        self.lib.MetricsParserManifestSummary.restype = c_char_p
        self.lib.MetricsParserManifestSummary.argtypes = [c_char_p]
        
        self.lib.MetricsParserParse.restype = c_void_p
        self.lib.MetricsParserParse.argtypes = [c_char_p, c_char_p,
                                                c_char_p, c_ulonglong, c_ulonglong, c_ulonglong]

        self.lib.MetricsParserResultHandleAdvance.restype = c_bool
        self.lib.MetricsParserResultHandleAdvance.argtypes = [c_void_p]

        self.lib.MetricsParserResultHandleFieldString.restype = c_char_p
        self.lib.MetricsParserResultHandleFieldString.argtypes = [c_void_p]

        self.lib.MetricsParserResultHandleMetricString.restype = c_char_p
        self.lib.MetricsParserResultHandleMetricString.argtypes = [c_void_p]

        self.lib.MetricsParserStringFree.restype = None
        self.lib.MetricsParserStringFree.argtypes = [c_char_p]

        self.lib.MetricsParserResultHandleSize.restype = c_ulonglong
        self.lib.MetricsParserResultHandleSize.argtypes = [c_void_p]

        self.lib.MetricsParserResultHandleCopyInto.restype = None
        self.lib.MetricsParserResultHandleCopyInto.argtypes = [c_void_p, c_ulonglong_p, c_double_p]
        
        self.lib.MetricsParserResultHandleFree.restype = None
        self.lib.MetricsParserResultHandleFree.argtypes = [c_void_p]        

        self.lib.MetricsParserBytesParse.restype = c_void_p
        self.lib.MetricsParserBytesParse.argtypes = [c_char_p, c_char_p,
                                                     c_char_p, c_ulonglong, c_ulonglong, c_ulonglong]

        self.lib.MetricsParserBytesResultHandleAdvance.restype = c_bool
        self.lib.MetricsParserBytesResultHandleAdvance.argtypes = [c_void_p]

        self.lib.MetricsParserBytesResultHandleSize.restype = c_ulonglong
        self.lib.MetricsParserBytesResultHandleSize.argtypes = [c_void_p]

        self.lib.MetricsParserBytesResultHandleBufferSize.restype = c_ulonglong
        self.lib.MetricsParserBytesResultHandleBufferSize.argtypes = [c_void_p, c_ulonglong]

        self.lib.MetricsParserBytesResultHandleCopyInto.restype = None
        self.lib.MetricsParserBytesResultHandleCopyInto.argtypes = [c_void_p, c_ulonglong_p, c_byte_p]
        
        self.lib.MetricsParserBytesResultHandleFree.restype = None
        self.lib.MetricsParserBytesResultHandleFree.argtypes = [c_void_p]

        self.lib.MetricsParserBytesResultHandleFieldString.restype = c_char_p
        self.lib.MetricsParserBytesResultHandleFieldString.argtypes = [c_void_p]

        self.lib.MetricsParserBytesResultHandleMetricString.restype = c_char_p
        self.lib.MetricsParserBytesResultHandleMetricString.argtypes = [c_void_p]

        self.metrics_file = metrics_file

    def ParseBytes(self, metric_regex, fields, min_timestamp=0,
                   max_timestamp=np.iinfo(np.uint64).max,
                   limiting_timestamp=0, blob_processor=None):
        result_handle = self.lib.MetricsParserBytesParse(self.metrics_file,
                                                         metric_regex,
                                                         fields,
                                                         min_timestamp,
                                                         max_timestamp,
                                                         limiting_timestamp)
        return_dict = {}
        while self.lib.MetricsParserBytesResultHandleAdvance(result_handle):
            blob_count = self.lib.MetricsParserBytesResultHandleSize(result_handle)
            # blob_count is the number of objects that the handle
            # holds. Each will have a potentially different size. Will
            # construct one big buffer that can hold all objects and
            # copy them back to back in it.
            blob_sizes = []
            for i in range(blob_count):
                blob_size = self.lib.MetricsParserBytesResultHandleBufferSize(result_handle, i)
                blob_sizes.append(blob_size)
            buf_type = c_byte * sum(blob_sizes)
            buf = buf_type()
            
            timestamps = np.zeros(blob_count, dtype=np.uint64)
            self.lib.MetricsParserBytesResultHandleCopyInto(result_handle,
                                                            timestamps.ctypes.data_as(c_ulonglong_p),
                                                            buf)

            str_p = self.lib.MetricsParserBytesResultHandleFieldString(result_handle)
            field_string = string_at(str_p)
            self.lib.MetricsParserStringFree(str_p)

            str_p = self.lib.MetricsParserBytesResultHandleMetricString(result_handle)
            metric_string = string_at(str_p)
            self.lib.MetricsParserStringFree(str_p)
            
            blobs = []
            offset = 0
            for i in range(blob_count):
                byte_buffer = buffer(buf, offset, blob_sizes[i])
                offset += blob_sizes[i]
                if blob_processor == None:
                    blobs.append(byte_buffer[:])
                else:
                    blobs.append(blob_processor(metric_string, field_string, byte_buffer))
            
            return_dict[(metric_string, field_string)] = (timestamps, blobs)
        return return_dict

    def ParseDistributions(self, metric_regex, fields, min_timestamp=0,
                           max_timestamp=np.iinfo(np.uint64).max,
                           limiting_timestamp=0):
        return self.ParseBytes(metric_regex, fields, min_timestamp,
                               max_timestamp, limiting_timestamp, ProcessDistributionBlob)
        

    def DumpManifest(self):
        str_p = self.lib.MetricsParserManifestSummary(self.metrics_file)
        return_string = string_at(str_p)
        return return_string
    
    def Parse(self, metric_regex, fields, min_timestamp=0,
              max_timestamp=np.iinfo(np.uint64).max,
              limiting_timestamp=0, deltas=False):
        result_handle = self.lib.MetricsParserParse(self.metrics_file,
                                                    metric_regex,
                                                    fields,
                                                    min_timestamp,
                                                    max_timestamp,
                                                    limiting_timestamp)
        return_dict = {}
        while self.lib.MetricsParserResultHandleAdvance(result_handle):
            str_p = self.lib.MetricsParserResultHandleFieldString(result_handle)
            field_string = string_at(str_p)
            self.lib.MetricsParserStringFree(str_p)

            str_p = self.lib.MetricsParserResultHandleMetricString(result_handle)
            metric_string = string_at(str_p)
            self.lib.MetricsParserStringFree(str_p)
            
            array_size = self.lib.MetricsParserResultHandleSize(result_handle)
            timestamps = np.zeros(array_size, dtype=np.uint64)
            values = np.zeros(array_size, dtype=np.double)
            self.lib.MetricsParserResultHandleCopyInto(result_handle,
                                                       timestamps.ctypes.data_as(c_ulonglong_p),
                                                       values.ctypes.data_as(c_double_p))
            if deltas:
                timestamps = timestamps[1:]
                values = values[1:] - values[:-1]

            return_dict[(metric_string, field_string)] = (timestamps, values)
        self.lib.MetricsParserResultHandleFree(result_handle)
        return return_dict

    def TimePlot(self, new_figure, metric_regex, fields,
                 min_timestamp=0,
                 max_timestamp=np.iinfo(np.uint64).max, max_series=10,
                 labels_override=[], plot_type=PLOT_LINE,
                 bin_size=1, **kwargs):
        # Plots a time vs metric plot. Up to max_series combinations
        # of metric fields will be plotted, each with a label derived
        # from label_extractor. The last argument will be passed to
        # plot.
        data = self.Parse(metric_regex, fields,
                          min_timestamp=min_timestamp, max_timestamp=max_timestamp)
        
        # Will first make sure all results belong to the same metric.
        metrics = set(i for i, _ in data.keys())
        assert(len(metrics) == 1)

        xlabel = DictExtract(kwargs, 'xlabel', 'x')
        ylabel = DictExtract(kwargs, 'ylabel', 'y')
        xscale = DictExtract(kwargs, 'xscale', 1.0)
        yscale = DictExtract(kwargs, 'yscale', 1.0)
        if new_figure:
            plt.figure()
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

        num_metrics = len(data)
        if labels_override:
            assert(len(labels_override) == num_metrics)

        i = 0
        for metric_and_label, values in data.items():
            x, y = values
            _, label = metric_and_label
            x = x.astype(np.float64)
            x *= xscale
            y *= yscale
            x, y = Bin(x, y, bin_size)

            if labels_override:
                label = labels_override[i]
            i += 1
            
            if plot_type == self.PLOT_SCATTER:
                plt.scatter(x, y, label=label, **kwargs)
            elif plot_type == self.PLOT_LINE:
                plt.plot(x, y, label=label, linewidth=2, **kwargs)
            elif plot_type == self.PLOT_STEP:
                plt.step(x, y, label=label, linewidth=2, **kwargs)
            else:
                print 'Unknown plot type', plot_type
                return
        plt.legend()

    def CDFPlot(self, new_figure, metric_regex, fields,
                min_timestamp=0, max_timestamp=np.iinfo(np.uint64).max,
                at=0, labels_override=[], plot_type=PLOT_LINE, **kwargs):
        # Plots a CDF of the values of all the matching fields in a
        # metric. The CDF will have a single value per metric -- the
        # one that is the closest in time to (but less than or equal
        # to) the 'at' argument. If 'at' is not provided will
        # plot all values, combining values from different fields.
        data = self.Parse(metric_regex, fields,
                          min_timestamp=min_timestamp, max_timestamp=max_timestamp,
                          limiting_timestamp=at)
        metrics = set(i for i, _ in data.keys())
        assert(len(metrics) == 1)

        xlabel = DictExtract(kwargs, 'xlabel', 'x')
        xscale = DictExtract(kwargs, 'xscale', 1.0)

        if new_figure:
            plt.figure()
            plt.xlabel(xlabel)

        num_metrics = len(data)
        if labels_override:
            assert(len(labels_override) == num_metrics)

        i = 0
        for metric_and_label, values in data.items():
            _, to_plot = values
            _, label = metric_and_label
            to_plot = to_plot.astype(np.float64)
            to_plot *= xscale

            if labels_override:
                label = labels_override[i]
            i += 1
                
            sorted=np.sort(to_plot)
            yvals=np.arange(len(sorted))/float(len(sorted))
            if plot_type == self.PLOT_SCATTER:
                plt.scatter(sorted, yvals, label=label, **kwargs)
            elif plot_type == self.PLOT_LINE:
                plt.plot(sorted, yvals, linewidth=2, label=label, **kwargs)
            else:
                print 'Invalid plot type for CDF', plot_type
                return
        plt.legend()        

def CDFPlot(metrics_files, new_figure, metric_regex, fields, 
            min_timestamp=0, max_timestamp=np.iinfo(np.uint64).max,
            at=0, labels_override=[], **kwargs):
    if new_figure:
        plt.figure()

    if labels_override:
        assert(len(labels_override) == len(metrics_files))
    for i, f in enumerate(metrics_files):
        parser = MetricsParser(f)
        label = f
        if labels_override:
            label = labels_override[i]

        parser.CDFPlot(False, metric_regex, fields, min_timestamp, 
                       max_timestamp, at, [label], **kwargs)

def Parse(metrics_files, metric_regex, fields, min_timestamp=0,
          max_timestamp=np.iinfo(np.uint64).max, at=0):
    return_list = []
    for f in metrics_files:
        parser = MetricsParser(f)
        res = parser.Parse(metric_regex, fields, min_timestamp,
                           max_timestamp, at)
        return_list.append(res)
    return return_list
