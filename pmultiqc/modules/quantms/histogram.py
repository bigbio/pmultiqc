from collections import OrderedDict

class Histogram:
    '''
    The histogram class categorizes instances by 'description' and 'plot_category', it sets up histogram bins 
    according to 'breaks' for accurate data statistics.

    ### __init__(self, description, plot_category, breaks = None)  
    Initialize parameters that include descriptions, data dictionaries, breaks, bins, cats, and nested 
    dictionaries generated from them.
    * self.description --> Declare the objects counted by the instance and use them to build self.cats.
    * self.plot_category --> Two kinds of histograms. One is to count the frequency of different data, and 
    the other is to count the frequency of different ranges of data.
    * self.cats --> A dictionary containing two keys needed by plotting functions, key 'name' representing 
    the name of bin and key 'description' representing the description of bin.
    * self.data --> A dictionary used to store statistics. The keys are bins and the values are statistics.
    * self.breaks --> The list of counted values or range intercepts, default to None.
    * self.bins --> The list of bins that represent the internal composition of the histogram.
    * self.dict --> A nested dictionary that the key 'data' corresponds to data statistics and the key 'cats' 
    corresponds to the dictionary used for plotting functions.

    ### addValue(self, value)  
    Update the value of the corresponding bin according to different plotting categorys.
    * self.threshold --> Any data greater than or equal to the threshold is counted into the last bin.

    ### to_dict(self, percentage = False, cats = None)
    Integrate statistical results and the dictionary 'cats' for ploting into a nested dictionary.
    '''
    def __init__(self, description, plot_category, breaks = None):
        '''
        Initialize parameters that include descriptions, data dictionaries, breaks, bins, cats, and nested 
        dictionaries generated from them.

        * description --> Declare the objects counted by the instance and use them to build A dictionary 
        containing two keys needed by plotting functions.

        * plot_category --> There are two ploting types, 'frequency' and 'range'. One is to count the frequency 
        of different data, and the other is to count the frequency of different ranges of data. 

        * breaks --> The list of counted values or range intercepts, default to None.
        '''
        # Initialise some parameters
        self.description = description
        self.data = dict()
        self.breaks = breaks
        self.bins = []
        self.cats = OrderedDict()
        self.dict = dict()
        self.dict['data'] = dict()

        if plot_category == 'frequency':
            self.plot_category = 1
        elif plot_category == 'range':
            self.plot_category = 2

        # Define the bins for the histogram, and initialize the data dictionary
        if breaks:
            breaks.sort()
            self.breaks = breaks
            self.threshold = breaks[-1]
        
            if self.plot_category == 1:
                self.bins.extend([str(i) for i in range(1, len(breaks))])
            elif self.plot_category == 2:
                self.bins.extend([(str(breaks[i]) + '-' + str(breaks[i+1]))
                                  for i in range(len(breaks) - 1)])

            self.bins.append('>= ' + str(breaks[-1]))
            for i in self.bins:
                self.data[i] = 0
        else: pass

    def addValue(self, value):
        '''
        Update the value of the corresponding bin according to different plotting methods.

        * value --> Statistics.
        '''

        if self.plot_category == 1:
            if self.breaks:
                if value < self.threshold:
                    self.data[str(value)] += 1
                else:
                    self.data[self.bins[-1]] += 1
            else:
                if value in self.bins:
                    self.data[value] += 1
                else:
                    self.data[value] = 1
                    self.bins.append(value)

        elif self.plot_category == 2:
            self.breaks.sort()
            if value < self.threshold:
                left, right = 0, len(self.breaks) - 1
                while left < right:
                    mid = (left + right + 1) // 2
                    if value >= self.breaks[mid]:
                        left = mid
                    else:
                        right = mid - 1
                # Each value in self.breaks is to the left of the range 
                # e.g. [left, right) in self.bins(same index)
                self.data[self.bins[left]] += 1
            else:
                self.data[self.bins[-1]] += 1


    def to_dict(self, percentage = False, cats = None):
        '''
        Integrate statistical results and the 'cats' dictionary for ploting into a nested dictionary.

        * percentage --> Whether to calculate the percentage of data representing each bins. Default is False, 
        and when True, dict['data'] is a nested dictionary of two keys, 'frequency' and 'percentage'. 
        According to the frequency of different bins in the histogram, the corresponding percentage is 
        generated and put into the nested dictionary dict['data']['percentage'].

        * cats --> A dictionary containing two keys needed by plotting functions, key 'name' representing the 
        name of bin and key 'description' representing the description of bin. Default is None.
        '''

        if percentage:
            total = sum(self.data.values())
            self.dict['data']['frequency'] = OrderedDict()
            self.dict['data']['percentage'] = OrderedDict()
            keys = sorted(self.data)  # sort
            for key in keys:
                self.dict['data']['frequency'][key] = OrderedDict()
                self.dict['data']['frequency'][key]['Frequency'] = self.data[key]
                self.dict['data']['percentage'][key] = OrderedDict()
                self.dict['data']['percentage'][key]['Percentage'] = round(self.data[key] /
                                                                           total * 100, 2)
        else:
            self.dict['data'] = self.data

        # Categorys of bins for ploting functions
        if cats:
            self.cats = cats
        else:
            for i in self.bins:
                self.cats[i] = dict()
                self.cats[i]['name'] = i
                if self.plot_category == 1:
                    self.cats[i]['description'] = self.description + ' ' + \
                        (self.bins[-1] if i == self.bins[-1] else 'is ' + i)
                elif self.plot_category == 2:
                    self.cats[i]['description'] = self.description + ' ' + \
                        (self.bins[-1] if i == self.bins[-1] else 'is between ' + \
                        i.split('-')[0] + ' and ' + i.split('-')[1])

        self.dict['cats'] = self.cats