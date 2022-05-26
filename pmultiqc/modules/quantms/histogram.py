from collections import OrderedDict

class Histogram:
    '''The histogram class categorizes instances by 'description' and 'plot_category', it sets up 
    histogram bins according to 'breaks' for accurate data statistics.

    :param description: Declare the objects counted by the instance and use them to build 'cats', 
        a dictionary needed by plotting functions
    :type description: str, optional
    :param plot_category: Two ploting types, 'frequency' refers to count the frequency of different 
        data, and 'range' refers to count the frequency of different ranges of data
    :type plot_category: str, optional
    :param breaks: The list of counted values or range intercepts, defaults to None
    :type breaks: list, optional
    '''
    
    def __init__(self, description, plot_category, breaks = None):
        '''Constructor method. This will also create two dictionaries named 'bins' and 'data', which 
        represent the internal composition of the histogram.
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
        '''Update the value of the corresponding bin according to different plotting methods.

        :param value: Data to be counted
        :type value: int, optional
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
        '''Integrate statistical results and the 'cats' dictionary for ploting into a nested dictionary.

        :param percentage: `True` if calculate the percentage of data representing each bins, defaults 
            to `False`
        :type percentage: bool, optional

        :param cats: A dictionary needed by plotting functions, it contains the name and the description 
            of bins, defaults to None
        :type cats: dict, optional
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
