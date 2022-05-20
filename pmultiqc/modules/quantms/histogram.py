from collections import OrderedDict, defaultdict

class Histogram:

    def classify(self, plot_cat, with_remains = False, data_cats = None):
        # with_remains -> Whether there are any remaining values outside the statistical range.
        # data_cats -> Make statistics and display according to different data categories.
        # plot_cat = 'frequency' -> Count the number of occurrences of different values.
        # plot_cat = 'range' -> Count the frequency with which a value occurs over several consecutive ranges.
        # plot_cat = 'record' -> Assigns values in order according to the different ranges provided.

        self.data = dict()
        self.data_cats = data_cats
        if plot_cat == 'frequency':
            if with_remains:
                self.plot_cat = 2
            else:
                self.plot_cat = 1
                self.frequency = dict()
        elif plot_cat == 'range':
            self.plot_cat = 3
        elif plot_cat == 'record':
            self.plot_cat = 4

    def setBreaks(self, breaks = None):

        if breaks:
            self.breaks = breaks
        else:
            self.breaks = []
            self.bins = []

        if self.plot_cat == 2 or self.plot_cat == 4:
            self.bins = [str(i) for i in range(1, len(breaks))]
            # if self.with_remains = True
            self.last_bin = '>= ' + str(self.breaks[-1])
            self.bins.append(self.last_bin)
        elif self.plot_cat == 3:
            self.bins = [(str(breaks[i]) + '-' + str(breaks[i+1])) for i in range(len(breaks) - 1)]
            # if self.with_remains = True
            self.last_bin = '>= ' + str(self.breaks[-1])
            self.bins.append(self.last_bin)
            
        if self.data_cats:
            for i in self.data_cats:
                self.data[i] = dict()
                for j in self.bins:
                    self.data[i][j] = None if self.plot_cat == 4 else 0

    def addValue(self, value, total, data_cat = None):
        self.total = total
        self.data_cat = data_cat
        
        if self.plot_cat == 3:
            self.add_compareValue(value)
        elif self.plot_cat == 1 or self.plot_cat == 2:
            self.add_countValue(value)
        else:
            for i in self.data[self.data_cat].keys():
                if self.data[self.data_cat][i] == None:
                    self.data[self.data_cat][i] = value
                    break


    def add_compareValue(self, value):
        if value >= self.breaks[-1]:
            bin = self.last_bin
        for i in range(len(self.breaks) - 1):
            left, right = self.breaks[i], self.breaks[i + 1]
            if left <= value < right:
                bin = str(left) + '-' + str(right)
        self.data[self.data_cat][bin] += 1

    def add_countValue(self, value):
        if self.plot_cat == 2:
            if value in self.breaks[0: -1]:
                bin = str(value)
            else:
                bin = self.last_bin
            self.data[self.data_cat][bin] += 1
            
        else:
            if value in self.breaks:
                self.frequency[value]['Frequency'] += 1
            else:
                self.bins.append(str(value))
                self.breaks.append(value)
                self.breaks.sort()
                self.frequency[value] = dict()
                self.frequency[value]['Frequency'] = 1

    def add_cats(self, cats = None):
        if cats:
            self.cats = cats
        else:
            self.cats = OrderedDict()
            for i in self.bins:
                self.cats[i] = dict()
                self.cats[i]['name'] = i
                self.cats[i]['description'] = ''

    def to_dict(self):
        self.dict = dict()
        if self.plot_cat == 1:
            self.dict['data'] = dict()
            self.dict['data']['frequency'] = OrderedDict()
            self.dict['data']['percentage'] = OrderedDict()
            keys = sorted(self.frequency)  # sort
            for key in keys:
                self.dict['data']['frequency'][key] = defaultdict(int)
                self.dict['data']['frequency'][key]['Frequency'] = self.frequency[key]['Frequency']
                self.dict['data']['percentage'][key] = defaultdict(float)
                self.dict['data']['percentage'][key]['Percentage'] = round(self.frequency[key]['Frequency'] / \
                                                                self.total * 100, 2)
        else:
            self.dict['data'] = self.data

        self.dict['bins'] = self.bins
        self.dict['cats'] = self.cats
