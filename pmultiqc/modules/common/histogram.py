from collections import OrderedDict

import numpy as np


class Histogram:
    """The histogram class categorizes instances by 'description' and 'plot_category', it sets up
    histogram bins according to 'breaks' for accurate data statistics.

    :param description: Declare the objects counted by the instance and use them to build 'cats',
        a dictionary needed by plotting functions
    :type description: str, optional
    :param plot_category: Two ploting types, 'frequency' refers to count the frequency of different
        data, and 'range' refers to count the frequency of different ranges of data
    :type plot_category: str, optional
    :param stacks: This list contains the types of stacks that apply to all bars
    :type stacks: list, optional
    :param breaks: The list of counted values or range intercepts, defaults to None
    :type breaks: list, optional
    """

    def __init__(self, description, plot_category, stacks=None, breaks=None):
        """Constructor method. This will also create two dictionaries named 'bins' and 'data', which
        represent the internal composition of the histogram.
        """
        # Initialise some parameters
        self.description = description
        self.data = OrderedDict()
        self.stacks = stacks
        self.breaks = breaks
        self.bins = []
        self.cats = OrderedDict()
        self.dict = dict()
        self.dict["data"] = dict()
        self.out_threshold = False

        if plot_category == "frequency":
            self.plot_category = 1
        elif plot_category == "range":
            self.plot_category = 2

        # Define the bins for the histogram, and initialize the data dictionary
        if breaks:
            breaks.sort()
            self.breaks = breaks
            self.threshold = breaks[-1]

            if self.plot_category == 1:
                self.bins.extend([str(i) for i in range(1, len(breaks))])
            elif self.plot_category == 2:
                self.bins.extend(
                    [(str(breaks[i]) + " ~ " + str(breaks[i + 1])) for i in range(len(breaks) - 1)]
                )

            # Whether to add the last bin to collect the remaining data
            if breaks[-1] != float("inf"):
                self.bins.append(">= " + str(breaks[-1]))

            # Initiate stacks in every bar
            if stacks and self.plot_category == 2:
                for i in self.bins:
                    self.data[i] = dict.fromkeys(stacks, 0)
            else:
                for i in self.bins:
                    self.data[i] = {"total": 0}
        else:
            pass

    def add_value(self, value, stack="total"):
        """Update the value of the corresponding bin according to different plotting methods.

        :param value: Data to be counted
        :type value: int, optional
        :param stack: The stack of bars to which the data belongs
        :type stack: str, optional
        """
        if value is None:
            return

        if self.plot_category == 1:
            if self.breaks:
                if value < self.threshold:
                    self.data[str(value)][stack] += 1
                else:
                    self.out_threshold = True if value > self.threshold else False
                    self.data[self.bins[-1]][stack] += 1

            else:
                if str(value) in self.bins:
                    self.data[str(value)][stack] += 1
                else:
                    self.bins.append(str(value))
                    if self.stacks:
                        self.data[str(value)] = dict.fromkeys(self.stacks, 0)
                        self.data[str(value)][stack] = 1
                    else:
                        self.data[str(value)] = {stack: 1}
                # sort
                if not isinstance(value, str) and not self.stacks:
                    data_keys = [type(value)(i) for i in self.data.keys()]
                    data_keys.sort()
                    self.data = OrderedDict(
                        {str(i): {"total": self.data[str(i)]["total"]} for i in data_keys}
                    )

        elif self.plot_category == 2:
            # Note: breaks are already sorted during __init__
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
                self.data[self.bins[left]][stack] += 1
            else:
                self.out_threshold = True if value > self.threshold else False
                self.data[self.bins[-1]][stack] += 1

    def add_values_batch(self, values, stack="total"):
        """Batch add multiple values to the histogram (optimized for performance).

        :param values: Array or Series of values to be counted
        :type values: array-like
        :param stack: The stack of bars to which the data belongs
        :type stack: str, optional
        """
        # Convert to numpy array and filter out NaN/None values
        if hasattr(values, 'values'):  # pandas Series
            arr = values.values
        else:
            arr = np.asarray(values)

        # Filter out NaN values
        if arr.dtype.kind == 'f':  # float array
            mask = ~np.isnan(arr)
            arr = arr[mask]
        elif arr.dtype.kind == 'O':  # object array
            mask = np.array([x is not None and (not isinstance(x, float) or not np.isnan(x)) for x in arr])
            arr = arr[mask]

        if len(arr) == 0:
            return

        if self.plot_category == 1:
            # Frequency mode
            if self.breaks:
                below_threshold = arr < self.threshold
                below_values = arr[below_threshold]
                above_count = np.sum(~below_threshold)

                # Count values below threshold
                for val in below_values:
                    val_str = str(int(val) if isinstance(val, (int, np.integer, float)) and float(val).is_integer() else val)
                    if val_str in self.data:
                        self.data[val_str][stack] += 1

                # Add above threshold count
                if above_count > 0:
                    self.out_threshold = True
                    self.data[self.bins[-1]][stack] += above_count
            else:
                # No breaks - use value_counts approach
                unique, counts = np.unique(arr, return_counts=True)
                for val, count in zip(unique, counts):
                    val_str = str(val)
                    if val_str in self.bins:
                        self.data[val_str][stack] += count
                    else:
                        self.bins.append(val_str)
                        if self.stacks:
                            self.data[val_str] = dict.fromkeys(self.stacks, 0)
                            self.data[val_str][stack] = count
                        else:
                            self.data[val_str] = {stack: count}

                # Sort if needed (only for non-string values)
                if len(arr) > 0 and not isinstance(arr[0], str) and not self.stacks:
                    try:
                        data_keys = [type(arr[0])(i) for i in self.data.keys()]
                        data_keys.sort()
                        self.data = OrderedDict(
                            {str(i): {"total": self.data[str(i)]["total"]} for i in data_keys}
                        )
                    except (ValueError, TypeError):
                        pass  # Keep original order if sorting fails

        elif self.plot_category == 2:
            # Range mode - use numpy searchsorted for vectorized bin assignment
            arr_float = arr.astype(float)
            below_threshold = arr_float < self.threshold
            below_values = arr_float[below_threshold]
            above_count = np.sum(~below_threshold)

            if len(below_values) > 0:
                # Use searchsorted for vectorized binary search
                # searchsorted returns index where value would be inserted to maintain order
                # We use side='right' and subtract 1 to get the correct bin
                bin_indices = np.searchsorted(self.breaks, below_values, side='right') - 1
                bin_indices = np.clip(bin_indices, 0, len(self.bins) - 2)  # Ensure valid indices

                # Count occurrences in each bin
                unique_bins, bin_counts = np.unique(bin_indices, return_counts=True)
                for bin_idx, count in zip(unique_bins, bin_counts):
                    self.data[self.bins[bin_idx]][stack] += count

            if above_count > 0:
                self.out_threshold = True
                self.data[self.bins[-1]][stack] += above_count

    def to_dict(self, percentage=False, cats=None):
        """Integrate statistical results and the 'cats' dictionary for ploting into a nested dictionary.

        :param percentage: `True` if calculate the percentage of data representing each bins, defaults
            to `False`
        :type percentage: bool, optional
        :param cats: A dictionary needed by plotting functions, it contains the name and the description
            of bins, defaults to None
        :type cats: dict, optional
        """
        if not self.stacks:
            for i in self.data.keys():
                self.data[i] = self.data[i]["total"]

        # Modify the last bin
        if not self.out_threshold and self.plot_category == 2:
            last_bin = self.bins[-1]
            self.bins[-1] = last_bin.strip(">")
            self.data[last_bin.strip(">")] = self.data.pop(last_bin)

        if percentage:
            total = sum(self.data.values())
            self.dict["data"]["frequency"] = OrderedDict()
            self.dict["data"]["percentage"] = OrderedDict()
            for key in self.data:
                self.dict["data"]["frequency"][key] = OrderedDict()
                self.dict["data"]["frequency"][key]["Frequency"] = self.data[key]
                self.dict["data"]["percentage"][key] = OrderedDict()
                self.dict["data"]["percentage"][key]["Percentage"] = round(
                    self.data[key] / total * 100, 2
                )
        else:
            self.dict["data"] = self.data

        # Categorys of bins for ploting functions
        if cats:
            self.cats = cats * len(self.bins) if self.stacks else cats
        else:
            for i in self.bins:
                self.cats[i] = dict()
                self.cats[i]["name"] = i

        self.dict["cats"] = self.cats
