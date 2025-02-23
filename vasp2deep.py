#!/usr/bin/env python3

import dpdata
import numpy as np

# load data from the OUTCAR file
data = dpdata.LabeledSystem('OUTCAR', fmt='vasp/outcar')
print("# the data contains %d frames" % len(data))

# Define indices specific to the GHz calculation
start_train = 20000
end_train = 30000
start_test = 30000
end_test = 40000

# Select training and testing data
data_training = data.sub_system(range(start_train, end_train))
data_testing = data.sub_system(range(start_test, end_test))


# all training data put into directory:"training_data"
data_training.to_deepmd_npy("training_data")

# all validation data put into directory:"validation_data"
data_testing.to_deepmd_npy("validation_data")

print("# the training data contains %d frames" % len(data_training))
print("# the validation data contains %d frames" % len(data_testing))
