# AM-DR Energy Model:

This is implementation of our paper "Quality-based and Energy-efficient Data Communication for the Internet of Things (IoT) Networks" that has been accepted *(to appear)* at IEEE Internet of Things (IoT) Journal, 2019. The preprint is available at:

### Citation

*Fathy, Y., Barnaghi, P. "Quality-based and Energy-efficient Data Communication for the Internet of Things (IoT) Networks", IEEE Internet of Things (IoT) Journal, 2019.*

Here are the two cases that are included in the paper. For any question, please e-mail: y.fathy@ucl.ac.uk 


### case_I:


Dataset 1: load_data.m
note: this m file is based on data.mat and its raw data (data.txt)

Dataset 2: load_energydata_complete.m
note: this m file is based on energydata_complete.mat and 
its raw data (energydata_complete.csv)

Dataset 3: AirQualityUCI.m
note: this m file is based on AirQualityUCI.csv file 

Note: all datasets are stored in "dataset" folder

### AM-DR:


The files (AMDR_dual_prediction_dataset1, AMDR_dual_prediction_dataset2,
AMDR_dual_prediction_dataset3) have the same implementation of our AM-DR, 
however, using different datasets. These files are only 
different in the way of presenting the data 
(e.g. axes range, load dataset file, axes labels).

AMDR_dual_prediction_dataset1 uses dataset 1
AMDR_dual_prediction_dataset2 uses dataset 2
AMDR_dual_prediction_dataset3 uses dataset 3



## baseline:

The files (dual_prediction_dataset1, dual_prediction_dataset2,
dual_prediction_dataset3) have the same implementation of the baseline 
approach, however, using different datasets. These files are only 
different in the way of presenting the data 
(e.g. axes range, load dataset file, axes labels).

dual_prediction_dataset1 uses dataset 1
dual_prediction_dataset2 uses dataset 2
dual_prediction_dataset3 uses dataset 3

The results of all these experiements are stored in
plot_results/case_I

------------------------------------------------------------

### case_II:

It includes:
 
-LEACH.m  => original implementation of LEACH routing protocol from
http://csr.bu.edu/sep/

-LEACH_Comments.m => same file of LEACH.m  with some comments

-Model_LEACH_AMDR_Integrate.m => integrate the proposed energy cost model 
of AM-DR with LEACH implementation




The results of all these experiements are stored in
plot_results/case_II

