In this repository are the parts of DeepJET that are not CMS specific and live outside CMS software.

=====================

The training has 2 steps outside CMSSW

0) Set up the environment. For details, please see environment/README

1) Convert the root tuples to the correct data format:

   cd convertFromRoot
   ./convertFromRoot.py -i /path/to/your/root/ntuples/train_val_samples.txt -o /directory/to/store/output/must/not/yet/exist -c data_format
   
   To store the output, it is recommended to use eos. The data format is defined by the classes TrainData_XX.py in the modules directory.
   For a first example, please use: TrainData_deepCSV_ST
   Ntuples to start with can be found at: /eos/cms/store/cmst3/group/dehep/DeepJet/NTuples
   
   The conversion will run for a while. Get yourself a coffee or two.
   If the conversion is interrupted, you can pick up where you left by:
   ./convertFromRoot.py -r /directory/you/stored/the/output/in/snapshot.dc -c TrainData_deepCSV_ST
   
   
   
   
3) Train the DNN:
   
   cd Train
   python DeepJetTrain_dense.py /directory/you/stored/the/output/in/dataCollection.dc /output/directory/must/not/exist/yet
   
   This will take a while. The output will not be big, so the output directory can be on the afs or somewhere else with low space
   
