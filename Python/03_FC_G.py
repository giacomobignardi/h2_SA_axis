#Author: Giacomo Bignardi
#Partially adapted from: Tutorial from Brainspace https://brainspace.readthedocs.io/en/latest/pages/getting_started.html
#Description: Reduce FCs computed on the tSchaefer 400 parcelation to the first 10 eigenvectors  
#Program: compute individual gradients ------------------------------------------------------------------------------------------------------------------------------
#load pyhton venv 
reticulate::py_config()

# import relevant packages
import pandas as pd
import numpy as np
import glob
from brainspace_modified.gradient import GradientMaps

# set Open Access working directories
wdOA_scripts = "/02_scripts"
wdOA_ImageOutput = "/05_Figures"

# set not Open Access working directories
wdNOA = os.getcwd()
wdNOA_Data = "/01_input"
wdNOA_output = "/03_outputs/processedData"
wdNOA_ImageOutput = "/05_Figures"
wdNOA_rawData = os.getcwd()[0:(len(os.getcwd())-len("/04_analysis_OA"))] + "/03_rawData/private"

# GROUP LEVEL####
# load mean template
FC_mean = pd.read_csv(wdNOA + wdNOA_output + "/02_FC_m.csv")
FC_gradients = pd.DataFrame()
# SUP: GRADIENT---------------------------------------
# calculate percentage of variance
GtFC = GradientMaps(n_components = 10,kernel ='normalized_angle', approach ='dm',random_state=42)
GtFC.fit(FC_mean.values[0:400])
GtFC.lambdas_[0]  / GtFC.lambdas_.sum()
(GtFC.lambdas_[0] + GtFC.lambdas_[1] + GtFC.lambdas_[2])   / GtFC.lambdas_.sum()
#---------------------------------------
# Save Gradient
pd.DataFrame(GtFC.gradients_).to_csv(wdNOA + wdNOA_output + "/03_FC_m_G.csv",index=False)

# INDIVIDUAL LEVEL####
# Pick ids for individual level gradients
ids = pd.read_csv(wdNOA + wdNOA_Data + "/subjects_for_giaco.csv",header=None)
# Load subjects IDs for sessions averaged within day
files_d1 = glob.glob(wdNOA_rawData  + "/rsfMRIavg/day1/*.csv")
files_d2 = glob.glob(wdNOA_rawData  + "/rsfMRIavg/day2/*.csv")
files = glob.glob(wdNOA_rawData  + "/rsfMRIavg/day_avg/*.csv")

# create a list of file to be extracted
sub_name_d1 = [] #day1
for i in list(range(len(files_d1))):
  sub_name_d1.append(files_d1[i].rsplit('/', 1)[-1].rsplit('.', 1)[0]) #split directory to the last /; select what last element; split directory to the last.; select the first element; append to an array
sub_name_d1.sort()

sub_name_d2 = [] #day2
for i in list(range(len(files_d2))):
  sub_name_d2.append(files_d2[i].rsplit('/', 1)[-1].rsplit('.', 1)[0]) #split directory to the last /; select what last element; split directory to the last.; select the first element; append to an array
sub_name_d2.sort()

sub_name = [] #average across days
for i in list(range(len(files))):
  sub_name.append(files[i].rsplit('/', 1)[-1].rsplit('.', 1)[0]) #split directory to the last /; select what last element; split directory to the last.; select the first element; append to an array
sub_name.sort()

sub_name_d1 = [int(x) for x in sub_name_d1]
sub_name_d2 = [int(x) for x in sub_name_d2]
sub_name = [int(x) for x in sub_name]

seti = list(ids[0])
sub_d1 = list(set(seti).intersection(list(map(int, sub_name_d1)))) 
sub_d1.sort()
sub_d2 = list(set(seti).intersection(list(map(int, sub_name_d2)))) 
sub_d2.sort()
sub = list(set(seti).intersection(list(map(int, sub_name)))) 
sub.sort()

# GRADIENT DAY 1---------------------------------------
# run diffusion map embedding for FC computed at day1
FCG_d1 = pd.DataFrame()
FCG_d1_er = pd.DataFrame()
# list(range(len(sub_name_1RL)))
for i in list(range(len(sub_d1))):
    #read subj FC matrix for different runs
  i_d1 = pd.DataFrame(); i_FC_d1 = pd.DataFrame();
  i_d1 = pd.read_csv(f"{wdNOA_rawData}/rsfMRIavg/day1/{sub_d1[i]}.csv", header=None)
  i_FC_d1 = i_d1.loc[0:399,0:399]
#######Compute indvidual gradientes and align to the template#######---------------------------------------  
  if (i_FC_d1.isnull().values.any() == False):
    G_FC_d1 = GradientMaps(n_components = 10,kernel ='normalized_angle', approach ='dm', alignment='procrustes',random_state=42)
    G_FC_d1.fit([FC_mean.values[0:400],i_FC_d1.values[0:400]]) 
    G_FC_d1_er = G_FC_d1.lambdas_[1][0]/ G_FC_d1.lambdas_[1].sum()
  else:
    G_FC_d1.aligned_ = (np.array([[np.nan]* 10]*400),np.array([[np.nan]* 10]*400))
    G_FC_d1_er = np.nan
  ##############---------------------------------------  
  FC_i_d1 = pd.DataFrame()
  FC_i_d1 = pd.DataFrame(G_FC_d1.aligned_[1])
  FC_i_d1 = pd.concat((pd.DataFrame([sub_d1[i]] * 400, columns = ['Sub']), pd.DataFrame([1.5] * 400, columns = ['Session']), pd.DataFrame(list(range(1,401,1)), columns = ['Parcel']),FC_i_d1), axis = 1)
  FCG_d1 = FCG_d1.append(FC_i_d1)
  FC_i_d1_er = pd.DataFrame()
  FC_i_d1_er = pd.DataFrame([G_FC_d1_er], columns = ['explanatory_ratio'])
  FC_i_d1_er = pd.concat((pd.DataFrame([sub_d1[i]], columns = ['Sub']),FC_i_d1_er), axis = 1)
  FCG_d1_er = FCG_d1_er.append(FC_i_d1_er)

# note  FCG_d1_er (explanatory ratio) might be needed for further follow up projects on cognitive-behavioral differences

# SAVE GRADIENTS Day1#### 
FCG_d1.iloc[:,0:6].to_csv(wdNOA + wdNOA_output + "/03_GFC_i_d1.csv",index=False)
FCG_d1_er.to_csv(wdNOA + wdNOA_output + "/03_GFC_i_er_d1.csv",index=False)

# GRADIENT DAY 2---------------------------------------
# run diffusion map embedding for FC computed at day2
FCG_d2 = pd.DataFrame()
FCG_d2_er = pd.DataFrame()
# list(range(len(sub_name_1RL)))
for i in list(range(len(sub_d2))):
    #read subj FC matrix for different runs
  i_d2 = pd.DataFrame(); i_FC_d2 = pd.DataFrame();
  i_d2 = pd.read_csv(f"{wdNOA_rawData}/rsfMRIavg/day2/{sub_d2[i]}.csv", header=None)
  i_FC_d2 = i_d2.loc[0:399,0:399]
#######Compute indvidual gradientes and align to the template#######---------------------------------------  
  if (i_FC_d2.isnull().values.any() == False):
    G_FC_d2 = GradientMaps(n_components = 10,kernel ='normalized_angle', approach ='dm', alignment='procrustes',random_state=42)
    G_FC_d2.fit([FC_mean.values[0:400],i_FC_d2.values[0:400]])
    G_FC_d2_er = G_FC_d2.lambdas_[1][0]/ G_FC_d2.lambdas_[1].sum()
  else:
    G_FC_d2.aligned_ = (np.array([[np.nan]* 10]*400),np.array([[np.nan]* 10]*400))
    G_FC_d2_er = np.nan
  ##############---------------------------------------  
  FC_i_d2 = pd.DataFrame()
  FC_i_d2 = pd.DataFrame(G_FC_d2.aligned_[1])
  FC_i_d2 = pd.concat((pd.DataFrame([sub_d2[i]] * 400, columns = ['Sub']), pd.DataFrame([2.5] * 400, columns = ['Session']), pd.DataFrame(list(range(1,401,1)), columns = ['Parcel']),FC_i_d2), axis = 1)
  FCG_d2 = FCG_d2.append(FC_i_d2)
  FC_i_d2_er = pd.DataFrame()
  FC_i_d2_er = pd.DataFrame([G_FC_d2_er], columns = ['explanatory_ratio'])
  FC_i_d2_er = pd.concat((pd.DataFrame([sub_d2[i]], columns = ['Sub']),FC_i_d2_er), axis = 1)
  FCG_d2_er = FCG_d2_er.append(FC_i_d2_er)

# SAVE GRADIENTS Day2####  
FCG_d2.iloc[:,0:6].to_csv(wdNOA + wdNOA_output + "/03_GFC_i_d2.csv",index=False)
FCG_d2_er.to_csv(wdNOA + wdNOA_output + "/03_GFC_i_er_d2.csv",index=False)

# GRADIENT average across days---------------------------------------
# run diffusion map embedding for FC averaged across days of scanning session
FCG = pd.DataFrame()
# list(range(len(sub_name_1RL)))
for i in list(range(len(sub))):
    #read subj FC matrix for different runs
  i_df = pd.DataFrame(); i_FC = pd.DataFrame();
  i_df = pd.read_csv(f"{wdNOA_rawData}/rsfMRIavg/day_avg/{sub[i]}.csv", header=None)
  i_FC = i_df.loc[0:399,0:399]
#######Compute indvidual gradientes and align to the template#######---------------------------------------  
  if (i_FC.isnull().values.any() == False):
    G_FC = GradientMaps(n_components = 10,kernel ='normalized_angle', approach ='dm', alignment='procrustes',random_state=42)
    G_FC.fit([FC_mean.values[0:400],i_FC.values[0:400]]) 
  else:
    G_FC.aligned_ = (np.array([[np.nan]* 10]*400),np.array([[np.nan]* 10]*400))
  ##############---------------------------------------  
  FC_i = pd.DataFrame()
  FC_i = pd.DataFrame(G_FC.aligned_[1])
  FC_i = pd.concat((pd.DataFrame([sub[i]] * 400, columns = ['Sub']), pd.DataFrame([0] * 400, columns = ['Session']), pd.DataFrame(list(range(1,401,1)), columns = ['Parcel']),FC_i), axis = 1)
  FCG = FCG.append(FC_i)

# SAVE GRADIENTS####  
FCG.iloc[:,0:6].to_csv(wdNOA + wdNOA_output + "/03_GFC_i.csv",index=False)
