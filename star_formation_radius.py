#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

fid = np.genfromtxt("../extracted_data/fid_sf_aux_quants.txt")[:,1]
cflow = np.genfromtxt("../extracted_data/cflow_sf_aux_quants.txt")[:,1]
tctff5 = np.genfromtxt("../extracted_data/tctff5_sf_aux_quants.txt")[:,1]
tctff20 = np.genfromtxt("../extracted_data/tctff20_sf_aux_quants.txt")[:,1]
linrot = np.genfromtxt("../extracted_data/linrot_sf_aux_quants.txt")[:,1]
norot = np.genfromtxt("../extracted_data/norot_sf_aux_quants.txt")[:,1]

cflow = np.where(cflow<0, 0, cflow) # there's one output w/ no new stars

means = np.array([np.mean(fid[40:81]), np.mean(cflow[40:81]), 
      np.mean(tctff5[40:81]), np.mean(tctff20[40:81]), 
      np.mean(linrot[40:81]), np.mean(norot[40:81])])

stds = np.array([np.std(fid[40:81]), np.std(cflow[40:81]), 
      np.std(tctff5[40:81]), np.std(tctff20[40:81]), 
      np.std(linrot[40:81]), np.std(norot[40:81])])

d = {'Mean':means, 'Std Dev':stds, 'Colors':['C0','C3','C2','C1','C4','C5']}

tab = pd.DataFrame(data=d, 
                   index=['Fid','Cflow','TcTff5','TcTff20','LinRot','NoRot'])

tab = tab.sort_values(by='Mean', ascending=False)
print(tab)


# # Significance Tests

# Level: $\alpha=0.05$. If True, you can reject the Null Hypothesis that the means are the same
print('----------------------------------')

print("Cflow larger than TcTff5:",
tab['Mean']['TcTff5'] - tab['Mean']['Cflow'] < -1.645*np.sqrt(tab['Std Dev']['TcTff5']**2/41 + tab['Std Dev']['Cflow']**2/41))

print("Cflow larger than TcTff20:",
tab['Mean']['TcTff20'] - tab['Mean']['Cflow'] < -1.645*np.sqrt(tab['Std Dev']['TcTff20']**2/41 + tab['Std Dev']['Cflow']**2/41))

print("Cflow larger than Fid:",
tab['Mean']['Fid'] - tab['Mean']['Cflow'] < -1.645*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['Cflow']**2/41))

print("Cflow larger than LinRot:",
tab['Mean']['LinRot'] - tab['Mean']['Cflow'] < -1.645*np.sqrt(tab['Std Dev']['LinRot']**2/41 + tab['Std Dev']['Cflow']**2/41))

print("Cflow larger than NoRot:",
tab['Mean']['NoRot'] - tab['Mean']['Cflow'] < -1.645*np.sqrt(tab['Std Dev']['NoRot']**2/41 + tab['Std Dev']['Cflow']**2/41))

print('----------------------------------')

print("TcTff5 larger than Fid:",
tab['Mean']['Fid'] - tab['Mean']['TcTff5'] < -1.645*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['TcTff5']**2/41))

print("TcTff5 larger than TcTff20:",
tab['Mean']['TcTff20'] - tab['Mean']['TcTff5'] < -1.645*np.sqrt(tab['Std Dev']['TcTff20']**2/41 + tab['Std Dev']['TcTff5']**2/41))

print("TcTff5 larger than LinRot:",
tab['Mean']['LinRot'] - tab['Mean']['TcTff5'] < -1.645*np.sqrt(tab['Std Dev']['LinRot']**2/41 + tab['Std Dev']['TcTff5']**2/41))

print("TcTff5 larger than NoRot:",
tab['Mean']['NoRot'] - tab['Mean']['TcTff5'] < -1.645*np.sqrt(tab['Std Dev']['NoRot']**2/41 + tab['Std Dev']['TcTff5']**2/41))

print("TcTff5 not equal to LinRot:",
np.abs(tab['Mean']['TcTff5'] - tab['Mean']['LinRot']) >= 1.96*np.sqrt(tab['Std Dev']['TcTff5']**2/41 + tab['Std Dev']['LinRot']**2/41))

print("TcTff5 not equal to NoRot:",
np.abs(tab['Mean']['TcTff5'] - tab['Mean']['NoRot']) >= 1.96*np.sqrt(tab['Std Dev']['TcTff5']**2/41 + tab['Std Dev']['NoRot']**2/41))

print('----------------------------------')

print("LinRot larger than TcTff20:",
tab['Mean']['TcTff20'] - tab['Mean']['LinRot'] < -1.645*np.sqrt(tab['Std Dev']['LinRot']**2/41 + tab['Std Dev']['TcTff20']**2/41))

print("LinRot larger than Fid:",
tab['Mean']['Fid'] - tab['Mean']['LinRot'] < -1.645*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['LinRot']**2/41))

print("LinRot larger than NoRot:",
tab['Mean']['NoRot'] - tab['Mean']['LinRot'] < -1.645*np.sqrt(tab['Std Dev']['LinRot']**2/41 + tab['Std Dev']['NoRot']**2/41))

print("LinRot not equal to Fid:",
np.abs(tab['Mean']['Fid'] - tab['Mean']['LinRot']) >= 1.96*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['LinRot']**2/41))

print("LinRot not equal to NoRot:",
np.abs(tab['Mean']['LinRot'] - tab['Mean']['NoRot']) >= 1.96*np.sqrt(tab['Std Dev']['LinRot']**2/41 + tab['Std Dev']['NoRot']**2/41))

print('----------------------------------')

print("NoRot larger than TcTff20:",
tab['Mean']['TcTff20'] - tab['Mean']['NoRot'] < -1.645*np.sqrt(tab['Std Dev']['NoRot']**2/41 + tab['Std Dev']['TcTff20']**2/41))

print("NoRot larger than Fid:",
tab['Mean']['Fid'] - tab['Mean']['NoRot'] < -1.645*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['NoRot']**2/41))

print("NoRot larger than LinRot:",
tab['Mean']['LinRot'] - tab['Mean']['NoRot'] < -1.645*np.sqrt(tab['Std Dev']['LinRot']**2/41 + tab['Std Dev']['NoRot']**2/41))

print("NoRot not equal to Fid:",
np.abs(tab['Mean']['Fid'] - tab['Mean']['NoRot']) >= 1.96*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['NoRot']**2/41))

print('----------------------------------')

print("Fid larger than TcTff20:",
tab['Mean']['TcTff20'] - tab['Mean']['Fid'] < -1.645*np.sqrt(tab['Std Dev']['Fid']**2/41 + tab['Std Dev']['TcTff20']**2/41))

print('----------------------------------')

plt.errorbar(tab.index, tab['Mean'], yerr=tab['Std Dev'], color='k',
             ecolor=tab['Colors'], marker='o', ls='none')
plt.minorticks_off()
plt.ylabel('Max Star Formation Radius  [kpc]')

plt.show()

