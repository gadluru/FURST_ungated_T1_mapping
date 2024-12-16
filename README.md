# Free-breathing Ungated Radial Simultaneous multi-slice cardiac T1 mapping

This repository contains code and example datasets for the paper 'Free-breathing Ungated Radial Simultaneous multi-slice cardiac T1 mapping'. 
Running Recon_UCAIR_ungated.m will reproduce some of the results demonstrated in the paper. This code was developed and tested on a RockLinux 8.6 operating system, 
with an AMD EPYC Milan 7543 32 2.8GHz 256 MB cache, 512 GB RAM, and Nvidia  A100 gpus. Systems with ~64 GB RAM may counter issues with insufficient memory 
that require code changes to reduce memory requirements during certain function calls. 

Le JV, Mendes JK, Sideris K, Bieging E, Carter S, Stehlik J, DiBella EVR, Adluru G. Free-Breathing Ungated Radial Simultaneous Multi-Slice Cardiac T1 Mapping. J Magn Reson Imaging. 2024 Dec 11. doi: 10.1002/jmri.29676. PMID: 39661447.

<br />
<br />

<p align="center" width="60%">
    <img width="60%" src="https://github.com/user-attachments/assets/a0f06f99-a8e8-4e02-bf6a-44adf27cfc34">
</p>
Figure 4. (A) Correlation and (B) Bland-Altman plots for myocardial pre-contrast T1 mapping with FURST and MOLLI in 13 canine subjects and 7 human subjects. The FURST and MOLLI pre-contrast myocardial T1 were significantly different (1212±128 ms and 1179±87 ms, respectively). (C) Correlation and (D) Bland-Altman plots for myocardial post-contrast T1 mapping with FURST and MOLLI. The FURST and MOLLI post-contrast myocardial T1 were significantly different (508±97 ms and 516±92 ms, respectively). (E) Correlation and (F) Bland-Altman plots for ECV mapping with FURST and MOLLI. There was no significant difference between FURST and MOLLI myocardial ECV (29±11 % and 28±11 %, respectively, p=0.05). The ECV mean difference was 0.48, with 95% CI:(6.0×10^(-4),0.96). Each dot (n=296) represents an averaged T1 value for an AHA segment of the myocardium. 

<br />
<br />

<p align="center" width="60%">
    <img width="60%" src="https://github.com/user-attachments/assets/57e96603-7711-4206-91b3-405e1393a189">
</p>
Figure 6. Pre-contrast T1 mapping results of FURST using (A) an all-data reconstruction to produce a cardiac phase averaged T1 map, (B) a systolic binning reconstruction to produce a systolic T1 map, and (C) a diastolic binning reconstruction to produce a diastolic T1 map. The FURST averaged phase myocardial T1 was 1260±31 ms. There was no significant difference between FURST and MOLLI segment-wise systolic pre-contrast myocardial T1 (1170±35 ms and 1205±46 ms, respectively, p=0.12). There were significant differences between FURST and MOLLI segment-wise diastolic pre-contrast myocardial T1 (1322±61 ms and 1219±41 ms, respectively).

<br />
<br />

Contact:

Johnathan Le 

le.johnv@outlook.com

Ganesh Adluru

gadluru@gmail.com
