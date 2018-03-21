# Project: Kinematics Pick & Place

## Kinematic Analysis

![DH.png](figures/DH.png)

Accordingly, we find the DH parameters of the KUKA arm to be as follows:

&alpha;   | a      | d    | &theta;
:--------:|:------:|:----:|:-------:
0         | 0      | 0.75 | q&#8321;
-&pi;/2   | 0.35   | 0    | q&#8322;-&pi;/2
0         | 1.25   | 0    | q&#8323;
-&pi;/2   | -0.054 | 1.50 | q&#8324;
&pi;/2    | 0      | 0    | q&#8325;
-&pi;/2   | 0      | 0    | q&#8326;
0         | 0      | 0.303| 0

These values are defined [here](kuka_kin.py#112).

## Project Implementation